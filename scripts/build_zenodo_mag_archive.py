#!/usr/bin/env python3
"""Build and validate the EHR02 MAG tarball for Zenodo.

The source index must contain the tab-separated columns ``genome_id``,
``filename`` and ``fasta_url``. Source files are downloaded concurrently and
reconciled within each assembly by exact genome length before being renamed to
``<genome_id>.fa.gz`` inside the archive, matching EHR01. EHR02 genome lengths
are unique within every assembly, so this also detects and repairs Core link
rotations deterministically. The final archive is streamed directly to disk so
the build needs only one archive's worth of free space.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


EXPECTED_MAG_COUNT = 11_538
READ_CHUNK = 1024 * 1024
USER_AGENT = "EHI-data-release-2-Zenodo-preparation/1.0"
# ERDA agrees with an earlier EHI metadata export for this record, while the
# current Core row contains a different MAG's derived statistics. The tuple is
# Total sequence bp. Keep it exact so further regressions still stop.
KNOWN_SIGNATURE_OVERRIDES = {
    "EHM021731": (2_117_230,),
}


@dataclass(frozen=True)
class Source:
    genome_id: str
    filename: str
    url: str

    @property
    def archive_name(self) -> str:
        return f"{self.genome_id}.fa.gz"


@dataclass(frozen=True)
class Download:
    source: Source
    payload: bytes
    fasta_size_bytes: int
    sequence_bp: int
    contigs: int
    n50: int
    md5: str
    sha256: str

    @property
    def signature(self) -> tuple[int]:
        return (self.sequence_bp,)


@dataclass(frozen=True)
class GenomeMetadata:
    assembly_id: str
    sequence_bp: int
    contigs: int
    n50: int

    @property
    def signature(self) -> tuple[int]:
        return (self.sequence_bp,)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_index", type=Path)
    parser.add_argument(
        "--genome-metadata",
        type=Path,
        default=Path("data/EHR02_genome.tsv"),
        help="TSV whose genome_id column defines the authoritative release set",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("zenodo/ehi_data_release_2_mags.tar.gz"),
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=6)
    return parser.parse_args()


def read_sources(path: Path) -> list[Source]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome_id", "filename", "fasta_url"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        sources = [
            Source(
                genome_id=(row["genome_id"] or "").strip(),
                filename=(row["filename"] or "").strip(),
                url=(row["fasta_url"] or "").strip(),
            )
            for row in reader
        ]

    if len(sources) != EXPECTED_MAG_COUNT:
        raise ValueError(
            f"source index has {len(sources):,} rows; expected {EXPECTED_MAG_COUNT:,}"
        )
    for source in sources:
        if not source.genome_id.startswith("EHM") or not source.genome_id[3:].isdigit():
            raise ValueError(f"invalid genome ID: {source.genome_id!r}")
        if not source.filename.endswith(".fa"):
            raise ValueError(f"unexpected source filename for {source.genome_id}: {source.filename!r}")
        if not source.url.startswith("https://") or not source.url.endswith(".fa.gz"):
            raise ValueError(f"unexpected FASTA URL for {source.genome_id}: {source.url!r}")
        if "/MAG//" in source.url:
            raise ValueError(f"missing batch directory in FASTA URL for {source.genome_id}")

    ids = [source.genome_id for source in sources]
    if len(ids) != len(set(ids)):
        raise ValueError("source index contains duplicate genome IDs")
    return sorted(sources, key=lambda source: source.genome_id)


def read_release_metadata(path: Path) -> dict[str, GenomeMetadata]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome_id", "assembly_id", "genome_size", "contigs", "genome_n50"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        rows = [
            (
                (row["genome_id"] or "").strip(),
                GenomeMetadata(
                    assembly_id=(row["assembly_id"] or "").strip(),
                    sequence_bp=int(float(row["genome_size"])),
                    contigs=int(float(row["contigs"])),
                    n50=int(float(row["genome_n50"])),
                ),
            )
            for row in reader
        ]
    ids = [genome_id for genome_id, _ in rows]
    if len(ids) != EXPECTED_MAG_COUNT:
        raise ValueError(
            f"genome metadata has {len(ids):,} rows; expected {EXPECTED_MAG_COUNT:,}"
        )
    if len(ids) != len(set(ids)):
        raise ValueError("genome metadata contains duplicate genome IDs")
    return dict(rows)


def validated_fasta_metrics(payload: bytes, genome_id: str) -> tuple[int, int, int, int]:
    if not payload.startswith(b"\x1f\x8b"):
        raise ValueError(f"{genome_id}: downloaded file is not gzip data")
    try:
        fasta = gzip.decompress(payload)
    except (EOFError, OSError) as exc:
        raise ValueError(f"{genome_id}: invalid gzip data: {exc}") from exc
    if not fasta.startswith(b">"):
        raise ValueError(f"{genome_id}: decompressed file is not FASTA")
    contig_lengths: list[int] = []
    current_length = 0
    for line in fasta.splitlines():
        if line.startswith(b">"):
            if current_length:
                contig_lengths.append(current_length)
            current_length = 0
        elif line:
            current_length += len(line.strip())
    if current_length:
        contig_lengths.append(current_length)
    if not contig_lengths:
        raise ValueError(f"{genome_id}: FASTA contains no sequences")

    sequence_bp = sum(contig_lengths)
    halfway = sequence_bp / 2
    cumulative = 0
    n50 = 0
    for length in sorted(contig_lengths, reverse=True):
        cumulative += length
        if cumulative >= halfway:
            n50 = length
            break
    return len(fasta), sequence_bp, len(contig_lengths), n50


def download(source: Source, retries: int) -> Download:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = Request(source.url, headers={"User-Agent": USER_AGENT})
            with urlopen(request, timeout=120) as response:
                payload = response.read()
                expected = response.headers.get("Content-Length")
            if expected is not None and len(payload) != int(expected):
                raise IOError(f"short download: got {len(payload)}, expected {expected}")
            fasta_size_bytes, sequence_bp, contigs, n50 = validated_fasta_metrics(
                payload, source.genome_id
            )
            return Download(
                source=source,
                payload=payload,
                fasta_size_bytes=fasta_size_bytes,
                sequence_bp=sequence_bp,
                contigs=contigs,
                n50=n50,
                md5=hashlib.md5(payload).hexdigest(),
                sha256=hashlib.sha256(payload).hexdigest(),
            )
        except (HTTPError, URLError, TimeoutError, OSError, ValueError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(min(30, 2 ** (attempt - 1)))
    raise RuntimeError(
        f"{source.genome_id}: download failed after {retries} attempts: {last_error}"
    )


def source_assembly(source: Source) -> str:
    assembly_id, marker, _ = source.filename.partition("_bin.")
    if not marker or not assembly_id.startswith("EHA"):
        raise ValueError(f"cannot derive assembly ID from {source.filename!r}")
    return assembly_id


def expected_signature(genome_id: str, metadata: GenomeMetadata) -> tuple[int]:
    return KNOWN_SIGNATURE_OVERRIDES.get(genome_id, metadata.signature)


def reconcile_assembly(
    assembly_id: str,
    targets: list[Source],
    downloads: list[Download],
    metadata: dict[str, GenomeMetadata],
) -> list[tuple[Source, Download]]:
    expected: dict[tuple[int], list[Source]] = defaultdict(list)
    actual: dict[tuple[int], list[Download]] = defaultdict(list)
    for target in targets:
        expected[expected_signature(target.genome_id, metadata[target.genome_id])].append(target)
    for item in downloads:
        actual[item.signature].append(item)

    if set(expected) != set(actual) or any(
        len(expected[signature]) != len(actual[signature])
        for signature in set(expected) | set(actual)
    ):
        missing = sorted(set(expected) - set(actual))[:5]
        extra = sorted(set(actual) - set(expected))[:5]
        raise ValueError(
            f"{assembly_id}: FASTA/metadata signatures cannot be reconciled; "
            f"missing={missing}, extra={extra}"
        )

    reconciled: list[tuple[Source, Download]] = []
    for signature, signature_targets in expected.items():
        signature_downloads = actual[signature]
        if len(signature_targets) == 1:
            reconciled.append((signature_targets[0], signature_downloads[0]))
            continue

        # Exact duplicate signatures are rare. Preserve an already matching
        # Core link only when every target can be paired unambiguously by ID.
        by_core_id = {item.source.genome_id: item for item in signature_downloads}
        if set(by_core_id) != {target.genome_id for target in signature_targets}:
            raise ValueError(
                f"{assembly_id}: ambiguous duplicate signature {signature} "
                f"for {[target.genome_id for target in signature_targets]}"
            )
        reconciled.extend(
            (target, by_core_id[target.genome_id]) for target in signature_targets
        )
    return sorted(reconciled, key=lambda pair: pair[0].genome_id)


def build_archive(
    output_path: Path,
    sources: list[Source],
    metadata: dict[str, GenomeMetadata],
    workers: int,
    retries: int,
) -> list[dict[str, object]]:
    partial_path = output_path.with_name(output_path.name + ".part")
    if partial_path.exists():
        partial_path.unlink()
    manifest: list[dict[str, object]] = []
    started = time.monotonic()
    source_bytes = 0
    processed = 0
    next_report = 100
    assembly_groups: dict[str, list[Source]] = {}
    for source in sources:
        assembly_id = source_assembly(source)
        if metadata[source.genome_id].assembly_id != assembly_id:
            raise ValueError(
                f"{source.genome_id}: source assembly {assembly_id} does not match "
                f"metadata assembly {metadata[source.genome_id].assembly_id}"
            )
        assembly_groups.setdefault(assembly_id, []).append(source)

    with partial_path.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=1,
            fileobj=raw_output,
            mtime=0,
        ) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT) as archive:
                with ThreadPoolExecutor(max_workers=workers) as executor:
                    for assembly_id, targets in assembly_groups.items():
                        downloads = list(
                            executor.map(
                                lambda source: download(source, retries),
                                targets,
                                buffersize=max(2, workers * 2),
                            )
                        )
                        reconciled = reconcile_assembly(
                            assembly_id, targets, downloads, metadata
                        )
                        remapped = [
                            (target.genome_id, item.source.genome_id)
                            for target, item in reconciled
                            if target.genome_id != item.source.genome_id
                        ]
                        if remapped:
                            print(
                                f"Reconciled {assembly_id}: remapped {len(remapped)} Core links",
                                file=sys.stderr,
                                flush=True,
                            )

                        for target, item in reconciled:
                            target_metadata = metadata[target.genome_id]
                            info = tarfile.TarInfo(target.archive_name)
                            info.size = len(item.payload)
                            info.mtime = 0
                            info.mode = 0o644
                            info.uid = info.gid = 0
                            info.uname = info.gname = ""
                            archive.addfile(info, io.BytesIO(item.payload))
                            source_bytes += len(item.payload)
                            processed += 1
                            manifest.append(
                                {
                                    "genome_id": target.genome_id,
                                    "archive_path": target.archive_name,
                                    "source_filename": item.source.filename,
                                    "core_genome_id_for_source": item.source.genome_id,
                                    "core_link_was_remapped": (
                                        target.genome_id != item.source.genome_id
                                    ),
                                    "compressed_size": len(item.payload),
                                    "fasta_size_bytes": item.fasta_size_bytes,
                                    "sequence_bp": item.sequence_bp,
                                    "metadata_sequence_bp": target_metadata.sequence_bp,
                                    "contigs": item.contigs,
                                    "metadata_contigs": target_metadata.contigs,
                                    "n50": item.n50,
                                    "metadata_n50": target_metadata.n50,
                                    "signature_matches_current_metadata": (
                                        item.sequence_bp == target_metadata.sequence_bp
                                        and item.contigs == target_metadata.contigs
                                        and item.n50 == target_metadata.n50
                                    ),
                                    "md5": item.md5,
                                    "sha256": item.sha256,
                                }
                            )

                        if processed >= next_report or processed == len(sources):
                            elapsed = max(time.monotonic() - started, 0.001)
                            rate = source_bytes / elapsed / 1_000_000
                            print(
                                f"Downloaded {processed:,}/{len(sources):,} MAGs; "
                                f"{source_bytes / 1_000_000_000:.2f} GB; {rate:.1f} MB/s",
                                flush=True,
                            )
                            next_report = ((processed // 100) + 1) * 100
    partial_path.replace(output_path)
    return manifest


def file_hashes(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * READ_CHUNK):
            md5.update(chunk)
            sha256.update(chunk)
    return md5.hexdigest(), sha256.hexdigest()


def validate_final_archive(output_path: Path, expected_names: list[str]) -> None:
    print(f"Reading back and validating {output_path}", flush=True)
    with tarfile.open(output_path, "r:gz") as archive:
        members = archive.getmembers()
    names = [member.name for member in members]
    if names != expected_names:
        raise ValueError("final archive member list does not match the EHR02 manifest")
    if any(not member.isfile() or member.size <= 0 for member in members):
        raise ValueError("final archive contains an empty or non-file member")


def write_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "genome_id",
        "archive_path",
        "source_filename",
        "core_genome_id_for_source",
        "core_link_was_remapped",
        "compressed_size",
        "fasta_size_bytes",
        "sequence_bp",
        "metadata_sequence_bp",
        "contigs",
        "metadata_contigs",
        "n50",
        "metadata_n50",
        "signature_matches_current_metadata",
        "md5",
        "sha256",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    if args.workers < 1 or args.retries < 1:
        raise ValueError("--workers and --retries must be positive")

    sources = read_sources(args.source_index)
    release_metadata = read_release_metadata(args.genome_metadata)
    source_ids = {source.genome_id for source in sources}
    release_ids = set(release_metadata)
    if source_ids != release_ids:
        missing = sorted(release_ids - source_ids)
        extra = sorted(source_ids - release_ids)
        raise ValueError(f"source/metadata ID mismatch; missing={missing[:5]}, extra={extra[:5]}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output.with_name(args.output.name.removesuffix(".tar.gz") + ".manifest.tsv")
    checksum_path = args.output.with_name(args.output.name + ".checksums.txt")

    manifest = build_archive(
        output_path=args.output,
        sources=sources,
        metadata=release_metadata,
        workers=args.workers,
        retries=args.retries,
    )
    expected_names = [str(row["archive_path"]) for row in manifest]
    if len(manifest) != len(sources):
        raise ValueError("staging archive is incomplete")
    if {row["genome_id"] for row in manifest} != release_ids:
        raise ValueError("archive genome-ID set does not match release metadata")

    validate_final_archive(args.output, expected_names)
    archive_md5, archive_sha256 = file_hashes(args.output)
    write_manifest(manifest_path, manifest)
    checksum_path.write_text(
        f"MD5 ({args.output.name}) = {archive_md5}\n"
        f"SHA256 ({args.output.name}) = {archive_sha256}\n"
    )
    summary = {
        "archive": str(args.output),
        "files": len(expected_names),
        "remapped_core_links": sum(
            bool(row["core_link_was_remapped"]) for row in manifest
        ),
        "metadata_exceptions": sum(
            not bool(row["signature_matches_current_metadata"]) for row in manifest
        ),
        "size_bytes": args.output.stat().st_size,
        "md5": archive_md5,
        "sha256": archive_sha256,
        "manifest": str(manifest_path),
        "checksums": str(checksum_path),
    }
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted; rerun the same command to restart the streaming build.", file=sys.stderr)
        raise SystemExit(130)

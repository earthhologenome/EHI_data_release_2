# EHR02 MAG source audit

The 11,538 `genome_id` values in `data/EHR02_genome.tsv` were matched exactly
against the production EHI Core `mags` table on 2026-09-28. Every matched Core
record had a FASTA filename and URL.

Eleven Core URLs had an empty assembly-batch path component and returned HTTP
404:

- `EHM053140`–`EHM053146` (`EHA03568`)
- `EHM078214`–`EHM078217` (`EHA05388`)

For the local Zenodo build index, those paths were completed with `ABB0422` and
`ABB0614`, respectively. These are the batch directories used by the remaining
MAGs from the same assemblies, and all 11 corrected URLs returned HTTP 200.
The production Core records should be corrected separately; no Core data were
changed during archive preparation.

Core also had an 18-record filename/URL rotation within assembly `EHA04646`
(`EHM064972`–`EHM064989`). The FASTA length and contig-count signatures form
an exact permutation of those 18 MAG metadata rows: `bin.19` belongs to the
first ID, followed by `bin.2`, `bin.20`–`bin.28`, and `bin.3`–`bin.9`. The
local build index was remapped accordingly. No sequence file was modified.

The complete assembly-level validation found the same class of rotated Core
link in 47 further assemblies. In total, 1,355 Core links across 48 assemblies
were reassigned locally by exact genome length (18 in the source index and
1,337 by the builder). Genome lengths are unique within every EHR02 assembly,
so all assignments were one-to-one and unambiguous. The local manifest records
the original Core-linked EHM ID, actual source filename, final EHM ID, file
size, and checksums for every archive member.

The archive builder validates each download as gzip-compressed FASTA and checks
its total sequence length against `genome_size` in `data/EHR02_genome.tsv`
before adding it to the tarball.

Nine linked FASTAs expose separate metadata-field regressions. `EHM021731` is
2,117,230 bp, 458 contigs, GC 0.373, and N50 5,637. Those values agree with an
earlier independent EHI MAG metadata export, but the current Core/EHR02 TSV
reports 4,571,350 bp, GC 0.417, and N50 87,867 for the same ID. The builder
allows only this exact, documented length exception and records both values in
the local per-file manifest. `EHM032091` agrees on length and contig count but
has N50 12,738 in FASTA versus 12,674 in metadata. `EHM053140`–`EHM053146`
agree on length and N50 but have shifted contig-count metadata. The Core
records and release metadata should be corrected separately; no ERDA FASTA was
altered.

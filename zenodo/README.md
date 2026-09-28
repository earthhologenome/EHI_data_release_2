# EHR02 Zenodo payload

This directory mirrors the first data release's Zenodo deposit:

- `ehi_data_release_2_mags.tar.gz` contains exactly 11,538 files at the archive
  root, named `<genome_id>.fa.gz`.
- Zenodo receives that exact archive as four bytewise parts named
  `ehi_data_release_2_mags.tar.gz.part-00` through
  `ehi_data_release_2_mags.tar.gz.part-03`, avoiding Zenodo's long-request
  gateway timeout. Concatenating them recovers the validated archive exactly.
- `readme` is the user-facing reconstruction and extraction note.
- `metadata.json` is a Zenodo deposition API payload matching record 16757580:
  Dataset, CC0, creator `Alberdi, Antton`, and EHI community.

Current unpublished Zenodo draft:

- Draft ID: `23011784`
- Reserved DOI: `10.5281/zenodo.23011784`
- Review URL: `https://zenodo.org/deposit/23011784`
- State verified after upload: `unsubmitted` (`submitted: false`)

The large archive, its per-file manifest, checksums, and the private source URL
index are intentionally ignored by Git.

Build and validate the archive with the streaming builder (only the final
archive occupies disk space):

```sh
python3 scripts/build_zenodo_mag_archive.py \
  zenodo/EHR02_mag_sources.tsv \
  --genome-metadata data/EHR02_genome.tsv
```

The builder verifies gzip/FASTA structure, reconciles every source file to an
EHM ID within its assembly by exact genome length, reads the finished tarball
back in full, and writes a per-file audit manifest plus MD5/SHA-256 checksums.

Validated build:

- Files: 11,538
- Size: 7,138,772,352 bytes
- MD5: `cb7da06d12469fcb786e5536cc69a4b8`
- SHA-256: `5004f20606b689c3a87ee152baf9371353401ebff776b5e563725d8950f82850`

Validated upload parts:

- `part-00`: 1,785,724,928 bytes; MD5
  `626b7b2404cf70f3aee07ab5c597579e`; SHA-256
  `a5dd4ac677601b3ebbdaaf90ee8170ac7fac91bd3749d158ab7c0e994e38b88b`
- `part-01`: 1,785,724,928 bytes; MD5
  `ede59fa494172c01aa890fae3a01f127`; SHA-256
  `b8d990c8ab0319999db97882228fc5f45bc5fc55888c4c2dec2b48b8cba9293d`
- `part-02`: 1,785,724,928 bytes; MD5
  `752ffb9281520fab7e66d838c048131d`; SHA-256
  `0d80c45c5754be2f2ca301c42daf6762298ab6b0b9b48c6de239d7aaac572487`
- `part-03`: 1,781,597,568 bytes; MD5
  `0bf18689f3e970337ddc3e9062e187ed`; SHA-256
  `2cb8c64cebeb7f7234043a6dbc05879e8b2839caf30c4ababe32b66951a3e55d`

Create an unpublished Zenodo draft, reserve its DOI, apply `metadata.json`, and
upload the readme and archive payload with:

```sh
export ZENODO_ACCESS_TOKEN=...  # token with deposit:write scope
scripts/create_zenodo_draft.sh
```

When the four archive parts exist, the script uploads them instead of attempting
the full archive in a single request.

If an upload is interrupted after draft creation, rerun with the printed draft
ID so a duplicate draft is not created:

```sh
scripts/create_zenodo_draft.sh DRAFT_ID
```

The script deliberately has no publish action.

Do not publish the Zenodo draft until its displayed file count, archive MD5,
title, creator, CC0 licence, and EHI community have been checked. Publication
is irreversible; reserving the DOI and saving an unpublished draft are not.

# ehi-core metadata assessment for EHR02

Assessment date: 2026-09-26.

## Outcome

The authenticated R collector successfully reads the live `ehi-core` grid API
and has the schema needed to generate the data tables, statistics, and
visualisations currently defined for the second EHI data release. However, the
live EHR02 rows are not yet complete enough to publish: two hologenomes marked
as requiring metagenomic preprocessing have no preprocessing record. Five
additional unpreprocessed hologenomes contain genomic-only data; their missing
metagenomic preprocessing is expected and does not count as an audit error.
Draft TSVs have been generated with unavailable values left blank; the
writer's default strict mode still refuses the two remaining errors for a
final publication export.

The live collection on 2026-09-26 found 988 EHR02 hologenomes, 874 samples,
375 successful individual assemblies, and 10,181 MAGs. Of those hologenomes,
49 are genomic-only and 939 require metagenomic preprocessing. Of the latter,
937 have preprocessing data and two do not.

Total sequenced gigabases are derived from `preprocessings.bases_pre_fastp`.
The core's `hologenomes.raw_data_gb` is the compressed-file size and is not the
same quantity as the historical `Total Data (GB)` field used by this analysis.

## Preprocessing classification and live row gaps

The following five records are marked `data_type = genome` and
`library_source = GENOMIC`; no preprocessing is required:

- `EHI03513`
- `EHI03514`
- `EHI03515`
- `EHI03516`
- `EHI03517`

The following two records still fail the audit:

- `EHI01165`
- `EHI01222`

`EHI01165` is marked `data_type = metagenome`; `EHI01222` has an `unknown`
data type but a `METAGENOMIC` library source. They need preprocessing records
with `bases_pre_fastp`, `host_bases`, and `metagenomic_bases`, or their core
classification must be corrected if they are actually genomic-only.

## Current release-link issue

The live public catalogue was checked on 2026-09-26. Catalogue version
`2026.09.25-1329` (core source snapshot `2026-09-25T13:22:55Z`) contains 997
EHR02 hologenomes but reports **0 EHR02 MAGs**. The authenticated live core now
contains 988 EHR02 hologenomes. The nine public-only records had their live
`release_id` cleared together on 2026-09-26 and have not yet disappeared from
the older public snapshot: `EHI00869`, `EHI00908`, `EHI00913`, `EHI00923`,
`EHI00937`, `EHI00943`, `EHI00985`, `EHI01025`, and `EHI01674`.

The essential relationship is not lost: an assembly reaches its release via
assembly sample -> preprocessing -> hologenome -> data release. The R collector
uses that chain and can therefore recover the EHR02 individual assemblies and
MAGs. However, `assemblies.release_code` is missing or inconsistent for those
assemblies, and the current public catalogue export relies on that field. It
must be populated (or the export must derive release membership through the
link chain) before EHR02 MAGs can appear correctly outside this analysis. None
of the 375 selected live individual assemblies currently has EHR02 in
`assemblies.release_code`; the R collector recovers the relationship through
the link chain and finds 10,181 MAGs.

The collector performs a row-level audit every time it runs. This matters
because having the necessary columns does not guarantee that every EHR02 row
has a value. Blocking omissions (for example a hologenome without a sample,
raw-data size or preprocessing composition) are labelled `error`; incomplete
but still renderable fields are labelled `warning`.

After excluding genomic-only records, the current non-blocking warnings are
168 missing/zero SingleM fractions, eight missing/zero Nonpareil completeness
estimates, and 59 MAGs without phylum or genus assignments. The affected
summaries and plots omit those values.

## Known schema gaps

The historical Airtable sample export contained three descriptors that are not
stored in `ehi-core`:

| Output column | Former Airtable field | Used downstream? | Effect |
|---|---|---:|---|
| `sample_origin` | `Origin` | No | Kept as `NA`; old sample TSV cannot be reproduced in full. |
| `freshness` | `Freshness` | No | Kept as `NA`; old sample TSV cannot be reproduced in full. |
| `time_to_freezing` | `Freezing` | No | Kept as `NA`; old sample TSV cannot be reproduced in full. |

These gaps do not block any calculation or visualisation in
`02_data_statistics.Rmd`, `03_host_tree.Rmd`, or `04_mag_statistics.Rmd`. If
the release must publish these descriptors as metadata, they need to be added
to the core's `samples` table and populated before publication.

## Required inputs outside ehi-core

The following visualization input is not core metadata and is currently absent
from `data/`:

- `EHR02_allhosts.nwk`, generated externally with TimeTree from the host list;

The four GTDB R214 tree and taxonomy files required by the MAG chapter have
been reused from the first-release repository. The chapter downloads the EHI
phylum-colour table from GitHub at render time, so a completely offline render
would still need a checked-in local copy. These are repository/render inputs
rather than missing fields in `ehi-core`.

## Fields that are available or derivable

| Release requirement | ehi-core source |
|---|---|
| EHR02 membership | `hologenomes.release_id` -> `data_releases.code` |
| Genomic-only versus metagenomic data | `hologenomes.data_type` and `hologenomes.library_source` |
| ENA accessions and compressed-file size | `hologenomes` |
| Total/host/metagenomic GB, SingleM and Nonpareil | `preprocessings` (base counts are converted to GB) |
| Sample -> capture -> specimen -> host | Core link columns across the four tables |
| Sampling place, date, environment and collaborators | `sampling_events` and participant lookup |
| Host taxonomy | `species` |
| Individual assemblies in EHR02 | Derived through assembly samples and preprocessings; batch type and success come from `assemblies` |
| MAG quality and taxonomy | `mags`, linked to release assemblies |

## Operational note

The existing grid API is sufficient but requires a core staff login. The R
helper uses the same Firebase email/password sign-in as the editor and never
stores credentials. For unattended rendering, credentials or a short-lived ID
token must be supplied through environment variables; they must not be added
to this repository.

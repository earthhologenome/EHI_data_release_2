source(testthat::test_path("..", "..", "R", "ehi_core.R"))

fixture_core <- function() {
  list(
    release = tibble::tibble(code = "EHR02"),
    hologenomes = tibble::tibble(
      code = "EHI00001", sample_id = "SAM01", raw_data_gb = 3,
      data_type = "metagenome", library_source = "METAGENOMIC",
      ena_experiment_accession = "ERX1", ena_run_accession = "ERR1",
      sequencing_batch_code = "SEB001"
    ),
    preprocessings = tibble::tibble(
      code = c("PR00001", "PR00002"), hologenome_id = c("EHI00001", "EHI00001"),
      bases_pre_fastp = c(5e9, 5e9),
      host_bases = c(1e9, 3e9), metagenomic_bases = c(2e9, 2e9),
      singlem_fraction = c(10, 20), nonpareil_c = c(0.8, 0.9)
    ),
    samples = tibble::tibble(
      code = "SAM01", capture_id = "SD0000101", sample_type = "Faecal",
      ena_sample_accession = "ERS1"
    ),
    captures = tibble::tibble(
      code = "SD0000101", specimen_id = "SD00001", sampling_event_id = "SE00001",
      length_mm = 40, weight_g = 8, development = "Adult"
    ),
    specimens = tibble::tibble(code = "SD00001", species = "Testudo example", sex = "Female"),
    sampling_events = tibble::tibble(
      code = "SE00001", participant_names = "Example Researcher", country = "Denmark",
      region = "Zealand", place = "Field", sampled_on = "2024-01-02",
      latitude = 55.7, longitude = 12.6, biome = "1000000 - Test biome",
      environment = "Woodland", weather = "Clear", temperature = 12
    ),
    species = tibble::tibble(
      code = "Testudo example", taxid = "123", tax_genus = "Testudo",
      tax_family = "Testudinidae", tax_order = "Testudines", tax_class = "Reptilia"
    ),
    assemblies = tibble::tibble(
      code = c("EHA00001", "EHA00002"), batch_type = c("Individual", "Coassembly"),
      successful = c(TRUE, TRUE), release_code = c("EHR02", "EHR02"),
      n50 = c(10000, 9000), l50 = c(20L, 25L), num_contigs = c(100L, 120L),
      num_bins = c(1L, 1L)
    ),
    assembly_samples = tibble::tibble(
      assembly_id = c("EHA00001", "EHA00002"),
      preprocessing_id = c("PR00001", "PR00002")
    ),
    mags = tibble::tibble(
      code = c("EHM000001", "EHM000002"), assembly_id = c("EHA00001", "EHA00002"),
      gtdb_release = "R214", gtdbtk_version = "2.3.2", tax_domain = "d__Bacteria", tax_phylum = "p__Bacillota",
      tax_class = "c__Bacilli", tax_order = "o__Test", tax_family = "f__Test",
      tax_genus = "g__Test", tax_species = "s__", completeness = 95,
      contamination = 1, size_bp = 2000000, gc = 45, n50 = 50000,
      coding_density = 0.9, contigs = 50L
    )
  )
}

as_metadata <- function(core = fixture_core()) {
  tables <- .build_ehi_release_tables(core, "EHR02")
  structure(
    c(tables, list(core = core, release = "EHR02", collected_at = Sys.time())),
    class = c("ehi_release_metadata", "list")
  )
}

test_that("core tables are converted to the release contracts", {
  metadata <- as_metadata()

  expect_identical(metadata$datasets, "EHI00001")
  expect_equal(metadata$rawdata_metadata$total_gb, 5)
  expect_equal(metadata$rawdata_metadata$host_gb, 2)
  expect_equal(metadata$rawdata_metadata$metagenomic_gb, 2)
  expect_equal(metadata$rawdata_metadata$lowqual_gb, 1)
  expect_equal(metadata$rawdata_metadata$microbial_fraction, 15)
  expect_equal(metadata$rawdata_metadata$nonpareil_completeness, 0.85)
  expect_true(metadata$rawdata_metadata$preprocessing_required)
  expect_identical(metadata$sample_metadata$rawdata_id, "EHI00001")
  expect_identical(metadata$capture_metadata$sample_id, "SAM01")
  expect_identical(metadata$specimen_metadata$capture_id, "SD0000101")
  expect_identical(metadata$sampling_metadata$capture_id, "SD0000101")
  expect_identical(metadata$host_metadata$host_order, "Testudines")
  expect_identical(metadata$assembly_metadata$assembly_id, "EHA00001")
  expect_identical(metadata$assembly_metadata$rawdata_id, "EHI00001")
  expect_identical(metadata$genome_metadata$genome_id, "EHM000001")
  expect_identical(metadata$genome_metadata$GTDB_version, "2.3.2")
  expect_identical(names(metadata$rawdata_metadata), c(
    "rawdata_id", "sample_id", "data_type", "library_source", "preprocessing_required",
    "accession_experiment", "accession_run", "sequencing_batch", "total_gb", "host_gb",
    "metagenomic_gb", "lowqual_gb", "microbial_fraction", "nonpareil_completeness"
  ))
  expect_identical(names(metadata$sample_metadata), c(
    "sample_id", "sample_type", "sample_origin", "freshness", "time_to_freezing",
    "rawdata_id", "capture_id", "accession"
  ))
  expect_identical(names(metadata$capture_metadata), c(
    "capture_id", "animal_length", "animal_weight", "animal_development", "specimen_id", "sample_id"
  ))
})

test_that("coassemblies are released only for single-individual host species", {
  metadata <- as_metadata()
  expect_false("EHA00002" %in% metadata$assembly_metadata$assembly_id)

  core <- fixture_core()
  core$specimens$species <- "Canis familiaris"
  core$species$code <- "Canis familiaris"
  metadata <- as_metadata(core)

  expect_setequal(metadata$assembly_metadata$assembly_id, c("EHA00001", "EHA00002"))
  expect_identical(
    metadata$assembly_metadata$assembly_type[metadata$assembly_metadata$assembly_id == "EHA00002"],
    "Coassembly"
  )
  expect_setequal(metadata$genome_metadata$genome_id, c("EHM000001", "EHM000002"))
  audit <- audit_ehi_core_release(metadata)
  expect_false("individual coassembly spans several specimens" %in% audit$check)

  core$hologenomes <- dplyr::bind_rows(
    core$hologenomes,
    dplyr::mutate(core$hologenomes, code = "EHI00002", sample_id = "SAM02")
  )
  core$preprocessings$hologenome_id[2] <- "EHI00002"
  core$samples <- dplyr::bind_rows(
    core$samples,
    dplyr::mutate(core$samples, code = "SAM02", capture_id = "SD0000201")
  )
  core$captures <- dplyr::bind_rows(
    core$captures,
    dplyr::mutate(core$captures, code = "SD0000201", specimen_id = "SD00002")
  )
  core$specimens <- dplyr::bind_rows(
    core$specimens,
    dplyr::mutate(core$specimens, code = "SD00002")
  )
  core$assembly_samples <- dplyr::bind_rows(
    core$assembly_samples,
    tibble::tibble(assembly_id = "EHA00002", preprocessing_id = "PR00001")
  )
  audit <- audit_ehi_core_release(as_metadata(core))
  expect_equal(audit$affected[audit$check == "individual coassembly spans several specimens"], 1L)
})

test_that("the legacy sample gaps are explicit but non-blocking", {
  metadata <- as_metadata()
  audit <- audit_ehi_core_release(metadata)

  expect_false(any(audit$severity == "error"))
  expect_setequal(
    audit$check[audit$severity == "schema_gap"],
    c("sample_origin", "freshness", "time_to_freezing")
  )
  expect_true(all(is.na(metadata$sample_metadata$sample_origin)))
})

test_that("blocking link and composition gaps are reported", {
  core <- fixture_core()
  core$hologenomes$sample_id <- NA_character_
  core$preprocessings$host_bases <- NA_real_
  metadata <- as_metadata(core)
  audit <- audit_ehi_core_release(metadata)

  expect_true("missing sample link" %in% audit$check[audit$severity == "error"])
  expect_true("missing required host GB" %in% audit$check[audit$severity == "error"])
})

test_that("genomic-only hologenomes do not require preprocessing", {
  core <- fixture_core()
  core$hologenomes$data_type <- "genome"
  core$hologenomes$library_source <- "GENOMIC"
  core$preprocessings <- core$preprocessings[0, ]
  metadata <- as_metadata(core)
  audit <- audit_ehi_core_release(metadata)

  expect_false(any(
    audit$severity == "error" &
      audit$check %in% c(
        "missing required total GB",
        "missing required host GB",
        "missing required metagenomic GB"
      )
  ))
  expect_equal(
    audit$affected[audit$check == "genomic-only without preprocessing"],
    1L
  )
})

test_that("grid rows retain scalar cells and encode link previews", {
  rows <- list(
    list(
      id = 7L,
      values = list(
        code = "EHI00001",
        missing = NULL,
        preprocessings = list(codes = list("PR00001"), total = 1L)
      )
    ),
    list(
      id = 8L,
      values = list(
        code = "EHI00002",
        missing = NULL,
        preprocessings = "No linked records"
      )
    )
  )
  result <- .grid_rows(rows)

  expect_identical(result$.row_id, c(7L, 8L))
  expect_identical(result$code, c("EHI00001", "EHI00002"))
  expect_true(all(is.na(result$missing)))
  expect_match(result$preprocessings[[1]], "PR00001", fixed = TRUE)
  expect_identical(result$preprocessings[[2]], "No linked records")
  expect_identical(class(result$preprocessings), "character")
})

test_that("complete credentials take precedence over a stale token", {
  original_login <- ehi_core_login
  on.exit(assign("ehi_core_login", original_login, envir = globalenv()), add = TRUE)
  assign("ehi_core_login", function(...) "fresh-token", envir = globalenv())

  expect_identical(
    ehi_core_token(token = "stale-token", email = "user@example.org", password = "secret"),
    "fresh-token"
  )
  expect_identical(
    ehi_core_token(token = "current-token", email = "", password = ""),
    "current-token"
  )
})

test_that("writer creates every release table and the audit", {
  directory <- tempfile("ehi-release-")
  metadata <- as_metadata()
  write_ehi_core_release(metadata, directory)

  expected <- c(
    paste0(
      "EHR02_",
      c("rawdata", "sample", "capture", "specimen", "sampling", "host", "assembly", "genome"),
      ".tsv"
    ),
    "EHR02_allhosts.txt", "EHR02_allhosts.tsv", "EHR02_metadata_audit.tsv"
  )
  expect_setequal(list.files(directory), expected)
})

test_that("writer preserves release tables when blocking checks fail", {
  directory <- tempfile("ehi-release-invalid-")
  dir.create(directory)
  marker <- file.path(directory, "EHR02_rawdata.tsv")
  writeLines("existing", marker)
  core <- fixture_core()
  core$preprocessings$host_bases <- NA_real_
  metadata <- as_metadata(core)

  expect_error(
    write_ehi_core_release(metadata, directory),
    "blocking metadata checks"
  )
  expect_identical(readLines(marker), "existing")
  expect_true(file.exists(file.path(directory, "EHR02_metadata_audit.tsv")))
})

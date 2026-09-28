# Helpers for collecting a reproducible data-release snapshot from ehi-core.
#
# ehi-core's grid API returns the same linked records the staff editor shows.
# This client signs in through the same Firebase email/password endpoint as the
# editor, pages through the API, and converts the linked core tables into the
# stable TSV contracts used by this book.

EHI_CORE_URL <- "https://ehi-core-backend-874005116816.europe-west1.run.app"

# This is Firebase's public web API key, not a credential. Access is enforced
# by Firebase authentication and ehi-core's own user/role checks.
EHI_FIREBASE_API_KEY <- "AIzaSyBKLbfOiolqSiasdkdyx_QDVPMoTJegSWU"

.require_ehi_packages <- function() {
  packages <- c("dplyr", "httr2", "jsonlite", "purrr", "readr", "tibble", "tidyr")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install the required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

.nonempty <- function(x) {
  length(x) == 1L && !is.na(x) && nzchar(trimws(x))
}

ehi_core_login <- function(
    email = Sys.getenv("EHI_CORE_EMAIL"),
    password = Sys.getenv("EHI_CORE_PASSWORD"),
    api_key = Sys.getenv("EHI_FIREBASE_API_KEY", EHI_FIREBASE_API_KEY)) {
  .require_ehi_packages()
  if (!.nonempty(email)) {
    stop("Set EHI_CORE_EMAIL or pass email= to ehi_core_login().", call. = FALSE)
  }
  if (!.nonempty(password)) {
    if (!interactive()) {
      stop(
        "Set EHI_CORE_PASSWORD when running non-interactively, or set EHI_CORE_TOKEN.",
        call. = FALSE
      )
    }
    if (!requireNamespace("askpass", quietly = TRUE)) {
      stop("Install the askpass package to enter the ehi-core password securely.", call. = FALSE)
    }
    password <- askpass::askpass("ehi-core password")
  }

  response <- httr2::request("https://identitytoolkit.googleapis.com/v1/accounts:signInWithPassword") |>
    httr2::req_url_query(key = api_key) |>
    httr2::req_body_json(list(
      email = trimws(email),
      password = password,
      returnSecureToken = TRUE
    )) |>
    httr2::req_timeout(30) |>
    httr2::req_perform()

  body <- httr2::resp_body_json(response, simplifyVector = TRUE)
  if (!.nonempty(body$idToken)) {
    stop("Firebase signed in but returned no ID token.", call. = FALSE)
  }
  body$idToken
}

ehi_core_token <- function(
    token = Sys.getenv("EHI_CORE_TOKEN"),
    email = Sys.getenv("EHI_CORE_EMAIL"),
    password = Sys.getenv("EHI_CORE_PASSWORD"),
    ...) {
  # Prefer a fresh sign-in when complete credentials are available. Firebase
  # ID tokens expire quickly, so an old EHI_CORE_TOKEN should not shadow the
  # durable non-interactive login configured alongside it.
  if (.nonempty(email) && .nonempty(password)) {
    return(ehi_core_login(email = email, password = password, ...))
  }
  if (.nonempty(token)) return(trimws(token))
  ehi_core_login(email = email, password = password, ...)
}

.ehi_core_get <- function(path, token, base_url, query = list()) {
  request <- httr2::request(paste0(sub("/+$", "", base_url), path)) |>
    httr2::req_auth_bearer_token(token) |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_retry(max_tries = 4) |>
    httr2::req_timeout(60)
  if (length(query)) {
    request <- do.call(httr2::req_url_query, c(list(request), query))
  }
  response <- httr2::req_perform(request)
  httr2::resp_body_json(response, simplifyVector = FALSE)
}

.grid_cell <- function(value) {
  if (is.null(value)) return(NA)
  if (is.atomic(value) && length(value) == 1L) return(value)
  # Link previews and file cells are not used by the release transformations,
  # but retaining them as JSON makes the generic table reader lossless.
  as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null"))
}

.grid_rows <- function(rows) {
  if (!length(rows)) return(tibble::tibble())
  purrr::map_dfr(rows, function(row) {
    values <- lapply(row$values, .grid_cell)
    tibble::as_tibble_row(c(list(.row_id = row$id), values))
  })
}

.core_table_columns <- list(
  preprocessings = c(
    "code", "hologenome_id", "bases_pre_fastp", "host_bases", "metagenomic_bases",
    "singlem_fraction", "nonpareil_c"
  ),
  samples = c("code", "capture_id", "sample_type", "ena_sample_accession"),
  captures = c(
    "code", "specimen_id", "sampling_event_id", "length_mm", "weight_g", "development"
  ),
  specimens = c("code", "species", "sex"),
  sampling_events = c(
    "code", "participant_names", "country", "region", "place", "sampled_on",
    "latitude", "longitude", "biome", "environment", "weather", "temperature"
  ),
  species = c("code", "taxid", "tax_genus", "tax_family", "tax_order", "tax_class"),
  assemblies = c(
    "code", "batch_type", "successful", "release_code", "n50", "l50", "num_contigs", "num_bins"
  ),
  assembly_samples = c("assembly_id", "preprocessing_id"),
  mags = c(
    "code", "assembly_id", "gtdb_release", "gtdbtk_version", "tax_domain", "tax_phylum", "tax_class",
    "tax_order", "tax_family", "tax_genus", "tax_species", "completeness",
    "contamination", "size_bp", "gc", "n50", "coding_density", "contigs"
  )
)

.empty_core_table <- function(table) {
  columns <- .core_table_columns[[table]]
  if (is.null(columns)) return(tibble::tibble())
  numeric_columns <- c(
    "bases_pre_fastp", "host_bases", "metagenomic_bases", "singlem_fraction", "nonpareil_c",
    "length_mm", "weight_g", "latitude", "longitude", "temperature", "n50", "l50",
    "num_contigs", "num_bins", "completeness", "contamination", "size_bp", "gc",
    "coding_density", "contigs"
  )
  values <- lapply(columns, function(column) {
    if (column == "successful") logical() else if (column %in% numeric_columns) numeric() else character()
  })
  tibble::as_tibble(stats::setNames(values, columns))
}

ehi_core_read_table <- function(
    table,
    token = ehi_core_token(),
    base_url = Sys.getenv("EHI_CORE_URL", EHI_CORE_URL),
    filters = NULL,
    page_size = 500L,
    progress = interactive()) {
  .require_ehi_packages()
  if (page_size < 1L || page_size > 500L) {
    stop("page_size must be between 1 and ehi-core's limit of 500.", call. = FALSE)
  }
  filter_json <- if (is.null(filters) || !length(filters)) {
    NULL
  } else {
    jsonlite::toJSON(
      list(match = "all", conditions = filters),
      auto_unbox = TRUE,
      null = "null"
    )
  }

  offset <- 0L
  total <- Inf
  pages <- list()
  while (offset < total) {
    query <- list(offset = offset, limit = as.integer(page_size))
    if (!is.null(filter_json)) query$filters <- filter_json
    page <- .ehi_core_get(
      paste0("/api/grid/tables/", table, "/rows"),
      token = token,
      base_url = base_url,
      query = query
    )
    total <- as.integer(page$total)
    pages[[length(pages) + 1L]] <- .grid_rows(page$rows)
    received <- length(page$rows)
    if (!received) break
    offset <- offset + received
    if (progress && total > page_size) {
      message(sprintf("ehi-core %-18s %d/%d", table, min(offset, total), total))
    }
  }
  dplyr::bind_rows(pages)
}

.read_codes <- function(read_table, table, column, values, chunk_size = 100L) {
  values <- sort(unique(stats::na.omit(as.character(values))))
  values <- values[nzchar(values)]
  if (!length(values)) return(.empty_core_table(table))
  chunks <- split(values, ceiling(seq_along(values) / chunk_size))
  result <- purrr::map_dfr(chunks, function(chunk) {
    read_table(
      table,
      filters = list(list(column = column, op = "is_any_of", values = unname(chunk)))
    )
  })
  if (!nrow(result) && !".row_id" %in% names(result)) return(.empty_core_table(table))
  dplyr::distinct(result, .data$.row_id, .keep_all = TRUE)
}

.mean_or_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

.collapse_codes <- function(x) {
  x <- sort(unique(stats::na.omit(as.character(x))))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(x, collapse = ", ")
}

.is_genomic_only <- function(data_type, library_source) {
  normalized_type <- tolower(trimws(dplyr::coalesce(as.character(data_type), "")))
  normalized_source <- toupper(trimws(dplyr::coalesce(as.character(library_source), "")))
  normalized_type == "genome" | normalized_source == "GENOMIC"
}

.missing_columns <- function(data, columns, table) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      "ehi-core table '", table, "' does not expose expected columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(data)
}

# Coassemblies are normally excluded because they pool samples from different
# animals. For these host species the coassemblies were built from repeated
# samples of a single individual, so their MAGs are released like those of
# individual assemblies.
EHI_INDIVIDUAL_COASSEMBLY_SPECIES <- c("Canis familiaris")

.build_ehi_release_tables <- function(
    core,
    release,
    individual_coassembly_species = EHI_INDIVIDUAL_COASSEMBLY_SPECIES) {
  hologenomes <- core$hologenomes
  preprocessings <- core$preprocessings
  samples <- core$samples
  captures <- core$captures
  specimens <- core$specimens
  sampling_events <- core$sampling_events
  species <- core$species
  assemblies <- core$assemblies
  assembly_samples <- core$assembly_samples
  mags <- core$mags

  .missing_columns(hologenomes, c(
    "code", "sample_id", "data_type", "library_source", "raw_data_gb",
    "ena_experiment_accession", "ena_run_accession", "sequencing_batch_code"
  ), "hologenomes")
  .missing_columns(preprocessings, c(
    "code", "hologenome_id", "bases_pre_fastp", "host_bases", "metagenomic_bases",
    "singlem_fraction", "nonpareil_c"
  ), "preprocessings")
  .missing_columns(samples, c("code", "capture_id", "sample_type", "ena_sample_accession"), "samples")
  .missing_columns(captures, c(
    "code", "specimen_id", "sampling_event_id", "length_mm", "weight_g", "development"
  ), "captures")
  .missing_columns(specimens, c("code", "species", "sex"), "specimens")
  .missing_columns(sampling_events, c(
    "code", "participant_names", "country", "region", "place", "sampled_on",
    "latitude", "longitude", "biome", "environment", "weather", "temperature"
  ), "sampling_events")
  .missing_columns(species, c(
    "code", "taxid", "tax_genus", "tax_family", "tax_order", "tax_class"
  ), "species")
  .missing_columns(assemblies, c(
    "code", "batch_type", "successful", "release_code", "n50", "l50", "num_contigs", "num_bins"
  ), "assemblies")
  .missing_columns(assembly_samples, c("assembly_id", "preprocessing_id"), "assembly_samples")
  .missing_columns(mags, c(
    "code", "assembly_id", "gtdb_release", "gtdbtk_version", "tax_domain", "tax_phylum", "tax_class",
    "tax_order", "tax_family", "tax_genus", "tax_species", "completeness",
    "contamination", "size_bp", "gc", "n50", "coding_density", "contigs"
  ), "mags")

  prep_summary <- preprocessings |>
    dplyr::group_by(.data$hologenome_id) |>
    dplyr::summarise(
      total_gb = .mean_or_na(.data$bases_pre_fastp) / 1e9,
      host_gb = .mean_or_na(.data$host_bases) / 1e9,
      metagenomic_gb = .mean_or_na(.data$metagenomic_bases) / 1e9,
      microbial_fraction = .mean_or_na(.data$singlem_fraction),
      nonpareil_completeness = .mean_or_na(.data$nonpareil_c),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      microbial_fraction = dplyr::na_if(.data$microbial_fraction, 0),
      nonpareil_completeness = dplyr::na_if(.data$nonpareil_completeness, 0)
    )

  rawdata_metadata <- hologenomes |>
    dplyr::transmute(
      rawdata_id = .data$code,
      sample_id = .data$sample_id,
      data_type = .data$data_type,
      library_source = .data$library_source,
      preprocessing_required = !.is_genomic_only(.data$data_type, .data$library_source),
      accession_experiment = .data$ena_experiment_accession,
      accession_run = .data$ena_run_accession,
      sequencing_batch = .data$sequencing_batch_code
    ) |>
    dplyr::left_join(prep_summary, by = c("rawdata_id" = "hologenome_id")) |>
    dplyr::mutate(lowqual_gb = .data$total_gb - .data$host_gb - .data$metagenomic_gb) |>
    dplyr::select(dplyr::all_of(c(
      "rawdata_id", "sample_id", "data_type", "library_source",
      "preprocessing_required", "accession_experiment", "accession_run",
      "sequencing_batch", "total_gb", "host_gb", "metagenomic_gb",
      "lowqual_gb", "microbial_fraction", "nonpareil_completeness"
    ))) |>
    dplyr::arrange(.data$rawdata_id)

  datasets <- rawdata_metadata$rawdata_id
  hologenomes_by_sample <- rawdata_metadata |>
    dplyr::filter(!is.na(.data$sample_id)) |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::summarise(rawdata_id = .collapse_codes(.data$rawdata_id), .groups = "drop")

  sample_metadata <- samples |>
    dplyr::transmute(
      sample_id = .data$code,
      sample_type = .data$sample_type,
      # These three legacy Airtable descriptors are not represented in ehi-core.
      sample_origin = NA_character_,
      freshness = NA_character_,
      time_to_freezing = NA_character_,
      capture_id = .data$capture_id,
      accession = .data$ena_sample_accession
    ) |>
    dplyr::left_join(hologenomes_by_sample, by = "sample_id") |>
    dplyr::select(dplyr::all_of(c(
      "sample_id", "sample_type", "sample_origin", "freshness",
      "time_to_freezing", "rawdata_id", "capture_id", "accession"
    ))) |>
    dplyr::arrange(.data$sample_id)

  samples_by_capture <- samples |>
    dplyr::filter(!is.na(.data$capture_id)) |>
    dplyr::group_by(.data$capture_id) |>
    dplyr::summarise(sample_id = .collapse_codes(.data$code), .groups = "drop")
  capture_metadata <- captures |>
    dplyr::transmute(
      capture_id = .data$code,
      animal_length = .data$length_mm,
      animal_weight = .data$weight_g,
      animal_development = .data$development,
      specimen_id = .data$specimen_id
    ) |>
    dplyr::left_join(samples_by_capture, by = "capture_id") |>
    dplyr::arrange(.data$capture_id)

  captures_by_specimen <- captures |>
    dplyr::filter(!is.na(.data$specimen_id)) |>
    dplyr::group_by(.data$specimen_id) |>
    dplyr::summarise(capture_id = .collapse_codes(.data$code), .groups = "drop")
  specimen_metadata <- specimens |>
    dplyr::transmute(
      specimen_id = .data$code,
      host_species = .data$species,
      sex = .data$sex
    ) |>
    dplyr::left_join(captures_by_specimen, by = "specimen_id") |>
    dplyr::arrange(.data$specimen_id)

  captures_by_event <- captures |>
    dplyr::filter(!is.na(.data$sampling_event_id)) |>
    dplyr::group_by(sampling_id = .data$sampling_event_id) |>
    dplyr::summarise(capture_id = .collapse_codes(.data$code), .groups = "drop")
  sampling_metadata <- sampling_events |>
    dplyr::transmute(
      sampling_id = .data$code,
      researcher = .data$participant_names,
      region = .data$region,
      date = .data$sampled_on,
      biome = .data$biome,
      environment = .data$environment,
      weather = .data$weather,
      temperature = .data$temperature,
      place = .data$place,
      latitude = .data$latitude,
      longitude = .data$longitude,
      country = .data$country
    ) |>
    dplyr::left_join(captures_by_event, by = "sampling_id") |>
    dplyr::arrange(.data$sampling_id)

  host_metadata <- species |>
    dplyr::transmute(
      host_species = .data$code,
      taxid = .data$taxid,
      host_genus = .data$tax_genus,
      host_family = .data$tax_family,
      host_order = .data$tax_order,
      host_class = .data$tax_class
    ) |>
    dplyr::arrange(.data$host_species)

  assembly_sources <- assembly_samples |>
    dplyr::inner_join(
      preprocessings |>
        dplyr::transmute(preprocessing_id = .data$code, rawdata_id = .data$hologenome_id),
      by = "preprocessing_id"
    ) |>
    dplyr::left_join(
      rawdata_metadata |> dplyr::select("rawdata_id", "sample_id"),
      by = "rawdata_id"
    ) |>
    dplyr::left_join(
      samples |> dplyr::select(sample_id = "code", "capture_id"),
      by = "sample_id"
    ) |>
    dplyr::left_join(
      captures |> dplyr::select(capture_id = "code", "specimen_id"),
      by = "capture_id"
    ) |>
    dplyr::left_join(
      specimens |> dplyr::select(specimen_id = "code", host_species = "species"),
      by = "specimen_id"
    ) |>
    dplyr::group_by(.data$assembly_id) |>
    dplyr::summarise(
      rawdata_id = .collapse_codes(.data$rawdata_id),
      individual_coassembly_host = !anyNA(.data$host_species) &&
        all(.data$host_species %in% individual_coassembly_species),
      .groups = "drop"
    )
  assembly_metadata <- assemblies |>
    dplyr::left_join(assembly_sources, by = c("code" = "assembly_id")) |>
    dplyr::filter(
      .data$code %in% assembly_sources$assembly_id,
      .data$batch_type == "Individual" |
        (.data$batch_type == "Coassembly" & .data$individual_coassembly_host %in% TRUE),
      .data$successful %in% TRUE,
      !is.na(.data$num_contigs),
      .data$num_contigs > 0
    ) |>
    dplyr::transmute(
      assembly_id = .data$code,
      assembly_type = .data$batch_type,
      assembly_n50 = .data$n50,
      assembly_l50 = .data$l50,
      assembly_contigs = .data$num_contigs,
      num_bins = .data$num_bins,
      rawdata_id = .data$rawdata_id
    ) |>
    dplyr::select(dplyr::all_of(c(
      "assembly_id", "assembly_type", "rawdata_id", "assembly_n50", "assembly_l50",
      "assembly_contigs", "num_bins"
    ))) |>
    dplyr::arrange(.data$assembly_id)

  genome_metadata <- mags |>
    dplyr::filter(.data$assembly_id %in% assembly_metadata$assembly_id) |>
    dplyr::transmute(
      genome_id = .data$code,
      assembly_id = .data$assembly_id,
      GTDB_version = .data$gtdbtk_version,
      domain = .data$tax_domain,
      phylum = .data$tax_phylum,
      class = .data$tax_class,
      order = .data$tax_order,
      family = .data$tax_family,
      genus = .data$tax_genus,
      species = .data$tax_species,
      completeness = .data$completeness,
      contamination = .data$contamination,
      genome_size = .data$size_bp,
      genome_gc = .data$gc,
      genome_n50 = .data$n50,
      coding_density = .data$coding_density,
      contigs = .data$contigs,
      Release = release
    ) |>
    dplyr::arrange(.data$genome_id)

  list(
    datasets = datasets,
    rawdata_metadata = rawdata_metadata,
    sample_metadata = sample_metadata,
    capture_metadata = capture_metadata,
    specimen_metadata = specimen_metadata,
    sampling_metadata = sampling_metadata,
    host_metadata = host_metadata,
    assembly_metadata = assembly_metadata,
    genome_metadata = genome_metadata
  )
}

collect_ehi_core_release <- function(
    release = "EHR02",
    token = ehi_core_token(),
    base_url = Sys.getenv("EHI_CORE_URL", EHI_CORE_URL),
    progress = interactive(),
    individual_coassembly_species = EHI_INDIVIDUAL_COASSEMBLY_SPECIES) {
  .require_ehi_packages()
  if (!grepl("^EHR[0-9]+$", release)) {
    stop("release must look like EHR02.", call. = FALSE)
  }
  read_table <- function(table, filters = NULL) {
    ehi_core_read_table(
      table,
      token = token,
      base_url = base_url,
      filters = filters,
      progress = progress
    )
  }

  release_row <- read_table(
    "data_releases",
    filters = list(list(column = "code", op = "is", value = release))
  )
  if (nrow(release_row) != 1L) {
    stop("ehi-core has no unique data release named ", release, ".", call. = FALSE)
  }
  hologenomes <- read_table(
    "hologenomes",
    filters = list(list(column = "release_id", op = "is", value = release))
  )
  if (!nrow(hologenomes)) {
    stop("ehi-core has no hologenomes linked to ", release, ".", call. = FALSE)
  }

  preprocessings <- .read_codes(read_table, "preprocessings", "hologenome_id", hologenomes$code)
  samples <- .read_codes(read_table, "samples", "code", hologenomes$sample_id)
  captures <- .read_codes(read_table, "captures", "code", samples$capture_id)
  specimens <- .read_codes(read_table, "specimens", "code", captures$specimen_id)
  sampling_events <- .read_codes(
    read_table, "sampling_events", "code", captures$sampling_event_id
  )
  species <- .read_codes(read_table, "species", "code", specimens$species)

  assembly_samples <- .read_codes(
    read_table, "assembly_samples", "preprocessing_id", preprocessings$code
  )
  assemblies <- .read_codes(
    read_table, "assemblies", "code", assembly_samples$assembly_id
  )
  mags <- .read_codes(read_table, "mags", "assembly_id", assemblies$code)

  raw <- list(
    release = release_row,
    hologenomes = hologenomes,
    preprocessings = preprocessings,
    samples = samples,
    captures = captures,
    specimens = specimens,
    sampling_events = sampling_events,
    species = species,
    assemblies = assemblies,
    assembly_samples = assembly_samples,
    mags = mags
  )
  tables <- .build_ehi_release_tables(raw, release, individual_coassembly_species)
  structure(
    c(tables, list(core = raw, release = release, collected_at = Sys.time())),
    class = c("ehi_release_metadata", "list")
  )
}

.audit_issue <- function(severity, area, check, affected, detail) {
  tibble::tibble(
    severity = severity,
    area = area,
    check = check,
    affected = as.integer(affected),
    detail = detail
  )
}

audit_ehi_core_release <- function(metadata) {
  stopifnot(inherits(metadata, "ehi_release_metadata"))
  issues <- list()
  add <- function(severity, area, check, affected, detail) {
    if (affected > 0L) {
      issues[[length(issues) + 1L]] <<- .audit_issue(severity, area, check, affected, detail)
    }
  }

  rawdata <- metadata$rawdata_metadata
  preprocessing_required <- rawdata$preprocessing_required %in% TRUE
  add("error", "raw data", "missing sample link", sum(is.na(rawdata$sample_id)),
      "A hologenome cannot be joined to its biological sample.")
  add("error", "raw data", "missing required total GB",
      sum(preprocessing_required & is.na(rawdata$total_gb)),
      "Metagenomic data statistics require preprocessings.bases_pre_fastp.")
  add("error", "raw data", "missing required host GB",
      sum(preprocessing_required & is.na(rawdata$host_gb)),
      "Metagenomic composition statistics require preprocessings.host_bases.")
  add("error", "raw data", "missing required metagenomic GB",
      sum(preprocessing_required & is.na(rawdata$metagenomic_gb)),
      "Metagenomic composition statistics require preprocessings.metagenomic_bases.")
  add("info", "raw data", "genomic-only without preprocessing",
      sum(!preprocessing_required & is.na(rawdata$total_gb)),
      "Expected: genomic-only records are excluded from metagenomic composition statistics.")
  add("error", "raw data", "negative low-quality GB",
      sum(rawdata$lowqual_gb < -1e-6, na.rm = TRUE),
      "Total GB is smaller than host plus metagenomic GB; inspect preprocessing units.")
  add("warning", "raw data", "missing SingleM fraction",
      sum(preprocessing_required & is.na(rawdata$microbial_fraction)),
      "The export remains usable, but microbial-fraction summaries will be incomplete.")
  add("warning", "raw data", "missing Nonpareil completeness",
      sum(preprocessing_required & is.na(rawdata$nonpareil_completeness)),
      "The export remains usable, but Nonpareil summaries will be incomplete.")

  samples <- metadata$sample_metadata
  add("error", "fieldwork", "sample absent from ehi-core",
      sum(!rawdata$sample_id %in% samples$sample_id, na.rm = TRUE),
      "Every release hologenome's linked sample must be present.")
  add("error", "fieldwork", "sample missing capture", sum(is.na(samples$capture_id)),
      "Samples without captures cannot reach specimen or sampling-event metadata.")
  add("error", "fieldwork", "capture missing specimen",
      sum(is.na(metadata$capture_metadata$specimen_id)),
      "Host identity cannot be assigned.")
  add("error", "fieldwork", "capture missing sampling event",
      sum(is.na(metadata$core$captures$sampling_event_id)),
      "Country and environmental metadata cannot be assigned.")
  add("error", "fieldwork", "specimen missing species",
      sum(is.na(metadata$specimen_metadata$host_species)),
      "Host taxonomy and host-based plots require a species.")
  add("warning", "fieldwork", "sample missing type", sum(is.na(samples$sample_type)),
      "Sample-type counts will include an unknown category.")
  add("warning", "host taxonomy", "species missing order",
      sum(is.na(metadata$host_metadata$host_order)),
      "Host-order summaries will include an unknown category.")
  add("warning", "host taxonomy", "species missing class",
      sum(is.na(metadata$host_metadata$host_class)),
      "Host-class summaries and facets will include an unknown category.")
  add("warning", "fieldwork", "sampling event missing country",
      sum(is.na(metadata$sampling_metadata$country)),
      "Country counts will omit these sampling events.")
  add("warning", "fieldwork", "sampling event missing biome",
      sum(is.na(metadata$sampling_metadata$biome)),
      "Environmental representation will include an unknown biome.")

  add("error", "assemblies", "no individual assemblies", nrow(metadata$assembly_metadata) == 0L,
      "Assembly and MAG sections cannot be generated for this release.")
  add("error", "assemblies", "assembly missing hologenome",
      sum(is.na(metadata$assembly_metadata$rawdata_id)),
      "The assembly cannot be joined back to its release sample.")
  coassembly_specimens <- metadata$assembly_metadata |>
    dplyr::filter(.data$assembly_type == "Coassembly") |>
    tidyr::separate_rows("rawdata_id", sep = ", ") |>
    dplyr::left_join(
      metadata$rawdata_metadata |> dplyr::select("rawdata_id", "sample_id"),
      by = "rawdata_id"
    ) |>
    dplyr::left_join(
      metadata$sample_metadata |> dplyr::select("sample_id", "capture_id"),
      by = "sample_id"
    ) |>
    dplyr::left_join(
      metadata$capture_metadata |> dplyr::select("capture_id", "specimen_id"),
      by = "capture_id"
    ) |>
    dplyr::group_by(.data$assembly_id) |>
    dplyr::summarise(specimens = dplyr::n_distinct(.data$specimen_id), .groups = "drop")
  add("warning", "assemblies", "individual coassembly spans several specimens",
      sum(coassembly_specimens$specimens > 1L),
      paste0(
        "These coassemblies are released as single-individual assemblies, but their samples ",
        "link to more than one specimen. Check the sample -> capture -> specimen links."
      ))
  selected_assemblies <- metadata$core$assemblies |>
    dplyr::filter(.data$code %in% metadata$assembly_metadata$assembly_id)
  add("warning", "assemblies", "missing or inconsistent explicit release",
      sum(is.na(selected_assemblies$release_code) |
            selected_assemblies$release_code != metadata$release),
      paste0(
        "Release membership was derived through assembly sample -> preprocessing -> hologenome. ",
        "Populate assemblies.release_code so other core exports can identify ", metadata$release, " MAGs."
      ))
  add("error", "MAGs", "no MAGs", nrow(metadata$genome_metadata) == 0L,
      "MAG statistics and visualisations cannot be generated for this release.")
  for (column in c("completeness", "contamination", "genome_size", "phylum", "genus")) {
    add("warning", "MAGs", paste("missing", column),
        sum(is.na(metadata$genome_metadata[[column]])),
        paste("MAG summaries or plots using", column, "will be incomplete."))
  }

  # Known schema gaps: retained explicitly so a successful audit cannot be
  # mistaken for full parity with the historical Airtable sample export.
  for (column in c("sample_origin", "freshness", "time_to_freezing")) {
    add("schema_gap", "sample metadata", column, nrow(samples),
        "Not stored in ehi-core. The current statistics and visualisations do not use it.")
  }

  if (!length(issues)) {
    return(.audit_issue("ok", "release", "all checks passed", 0L, "No issues found."))
  }
  dplyr::bind_rows(issues)
}

write_ehi_core_release <- function(
    metadata,
    directory = "data",
    write_audit = TRUE,
    strict = TRUE) {
  stopifnot(inherits(metadata, "ehi_release_metadata"))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  release <- metadata$release
  audit <- audit_ehi_core_release(metadata)
  if (write_audit) {
    readr::write_tsv(
      audit,
      file.path(directory, paste0(release, "_metadata_audit.tsv"))
    )
  }
  errors <- audit[audit$severity == "error", , drop = FALSE]
  if (strict && nrow(errors)) {
    stop(
      release, " has ", nrow(errors),
      " blocking metadata checks. The audit was written, but release tables were not replaced. ",
      "Resolve the errors or call write_ehi_core_release(..., strict = FALSE) for inspection only.",
      call. = FALSE
    )
  }
  tables <- c(
    rawdata = "rawdata_metadata",
    sample = "sample_metadata",
    capture = "capture_metadata",
    specimen = "specimen_metadata",
    sampling = "sampling_metadata",
    host = "host_metadata",
    assembly = "assembly_metadata",
    genome = "genome_metadata"
  )
  purrr::iwalk(tables, function(object, suffix) {
    readr::write_tsv(metadata[[object]], file.path(directory, paste0(release, "_", suffix, ".tsv")))
  })

  readr::write_lines(
    sort(unique(stats::na.omit(metadata$host_metadata$host_species))),
    file.path(directory, paste0(release, "_allhosts.txt"))
  )
  readr::write_tsv(
    metadata$host_metadata,
    file.path(directory, paste0(release, "_allhosts.tsv"))
  )
  invisible(metadata)
}

print.ehi_release_metadata <- function(x, ...) {
  cat(
    x$release, " collected from ehi-core at ", format(x$collected_at, tz = "UTC"), " UTC\n",
    "  hologenomes: ", nrow(x$rawdata_metadata), "\n",
    "  samples:     ", nrow(x$sample_metadata), "\n",
    "  assemblies:  ", nrow(x$assembly_metadata), "\n",
    "  MAGs:        ", nrow(x$genome_metadata), "\n",
    sep = ""
  )
  invisible(x)
}

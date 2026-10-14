#!/usr/bin/env Rscript

# fetch_SRA Pipeline
# Author: Moises A. Rojas
# Batch version of the original R Markdown workflow.

suppressMessages(library(tidyverse))
suppressMessages(library(purrr))
suppressMessages(library(GEOquery))
suppressMessages(library(rentrez))
suppressMessages(library(recount3))
suppressMessages(library(edgeR))

# Input studies: `accession` (GSE ID) is required; `disease` is retained as the target cohort label.
input_studies <- read.csv(
  "config/input_studies.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("accession", "disease") %in% names(input_studies))) {
  stop("config/input_studies.csv must contain accession and disease columns.", call. = FALSE)
}

if (!all(grepl("^GSE[0-9]+$", input_studies$accession))) {
  stop("All accession values must be valid GSE accessions.", call. = FALSE)
}

# Inspect one GSE without downloading count matrices.
inspect_gse <- function(gse_id) {
  stopifnot(grepl("^GSE[0-9]+$", gse_id))

  geo <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = tempdir())
  metadata <- map_dfr(geo, ~ Biobase::pData(.x) %>%
    tibble::rownames_to_column("sample_id"))

  supp <- tryCatch(
    getGEOSuppFiles(gse_id, fetch_files = FALSE, baseDir = tempdir()),
    error = function(e) NULL
  )

  count_files <- if (!is.null(supp)) {
    supp[grepl("raw_counts|counts|count_matrix", rownames(supp), ignore.case = TRUE), , drop = FALSE]
  } else {
    NULL
  }

  list(
    accession = gse_id,
    metadata = metadata,
    n_samples = nrow(metadata),
    supplementary = supp,
    count_files = count_files
  )
}

gse_info <- map(input_studies$accession, inspect_gse)

# Download NCBI-computed raw gene counts.
download_gse_counts <- function(
    gse_id,
    counts_dir = "data_ncbi/counts",
    annotation_dir = "data_ncbi/annotation",
    overwrite = FALSE
) {
  gse_id <- toupper(trimws(gse_id))
  stopifnot(grepl("^GSE[0-9]+$", gse_id))

  dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)

  counts_file <- paste0(gse_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  annotation_file <- "Human.GRCh38.p13.annot.tsv.gz"

  make_url <- function(filename, accession = NULL) {
    base_url <- "https://www.ncbi.nlm.nih.gov/geo/download/"
    params <- c(
      "format=file",
      "type=rnaseq_counts",
      if (!is.null(accession)) paste0("acc=", accession),
      paste0("file=", filename)
    )
    paste0(base_url, "?", paste(params, collapse = "&"))
  }

  counts_url <- make_url(counts_file, accession = gse_id)
  annotation_url <- make_url(annotation_file)
  cat("Counts:", counts_url, "\n")
  cat("Annotation:", annotation_url, "\n")

  counts_path <- file.path(counts_dir, counts_file)
  annotation_path <- file.path(annotation_dir, annotation_file)

  valid_gzip <- function(path) {
    if (!file.exists(path) || file.info(path)$size < 100) return(FALSE)
    con <- file(path, open = "rb")
    magic <- readBin(con, "raw", n = 2)
    close(con)
    if (!identical(magic, as.raw(c(0x1f, 0x8b)))) return(FALSE)
    tryCatch({
      con <- gzfile(path, "rt")
      on.exit(close(con), add = TRUE)
      length(readLines(con, n = 1, warn = FALSE)) == 1
    }, error = function(e) FALSE)
  }

  fetch_file <- function(url, dest) {
    if (!overwrite && valid_gzip(dest)) {
      message("Already downloaded: ", basename(dest))
      return(dest)
    }
    tmp <- paste0(dest, ".part")
    if (file.exists(tmp)) unlink(tmp)
    message("Downloading: ", basename(dest))
    tryCatch({
      download.file(url, destfile = tmp, method = "libcurl", mode = "wb", quiet = TRUE, timeout = 600)
      if (!valid_gzip(tmp)) stop("Invalid gzip file or unavailable resource")
      if (file.exists(dest)) unlink(dest)
      if (!file.rename(tmp, dest)) stop("Could not move downloaded file")
      dest
    }, error = function(e) {
      if (file.exists(tmp)) unlink(tmp)
      warning("Failed: ", basename(dest), " — ", conditionMessage(e))
      NA_character_
    })
  }

  counts_result <- fetch_file(counts_url, counts_path)
  # annotation_result <- fetch_file(annotation_url, annotation_path)
  tibble(
    accession = gse_id,
    counts_path = counts_result,
    counts_downloaded = !is.na(counts_result)
  )
}

download_result <- map(input_studies$accession, download_gse_counts)

# Resolve a GEO series to its SRA project.
resolve_srp_entrez <- function(gse_id) {
  gds <- entrez_search(db = "gds", term = paste0(gse_id, "[ACCN] AND gse[ETYP]"))
  if (length(gds$ids) == 0L) return(NA_character_)

  links <- entrez_link(dbfrom = "gds", db = "sra", id = gds$ids, cmd = "neighbor")
  sra_ids <- links$links$gds_sra
  if (is.null(sra_ids) || length(sra_ids) == 0L) return(NA_character_)

  csv_text <- entrez_fetch(db = "sra", id = sra_ids, rettype = "runinfo", retmode = "text")
  runinfo <- read_csv(I(csv_text), col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  if (!"SRAStudy" %in% names(runinfo)) return(NA_character_)

  srp <- unique(runinfo$SRAStudy)
  srp <- srp[!is.na(srp) & grepl("^SRP[0-9]+$", srp)]
  if (length(srp) == 1L) return(srp)
  if (length(srp) > 1L) warning("Multiple SRP accessions found for ", gse_id, ": ", paste(srp, collapse = ", "))
  NA_character_
}

resolve_sra_accession <- function(accession, strict = FALSE) {
  accession <- toupper(trimws(accession))
  fail <- function(message) {
    if (strict) stop(message, call. = FALSE)
    warning(message, call. = FALSE)
    NA_character_
  }
  if (length(accession) != 1L || is.na(accession) || !grepl("^(GSE|SRP)[0-9]+$", accession)) {
    return(fail("Invalid GSE/SRP accession format"))
  }
  if (grepl("^SRP[0-9]+$", accession)) return(accession)

  message("Resolving GSE to SRP: ", accession)
  srp <- tryCatch(resolve_srp_entrez(accession), error = function(e) {
    warning("Entrez lookup failed: ", conditionMessage(e))
    NA_character_
  })
  if (length(srp) == 1L && !is.na(srp)) {
    message(accession, " -> ", srp)
    return(srp)
  }

  message("Entrez lookup unsuccessful. Trying GEOquery...")
  geo <- tryCatch(
    getGEO(accession, GSEMatrix = FALSE, getGPL = FALSE, destdir = tempdir()),
    error = function(e) { message("GEOquery lookup failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(geo)) return(fail(paste("Could not retrieve GEO accession:", accession)))
  relations <- GEOquery::Meta(geo)$relation
  if (is.null(relations)) return(fail(paste("No SRA relation found for:", accession)))
  relations <- unlist(relations, use.names = FALSE)
  srp <- unique(unlist(regmatches(relations, gregexpr("SRP[0-9]+", relations)), use.names = FALSE))
  if (length(srp) == 0L) return(fail(paste("No linked SRP found for:", accession)))
  if (length(srp) > 1L) return(fail(paste("Multiple SRP accessions found for:", accession, paste(srp, collapse = ", "))))
  message(accession, " -> ", srp)
  srp
}

download_sra_metadata <- function(accession, outdir = "data_ncbi/metadata", overwrite = FALSE) {
  accession <- toupper(trimws(accession))
  stopifnot(grepl("^(GSE|SRP)[0-9]+$", accession))
  study_dir <- file.path(outdir, accession)
  dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(study_dir, "SraRunInfo.csv")
  if (file.exists(dest) && !overwrite) {
    message("Already downloaded: ", dest)
    return(dest)
  }
  search <- entrez_search(db = "sra", term = paste0(accession, "[All Fields]"), use_history = TRUE, retmax = 0)
  if (search$count == 0) stop("No SRA records found: ", accession)
  csv_text <- entrez_fetch(db = "sra", web_history = search$web_history, rettype = "runinfo", retmode = "text")
  metadata <- read_csv(I(csv_text), col_types = cols(.default = "c"), show_col_types = FALSE)
  if (!"Run" %in% names(metadata) || nrow(metadata) == 0) stop("Invalid RunInfo response: ", accession)
  write_csv(metadata, dest)
  message("Saved ", nrow(metadata), " sequencing runs: ", dest)
  dest
}

download_sra_run_table <- function(accession, outdir = "data_ncbi/metadata", batch_size = 100L, overwrite = FALSE) {
  accession <- toupper(trimws(accession))
  stopifnot(grepl("^(GSE|SRP)[0-9]+$", accession))
  study_dir <- file.path(outdir, accession)
  dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(study_dir, "SraRunTable.csv")
  if (file.exists(dest) && !overwrite) {
    message("Already downloaded: ", dest)
    return(dest)
  }
  search <- entrez_search(db = "sra", term = paste0(accession, "[All Fields]"), retmax = 0, use_history = TRUE)
  if (search$count == 0) stop("No SRA records found for ", accession)

  parse_package <- function(pkg) {
    experiment <- xml2::xml_attr(xml2::xml_find_first(pkg, "./EXPERIMENT"), "accession")
    runs <- xml2::xml_attr(xml2::xml_find_all(pkg, "./RUN_SET/RUN"), "accession")
    if (length(runs) == 0) return(NULL)
    nodes <- xml2::xml_find_all(pkg, "./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")
    tags <- xml2::xml_text(xml2::xml_find_all(nodes, "./TAG"))
    values <- xml2::xml_text(xml2::xml_find_all(nodes, "./VALUE"))
    attrs <- tibble(name = tags, value = values) %>%
      dplyr::filter(!is.na(name), nzchar(name)) %>%
      dplyr::mutate(name = trimws(name), value = trimws(value)) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(value = paste(unique(value), collapse = "; "), .groups = "drop")
    wide <- if (nrow(attrs) > 0) {
      pivot_wider(attrs, names_from = name, values_from = value, names_prefix = "attr_")
    } else {
      tibble(.placeholder = NA_character_) %>% dplyr::select(-.placeholder)
    }
    tibble(Run = runs, Experiment = experiment) %>% bind_cols(wide[rep(1, length(runs)), , drop = FALSE])
  }

  starts <- seq(0, search$count - 1, by = batch_size)
  results <- map_dfr(starts, function(start) {
    message("Fetching records ", start + 1, "–", min(start + batch_size, search$count))
    xml_text <- rentrez::entrez_fetch(db = "sra", web_history = search$web_history, rettype = "xml", retmode = "xml", retstart = start, retmax = batch_size)
    doc <- xml2::read_xml(xml_text)
    packages <- xml2::xml_find_all(doc, ".//EXPERIMENT_PACKAGE")
    map_dfr(packages, parse_package)
  })
  if (nrow(results) == 0) stop("No run metadata extracted")
  runinfo_path <- download_sra_metadata(accession, outdir = outdir)
  runinfo <- read_csv(runinfo_path, col_types = cols(.default = "c"))
  results <- runinfo %>% left_join(results %>% dplyr::select(-Experiment), by = "Run")
  write_csv(results, dest)
  message("Saved ", nrow(results), " runs: ", dest)
  dest
}

input_studies <- input_studies %>%
  dplyr::mutate(srp_id = map_chr(accession, resolve_sra_accession)) %>%
  dplyr::mutate(metadata_path = purrr::map_chr(srp_id, \(x) {
    if (is.na(x)) NA_character_ else download_sra_run_table(x)
  }))

# Read and align every study before constructing the merged edgeR object.
read_study_data <- function(accession, disease, srp_id, metadata_path, counts_dir = "data_ncbi/counts") {
  counts_path <- file.path(counts_dir, paste0(accession, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
  if (!file.exists(counts_path)) stop("Count matrix not found: ", counts_path)
  counts_df <- read_tsv(counts_path, show_col_types = FALSE)
  gene_ids <- as.character(counts_df[[1]])
  if (anyNA(gene_ids) || anyDuplicated(gene_ids)) stop("Missing or duplicated gene ids: ", accession)
  count_data <- counts_df[-1]
  if (!all(vapply(count_data, is.numeric, logical(1)))) stop("Non-numeric count columns: ", accession)
  counts <- as.matrix(count_data)
  rownames(counts) <- gene_ids
  stopifnot(!anyNA(counts), all(is.finite(counts)), all(counts >= 0), all(counts == round(counts)), !anyDuplicated(colnames(counts)))

  if (!file.exists(metadata_path)) stop("Metadata not found: ", metadata_path)
  metadata <- read_csv(metadata_path, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  if (!"SampleName" %in% names(metadata)) stop("SampleName column missing: ", accession)
  if (anyNA(metadata$SampleName) || anyDuplicated(metadata$SampleName)) {
    stop("Missing or duplicated SampleName values in ", accession, ". Resolve run-to-sample relationships first.")
  }
  idx <- match(colnames(counts), metadata$SampleName)
  if (anyNA(idx)) stop("Samples without matching metadata in ", accession, ": ", paste(colnames(counts)[is.na(idx)], collapse = ", "))
  metadata <- metadata[idx, , drop = FALSE]
  stopifnot(identical(colnames(counts), metadata$SampleName))
  metadata$Study <- accession
  metadata$SRAStudy <- srp_id
  metadata$target_disease <- disease
  list(counts = counts, metadata = metadata)
}

study_data <- pmap(input_studies %>% dplyr::select(accession, disease, srp_id, metadata_path), read_study_data)
names(study_data) <- input_studies$accession

gene_lists <- map(study_data, ~ rownames(.x$counts))
gene_lists <- Reduce(intersect, gene_lists)
message("Common genes across studies: ", length(gene_lists))
if (length(gene_lists) == 0L) stop("No common gene ids across studies.")

counts_merged <- map(study_data, ~ .x$counts[gene_lists, , drop = FALSE])
counts_merged <- do.call(cbind, counts_merged)
stopifnot(!anyDuplicated(colnames(counts_merged)), !anyDuplicated(rownames(counts_merged)), all(is.finite(counts_merged)), all(counts_merged >= 0))

metadata_merged <- study_data %>% map("metadata") %>% bind_rows()
metadata_merged <- metadata_merged[match(colnames(counts_merged), metadata_merged$SampleName), , drop = FALSE]

if ("attr_diagnosis" %in% names(metadata_merged)) {
  print(table(metadata_merged$Study, metadata_merged$attr_diagnosis, useNA = "ifany"))
}
stopifnot(!anyNA(metadata_merged$SampleName), identical(colnames(counts_merged), metadata_merged$SampleName))
metadata_merged <- column_to_rownames(metadata_merged, "SampleName")

annotation_path <- file.path("data_ncbi", "annotation", "Human.GRCh38.p13.annot.tsv.gz")
if (!file.exists(annotation_path)) stop("Annotation file not found: ", annotation_path)
annotation <- read_tsv(annotation_path, show_col_types = FALSE)
annotation_ids <- as.character(annotation[["GeneID"]])
if (anyDuplicated(annotation_ids)) stop("Duplicated gene ids in annotation.")
idx <- match(rownames(counts_merged), annotation_ids)
message("Annotated genes: ", sum(!is.na(idx)), " / ", nrow(counts_merged))
if (all(is.na(idx))) stop("No gene ids matched the annotation.")
annotation <- as.data.frame(annotation[idx, , drop = FALSE])
rownames(annotation) <- rownames(counts_merged)
annotation <- annotation %>% dplyr::mutate(GeneID = as.character(GeneID))
stopifnot(identical(colnames(counts_merged), rownames(metadata_merged)), identical(rownames(counts_merged), rownames(annotation)))

dge_raw <- edgeR::DGEList(counts = counts_merged, samples = metadata_merged, genes = annotation)

outdir <- "out_edgeR"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
saveRDS(dge_raw, file.path(outdir, "dge_raw.rds"), compress = "bzip2")
message("Saved merged DGE object: ", nrow(dge_raw), " genes and ", ncol(dge_raw), " samples from ", length(unique(dge_raw$samples$Study)), " studies.")

# `dge_raw` is one merged object; summarize its samples by study rather than
# referring to the undefined per-study object list from the R Markdown draft.
dge_summary <- dge_raw$samples %>%
  rownames_to_column("SampleName") %>%
  dplyr::group_by(Study, SRAStudy, target_disease) %>%
  dplyr::summarise(
    n_samples = n(),
    n_genes = nrow(dge_raw),
    total_counts = sum(lib.size),
    median_libsize = median(lib.size),
    .groups = "drop"
  ) %>%
  dplyr::rename(accession = Study, srp_id = SRAStudy)
print(dge_summary)
write_csv(dge_summary, file.path(outdir, "dge_raw_summary.csv"))

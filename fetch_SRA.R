#!/usr/bin/env Rscript

# fetch_SRA Pipeline
# Author: Moises A. Rojas
# Batch version of the original R Markdown workflow.
# Run: Rscript fetch_SRA.R config/input_studies.csv

suppressMessages(library(tidyverse))
suppressMessages(library(GEOquery))
suppressMessages(library(rentrez))
suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args)) args[[1]] else "config/input_studies.csv"
if (!file.exists(config_path)) stop("Configuration not found: ", config_path, call. = FALSE)
input_studies <- read_csv(config_path, show_col_types = FALSE)
if (!"accession" %in% names(input_studies)) stop("Configuration needs an accession column.", call. = FALSE)
if (!"disease" %in% names(input_studies)) input_studies$disease <- NA_character_
input_studies <- input_studies %>% mutate(accession = toupper(trimws(accession)))
if (!all(grepl("^GSE[0-9]+$", input_studies$accession))) stop("All accession values must be GSE accessions.", call. = FALSE)

inspect_gse <- function(gse_id) {
  geo <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = tempdir())
  metadata <- map_dfr(geo, ~ Biobase::pData(.x) %>% tibble::rownames_to_column("sample_id"))
  geo_samples <- Biobase::pData(geo[[1]]) %>% tibble::rownames_to_column("GSM") %>% select(GSM, title)
  supp <- tryCatch(getGEOSuppFiles(gse_id, fetch_files = FALSE, baseDir = tempdir()), error = function(e) NULL)
  list(accession = gse_id, metadata = metadata, n_samples = nrow(metadata), geo_samples = geo_samples,
       supplementary = supp, count_files = if (!is.null(supp)) supp[grepl("raw_counts|counts|count_matrix", rownames(supp), ignore.case = TRUE), , drop = FALSE] else NULL)
}

download_gse_counts <- function(gse_id, counts_dir = "data_ncbi/counts", annotation_dir = "data_ncbi/annotation", overwrite = FALSE) {
  gse_id <- toupper(trimws(gse_id)); stopifnot(grepl("^GSE[0-9]+$", gse_id))
  dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE); dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
  make_url <- function(filename, accession = NULL) paste0("https://www.ncbi.nlm.nih.gov/geo/download/?", paste(c("format=file", "type=rnaseq_counts", if (!is.null(accession)) paste0("acc=", accession), paste0("file=", filename)), collapse = "&"))
  valid_gzip <- function(path) {
    if (!file.exists(path) || file.info(path)$size < 100) return(FALSE)
    con <- file(path, "rb"); on.exit(close(con), add = TRUE)
    if (!identical(readBin(con, "raw", n = 2), as.raw(c(0x1f, 0x8b)))) return(FALSE)
    tryCatch({ z <- gzfile(path, "rt"); on.exit(close(z), add = TRUE); length(readLines(z, n = 1, warn = FALSE)) == 1 }, error = function(e) FALSE)
  }
  fetch_file <- function(url, dest) {
    if (!overwrite && valid_gzip(dest)) return(dest)
    tmp <- paste0(dest, ".part"); if (file.exists(tmp)) unlink(tmp)
    tryCatch({ download.file(url, tmp, method = "libcurl", mode = "wb", quiet = TRUE, timeout = 600); if (!valid_gzip(tmp)) stop("Invalid gzip file or unavailable resource"); if (file.exists(dest)) unlink(dest); if (!file.rename(tmp, dest)) stop("Could not move downloaded file"); dest }, error = function(e) { if (file.exists(tmp)) unlink(tmp); warning("Failed: ", basename(dest), " — ", conditionMessage(e)); NA_character_ })
  }
  counts_file <- paste0(gse_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  annotation_file <- "Human.GRCh38.p13.annot.tsv.gz"
  tibble(accession = gse_id,
    counts_path = fetch_file(make_url(counts_file, gse_id), file.path(counts_dir, counts_file)),
    annotation_path = fetch_file(make_url(annotation_file), file.path(annotation_dir, annotation_file)))
}

resolve_srp_entrez <- function(gse_id) {
  gds <- entrez_search(db = "gds", term = paste0(gse_id, "[ACCN] AND gse[ETYP]"))
  if (!length(gds$ids)) return(NA_character_)
  ids <- entrez_link(dbfrom = "gds", db = "sra", id = gds$ids, cmd = "neighbor")$links$gds_sra
  if (is.null(ids) || !length(ids)) return(NA_character_)
  runinfo <- read_csv(I(entrez_fetch(db = "sra", id = ids, rettype = "runinfo", retmode = "text")), col_types = cols(.default = "c"), show_col_types = FALSE)
  srp <- unique(runinfo$SRAStudy[!is.na(runinfo$SRAStudy) & grepl("^SRP[0-9]+$", runinfo$SRAStudy)])
  if (length(srp) == 1L) srp else NA_character_
}

resolve_sra_accession <- function(accession, strict = FALSE) {
  accession <- toupper(trimws(accession)); fail <- function(x) { if (strict) stop(x, call. = FALSE); warning(x, call. = FALSE); NA_character_ }
  if (length(accession) != 1L || is.na(accession) || !grepl("^(GSE|SRP)[0-9]+$", accession)) return(fail("Invalid GSE/SRP accession format"))
  if (grepl("^SRP", accession)) return(accession)
  srp <- tryCatch(resolve_srp_entrez(accession), error = function(e) NA_character_)
  if (!is.na(srp)) return(srp)
  geo <- tryCatch(getGEO(accession, GSEMatrix = FALSE, getGPL = FALSE, destdir = tempdir()), error = function(e) NULL)
  relations <- if (is.null(geo)) character() else unlist(GEOquery::Meta(geo)$relation)
  srp <- unique(unlist(regmatches(relations, gregexpr("SRP[0-9]+", relations))))
  if (length(srp) != 1L) return(fail(paste("Could not resolve one linked SRP for", accession)))
  srp
}

download_sra_metadata <- function(accession, outdir = "data_ncbi/metadata", overwrite = FALSE) {
  dest <- file.path(outdir, accession, "SraRunInfo.csv"); dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && !overwrite) return(dest)
  search <- entrez_search(db = "sra", term = paste0(accession, "[All Fields]"), use_history = TRUE, retmax = 0)
  if (!search$count) stop("No SRA records found: ", accession)
  x <- read_csv(I(entrez_fetch(db = "sra", web_history = search$web_history, rettype = "runinfo", retmode = "text")), col_types = cols(.default = "c"), show_col_types = FALSE)
  if (!"Run" %in% names(x) || !nrow(x)) stop("Invalid RunInfo response: ", accession)
  write_csv(x, dest); dest
}

download_sra_run_table <- function(accession, outdir = "data_ncbi/metadata", batch_size = 100L, overwrite = FALSE) {
  dest <- file.path(outdir, accession, "SraRunTable.csv"); dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && !overwrite) return(dest)
  search <- entrez_search(db = "sra", term = paste0(accession, "[All Fields]"), retmax = 0, use_history = TRUE)
  if (!search$count) stop("No SRA records found for ", accession)
  parse_pkg <- function(pkg) {
    runs <- xml2::xml_attr(xml2::xml_find_all(pkg, "./RUN_SET/RUN"), "accession"); if (!length(runs)) return(NULL)
    nodes <- xml2::xml_find_all(pkg, "./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")
    attrs <- tibble(name = xml2::xml_text(xml2::xml_find_all(nodes, "./TAG")), value = xml2::xml_text(xml2::xml_find_all(nodes, "./VALUE"))) %>% filter(!is.na(name), nzchar(name)) %>% mutate(name = trimws(name), value = trimws(value)) %>% group_by(name) %>% summarise(value = paste(unique(value), collapse = "; "), .groups = "drop")
    wide <- if (nrow(attrs)) pivot_wider(attrs, names_from = name, values_from = value, names_prefix = "attr_") else tibble()
    tibble(Run = runs) %>% bind_cols(wide[rep(1, length(runs)), , drop = FALSE])
  }
  attrs <- map_dfr(seq(0, search$count - 1, by = batch_size), function(start) {
    xml <- entrez_fetch(db = "sra", web_history = search$web_history, rettype = "xml", retmode = "xml", retstart = start, retmax = batch_size)
    map_dfr(xml2::xml_find_all(xml2::read_xml(xml), ".//EXPERIMENT_PACKAGE"), parse_pkg)
  })
  if (!nrow(attrs)) stop("No run metadata extracted")
  out <- left_join(read_csv(download_sra_metadata(accession, outdir), col_types = cols(.default = "c"), show_col_types = FALSE), attrs, by = "Run")
  write_csv(out, dest); dest
}

collapse_sra_metadata <- function(metadata) {
  stopifnot("SampleName" %in% names(metadata))
  if (anyNA(metadata$SampleName) || any(!nzchar(trimws(metadata$SampleName)))) stop("Missing or empty SampleName values.")
  for (col in c(intersect("BioSample", names(metadata)), grep("^attr_", names(metadata), value = TRUE))) {
    bad <- metadata %>% group_by(SampleName) %>% summarise(n = n_distinct(.data[[col]][!is.na(.data[[col]]) & nzchar(trimws(.data[[col]]))]), .groups = "drop") %>% filter(n > 1)
    if (nrow(bad)) stop("Conflicting values in ", col, " for SampleName: ", paste(bad$SampleName, collapse = ", "))
  }
  metadata %>% group_by(SampleName) %>% summarise(across(everything(), ~ { x <- unique(.x[!is.na(.x) & nzchar(trimws(.x))]); if (!length(x)) NA_character_ else paste(x, collapse = "; ") }), .groups = "drop")
}

read_study_data <- function(accession, disease, srp_id, metadata_path, geo_samples, counts_dir = "data_ncbi/counts") {
  p <- file.path(counts_dir, paste0(accession, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")); if (!file.exists(p)) stop("Count matrix not found: ", p)
  x <- read_tsv(p, show_col_types = FALSE); genes <- as.character(x[[1]]); if (anyNA(genes) || anyDuplicated(genes)) stop("Missing or duplicated gene ids: ", accession)
  counts <- as.matrix(x[-1]); rownames(counts) <- genes
  if (!all(vapply(x[-1], is.numeric, logical(1))) || anyNA(counts) || any(counts < 0) || any(counts != round(counts)) || anyDuplicated(colnames(counts))) stop("Invalid raw count matrix: ", accession)
  metadata <- read_csv(metadata_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  if (!"SampleName" %in% names(metadata)) stop("SampleName column missing: ", accession)
  if (nrow(metadata) > n_distinct(metadata$SampleName)) metadata <- collapse_sra_metadata(metadata)
  stopifnot(!anyNA(metadata$SampleName), !anyDuplicated(metadata$SampleName), all(c("GSM", "title") %in% names(geo_samples)), !anyDuplicated(geo_samples$GSM))
  idx_geo <- match(metadata$SampleName, geo_samples$GSM); metadata$attr_geo_title <- geo_samples$title[idx_geo]
  message(accession, ": ", sum(!is.na(idx_geo)), " / ", nrow(metadata), " samples matched to GEO titles.")
  idx <- match(colnames(counts), metadata$SampleName); if (anyNA(idx)) stop("Samples without matching metadata in ", accession, ": ", paste(colnames(counts)[is.na(idx)], collapse = ", "))
  metadata <- metadata[idx, , drop = FALSE]; metadata$Study <- accession; metadata$SRAStudy <- srp_id; metadata$target_disease <- disease
  list(counts = counts, metadata = metadata)
}

gse_info <- map(input_studies$accession, inspect_gse); names(gse_info) <- input_studies$accession
input_studies <- input_studies %>% mutate(geo_samples = map(accession, ~ gse_info[[.x]]$geo_samples), srp_id = map_chr(accession, resolve_sra_accession), metadata_path = map_chr(srp_id, ~ if (is.na(.x)) NA_character_ else download_sra_run_table(.x)))
download_result <- map_dfr(input_studies$accession, download_gse_counts)
study_data <- pmap(input_studies %>% select(accession, disease, srp_id, metadata_path, geo_samples), read_study_data)
gene_ids <- Reduce(intersect, map(study_data, ~ rownames(.x$counts))); if (!length(gene_ids)) stop("No common gene ids across studies.")
counts_merged <- do.call(cbind, map(study_data, ~ .x$counts[gene_ids, , drop = FALSE]))
metadata_merged <- bind_rows(map(study_data, "metadata")); metadata_merged <- metadata_merged[match(colnames(counts_merged), metadata_merged$SampleName), , drop = FALSE]
stopifnot(!anyDuplicated(colnames(counts_merged)), identical(colnames(counts_merged), metadata_merged$SampleName))
metadata_merged <- tibble::column_to_rownames(metadata_merged, "SampleName")
annotation_path <- "data_ncbi/annotation/Human.GRCh38.p13.annot.tsv.gz"; if (!file.exists(annotation_path)) stop("Annotation not found: ", annotation_path)
annotation <- read_tsv(annotation_path, show_col_types = FALSE); annotation_ids <- as.character(annotation$GeneID); if (anyDuplicated(annotation_ids)) stop("Duplicated gene ids in annotation.")
idx <- match(rownames(counts_merged), annotation_ids); if (all(is.na(idx))) stop("No gene ids matched the annotation.")
annotation <- as.data.frame(annotation[idx, , drop = FALSE]); rownames(annotation) <- rownames(counts_merged); annotation$GeneID <- as.character(annotation$GeneID)
dge_raw <- edgeR::DGEList(counts = counts_merged, samples = metadata_merged, genes = annotation)
outdir <- "out_edgeR"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE); saveRDS(dge_raw, file.path(outdir, "dge_raw.rds"), compress = "bzip2")
dge_summary <- dge_raw$samples %>% tibble::rownames_to_column("SampleName") %>% group_by(Study, SRAStudy, target_disease) %>% summarise(n_samples = n(), n_genes = nrow(dge_raw), total_counts = sum(lib.size), median_libsize = median(lib.size), .groups = "drop") %>% rename(accession = Study, srp_id = SRAStudy)
write_csv(dge_summary, file.path(outdir, "dge_raw_summary.csv"))
message("Saved merged DGE object: ", nrow(dge_raw), " genes and ", ncol(dge_raw), " samples.")

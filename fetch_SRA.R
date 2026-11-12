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
  supp_names <- if (is.null(supp)) character() else if ("fname" %in% names(supp)) as.character(supp$fname) else basename(rownames(supp))
  list(accession = gse_id, metadata = metadata, n_samples = nrow(metadata), geo_samples = geo_samples,
       supplementary = supp, count_files = if (!is.null(supp)) supp[grepl("raw_counts|counts|count_matrix", supp_names, ignore.case = TRUE), , drop = FALSE] else NULL)
}

download_gse_counts <- function(gse_id, counts_dir = "data_ncbi/counts", annotation_dir = "data_ncbi/annotation", overwrite = FALSE) {
  gse_id <- toupper(trimws(gse_id)); stopifnot(grepl("^GSE[0-9]+$", gse_id))
  dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE); dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
  make_url <- function(filename, accession = NULL) paste0("https://www.ncbi.nlm.nih.gov/geo/download/?", paste(c("format=file", "type=rnaseq_counts", if (!is.null(accession)) paste0("acc=", accession), paste0("file=", filename)), collapse = "&"))
  valid_file <- function(path, compressed = grepl("\\.gz$", path, ignore.case = TRUE)) {
    if (!file.exists(path) || file.info(path)$size < 100) return(FALSE)
    if (!compressed) return(tryCatch(length(readLines(path, n = 1, warn = FALSE)) == 1, error = function(e) FALSE))
    con <- file(path, "rb"); on.exit(close(con), add = TRUE)
    if (!identical(readBin(con, "raw", n = 2), as.raw(c(0x1f, 0x8b)))) return(FALSE)
    tryCatch({ z <- gzfile(path, "rt"); on.exit(close(z), add = TRUE); length(readLines(z, n = 1, warn = FALSE)) == 1 }, error = function(e) FALSE)
  }
  fetch_file <- function(url, dest) {
    if (!overwrite && valid_file(dest)) return(dest)
    tmp <- paste0(dest, ".part"); if (file.exists(tmp)) unlink(tmp)
    tryCatch({ download.file(url, tmp, method = "libcurl", mode = "wb", quiet = TRUE, timeout = 600); if (!valid_file(tmp, grepl("\\.gz$", dest, ignore.case = TRUE))) stop("Invalid or unavailable resource"); if (file.exists(dest)) unlink(dest); if (!file.rename(tmp, dest)) stop("Could not move downloaded file"); dest }, error = function(e) { if (file.exists(tmp)) unlink(tmp); warning("Failed: ", basename(dest), " — ", conditionMessage(e)); NA_character_ })
  }
  fetch_supplementary_counts <- function() {
    supp <- tryCatch(getGEOSuppFiles(gse_id, fetch_files = FALSE, baseDir = tempdir()), error = function(e) NULL)
    if (is.null(supp) || !nrow(supp)) return(tibble(counts_path = NA_character_, counts_source = NA_character_, counts_file = NA_character_))
    filenames <- if ("fname" %in% names(supp)) as.character(supp$fname) else basename(rownames(supp))
    is_delimited <- grepl("\\.(csv|tsv|txt)(\\.gz)?$", filenames, ignore.case = TRUE)
    is_count <- grepl("(raw[_. -]*)?counts?|count[_. -]*matrix", filenames, ignore.case = TRUE)
    is_processed <- grepl("(tpm|fpkm|rpkm|cpm|normaliz|log[0-9]*|vst|rld)", filenames, ignore.case = TRUE)
    candidates <- filenames[is_delimited & is_count & !is_processed]
    if (!length(candidates)) return(tibble(counts_path = NA_character_, counts_source = NA_character_, counts_file = NA_character_))
    score <- ifelse(grepl("raw[_. -]*counts?", candidates, ignore.case = TRUE), 3L,
      ifelse(grepl("count[_. -]*matrix", candidates, ignore.case = TRUE), 2L, 1L))
    candidates <- candidates[score == max(score)]
    if (length(candidates) != 1L) {
      warning("Ambiguous supplementary count files for ", gse_id, ": ", paste(candidates, collapse = ", "), call. = FALSE)
      return(tibble(counts_path = NA_character_, counts_source = NA_character_, counts_file = NA_character_))
    }
    filename <- candidates[[1]]
    file_index <- match(filename, filenames)
    group <- sub("[0-9]{3}$", "nnn", gse_id)
    url <- if ("url" %in% names(supp) && !is.na(supp$url[[file_index]])) as.character(supp$url[[file_index]]) else
      paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", group, "/", gse_id, "/suppl/", utils::URLencode(filename, reserved = TRUE))
    dest <- file.path(counts_dir, filename)
    path <- fetch_file(url, dest)
    if (is.na(path)) return(tibble(counts_path = NA_character_, counts_source = NA_character_, counts_file = NA_character_))
    tibble(counts_path = path, counts_source = "geo_supplementary", counts_file = filename)
  }
  counts_file <- paste0(gse_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  annotation_file <- "Human.GRCh38.p13.annot.tsv.gz"
  ncbi_path <- fetch_file(make_url(counts_file, gse_id), file.path(counts_dir, counts_file))
  counts <- if (!is.na(ncbi_path)) tibble(counts_path = ncbi_path, counts_source = "ncbi_rnaseq_counts", counts_file = counts_file) else fetch_supplementary_counts()
  bind_cols(tibble(accession = gse_id), counts,
    tibble(annotation_path = fetch_file(make_url(annotation_file), file.path(annotation_dir, annotation_file))))
}

download_geo_supp_metadata <- function(gse_id, outdir = "data_ncbi/metadata", overwrite = FALSE) {
  supp <- tryCatch(getGEOSuppFiles(gse_id, fetch_files = FALSE, baseDir = tempdir()), error = function(e) NULL)
  if (is.null(supp) || !nrow(supp)) return(NA_character_)
  filenames <- if ("fname" %in% names(supp)) as.character(supp$fname) else basename(rownames(supp))
  candidate <- grepl("(metadata|sample[_ .-]*(information|info)?|phenotype)", filenames, ignore.case = TRUE) &
    grepl("\\.(csv|tsv|txt)(\\.gz)?$", filenames, ignore.case = TRUE)
  filenames <- filenames[candidate]
  if (length(filenames) != 1L) return(NA_character_)
  filename <- filenames[[1]]
  source_index <- match(filename, if ("fname" %in% names(supp)) as.character(supp$fname) else basename(rownames(supp)))
  url <- if ("url" %in% names(supp)) as.character(supp$url[[source_index]]) else NA_character_
  if (is.na(url) || !nzchar(url)) return(NA_character_)
  dest <- file.path(outdir, gse_id, filename); dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && !overwrite && file.info(dest)$size > 100) return(dest)
  tmp <- paste0(dest, ".part"); if (file.exists(tmp)) unlink(tmp)
  tryCatch({
    download.file(url, tmp, method = "libcurl", mode = "wb", quiet = TRUE, timeout = 600)
    if (!file.exists(tmp) || file.info(tmp)$size <= 100) stop("Downloaded file is empty")
    if (file.exists(dest)) unlink(dest)
    if (!file.rename(tmp, dest)) stop("Could not move downloaded metadata")
    dest
  }, error = function(e) { if (file.exists(tmp)) unlink(tmp); warning("Failed GEO metadata for ", gse_id, ": ", conditionMessage(e), call. = FALSE); NA_character_ })
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

read_count_matrix <- function(path, accession, source) {
  if (!file.exists(path)) stop("Count matrix not found: ", path)
  delim <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
  preview <- readLines(path, n = 2L, warn = FALSE)
  if (length(preview) < 2L) stop("Count matrix has no header and data row: ", path)
  n_fields <- function(line) length(strsplit(line, delim, fixed = TRUE)[[1]])
  # Some author matrices omit the gene-column label, leaving one fewer field
  # in the header than in each data row (e.g. GSE318560).
  if (n_fields(preview[[1]]) == n_fields(preview[[2]]) - 1L) {
    header <- read_delim(I(preview[[1]]), delim = delim, col_names = FALSE, show_col_types = FALSE) %>% unlist(use.names = FALSE) %>% as.character()
    x <- read_delim(path, delim = delim, skip = 1L, col_names = FALSE, show_col_types = FALSE, name_repair = "minimal")
    names(x) <- c("gene_id", header)
  } else {
    x <- read_delim(path, delim = delim, show_col_types = FALSE, name_repair = "minimal")
  }
  if (ncol(x) < 2L || !nrow(x)) stop("Count matrix has no gene column and sample columns: ", path)
  genes <- trimws(as.character(x[[1]]))
  is_ensembl <- mean(grepl("^ENS[A-Z]*[GTP][0-9]+(\\.[0-9]+)?$", genes)) > 0.9
  if (is_ensembl) genes <- sub("\\.[0-9]+$", "", genes)
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) stop("Missing or duplicated gene ids", if (is_ensembl) " after removing Ensembl version suffixes" else "", ": ", accession)
  counts <- as.matrix(x[-1]); rownames(counts) <- genes
  if (!all(vapply(x[-1], is.numeric, logical(1))) || anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) || anyDuplicated(colnames(counts))) stop("Invalid count matrix: ", accession)
  if (any(counts != round(counts))) {
    if (identical(source, "ncbi_rnaseq_counts")) stop("NCBI count matrix contains non-integer values: ", accession)
    warning(accession, ": GEO supplementary matrix has non-integer values; retaining them as author-supplied estimated counts.", call. = FALSE)
  }
  counts
}

harmonize_gene_ids <- function(study, annotation) {
  if (!"GeneID" %in% names(annotation)) stop("Annotation is missing its GeneID column.")
  annotation_ids <- trimws(as.character(annotation$GeneID))
  if (anyNA(annotation_ids) || any(!nzchar(annotation_ids)) || anyDuplicated(annotation_ids)) stop("Annotation has missing or duplicated GeneID values.")
  genes <- rownames(study$counts)
  unique_lookup <- function(keys, query, normalize = identity) {
    keys <- normalize(trimws(as.character(keys))); query <- normalize(trimws(as.character(query)))
    usable <- !is.na(keys) & nzchar(keys) & !duplicated(keys) & !duplicated(keys, fromLast = TRUE)
    annotation_ids[match(query, ifelse(usable, keys, NA_character_))]
  }
  # GeneID is the canonical output ID.  The other fields are only lookups.
  mapped_ids <- unique_lookup(annotation_ids, genes)
  mapping_label <- "GeneID"
  ensembl_column <- names(annotation)[tolower(names(annotation)) == "ensemblgeneid"]
  if (length(ensembl_column)) {
    needs_mapping <- is.na(mapped_ids)
    mapped_ids[needs_mapping] <- unique_lookup(annotation[[ensembl_column[[1]]]], sub("\\.[0-9]+$", "", genes[needs_mapping]))
    mapping_label <- paste0(mapping_label, ", EnsemblGeneID")
  }
  symbol_candidates <- names(annotation)[gsub("[^a-z0-9]", "", tolower(names(annotation))) %in% c("symbol", "genesymbol", "hgncsymbol")]
  if (length(symbol_candidates)) {
    needs_mapping <- is.na(mapped_ids)
    mapped_ids[needs_mapping] <- unique_lookup(annotation[[symbol_candidates[[1]]]], genes[needs_mapping], toupper)
    mapping_label <- paste0(mapping_label, ", ", symbol_candidates[[1]])
  }
  # Some NCBI annotations leave Symbol blank but begin Description with it.
  if ("Description" %in% names(annotation)) {
    needs_mapping <- is.na(mapped_ids)
    description_symbol <- sub("\\s+.*$", "", trimws(as.character(annotation$Description)))
    mapped_ids[needs_mapping] <- unique_lookup(description_symbol, genes[needs_mapping], toupper)
    mapping_label <- paste0(mapping_label, ", Description prefix")
  }
  keep <- !is.na(mapped_ids)
  if (!any(keep)) stop(study$metadata$Study[[1]], ": no author gene identifiers matched ", mapping_label, ".")
  if (sum(keep) < length(keep)) message(study$metadata$Study[[1]], ": retaining ", sum(keep), " / ", length(keep), " genes after mapping through ", mapping_label, ".")
  study$counts <- study$counts[keep, , drop = FALSE]
  rownames(study$counts) <- mapped_ids[keep]
  if (anyDuplicated(rownames(study$counts))) {
    before <- nrow(study$counts)
    study$counts <- rowsum(study$counts, group = rownames(study$counts), reorder = FALSE)
    message(study$metadata$Study[[1]], ": collapsed ", before - nrow(study$counts), " duplicate rows after GeneID mapping.")
  }
  study
}

match_count_columns <- function(count_names, metadata, geo_samples, accession) {
  canonical_sample_id <- function(x) {
    x <- toupper(gsub("[^A-Za-z0-9]", "", trimws(x)))
    sub("NX$", "EC", x)
  }
  idx <- match(count_names, metadata$SampleName)
  if (!anyNA(idx)) return(idx)
  alias <- bind_rows(
    transmute(geo_samples, alias = GSM, GSM),
    transmute(geo_samples, alias = title, GSM),
    transmute(geo_samples, alias = trimws(sub("[,;].*$", "", title)), GSM)
  ) %>% filter(!is.na(alias), nzchar(alias)) %>% group_by(alias) %>%
    summarise(GSM = if (n_distinct(GSM) == 1L) first(GSM) else NA_character_, .groups = "drop")
  unresolved <- which(is.na(idx))
  gsm <- alias$GSM[match(count_names[unresolved], alias$alias)]
  idx[unresolved] <- match(gsm, metadata$SampleName)
  # Authors often add punctuation or use NX where GEO titles use EC.  Apply
  # that normalization only when it identifies a unique GEO title.
  if (anyNA(idx)) {
    title_key <- canonical_sample_id(geo_samples$title)
    count_key <- canonical_sample_id(count_names[is.na(idx)])
    unique_titles <- !duplicated(title_key) & !duplicated(title_key, fromLast = TRUE)
    gsm <- geo_samples$GSM[match(count_key, ifelse(unique_titles, title_key, NA_character_))]
    idx[is.na(idx)] <- match(gsm, metadata$SampleName)
  }
  # A final one-character rescue is permitted only for an otherwise unique,
  # unused GEO title (e.g. author label 451.CA1 versus GEO title 45CA1).
  if (anyNA(idx)) {
    missing <- which(is.na(idx)); available <- setdiff(seq_len(nrow(metadata)), idx[!is.na(idx)])
    title_key <- canonical_sample_id(geo_samples$title)
    for (i in missing) {
      candidate_gsm <- geo_samples$GSM[adist(canonical_sample_id(count_names[[i]]), title_key) == 1L]
      candidate <- match(candidate_gsm, metadata$SampleName)
      candidate <- intersect(candidate, available)
      if (length(candidate) == 1L) idx[[i]] <- candidate
    }
  }
  if (anyNA(idx)) stop("Samples without matching metadata in ", accession, ": ", paste(count_names[is.na(idx)], collapse = ", "))
  if (anyDuplicated(idx)) stop("Multiple count columns map to the same metadata sample in ", accession)
  idx
}

read_study_data <- function(accession, disease, srp_id, metadata_path, geo_samples, counts_path, counts_source) {
  counts <- read_count_matrix(counts_path, accession, counts_source)
  metadata <- read_csv(metadata_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  if (!"SampleName" %in% names(metadata)) {
    normalized_names <- gsub("[^a-z0-9]", "", tolower(names(metadata)))
    sample_column <- names(metadata)[normalized_names == "samplename"]
    if (!length(sample_column)) sample_column <- names(metadata)[normalized_names == "sampleid"]
    if (length(sample_column) == 1L) names(metadata)[match(sample_column, names(metadata))] <- "SampleName"
  }
  if (!"SampleName" %in% names(metadata)) stop("SampleName column missing: ", accession)
  if (nrow(metadata) > n_distinct(metadata$SampleName)) metadata <- collapse_sra_metadata(metadata)
  stopifnot(!anyNA(metadata$SampleName), !anyDuplicated(metadata$SampleName), all(c("GSM", "title") %in% names(geo_samples)), !anyDuplicated(geo_samples$GSM))
  idx_geo <- match(metadata$SampleName, geo_samples$GSM)
  idx_title <- match(metadata$SampleName, geo_samples$title)
  idx_geo[is.na(idx_geo)] <- idx_title[is.na(idx_geo)]
  metadata$attr_geo_title <- geo_samples$title[idx_geo]
  message(accession, ": ", sum(!is.na(idx_geo)), " / ", nrow(metadata), " samples matched to GEO GSM IDs or titles.")
  idx <- match_count_columns(colnames(counts), metadata, geo_samples, accession)
  metadata <- metadata[idx, , drop = FALSE]
  # The author matrix can use a GEO-title alias (e.g. BBN IDs) while SRA
  # metadata uses GSM IDs.  Keep the latter, but use the count-column ID as
  # the merge key so metadata and counts stay in exactly the same order.
  metadata$SRA_SampleName <- metadata$SampleName
  metadata$SampleName <- colnames(counts)
  metadata$Study <- accession; metadata$SRAStudy <- srp_id; metadata$target_disease <- disease
  metadata$count_source <- counts_source
  list(counts = counts, metadata = metadata)
}

gse_info <- map(input_studies$accession, inspect_gse); names(gse_info) <- input_studies$accession
input_studies <- input_studies %>% mutate(geo_samples = map(accession, ~ gse_info[[.x]]$geo_samples), srp_id = map_chr(accession, resolve_sra_accession))
input_studies <- input_studies %>% mutate(
  geo_metadata_path = map2_chr(accession, srp_id, ~ if (is.na(.y)) download_geo_supp_metadata(.x) else NA_character_),
  metadata_path = map2_chr(srp_id, geo_metadata_path, ~ if (!is.na(.x)) download_sra_run_table(.x) else .y)
)
download_result <- map_dfr(input_studies$accession, download_gse_counts)
input_studies <- left_join(input_studies, download_result %>% select(accession, counts_path, counts_source, counts_file), by = "accession")
outdir <- "out_edgeR"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
skipped_studies <- input_studies %>% mutate(skip_reason = case_when(
  is.na(counts_path) ~ "No usable NCBI or unambiguous GEO supplementary delimited count matrix",
  is.na(metadata_path) ~ "SRA run metadata unavailable",
  TRUE ~ NA_character_
)) %>% filter(!is.na(skip_reason)) %>% select(accession, disease, skip_reason)
if (nrow(skipped_studies)) write_csv(skipped_studies, file.path(outdir, "skipped_studies.csv"))
input_studies <- input_studies %>% filter(!accession %in% skipped_studies$accession)
if (!nrow(input_studies)) stop("No studies with usable counts and metadata; see out_edgeR/skipped_studies.csv", call. = FALSE)
study_data <- pmap(input_studies %>% select(accession, disease, srp_id, metadata_path, geo_samples, counts_path, counts_source), read_study_data)
annotation_path <- "data_ncbi/annotation/Human.GRCh38.p13.annot.tsv.gz"; if (!file.exists(annotation_path)) stop("Annotation not found: ", annotation_path)
annotation <- read_tsv(annotation_path, show_col_types = FALSE)
study_data <- map(study_data, harmonize_gene_ids, annotation = annotation)
gene_ids <- Reduce(intersect, map(study_data, ~ rownames(.x$counts))); if (!length(gene_ids)) stop("No common gene ids across studies.")
counts_merged <- do.call(cbind, map(study_data, ~ .x$counts[gene_ids, , drop = FALSE]))
metadata_merged <- bind_rows(map(study_data, "metadata")); metadata_merged <- metadata_merged[match(colnames(counts_merged), metadata_merged$SampleName), , drop = FALSE]
stopifnot(!anyDuplicated(colnames(counts_merged)), identical(colnames(counts_merged), metadata_merged$SampleName))
metadata_merged <- tibble::column_to_rownames(metadata_merged, "SampleName")
annotation_ids <- trimws(as.character(annotation$GeneID)); if (anyDuplicated(annotation_ids)) stop("Duplicated GeneID values in annotation.")
idx <- match(rownames(counts_merged), annotation_ids); if (anyNA(idx)) stop("Merged count matrix contains gene IDs absent from annotation.")
annotation <- as.data.frame(annotation[idx, , drop = FALSE]); rownames(annotation) <- rownames(counts_merged); annotation$GeneID <- rownames(counts_merged)
dge_raw <- edgeR::DGEList(counts = counts_merged, samples = metadata_merged, genes = annotation)
saveRDS(dge_raw, file.path(outdir, "dge_raw.rds"), compress = "bzip2")
dge_summary <- dge_raw$samples %>% tibble::rownames_to_column("SampleName") %>% group_by(Study, SRAStudy, target_disease) %>% summarise(n_samples = n(), n_genes = nrow(dge_raw), total_counts = sum(lib.size), median_libsize = median(lib.size), .groups = "drop") %>% rename(accession = Study, srp_id = SRAStudy)
write_csv(dge_summary, file.path(outdir, "dge_raw_summary.csv"))
message("Saved merged DGE object: ", nrow(dge_raw), " genes and ", ncol(dge_raw), " samples.")

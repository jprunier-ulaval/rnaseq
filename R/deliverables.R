#' Produce deliverables
#'
#' counts.csv: counts at the gene level and at the transcript level
#' de_*.csv: DE analysis at the gene level and at the transcript level for
#'           specified contrasts
#'
#' @param dir_kallisto Directory with Kallisto quantifications
#' @param anno The version of the annotation to use. Default: "Hs.Ensembl91"
#'             Currently available:
#'             * Hs.Gencode19
#'             * Hs.Gencode27
#'             * Hs.Gencode32
#'             * Hs.Gencode35
#'             * Hs.Gencode37
#'             * Hs.Ensembl79
#'             * Hs.Ensembl91
#'             * Hs.Ensembl95
#'             * Hs.Ensembl97
#'             * Hs.Ensembl98
#'             * Hs.Ensembl100
#'             * Hs.Ensembl101
#'             * Mm.Ensembl91
#'             * Mm.Ensembl92
#'             * Mm.Ensembl94
#'             * Mm.Ensembl97
#'             * Mm.Ensembl99
#'             * Mm.Ensembl100
#'             * Mm.Ensembl101
#'             * Rn.Ensembl76
#'             * Rn.Ensembl79
#'             * Rn.Ensembl92
#'             * Rn.Ensembl98
#'             * Bt.Ensembl99
#'             * Mmu.Ensembl101
#'             * Mmu.Ensembl103
#'             * peaux_colonisees
#' @param design The experimental design
#'               (1st column: sample; 2nd column: group)
#' @param contrasts List of contrasts
#' @param dir_output Directory where to write outputs
#' @param file_type Abundance file format to use (h5 or tsv).
#' @param digits Integer indicating the number of decimal places
#' @param ignoreTxVersion Should the version number in the IDs be ignored
#'        during import. Default: \code{FALSE}.
#' @param use_ruv Use RUVg normalization. Default: \code{FALSE}.
#' @param housekeeping_genes A \code{vector} of gene symbols
#' @param ncores Number of cores to use for de analysis. Default \code{1}.
#'
#' @return Invisibly returns a list with txi_tx, txi_genes, df_tx (PCA),
#'         df_genes (PCA), design, contrasts, counts and de values.
#'
#' @examples
#' dir_kallisto <- get_demo_kallisto_dir()
#' contrasts <- list(comp = c("group", "A", "B"))
#' design <- get_demo_design()
#' file_anno <- get_demo_anno_file()
#' produce_deliverables(dir_kallisto,
#'                      anno = file_anno,
#'                      design = design,
#'                      contrast = contrasts,
#'                      dir_output = ".",
#'                      file_type = "tsv",
#'                      ignoreTxVersion = FALSE)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom purrr map
#' @importFrom purrr iwalk
#' @importFrom readr write_csv
#'
#' @export
produce_deliverables <- function (dir_kallisto, anno, design, contrasts,
                                  dir_output, file_type = "h5", digits = 4,
                                  ignoreTxVersion = FALSE, use_ruv = FALSE,
                                  housekeeping_genes = get_human_hsk(),
                                  ncores = 1) {

    stopifnot(dir.exists(dir_kallisto))
    stopifnot(dir.exists(dir_output))
    stopifnot(file_type %in% c("tsv", "h5"))
    stopifnot(colnames(design) == c("sample", "group"))
    if (!file.exists(anno)) {
        validate_anno(anno)
    }
    stopifnot(is.logical(ignoreTxVersion))
    stopifnot(is.logical(use_ruv))
    stopifnot(is.numeric(ncores))
    stopifnot(as.integer(ncores) == ncores)
    stopifnot(ncores > 0)

    # Import quantifications
    files <- get_filenames(dir_kallisto, file_type)
    txi <- produce_txi(files = files, anno = anno, ignoreTxVersion = ignoreTxVersion)

    # PCA
    pdf(file.path(dir_output, "PCA_genes.pdf"))
    df_genes <- produce_pca(txi$genes, use_ruv = use_ruv)
    dev.off()
    pdf(file.path(dir_output, "PCA_tx.pdf"))
    df_tx <- produce_pca(txi$tx, use_ruv = use_ruv)
    dev.off()

    # Produce counts
    counts_genes <- format_counts(txi$genes, digits = digits)
    counts_tx <- format_counts(txi$tx, digits = digits)
    counts <- rbind(counts_tx, counts_genes)
    readr::write_csv(counts, file.path(dir_output, "counts.csv"))

    # Produce DE
    de <- produce_de(txi, design,  ~ group, use_ruv = use_ruv, ncores = ncores)
    purrr::iwalk(de, ~ readr::write_csv(.x, file.path(dir_output, paste0(.y, ".csv"))))

    invisible(list(txi_tx = txi$tx,
                   txi_genes = txi$genes,
                   df_tx = df_tx,
                   df_genes = df_genes,
                   design = design,
                   contrasts = contrasts,
                   counts = counts,
                   de = de,
                   dds_genes = dds$genes,
                   dds_tx = dds$tx))
}

# TODO: params
# TODO: example
# TODO: stopifnot
#' Produce DE table
#'
#' The differential expression (DE) table will be produced by merging results
#' at the gene and transcripts level.
#'
#' @return A data.frame with the gene and tx DE merged in a single table
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom purrr map
#' @importFrom parallel mclapply
#'
#' @export
produce_de <- function(txi, design, formula = ~ group, use_ruv = FALSE, ncores = 1) {
    # Make sure design is in correct order
    stopifnot(all(as.character(design$sample) %in% colnames(txi$genes$counts)))
    design <- dplyr::mutate(design, sample = factor(sample, levels = colnames(txi$genes$counts))) %>%
                     dplyr::arrange(sample)
    stopifnot(identical(colnames(txi$genes$counts), as.character(design$sample)))

    dds <- list()
    if (ncores == 1) {
        dds[["genes"]] <- deseq2_analysis(txi$genes, design, formula, use_ruv = use_ruv)
        dds[["tx"]] <- deseq2_analysis(txi$tx, design, formula, use_ruv = use_ruv)
    } else {
        dds <- parallel::mclapply(list(txi$genes, txi$tx), function(x) deseq2_analysis(x, design, formula, use_ruv = use_ruv), mc.cores = 2)
        names(dds) <- c("genes", "tx")
    }

    if (ncores == 1) {
        de_genes <- purrr::map(contrasts, ~ format_de(dds$genes, txiDgenes, .x, ignoreTxVersion, digits = digits))
        de_tx <- purrr::map(contrasts, ~ format_de(dds$tx, txi$tx, .x, ignoreTxVersion, digits = digits))
    } else {
        de_genes <- parallel::mclapply(contrasts, function(x) format_de(dds$genes, txi$genes, x, ignoreTxVersion, digits = digits), mc.cores = ncores)
        de_tx <- parallel::mclapply(contrasts, function(x) format_de(dds$tx, txi$tx, x, ignoreTxVersion, digits = digits), mc.cores = ncores)
    }

    rbind_df <- function(n) {
        tx <- de_tx[[n]]
        genes <- de_genes[[n]]
        rbind(tx, genes)
    }
    de <- purrr::map(names(de_tx), rbind_df)
    names(de) <- names(de_tx)
    de
}

# TODO: params
# TODO: example
# TODO: stopifnot
#' Produce txi objects
#'
#' Produce the txi objects at the gene and at the transcript level.
#'
#' @return A list with the txi objects at the gene and at the transcript level
#'
#' @export
produce_txi <- function(files, anno, ignoreTxVersion = TRUE, use_ruv = FALSE) {
    txi <- list()
    txi$tx <- import_kallisto(files, anno = anno, txOut = TRUE, ignoreTxVersion = ignoreTxVersion)
    txi$genes <- summarize_to_gene(txi$tx, anno = anno, ignoreTxVersion = ignoreTxVersion)

    # RUV
    if (use_ruv) {
        txi$tx <- ruvg_normalization(txi$tx,
                                     housekeeping_genes = housekeeping_genes,
                                     ignoreTxVersion = ignoreTxVersion)
        txi$genes <- ruvg_normalization(txi$genes,
                                        housekeeping_genes = housekeeping_genes,
                                        ignoreTxVersion = ignoreTxVersion)
    }
    txi
}

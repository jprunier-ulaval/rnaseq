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
#' @import purrr
#' @import dplyr
#'
#' @export
produce_deliverables <- function (dir_kallisto, anno, design, contrasts,
                                  dir_output, file_type = "h5", digits = 4,
                                  ignoreTxVersion = FALSE, use_ruv = FALSE,
                                  housekeeping_genes = get_human_hsk()) {

    stopifnot(dir.exists(dir_kallisto))
    stopifnot(dir.exists(dir_output))
    stopifnot(file_type %in% c("tsv", "h5"))
    stopifnot(colnames(design) == c("sample", "group"))
    if (!file.exists(anno)) {
        validate_anno(anno)
    }
    stopifnot(is.logical(ignoreTxVersion))
    stopifnot(is.logical(use_ruv))

    file_type <-paste0(file_type, "$")
    files <- dir(dir_kallisto, pattern = file_type, recursive = TRUE, full.names = TRUE)
    names(files) <- basename(dirname(files))

    # Import quantifications
    txi_tx <- import_kallisto(files, anno = anno, txOut = TRUE, ignoreTxVersion = ignoreTxVersion)
    txi_genes <- summarize_to_gene(txi_tx, anno = anno, ignoreTxVersion = ignoreTxVersion)

    # RUV
    if (use_ruv) {
        txi_tx <- ruvg_normalization(txi_tx,
                                     housekeeping_genes = housekeeping_genes,
                                     ignoreTxVersion = TRUE)
        txi_genes <- ruvg_normalization(txi_genes,
                                        housekeeping_genes = housekeeping_genes,
                                        ignoreTxVersion = TRUE)
    }

    # Make sure design is in correct order
    stopifnot(all(as.character(design$sample) %in% colnames(txi_genes$counts)))
    design <- mutate(design, sample = factor(sample, levels = colnames(txi_genes$counts))) %>%
                     arrange(sample)
    stopifnot(identical(colnames(txi_genes$counts), as.character(design$sample)))

    # PCA
    pdf(file.path(dir_output, "PCA_genes.pdf"))
    df_genes <- produce_pca(txi_genes, use_ruv = use_ruv)
    dev.off()
    pdf(file.path(dir_output, "PCA_tx.pdf"))
    df_tx <- produce_pca(txi_tx, use_ruv = use_ruv)
    dev.off()

    # Produce counts
    counts_genes <- format_counts(txi_genes, digits = digits)
    counts_tx <- format_counts(txi_tx, digits = digits)
    counts <- rbind(counts_tx, counts_genes)
    write_csv(counts, file.path(dir_output, "counts.csv"))

    # Produce DE
    dds_genes <- deseq2_analysis(txi_genes, design, ~ group, use_ruv = use_ruv)
    dds_tx <- deseq2_analysis(txi_tx, design, ~ group, use_ruv = use_ruv)

    de_genes <- map(contrasts, ~ format_de(dds_genes, txi_genes, .x, ignoreTxVersion, digits = digits))
    de_tx <- map(contrasts, ~ format_de(dds_tx, txi_tx, .x, ignoreTxVersion, digits = digits))

    rbind_df <- function(n) {
        tx <- de_tx[[n]]
        genes <- de_genes[[n]]
        rbind(tx, genes)
    }
    de <- map(names(de_tx), rbind_df)
    names(de) <- names(de_tx)

    iwalk(de, ~ write_csv(.x, file.path(dir_output, paste0(.y, ".csv"))))

    invisible(list(txi_tx = txi_tx,
                   txi_genes = txi_genes,
                   df_tx = df_tx,
                   df_genes = df_genes,
                   design = design,
                   contrasts = contrasts,
                   counts = counts,
                   de = de,
                   dds_genes = dds_genes,
                   dds_tx = dds_tx))
}

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
#'             * Hs.Ensembl79
#'             * Hs.Ensembl91
#'             * Hs.Ensembl95
#'             * Hs.Ensembl97
#'             * Hs.Ensembl98
#'             * Mm.Ensembl91
#'             * Mm.Ensembl92
#'             * Mm.Ensembl94
#'             * Mm.Ensembl97
#'             * Rn.Ensembl76
#'             * Rn.Ensembl79
#'             * Rn.Ensembl92
#'             * Rn.Ensembl98
#'             * peaux_colonisees
#' @param design The experimental design
#'               (1st column: sample; 2nd column: group)
#' @param contrasts List of contrasts
#' @param dir_output Directory where to write outputs
#' @param file_type Abundance file format to use (h5 or tsv).
#' @param digits Integer indicating the number of decimal places
#'
#' @return Invisibly returns a list with txi_tx, txi_genes, df_tx (PCA),
#'         df_genes (PCA), design, contrasts, counts and de values.
#'
#' @examples
#' dir_kallisto <- get_demo_kallisto_dir()
#' contrasts <- list(comp = c("group", "A", "B"))
#' design <- get_demo_design()
#' produce_deliverables(dir_kallisto,
#'                      anno = "Hs.Ensembl79",
#'                      design = design,
#'                      contrast = contrasts,
#'                      dir_output = ".",
#'                      file_type = "tsv")
#'
#' @import purrr
#'
#' @export
produce_deliverables <- function (dir_kallisto, anno, design, contrasts, dir_output, file_type = "h5", digits = 4) {
    stopifnot(dir.exists(dir_kallisto))
    stopifnot(dir.exists(dir_output))
    stopifnot(file_type %in% c("tsv", "h5"))
    stopifnot(colnames(design) == c("sample", "group"))
    validate_anno(anno)

    file_type <-paste0(file_type, "$")
    files <- dir(dir_kallisto, pattern = file_type, recursive = TRUE, full.names = TRUE)
    names(files) <- basename(dirname(files))

    # Import quantifications
    txi_genes <- import_kallisto(files, anno = anno, txOut = FALSE)
    txi_tx <- import_kallisto(files, anno = anno, txOut = TRUE)
    stopifnot(identical(rownames(txi_genes$counts), design$sample))

    # PCA
    pdf(file.path(dir_output, "PCA_genes.pdf"))
    df_genes <- produce_pca(txi_genes)
    dev.off()
    pdf(file.path(dir_output, "PCA_tx.pdf"))
    df_tx <- produce_pca(txi_tx)
    dev.off()

    # Produce counts
    counts_genes <- format_counts(txi_genes, digits = digits)
    counts_tx <- format_counts(txi_tx, digits = digits)
    counts <- rbind(counts_tx, counts_genes)
    write_csv(counts, file.path(dir_output, "counts.csv"))

    # Produce DE
    dds_genes <- deseq2_analysis(txi_genes, design, ~ group)
    dds_tx <- deseq2_analysis(txi_tx, design, ~ group)

    de_genes <- map(contrasts, ~ format_de(dds_genes, txi_genes, .x, digits = digits))
    de_tx <- map(contrasts, ~ format_de(dds_tx, txi_tx, .x, digits = digits))

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

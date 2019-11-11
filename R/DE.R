#' DESeq2 analysis
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param design The experimental design (see ?DESeqDataSetFromTximport).
#' @param formula The design formula in data.frame format (see
#'                ?DESeqDataSetFromTximport).
#' @param filter The minimum number of reads detected for a feature across all
#'               samples. Default: 2
#'
#' @return A DESeqDataSet object.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- DESeq2::results(dds, contrast = c("group", "A", "B"))
#'
#' @import DESeq2
#'
#' @export
deseq2_analysis <- function(txi, design, formula, filter = 2) {
    stopifnot(ncol(design) == 2)
    stopifnot(colnames(design) != c("sample", "group"))
    dds <- DESeqDataSetFromTximport(txi, design, formula)
    dds <- dds[rowSums(counts(dds)) >= filter]
    dds <- DESeq(dds)
}

#' Prepare formated DE table.
#'
#' The table contains annotation and the DE results.
#'
#' @param dds The DESeqDataSet object returned by deseq2_analysis.
#' @param txi The txi object returned by the import_kallisto function.
#' @param contrast The contrast for the comparison (see ?DESeq2::results).
#' @param digits Integer indicating the number of decimal places
#'
#' @return A data.frame with the anno and the merged counts values.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- format_de(dds, txi, c("group", "A", "B"))
#'
#' @import DESeq2
#' @import dplyr
#' @import tibble
#'
#' @export
format_de <- function(dds, txi, contrast, digits = 4) {
    res <- results(dds, contrast = contrast) %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        left_join(txi$anno, by = "id") %>%
        mutate(mean_TPM_grp1 = get_mean_tpm(dds, txi, contrast[2]),
               mean_TPM_grp2 = get_mean_tpm(dds, txi, contrast[3]),
               fold_change = 2^log2FoldChange) %>%
        splicing_analysis(txi)

    res <- dplyr::select(res, id, ensembl_gene, symbol, entrez_id, transcript_type,
                  mean_TPM_grp1, mean_TPM_grp2, pV = pvalue, qV = padj,
                  percent_grp1, percent_grp2, main_isoform_grp1,
                  main_isoform_grp2, baseMean, lfcSE, fold_change,
                  log2FoldChange, stat)

    res <- mutate(res,
           mean_TPM_grp1 = round(mean_TPM_grp1, digits) %>% format(scientific = FALSE),
           mean_TPM_grp2 = round(mean_TPM_grp2, digits) %>% format(scientific = FALSE),
           pV = round(pV, digits) %>% format(scientific = FALSE),
           qV = round(qV, digits) %>% format(scientific = FALSE),
           percent_grp1 = round(percent_grp1, digits) %>% format(scientific = FALSE),
           percent_grp2 = round(percent_grp2, digits) %>% format(scientific = FALSE),
           baseMean = round(baseMean, digits) %>% format(scientific = FALSE),
           lfcSE = round(lfcSE, digits) %>% format(scientific = FALSE),
           fold_change = round(fold_change, digits) %>% format(scientific = FALSE),
           log2FoldChange = round(log2FoldChange, digits) %>% format(scientific = FALSE),
           stat = round(stat, digits) %>% format(scientific = FALSE))
    as.data.frame(res)
}

get_mean_tpm <- function(dds, txi, group) {
    samples <- dds@colData[,"sample", drop = TRUE]
    samples <- samples[dds@colData[,"group", drop = TRUE] == group]
    mean_tpm <- rowMeans(txi$abundance[,samples])
    mean_tpm[names(dds)]
}

splicing_analysis <- function(res, txi) {
    is.max <- function(x) seq_along(x) == which.max(x)
    if (txi$txOut) {
        sums <- group_by(res, ensembl_gene) %>%
            summarize(sum_g1 = sum(mean_TPM_grp1),
            sum_g2 = sum(mean_TPM_grp2))
        left_join(res, sums, by = "ensembl_gene") %>%
            group_by(ensembl_gene) %>%
            mutate(percent_grp1 = mean_TPM_grp1 / sum_g1,
                   percent_grp2 = mean_TPM_grp2 / sum_g2) %>%
            mutate(percent_grp1 = if_else(is.nan(percent_grp1),
                                          0, percent_grp1),
                   percent_grp2 = if_else(is.nan(percent_grp2),
                                          0, percent_grp2)) %>%
            mutate(main_isoform_grp1 = is.max(mean_TPM_grp1),
                   main_isoform_grp2 = is.max(mean_TPM_grp2)) %>%
            dplyr::select(-sum_g1, -sum_g2)
    } else {
        mutate(res,
               percent_grp1 = NA,
               percent_grp2 = NA,
               main_isoform_grp1 = NA,
               main_isoform_grp2 = NA)
    }
}

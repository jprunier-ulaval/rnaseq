#' DESeq2 analysis
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param design The experimental design (see ?DESeqDataSetFromTximport).
#' @param formula The design formula in data.frame format (see
#'                ?DESeqDataSetFromTximport).
#' @param filter The minimum number of reads detected for a feature across all
#'               samples. Default: 2
#' @param use_ruv Use RUVg normalization? Needs to be pre-computed using the
#'                \code{ruvg_normalization} function.
#' @param ... Extract param for the DESeq2::DESeq function
#'
#' @return A DESeqDataSet object.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- DESeq2::results(dds, contrast = c("group", "A", "B"))
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom DESeq2 DESeq
#'
#' @export
deseq2_analysis <- function(txi, design, formula, filter = 2, use_ruv = FALSE, ...) {
    stopifnot(ncol(design) == 2)
    stopifnot(all(colnames(design) == c("sample", "group")))
    stopifnot(identical(colnames(txi$counts), as.character(design$sample)))
    stopifnot(is(use_ruv, "logical"))

    if (use_ruv) {
        stopifnot("ruvg_counts" %in% names(txi))
        stopifnot(is(txi$ruvg_counts, "matrix"))
        dds <- DESeq2::DESeqDataSetFromMatrix(txi$ruvg_counts, design, formula)
    } else {
        dds <- DESeq2::DESeqDataSetFromTximport(txi, design, formula)
    }
    dds <- dds[rowSums(counts(dds)) >= filter]
    dds <- DESeq2::DESeq(dds, ...)
    dds
}

#' Prepare formated DE table.
#'
#' The table contains annotation and the DE results.
#'
#' @param dds The DESeqDataSet object returned by deseq2_analysis.
#' @param txi The txi object returned by the import_kallisto function.
#' @param contrast The contrast for the comparison (see ?DESeq2::results).
#' @param ignoreTxVersion Should the transcript version be ignored for anno
#' mapping. Default: \code{FALSE}.
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
#' @importFrom magrittr %>%
#' @importFrom DESeq2 results
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom stringr str_replace
#'
#' @export
format_de <- function(dds, txi, contrast, ignoreTxVersion = FALSE, digits = 4) {
    res <- DESeq2::results(dds, contrast = contrast) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("id")
    if (!ignoreTxVersion) {
        res <- dplyr::left_join(res, txi$anno, by = "id")
    } else {
        res$id <- stringr::str_replace(res$id, "\\..*$", "")
        res <- dplyr::left_join(res,
                                dplyr::mutate(txi$anno,
                                              id = str_replace(id, "\\..*$", "")),
                                by = "id")
    }
    res <- dplyr::mutate(res,
               mean_TPM_grp1 = get_mean_tpm(dds, txi, contrast[2]),
               mean_TPM_grp2 = get_mean_tpm(dds, txi, contrast[3]),
               ratio = 2^log2FoldChange,
               fold_change = if_else(ratio < 1, -1 * (1/ratio), ratio)) %>%
               splicing_analysis(txi)

    res <- dplyr::select(res, id, ensembl_gene, symbol, entrez_id, transcript_type,
                  mean_TPM_grp1, mean_TPM_grp2, pV = pvalue, qV = padj,
                  percent_grp1, percent_grp2, main_isoform_grp1,
                  main_isoform_grp2, baseMean, lfcSE, log2FoldChange,
                  fold_change, ratio, stat)

    res <- dplyr::mutate(res,
           mean_TPM_grp1 = round_values(mean_TPM_grp1, digits),
           mean_TPM_grp2 = round_values(mean_TPM_grp2, digits),
           pV = round_values(pV, digits),
           qV = round_values(qV, digits),
           percent_grp1 = round_values(percent_grp1, digits),
           percent_grp2 = round_values(percent_grp2, digits),
           baseMean = round_values(baseMean, digits),
           lfcSE = round_values(lfcSE, digits),
           fold_change = round_values(fold_change, digits),
           log2FoldChange = round_values(log2FoldChange, digits),
           stat = round_values(stat, digits))
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
        sums <- dplyr::group_by(res, ensembl_gene) %>%
            dplyr::summarize(sum_g1 = sum(mean_TPM_grp1),
            sum_g2 = sum(mean_TPM_grp2))
        dplyr::left_join(res, sums, by = "ensembl_gene") %>%
            dplyr::group_by(ensembl_gene) %>%
            dplyr::mutate(percent_grp1 = mean_TPM_grp1 / sum_g1,
                   percent_grp2 = mean_TPM_grp2 / sum_g2) %>%
            dplyr::mutate(percent_grp1 = if_else(is.nan(percent_grp1),
                                          0, percent_grp1),
                   percent_grp2 = if_else(is.nan(percent_grp2),
                                          0, percent_grp2)) %>%
            dplyr::mutate(main_isoform_grp1 = is.max(mean_TPM_grp1),
                   main_isoform_grp2 = is.max(mean_TPM_grp2)) %>%
            dplyr::select(-sum_g1, -sum_g2)
    } else {
        dplyr::mutate(res,
               percent_grp1 = NA,
               percent_grp2 = NA,
               main_isoform_grp1 = NA,
               main_isoform_grp2 = NA)
    }
}

round_values <- function(values, digits = 4) {
    stopifnot(is.numeric(values))
    i <- is.na(values)
    new_values <- round(values, digits) %>% format(scientific = FALSE)
    new_values[i] <- NA
    as.numeric(new_values)
}

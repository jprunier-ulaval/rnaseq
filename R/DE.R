#' DESeq2 analysis
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param design The experimental design (see ?DESeqDataSetFromTximport).
#' @param formula The design formula in data.frame format (see
#'                ?DESeqDataSetFromTximport).
#' @param contrast The contrast for the comparison (see ?DESeq2::results).
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
    dds <- DESeqDataSetFromTximport(txi, design, formula)
    dds <- dds[rowSums(counts(dds)) >= filter]
    dds <- DESeq(dds)
}

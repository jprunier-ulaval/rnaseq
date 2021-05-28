#' Prepare formated table counts.
#'
#' The table contains annotation and the counts (raw_counts/tpm/fpkm).
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param digits Integer indicating the number of decimal places
#' @param use_ruv Use RUVg normalization? Needs to be pre-computed using the
#'                \code{ruvg_normalization} function. Default: \code{FALSE}.
#'
#' @return A data.frame with the anno and the merged counts values.
#'
#' @examples
#' txi <- get_demo_txi()
#' counts <- format_counts(txi)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr everything
#'
#' @export
format_counts <- function(txi, digits = 4, use_ruv = FALSE) {
    names_txi <- c("abundance", "counts", "length", "countsFromAbundance",
                   "fpkm", "anno", "txOut")
    stopifnot(all(names_txi %in% names(txi)))
    stopifnot(identical(rownames(txi$counts), rownames(txi$abundance)))
    stopifnot(identical(colnames(txi$counts), colnames(txi$abundance)))
    stopifnot(identical(rownames(txi$counts), rownames(txi$fpkm)))
    stopifnot(identical(colnames(txi$counts), colnames(txi$fpkm)))
    if (use_ruv) {
        stopifnot(identical(rownames(txi$counts), rownames(txi$ruvg_counts)))
        stopifnot(identical(colnames(txi$counts), colnames(txi$ruvg_counts)))
    }

    # Extract values
    if (!use_ruv) {
        raw_counts <- round(txi$counts, digits)
    } else {
        raw_counts <- round(txi$ruvg_counts, digits)
    }
    tpm <- round(txi$abundance, digits)
    fpkm <- round(txi$fpkm, digits)

    # Merge values
    res <- paste(paste(raw_counts, tpm, sep = "/"), fpkm, sep = "/") %>%
        matrix(ncol = ncol(raw_counts)) %>%
        as.data.frame %>%
        setNames(colnames(raw_counts)) %>%
        dplyr::mutate(id = rownames(raw_counts))

    # Add anno
    dplyr::left_join(res, txi$anno, by = "id") %>%
        dplyr::select(id:transcript_type, dplyr::everything())
}

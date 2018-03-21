#' Prepare formated table counts.
#'
#' The table contains annotation and the counts (raw_counts/tpm/fpkm).
#'
#' @param txi The txi object returned by the import_kallisto function.
#'
#' @return A data.frame with the anno and the merged counts values.
#'
#' @examples
#' txi <- get_demo_txi()
#' counts <- format_counts(txi)
#'
#' @import dplyr
#'
#' @export
format_counts <- function(txi) {
    names_txi <- c("abundance", "counts", "length", "countsFromAbundance",
                   "fpkm", "anno", "txOut")
    stopifnot(all(names_txi %in% names(txi)))
    stopifnot(identical(rownames(txi$counts), rownames(txi$abundance)))
    stopifnot(identical(colnames(txi$counts), colnames(txi$abundance)))
    stopifnot(identical(rownames(txi$counts), rownames(txi$fpkm)))
    stopifnot(identical(colnames(txi$counts), colnames(txi$fpkm)))

    # Extract values
    raw_counts <- round(txi$counts, 4)
    tpm <- round(txi$abundance, 4)
    fpkm <- round(txi$fpkm, 4)

    # Merge values
    res <- paste(paste(raw_counts, tpm, sep = "/"), fpkm, sep = "/") %>%
        matrix(ncol = ncol(raw_counts)) %>%
        as.data.frame %>%
        setNames(colnames(raw_counts)) %>%
        mutate(id = rownames(raw_counts))

    # Add anno
    left_join(res, txi$anno, by = "id") %>%
        dplyr::select(id:transcript_type, everything())
}

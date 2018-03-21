#' Prepare formated table counts.
#'
#' The table contains annotation and the counts (raw_counts/tpm/fpkm).
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param anno The anno used when creating the txi object.
#' @param level The level of the analysis (gene or transcript)
#'
#' @return A data.frame with the anno and the merged counts values.
#'
#' @examples
#' txi <- get_demo_txi()
#' counts <- format_counts(txi, "Hs.Ensembl79", "gene")
#'
#' @import dplyr
#'
#' @export
format_counts <- function(txi, anno, level) {
    stopifnot(level %in% c("gene", "transcript"))

    # Extract values
    raw_counts <- round(txi$counts, 4)
    tpm <- round(txi$abundance, 4)
    fpkm <- round(get_fpkm(txi), 4)

    stopifnot(identical(rownames(raw_counts), rownames(tpm)))
    stopifnot(identical(colnames(raw_counts), colnames(tpm)))

    # Merge values
    res <- paste(paste(raw_counts, tpm, sep = "/"), fpkm, sep = "/") %>%
        matrix(ncol = ncol(raw_counts)) %>%
        as.data.frame %>%
        setNames(colnames(raw_counts)) %>%
        mutate(id = rownames(raw_counts))

    # Add anno
    anno <- get_anno(anno, level = level)
    left_join(res, anno, by = "id") %>%
        dplyr::select(id:transcript_type, everything())
}

get_fpkm <- function(txi) {
    (txi$counts * 10^6) / (colSums(txi$counts) * (txi$length))
}

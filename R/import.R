#' Import quantifications from Kallisto
#'
#' @param filenames Paths to the abundance files.
#' @param anno The version of the annotation to use. Default: "Hs.Ensembl91"
#'             Currently available:
#'             * Hs.Ensembl91
#'             * Hs.Ensembl79
#' @param txOut Return counts and abundance at the transcript level. Default:
#'              FALSE
#' @param ignoreTxVersion Ignore version of tx. Default = FALSE
#'
#' @return A txi object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' txi <- import_kallisto(abundances, anno = "Hs.Ensembl79")
#'
#' @import readr
#' @import dplyr
#' @import stringr
#' @import tximport
#'
#' @export
import_kallisto <- function(filenames, anno = "Hs.Ensembl91", txOut = FALSE,
                            ignoreTxVersion = FALSE) {
    stopifnot(all(file.exists(filenames)))
    if (stringr::str_detect(anno, "Ensembl")) {
        anno <- get_anno(anno)
        tx2gene <- dplyr::select(anno, TXNAME = id, GENEID = ensembl_gene)
    }
    if (txOut == TRUE) {
        tximport(filenames, type = "kallisto", tx2gene = tx2gene, txOut = TRUE,
                 ignoreTxVersion = ignoreTxVersion)
    } else {
        tximport(filenames, type = "kallisto", tx2gene = tx2gene,
                 ignoreTxVersion = ignoreTxVersion)
    }
}

get_anno <- function(anno) {
    valid_anno <- c("Hs.Ensembl91", "Hs.Ensembl79")
    stopifnot(anno %in% valid_anno)
    get(anno)
}

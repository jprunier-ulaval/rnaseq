#' Import quantifications from Kallisto
#'
#' @param filenames Paths to the abundance files.
#' @param anno The version of the annotation to use. Default: "Hs.Ensembl91"
#'             Currently available:
#'             * Hs.Gencode27
#'             * Hs.Ensembl91
#'             * Hs.Ensembl79
#'             * Mm.Ensembl91
#'             * Mm.Ensembl92
#'             * Rn.Ensembl92
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
    stopifnot(txOut %in% c(TRUE, FALSE))
    stopifnot(ignoreTxVersion %in% c(TRUE, FALSE))

    tx2gene <- get(anno) %>%
        dplyr::select(TXNAME = id, GENEID = ensembl_gene)
    if (txOut == TRUE) {
        txi <- tximport(filenames, type = "kallisto", tx2gene = tx2gene, txOut = TRUE,
                 ignoreTxVersion = ignoreTxVersion)
    } else {
        txi <- tximport(filenames, type = "kallisto", tx2gene = tx2gene,
                 ignoreTxVersion = ignoreTxVersion)
    }
    txi$fpkm <- get_fpkm(txi)
    txi$anno <- get_anno(anno, txOut)
    txi$txOut <- txOut
    stopifnot(all(rownames(txi$fpkm) %in% txi$anno$id))
    txi
}

get_anno <- function(anno, txOut) {
    validate_anno(anno)
    anno <- get(anno)
    if (!txOut) {
        anno <- mutate(anno, id = ensembl_gene) %>%
            filter(!duplicated(ensembl_gene))
    }
    anno
}

validate_anno <- function(anno) {
    valid_anno <- c("Hs.Gencode27", "Hs.Ensembl91", "Hs.Ensembl79",
		    "Mm.Ensembl91", "Mm.Ensembl92", "Rn.Ensembl92")
    stopifnot(anno %in% valid_anno)
}

get_fpkm <- function(txi) {
    (txi$counts * 10^6) / (colSums(txi$counts) * (txi$length))
}

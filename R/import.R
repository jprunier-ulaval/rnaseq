#' Import quantifications from Kallisto
#'
#' @param filenames Paths to the abundance files.
#' @param anno The version of the annotation to use. Default: "Hs.Ensembl91"
#'             Currently available:
#'             * Hs.Ensembl91
#'             * Hs.Ensembl79
#' @param txOut Return counts and abundance at the transcript level. Default:
#'              FALSE
#'
#' @return A txi object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' txi <- import_kallisto(abundances, anno = "Hs.Ensembl91")
#'
#' @import readr
#' @import dplyr
#' @import stringr
#' @import tximport
#'
#' @export
import_kallisto <- function(filenames, anno = "Hs.Ensembl91", txOut = FALSE) {
    stopifnot(all(file.exists(filenames)))
    valid_anno <- c("Hs.Ensembl91", "Hs.Ensembl79")
    stopifnot(anno %in% valid_anno)
    if (stringr::str_detect(anno, "Ensembl")) {
        if (anno == "Hs.Ensembl91") {
            anno <- system.file("extdata/Homo_sapiens.GRCh38.91.all.anno.csv.gz", package="rnaseq")
        }
        if (anno == "Hs.Ensembl79") {
            anno <- system.file("extdata/Homo_sapiens.GRCh38.79.all.anno.csv.gz", package="rnaseq")
        }
        anno <- readr::read_csv(anno, col_types = "ciiicccnicccccccccccccccccc")
        tx2gene <- dplyr::select(anno, TXNAME = transcript_id, GENEID = gene_id)
    }
    if (txOut == TRUE) {
        tximport(filenames, type = "kallisto", tx2gene = tx2gene, txOut = TRUE)
    } else {
        tximport(filenames, type = "kallisto", tx2gene = tx2gene)
    }
}

#' Import quantifications from Kallisto
#'
#' @param filenames Paths to the abundance files.
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
#'             * Hs.Ensembl100
#'             * Hs.Ensembl101
#'             * Mm.Ensembl91
#'             * Mm.Ensembl92
#'             * Mm.Ensembl94
#'             * Mm.Ensembl97
#'             * Mm.Ensembl99
#'             * Mm.Ensembl100
#'             * Rn.Ensembl76
#'             * Rn.Ensembl79
#'             * Rn.Ensembl92
#'             * Rn.Ensembl98
#'             * Bt.Ensembl99
#'             * peaux_colonisees
#' @param txOut Return counts and abundance at the transcript level. Default:
#'              FALSE
#' @param ignoreTxVersion Ignore version of tx. Default = FALSE
#' @param ercc92 Include ERCC92 annotation when importing. Default = FALSE
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
                            ignoreTxVersion = FALSE, ercc92 = FALSE) {
    stopifnot(all(file.exists(filenames)))
    stopifnot(txOut %in% c(TRUE, FALSE))
    stopifnot(ignoreTxVersion %in% c(TRUE, FALSE))

    tx2gene <- get(anno) %>%
        dplyr::select(TXNAME = id, GENEID = ensembl_gene)
    if (ercc92 == TRUE) {
        tx2gene_ercc92 <- dplyr::select(ERCC92, TXNAME = id, GENEID = ensembl_gene) 
        tx2gene <- rbind(tx2gene, tx2gene_ercc92)
    }
    if (txOut == TRUE) {
        txi <- tximport(filenames, type = "kallisto", tx2gene = tx2gene, txOut = TRUE,
                 ignoreTxVersion = ignoreTxVersion)
    } else {
        txi <- tximport(filenames, type = "kallisto", tx2gene = tx2gene,
                 ignoreTxVersion = ignoreTxVersion)
    }
    txi$fpkm <- get_fpkm(txi)
    txi$anno <- get_anno(anno, txOut)
    if (ercc92 == TRUE) {
        txi$anno <- rbind(txi$anno, ERCC92)
    }
    txi$txOut <- txOut
    if (!ignoreTxVersion) {
        stopifnot(all(rownames(txi$fpkm) %in% txi$anno$id))
    } else {
        row_names_txi_fpkm <- str_replace(rownames(txi$fpkm), "\\..*$", "")
        txi_anno_id <- str_replace(txi$anno$id, "\\..*$", "")
        stopifnot(all(row_names_txi_fpkm %in% txi_anno_id))
    }
    txi
}

summarize_to_gene <- function(txi_tx, anno, ignoreTxVersion = FALSE) {
    stopifnot(ignoreTxVersion %in% c(TRUE, FALSE))

    tx2gene <- get_anno(anno) %>%
        dplyr::select(TXNAME = id, GENEID = ensembl_gene)

    txi <- summarizeToGene(txi_tx, tx2gene = tx2gene, ignoreTxVersion = ignoreTxVersion)
    txi$fpkm <- get_fpkm(txi)
    txi$anno <- get_anno(anno, txOut = FALSE)
    txi$txOut <- FALSE
    stopifnot(all(rownames(txi$fpkm) %in% txi$anno$id))
    txi
}

get_anno <- function(anno, txOut = TRUE) {
    validate_anno(anno)
    anno <- get(anno)
    if (!txOut) {
        anno <- mutate(anno, id = ensembl_gene) %>%
            filter(!duplicated(ensembl_gene))
    }
    anno
}

validate_anno <- function(anno) {
    valid_anno <- c("Hs.Gencode19", "Hs.Gencode27", "Hs.Ensembl79",
                    "Hs.Ensembl91", "Hs.Ensembl95", "Hs.Ensembl97",
                    "Hs.Ensembl98", "Hs.Ensembl100", "Hs.Ensembl101",
                    "Mm.Ensembl91", "Mm.Ensembl92", "Mm.Ensembl94",
                    "Mm.Ensembl99", "Mm.Ensembl100", "Rn.Ensembl76",
                    "Rn.Ensembl79", "Rn.Ensembl92", "Rn.Ensembl98",
                    "Bt.Ensembl99", "peaux_colonisees")
    stopifnot(anno %in% valid_anno)
}

get_fpkm <- function(txi) {
    (txi$counts * 10^6) / (colSums(txi$counts) * (txi$length/1000))
}

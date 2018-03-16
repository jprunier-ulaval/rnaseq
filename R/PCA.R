#' Produce a PCA plot
#'
#' @param txi The txi object returned by the import_kallisto function.
#'
#' @return Produce the PCA and silently returns the data.frame used.
#'
#' @examples
#' txi <- get_demo_txi()
#' \dontrun{
#' df <- produce_pca(txi)
#' }
#'
#' @import FactoMineR
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ggrepel
#'
#' @export
produce_pca <- function(txi) {
    tpm <- as.data.frame(txi$abundance) %>%
        mutate(ensembl_gene = rownames(txi$abundance)) %>%
        tidyr::gather(sample, tpm, -ensembl_gene) %>%
        tbl_df

    min_tpm <- group_by(tpm, ensembl_gene) %>%
        summarize(tpm = sum(tpm)) %>%
        filter(tpm >= 5) %>%
        pull(ensembl_gene)

    tpm_filter <- filter(tpm, ensembl_gene %in% min_tpm)

    df <- spread(tpm_filter, sample, tpm) %>%
        as.data.frame
    rownames(df) <- df$ensembl_gene
    m <- dplyr::select(df, -ensembl_gene) %>%
        as.matrix %>%
        t

    pca <- PCA(m, graph = FALSE)
    coord <- pca$ind$coord
    df <- data.frame(Dim1 = coord[,1], Dim2 = coord[,2], Dim3 = coord[,3])
    df <- df %>%
        mutate(sample = rownames(df)) %>%
        tbl_df

    p <- ggplot(df, aes(x = Dim1, y = Dim2)) +
        geom_point(size = 3) +
        geom_text_repel(aes(label = sample), color = "black", force = 10) +
        theme_bw() +
        xlab(paste0("Dim1 (", pca$eig[1,2] %>% round(2), "%)")) +
        ylab(paste0("Dim2 (", pca$eig[2,2] %>% round(2), "%)"))
    print(p)
    invisible(df)
}


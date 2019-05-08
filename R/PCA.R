#' Produce a PCA plot from txi
#'
#' @param txi The \code{txi} object returned by the \code{import_kallisto}
#' function.
#' @param txi Produce the graph. \code{TRUE} or \code{FALSE}. Default:
#' \code{TRUE}.
#'
#' @return Produce the PCA and silently returns the \code{data.frame} and the
#' result from the \code{PCA} function as a \code{list} with 2 elements
#' (\code{df} and \code{pca}).
#'
#' @examples
#' txi <- get_demo_txi()
#' df <- produce_pca(txi, graph = FALSE)
#'
#' @import FactoMineR
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ggrepel
#'
#' @export
produce_pca <- function(txi, graph = TRUE) {
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
    p <- plot_pca(df)
    print(p)

    invisible(list(df = df, pca = pca))
}

#' Produce a PCA plot from produce_pca results
#'
#' @param txi The txi object returned by the import_kallisto function.
#'
#' @return Returns the \code{ggplot} object
#'
#' @examples
#' txi <- get_demo_txi()
#' res_pca <- produce_pca(txi, graph = FALSE)
#' p <- plot_pca(res_pca$df, res_pca$pca)
#'
#' @import tidyr
#' @import ggplot2
#' @import ggrepel
#'
#' @export
plot_pca <- function(df, pca) {
    ggplot(df, aes(x = Dim1, y = Dim2)) +
        geom_point(size = 3) +
        geom_text_repel(aes(label = sample), color = "black", force = 10) +
        theme_bw() +
        xlab(paste0("Dim1 (", pca$eig[1,2] %>% round(2), "%)")) +
        ylab(paste0("Dim2 (", pca$eig[2,2] %>% round(2), "%)"))
}

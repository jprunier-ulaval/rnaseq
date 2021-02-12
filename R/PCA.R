#' Produce a PCA plot from txi
#'
#' @param txi The \code{txi} object returned by the \code{import_kallisto}
#' function.
#' @param graph Produce the graph. \code{TRUE} or \code{FALSE}. Default:
#' \code{TRUE}.
#' @param use_ruv Use RUVg normalization? Needs to be pre-computed using the
#'                \code{ruvg_normalization} function. Default: \code{FALSE}.
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
produce_pca <- function(txi, graph = TRUE, use_ruv = FALSE) {
    if (!use_ruv) {
        tpm <- as.data.frame(txi$abundance)
    } else {
        stopifnot("ruvg_counts" %in% names(txi))
        stopifnot(is(txi$ruvg_counts, "matrix"))
        tpm <- as.data.frame(txi$ruvg_counts)
    }
    tpm <- tpm %>%
            mutate(ensembl_gene = rownames(txi$abundance)) %>%
            tidyr::gather(sample, tpm, -ensembl_gene) %>%
            as_tibble

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
    df <- data.frame(Dim1 = coord[,1],
                     Dim2 = coord[,2],
                     Dim3 = coord[,3],
                     Dim4 = coord[,4],
                     Dim5 = coord[,5],
                     Dim6 = coord[,6])
    df <- df %>%
        mutate(sample = rownames(df)) %>%
        as_tibble
    res <- list(df = df, pca = pca)

    p <- plot_pca(res)
    if (graph) print(p)

    invisible(res)
}

#' Produce a PCA plot from produce_pca results
#'
#' @param res_pca The \code{list} returned from the produce_pca function
#' @param color The name of the column in \code{res_pca$df} to use to color the
#' points in the PCA. If \code{NULL}, won't add color. Default: \code{NULL}.
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
plot_pca <- function(res_pca, color = NULL) {
    if (is.null(color)) {
        p <- ggplot(res_pca$df, aes(x = Dim1, y = Dim2))
    } else {
        p <- ggplot(res_pca$df, aes_string(x = "Dim1", y = "Dim2", color = color))
    }
    p + geom_point(size = 3) +
        geom_text_repel(aes(label = sample), color = "black", force = 10) +
        theme_bw() +
        xlab(paste0("Dim1 (", res_pca$pca$eig[1,2] %>% round(2), "%)")) +
        ylab(paste0("Dim2 (", res_pca$pca$eig[2,2] %>% round(2), "%)"))
}

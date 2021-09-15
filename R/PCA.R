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
#' @importFrom magrittr %>%
#' @importFrom FactoMineR PCA
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom stringr str_extract
#' @importFrom utils tail
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
            dplyr::mutate(ensembl_gene = rownames(tpm)) %>%
            tidyr::gather(sample, tpm, -ensembl_gene)

    min_tpm <- dplyr::group_by(tpm, ensembl_gene) %>%
        dplyr::summarize(tpm = sum(tpm)) %>%
        dplyr::filter(tpm >= 5) %>%
        dplyr::pull(ensembl_gene)

    tpm_filter <- dplyr::filter(tpm, ensembl_gene %in% min_tpm)

    df <- tidyr::spread(tpm_filter, sample, tpm) %>%
        as.data.frame
    rownames(df) <- df$ensembl_gene
    m <- dplyr::select(df, -ensembl_gene) %>%
        as.matrix %>%
        t

    pca <- FactoMineR::PCA(m, graph = FALSE)
    coord <- pca$ind$coord
    dims <- stringr::str_extract(colnames(coord), "[0-9]*$") %>%
        as.numeric
    max_dim <- min(tail(dims, 1), 5)
    df <- data.frame(Dim1 = coord[,1])
    for (i in seq(2, max_dim)) {
        df[[paste0("Dim", i)]] <- coord[,i]
    }
    df <- df %>%
        dplyr::mutate(sample = rownames(df))
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
#' p <- plot_pca(res_pca)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggrepel geom_text_repel
#'
#' @export
plot_pca <- function(res_pca, color = NULL) {
    if (is.null(color)) {
        p <- ggplot2::ggplot(res_pca$df, ggplot2::aes(x = Dim1, y = Dim2))
    } else {
        p <- ggplot2::ggplot(res_pca$df,
                             ggplot2::aes_string(x = "Dim1",
                                                 y = "Dim2",
                                                 color = color))
    }
    p + ggplot2::geom_point(size = 3) +
        ggrepel::geom_text_repel(ggplot2::aes(label = sample),
                                 color = "black", force = 10) +
        ggplot2::theme_bw() +
        ggplot2::xlab(paste0("Dim1 (",
                             res_pca$pca$eig[1,2] %>% round(2), "%)")) +
        ggplot2::ylab(paste0("Dim2 (",
                             res_pca$pca$eig[2,2] %>% round(2), "%)"))
}

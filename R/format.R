#' Extract raw counts with annotation data.frame from txi
#'
#' @param txi: The txi object returned from the `import_kallisto` function
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#' raw_counts <- get_raw_count_anno_df(txi)
#'
#' @import tibble
#' @import dplyr
#'
#' @export
get_raw_count_anno_df <- function(txi) {
    as.data.frame(txi$counts) %>%
        rownames_to_column("id") %>%
        left_join(txi$anno, by = "id") %>%
        dplyr::select(one_of(colnames(txi$anno)), everything())
}

#' Extract TPM with annotation data.frame from txi
#'
#' @param txi: The txi object returned from the `import_kallisto` function
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#' tpm <- get_tpm_anno_df(txi)
#'
#' @import tibble
#' @import dplyr
#'
#' @export
get_tpm_anno_df <- function(txi) {
    as.data.frame(txi$abundance) %>%
        rownames_to_column("id") %>%
        left_join(txi$anno, by = "id") %>%
        dplyr::select(one_of(colnames(txi$anno)), everything())
}

#' Extract RUVg counts with annotation data.frame from txi
#'
#' @param txi: The txi object returned from the `import_kallisto` function and
#'             the `ruvg_normalisation` function
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#' ruvg <- get_ruvg_anno_df(txi)
#'
#' @import tibble
#' @import dplyr
#'
#' @export
get_ruvg_anno_df <- function(txi) {
    as.data.frame(txi$ruvg_counts) %>%
        rownames_to_column("id") %>%
        left_join(txi$anno, by = "id") %>%
        dplyr::select(one_of(colnames(txi$anno)), everything())
}

#' Extract counts with annotation data.frame from txi
#'
#' @param txi: The txi object returned from the `import_kallisto` function
#' @param colname: The txi object returned from the `import_kallisto` function
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#' raw_counts <- get_anno_df(txi, "raw_count")
#' tpm <- get_anno_df(txi, "tpm")
#' fpkm <- get_anno_df(txi, "fpkm")
#' ruvg <- get_anno_df(txi, "ruvg")
#'
#' @import tibble
#' @import dplyr
#'
#' @export
get_anno_df <- function(txi, col_name) {
    as.data.frame(txi$fpkm) %>%
        rownames_to_column("id") %>%
        left_join(txi$anno, by = "id") %>%
        dplyr::select(one_of(colnames(txi$anno)), everything())
}

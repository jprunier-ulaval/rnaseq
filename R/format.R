#' Extract raw counts with annotation data.frame from txi 
#'
#' @param txi: The txi object returned from the `import_kallisto` function
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' txi <- import_kallisto(abundances, anno = "Hs.Ensembl79")
#' raw_counts <- get_raw_count_anno_df(txi)
#'
#' @import tibble
#' @import dplyr
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
#' txi <- import_kallisto(abundances, anno = "Hs.Ensembl79")
#' tpm <- get_tpm_anno_df(txi)
#'
#' @import tibble
#' @import dplyr
get_tpm_anno_df <- function(txi) {
    as.data.frame(txi$abundance) %>%
        rownames_to_column("id") %>%
        left_join(txi$anno, by = "id") %>%
        dplyr::select(one_of(colnames(txi$anno)), everything())
}

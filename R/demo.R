#' Get demo kallisto abundance files
#' 
#' @return A vector of kallisto abundance filenames
#' 
#' @examples
#' abundances <- get_demo_abundance_files()
#'
#' @export
get_demo_abundance_files <- function() {
    c(system.file("extdata/a/abundance.tsv", package="rnaseq"),
      system.file("extdata/b/abundance.tsv", package="rnaseq"),
      system.file("extdata/c/abundance.tsv", package="rnaseq"),
      system.file("extdata/d/abundance.tsv", package="rnaseq"))
}

#' Get demo txi file
#' 
#' @return A txi object
#' 
#' @examples
#' txi <- get_demo_txi()
#'
#' @export
get_demo_txi <- function() {
    abundances <- get_demo_abundance_files()
    import_kallisto(abundances, anno = "Hs.Ensembl79")
}

#' Get demo design
#' 
#' @return A data.frame corresponding to the demo txi object. See
#'         ?get_demo_txi.
#' 
#' @examples
#' design <- get_demo_design()
#'
#' @export
get_demo_design <- function() {
    data.frame(group = c("A", "A", "B", "B"), sample = letters[1:4])
}

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
      system.file("extdata/b/abundance.tsv", package="rnaseq"))
}

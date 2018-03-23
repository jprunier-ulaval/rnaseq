library(tidyverse)
library(rtracklayer)
source("scripts/add_refseq.R")

argv <- commandArgs(trailingOnly = TRUE)
input <- argv[1]
release <- argv[2]
var <- paste0("Hs.Gencode", release)
output <- paste0(var, ".RData")

stopifnot(file.exists(input))

fa <- readDNAStringSet(input)
x <- names(fa) %>% str_extract("^[^|]*")
fa <- fa[!str_detect(x, "PAR_Y")]

col_names <- c("id",
               "ensembl_gene",
               "havana_gene",
               "havana_transcript",
               "transcript_name",
               "symbol",
               "length",
               "transcript_type",
               "filler")

df <- tibble(full_name = names(fa)) %>%
    separate(full_name, into = col_names, sep = "\\|")

df <- add_refseq(df, "ensembl_gene") %>%
    dplyr::select(id, ensembl_gene, symbol, entrez_id = entrezgene,
                  transcript_type)

# Save results
assign(var, df)
save_cmd <- paste0("save(", var, ', file = "', output, '", compress = "xz")')
eval(parse(text = save_cmd))

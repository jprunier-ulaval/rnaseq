require(biomaRt)
add_refseq <- function(x, y) {
    ensembl <- useMart("ensembl")
    ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
    attr <- c("ensembl_gene_id", "entrezgene")
    all.entrezgene <- getBM(attributes = attr, values = "*", mart = ensembl) %>%
        unique %>%
        filter(!duplicated(ensembl_gene_id))
    left_join(x, all.entrezgene, by = setNames("ensembl_gene_id", y))
}

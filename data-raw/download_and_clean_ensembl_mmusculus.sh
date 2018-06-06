#!/bin/bash

# Note: does not work with Ensembl release < 82

RELEASE=$1
URL_FTP=ftp://ftp.ensembl.org/pub/release-${RELEASE}
FILE_GTF=Mus_musculus.GRCm38.${RELEASE}.chr.gtf.gz

wget ${URL_FTP}/gtf/mus_musculus/${FILE_GTF}
Rscript scripts/prepare_anno_ensembl.R Mm ${FILE_GTF} ${RELEASE}

rm -rf ${FILE_GTF}

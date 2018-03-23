#!/bin/bash

RELEASE=$1
URL_FTP=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${RELEASE}
FILE_FA=gencode.v${RELEASE}.transcripts.fa.gz

wget ${URL_FTP}/${FILE_FA}

echo Rscript scripts/prepare_anno_gencode.R ${FILE_GTF} ${RELEASE}

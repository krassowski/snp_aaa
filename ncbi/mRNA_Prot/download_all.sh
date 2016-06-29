#!/usr/bin/env bash
wget -r --no-parent -A 'human.*.rna.gbff.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/
mv ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/* .
rm -r ftp.ncbi.nlm.nih.gov


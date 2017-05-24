#!/usr/bin/env bash
mkdir -p ncbi
cd ncbi
mkdir -p dbsnp_150-grch37p13
cd dbsnp_150-grch37p13
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz.tbi

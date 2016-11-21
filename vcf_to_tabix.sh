#!/bin/bash
# bgzip and tabix comes from htslib library
filename=$1
gunzip -c $filename | bgzip -c > $filename.bgz
tabix -p vcf $filename.bgz

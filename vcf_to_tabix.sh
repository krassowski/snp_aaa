#!/bin/bash
# bgzip and tabix comes from htslib library
filename=$1
bgzip -c $filename > $filename.gz
tabix -p vcf $filename.gz

#!/bin/bash
# bgzip and tabix comes from htslib library
filename=$1

(gunzip -c $filename | grep ^#; gunzip -c $filename | grep -v ^# | sort -k1,1d -k2,2n) | bgzip -c > $filename.bgz
tabix -p vcf $filename.bgz

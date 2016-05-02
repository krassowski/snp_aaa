#!/bin/bash
filename=$1
bgzip -c $filename > $filename.gz
tabix -p vcf $filename.gz

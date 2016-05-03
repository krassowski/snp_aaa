#!/bin/bash
(
    (
        head -n 1 $1 |  # take care of header
        sed 's/Chromosome:G_Start\.\.G_Stop/Chromosome\tG_Start\tG_Stop/'   # split column names for chromosomal position
    ) &&
    (
        tail -n +2 $1 | # take care for all lines after header
        sed -E 's/(.):([0-9]+)\.\./\1\t\2\t/' | # replace Y:12..15 to Y  12   15
        sort -t $'\t' -k20,20 -k21,21n -k22,22n # sort the stuff: use tab as sep (crucial!), 20 - chrom name, 21 - start, 22 - end
    )
) | bgzip > $1.tabix.gz;    # block gzip with a tool from htslib

tabix -s 20 -b 21 -e 22 -f -S 1 $1.tabix.gz;    # create index with tabix from htslib: -S 1 = skip first line (header)
                                                # -s 20 = chrom, -b 21 = start, -e 22 = end, -f = update old file if exists

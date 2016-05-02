#!/bin/bash
((head -n 1 $1 | sed 's/Chromosome:G_Start\.\.G_Stop/Chromosome\tG_Start\tG_Stop/') && (tail -n +2 $1 | sed -E 's/(.):([0-9]+)\.\./\1\t\2\t/' | sort -k20,20 -k21,21n -k22,22n)) | bgzip > $1.tabix.gz;
tabix -s 20 -b 21 -e 22 -f -S 1 $1.tabix.gz;
# -S 1 = skip first line (header)
# -s 20 = chromosome column
# -f = overwrite existing file if any

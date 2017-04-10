#Note: hgsid may require an update each time when downloading data. just open the website and gra it from cookies/request
wget http://genome.ucsc.edu/cgi-bin/hgTables --post-data "hgsid=587035373_LgpTD7UIYEEaePMt5dz5SQF74Oie&jsh_pageVertPos=0&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=refGene&hgta_table=refGene&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=primaryTable&boolshad.sendToGalaxy=0&boolshad.sendToGreat=0&boolshad.sendToGenomeSpace=0&hgta_outFileName=hg19_refGene_table&hgta_compressType=none&hgta_doTopSubmit=get+output"


wget http://genome.ucsc.edu/cgi-bin/hgTables --post-data "hgsid=587035373_LgpTD7UIYEEaePMt5dz5SQF74Oie&jsh_pageVertPos=0&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=refGene&hgta_table=refGene&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=sequence&boolshad.sendToGalaxy=0&boolshad.sendToGreat=0&boolshad.sendToGenomeSpace=0&hgta_outFileName=hg19_refGene_table&hgta_compressType=none&hgta_doTopSubmit=get+output" -O temp
wget "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=587035373_LgpTD7UIYEEaePMt5dz5SQF74Oie&hgta_geneSeqType=mRNA&hgta_doGenePredSequence=submit" -O sequence.fasta
rm temp

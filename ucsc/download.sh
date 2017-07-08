#!/usr/bin/env bash

wget https://gist.githubusercontent.com/krassowski/8c9710fa20ac944ec8d47ac4a0ac5b4a/raw/444fcc584bc10b5e504c05a6063a281cee808c9c/ucsc_download.sh
source ucsc_download.sh


get_from_ucsc /dev/null <<-QUERY
	hgta_database:hg19
	hgta_table:knownToRefSeq
	hgta_fs.linked.hg19.knownToEnsembl:on
	hgta_outFileName:output
	hgta_compressType:gzip
	hgta_doSelectFieldsMore:allow selection from checked tables
QUERY


get_from_ucsc refseq_to_ensembl.tsv.gz <<-QUERY
	hgta_database:hg19
	hgta_table:knownToRefSeq
	hgta_fs.check.hg19.knownToRefSeq.value:on
	hgta_fs.check.hg19.knownToEnsembl.value:on
	hgta_doPrintSelectedFields:get output
QUERY

unset sid

get_from_ucsc ref_gene.tsv.gz <<-QUERY
	clade:mammal
	org:Human
	db:hg19
	hgta_group:genes
	hgta_track:refGene
	hgta_table:refGene
	hgta_regionType:genome
	hgta_outputType:primaryTable
	hgta_outFileName:output
	hgta_compressType:gzip
	hgta_doTopSubmit:get output
QUERY

unset sid

get_from_ucsc /dev/null <<-QUERY
	clade:mammal
	org:Human
	db:hg19
	hgta_group:genes
	hgta_track:refGene
	hgta_table:refGene
	hgta_regionType:genome
	hgta_outputType:sequence
	hgta_outFileName:hg19_refGene_table
	hgta_compressType:gzip
	hgta_doTopSubmit:get output
QUERY


get_from_ucsc sequences.fasta.gz <<-QUERY
	hgta_geneSeqType:mRNA
	hgta_doGenePredSequence:submit
QUERY


gunzip sequences.fasta.gz
rm ucsc_download.sh
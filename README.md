# Pipeline for SNP analysis with regard to influence of poly(A) tracks on protein creation in Human #

The repository contains multiple Python modules, useful when analyzing data from Ensembl, Cosmic, NCBI and other services while searching for information about effect of particular SNPs. Most of the modules are written in Python 2.7.4, some uses Python 3.4. The major part of the analysis revolves around data from [PATACSDB](https://peerj.com/articles/cs-45/) and expands on the research published in Science Advances: ["Translational control by lysine-encoding A-rich sequences"](http://advances.sciencemag.org/content/1/6/e1500154). The core module for analysis of somatic SNPs (coming mostly from COSMIC) is called `snp_parser`. Other, started but unfinished module `gte_parser` was intended to perform analysis on all mutations in protein-coding DNA sequences with respect to expression level taken from [GTEx portal](http://www.gtexportal.org/home/).


## What is the workflow of particular modules? ##

### Somatic SNPs analysis  ###

The primary goal of this analysis was to create an association table showing relation between:

* expression of genes which might be affected by particular SNPs, and
* introduction, removal or change in length of poly(A) tracks in DNA sequence, caused by these SNPs.

The polyadenylate track [poly(A)] was defined (the same as in PATACSDB) as 12A-1 (sequence of twelve subsequent adenines with one mismatch allowed). The code responsible for poly(A) detection is located in `poly_a` module.

#### Data retrieval

* Variants' names were obtained from Ensembl with use of Biomart API, with restrict to hsapiens_snp_som dataset (somatic mutations only),
* the variant's alleles data were taken from VCF files from three sources: [Ensembl](ftp://ftp.ensembl.org/pub/release-87/variation/vcf/homo_sapiens/), [COSMIC](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/VCF/CosmicCodingMuts.vcf.gz) and [NCBI](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh38p2/VCF/00-All.vcf.gz)  [tabix](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh38p2/VCF/00-All.vcf.gz.tbi),
* gene sequences (also from Ensembl) was downloaded as [chromosome-grouped DNA](ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/), [CDS](ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cds/) and [cDNA](ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cdna/) FASTA files.
* from COSMIC following data were downloaded:
    * [gene expression level 3 data](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/CosmicCompleteGeneExpression.tsv.gz) - from the TCGA portal
    * [all copy number abberations](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/CosmicCompleteCNA.tsv.gz) - used to map gene expression from samples to gene transcripts
    * [all COSMIC genes cDNA sequences](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/All_COSMIC_Genes.fasta.gz) - names mappings
* From [PATACSDB](http://sysbio.ibb.waw.pl/patacsdb) list of gene transcript identifiers having poly(A) tracks was retrieved (in CSV format)

The analysed variants set has been limited to: synonymous, stop gained, coding sequence and missense variants.  Only SNPs from genes present in PATACSDB (Homo sapiens dataset) were downloaded (so we had guarantee of at least one poly(A) track presence). I was not able to get the alleles from Ensembl's biomart as it is not providing the alleles of variants from COSMIC database at the time of writing (and those were crucial for the analysis!). From COSMIC genes cDNA sequences only headers were used in order to map the names of transcripts between Ensembl and COSMIC databases (worth noting, the distribution of data between transcripts from COSMIC were later in the analysis treated as uninformative, due to purported lack of proper curation i.e. there are serious premises multiple transcript-specific entries were nonetheless assigned to canonical transcript when incorporated to the COSMIC database).

Versions: all the data come from GRCh38 and COSMIC v79, Ensembl 87.

#### General workflow

Variants retrieved from biomart and populated with additional data from VCF files are sorted by genes. Later in surroundings of each variant the search for poly(A) is performed - both in original and mutated sequence.

#### Future development ideas

It would be beneficial to split the final association table to separate tables for different variant consequences (synonymous, stop gained, coding sequence and missense) and even to include less restrictive filters on consequences (i.e. take all variants) but then have multiple tables showing associations with respect to consequences of variants. Some types of variant's consequences might be used as control - for example we are expecting to find no correlation in UTR related variants.

### mRNA expression levels analysis - GTEx ###

#### The workflow idea

1. Restrict analyse data to all genes where poly(A) is guaranteed to occur (PATACSDB, human)
2. Fetch all SNP variants (hsapiens_snp) from Ensembl's biomart, restrict to variants from dbSNP (as other variants are not present in GTEx)
3. With use of tuples: (snp, ensemble_gene_id) mine the GTEx database saving average level of expression associated with particular tuple GTEx_Analysis_v6_eQTLs.tar.gz, spiting the analysis between distinct tissues (as reported by GTEx)
4. Prepare and analyse accumulated data

### What already has been tried

Using gbff files from NCBI in order to gather variants does not look as the best idea - there are a lot of inconsistencies and edge cases to think about.
Parallel processing of GTEx data might be useful when analyzing the GTEx data (divided by tissues), some related code is present in the `gte_parser` file.
The file with all required data is called `GTEx_Analysis_V6_eQTLs.tar.gz` and is available to download from GTEx portal after [registration](http://www.gtexportal.org/home/register).


##  Summary of set up

### Dependencies

To use biomart API, the parser uses [fork of pip-distributed biomart](https://github.com/krassowski/biomart) package (which is 2 features ahead of the original repository, pull request has been already sent).
To download it, type following from snp_aaa directory:
```
git clone https://github.com/krassowski/biomart
```

Other dependencies (vcf and tabix parser) might be installed easily with use of pip:

```
pip install -r requirements.txt
```

You will also need to have `samtools` or at least `htslib` installed.

### Databases download

#### Directories to create

```
mkdir cosmic ensembl
```

#### COSMIC
To [download exported data](http://cancer.sanger.ac.uk/cosmic/download) from COSMIC databases [registration](https://cancer.sanger.ac.uk/cosmic/register) is required. Links to particular files are included in _data retrieval_ section. All the COSMIC files should go to the `cosmic` directory.

VCF files from Cosmic will need an additional processing - creation of tabix file will be needed. Make sure that you have [htslib](http://www.htslib.org) installed and run:
```
./create_tabix.sh filename
```

#### Ensembl

All files from Ensembl should be placed in `ensembl` directory. Those files could be downloaded using following script:
```
version=75
assembly=37

cd ensembl
if [ -d "$version" ]
	echo "This version of ensembl data is already installed"
else
	mkdir $version

	if (( $assembly < 76 ))
	then
		assembly="$assembly.$version"
	fi

	wget ftp://ftp.ensembl.org/pub/release-$ensembl_version/fasta/homo_sapiens/cds/Homo_sapiens.GRCh$assembly.cds.all.fa.gz
	wget ftp://ftp.ensembl.org/pub/release-$ensembl_version/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh$assembly.cdna.all.fa.gz
	wget ftp://ftp.ensembl.org/pub/release-$ensembl_version/fasta/homo_sapiens/dna/Homo_sapiens.GRCh$assembly.dna.chromosome.*.fa.gz
	wget ftp://ftp.ensembl.org/pub/release-$ensembl_version/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz

	if (( $assembly > 75 ))
	then
		wget ftp://ftp.ensembl.org/pub/release-$ensembl_version/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz.tbi
	else
		#./create_tabix.sh Homo_sapiens.vcf.gz
		echo "TODO: TABIX CREATION"
	fi
fi
```

When analyzing only somatic mutations one might want to use `Homo_sapiens_somatic.vcf.gz` instead of `Homo_sapiens.vcf.gz` since it is an order of magnitude smaller in size.

### Running somatic SNPs analysis

The script is written in Python 2.7 and has some basic command line interface. To start full analysis just run it without any additional options (not recommended if you have slow computer):

```
./snp_parser.py
```

You can adjust number of analysed genes with `-n` parameter:

```
./snp_parser.py -n 5
```

to get full list of available options, use `./snp_parser.py -h`.

Probably the most useful option is caching: `./snp_parser.py cache save` and `./snp_parser.py cache load`. I strongly recommend to use it always when experimenting with code or performing bigger analyses.

### Output

You will find the output (several text files) in `reports` directory.

The most interesting file for somatic SNPs analysis will be called 'Poly_A_and_expression_table.txt' and it would be in the following format:

```
# Poly A and expression table
#Gene	CNV+	CNV-	AAA+	AAA-
ENSG00000130066    4    1    2    1
```

Where:

* CNV+/- represents count of samples where the increase or decrease of expression was called in COSMIC database (compare the example data with corresponding COSMIC entry: [SAT1 ENSG00000130066](https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=SAT1#cnv_t), `CN type` column) 
* AAA+/- shows number of variants that create/elongate/shorten to poly(A) track.


## How to run tests

You can run basic tests of poly_aaa module simply running it from command line (as opposed to importing from another python package).

## About

The code in the repository as available before 2016-06-30 was written during 60-hours internship in Institute of Biochemistry and Biophysics Polish Academy of Sciences, under supervision and guidance of Paweł Szczęsny.


## Speed startup

```bash
python2 -OO snp_parser.py -g # generates GTEx key-value database from raw GTEx data, use just once
python2 -OO snp_parser.py -r --dataset hsapiens_snp -s 5 # downloads all mutations, use every time when changing datasets. Takes about 1 day
python2 -OO snp_parser.py --report list_poly_aaa_variants copy_number_expression poly_aaa_vs_expression --cache save # Two things here: 1. parse all mutations (so poly_aaa variants will be discovered here) (about 1 day) 2. report specified reports
# once mutations have been analyse, use 'load' instead of save (we keep results in cache) to re-run desired reports, e.g.:
python2 -OO snp_parser.py --report  poly_aaa_vs_expression --cache load
```

cd ncbi
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20161122.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20161122.vcf.gz.tbi
mv All_20161122.vcf.gz 00-All.vcf.gz
mv All_20161122.vcf.gz.tbi 00-All.vcf.gz.tbi

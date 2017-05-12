# Pipeline for SNP analysis with emphasis on influence of poly(A) tracks on protein expression in Human #

The repository contains multiple Python modules, useful when analyzing data from Ensembl, Cosmic, NCBI and other services while searching for information about effect of particular SNPs.
The pipeline is written in Python 2.7.
The major part of the analysis revolves around data from [PATACSDB](https://peerj.com/articles/cs-45/) and expands the research published in Science Advances: ["Translational control by lysine-encoding A-rich sequences"](http://advances.sciencemag.org/content/1/6/e1500154).

### Analyses

#### Poly(A)

Aim: Select only poly(A) related variants (such that make the track longer, shorter or does not change its length).

Data source: Ensembl's biomart (which incorporates data from NCBI and Cosmic too).

#### CNV vs poly(A)

Aim: Compare aggregations of mutations causing lengthening/shortening of poly(A) tracks and net effect of such groups with copy number expression of relevant region as reported by Cosmic.

Data sources: Cosmic v79

Genome: GRCh38 (Ensembl 87)

#### Expression vs poly(A)


Aim: Compare length / lengthening of poly(A) tracks by a variant with variant's effect sizes as defined in GTEx.

Data source: GTEx

Genome: GRCh37 (Ensembl 75)

Takeaway: not enough poly(A) related variants in GTEx as for 2017 to perform this analysis.

#### SPIDEX vs poly(A)

Aim: Compare length / lengthening of poly(A) tracks by a variant with its predicted impact on splicing inclusion index.

Data source: SPIDEX

Full genome:
./snp_parser.py -v ensembl --report poly_aaa_vs_spidex
Genome: GRCh37 (Ensembl 88)

Only genes known to have poly_aaa
./snp_parser.py -v biomart --report poly_aaa_vs_spidex
Genome: GRCh37 (Ensembl 75)


##### ZCRB1:c.411G>A

Aim: As above, but taking into account well tested mutation with known poly(A) track shortening effect.

Data source: Genomic coordinates from Cosmic

Genome: GRCh37 (Ensembl 75)

#### GTEx vs SPIDEX

Aim: Verification of an assumption about extended predictive capabilities of SPIDEX database.

./snp_parser.py -n --report gtex_vs_spidex
Genome: GRCh37 (Ensembl 75)

##### Motifs discovery

Preprare motifs analysis of mutations which are changed in the same direction in both spidex and GTEX
./snp_parser.py -n --report gtex_vs_spidex_for_motifs

## What is the workflow of particular modules? ##

### Somatic SNPs analysis  ###

The primary goal of this analysis was to create an association table showing relation between:

* expression of genes which might be affected by particular SNPs, and
* introduction, removal or change in length of poly(A) tracks in DNA sequence, caused by these SNPs.


#### Data retrieval

* Variants' names were obtained from Ensembl with use of Biomart API, with restrict to hsapiens_snp_som dataset (somatic mutations only),
* the variant's alleles data were taken from VCF files from three sources: [Ensembl](ftp://ftp.ensembl.org/pub/release-87/variation/vcf/homo_sapiens/), [COSMIC](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/VCF/CosmicCodingMuts.vcf.gz) and [NCBI](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh38p2/VCF/00-All.vcf.gz)  [tabix](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh38p2/VCF/00-All.vcf.gz.tbi),
* gene sequences (also from Ensembl) were downloaded as [chromosome-grouped DNA](ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/), [CDS](ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cds/) and [cDNA](ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cdna/) FASTA files.
* from COSMIC following data were downloaded:
    * [gene expression level 3 data](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/CosmicCompleteGeneExpression.tsv.gz) - from the TCGA portal
    * [all copy number abberations](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/CosmicCompleteCNA.tsv.gz) - used to map gene expression from samples to gene transcripts
    * [all COSMIC genes cDNA sequences](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/All_COSMIC_Genes.fasta.gz) - names mappings
* From [PATACSDB](http://sysbio.ibb.waw.pl/patacsdb) list of gene transcript identifiers having poly(A) tracks was retrieved (in CSV format)

The analysed variants set has been limited to: synonymous, stop gained, coding sequence and missense variants.  Only SNPs from genes present in PATACSDB (Homo sapiens dataset) were downloaded (so we had guarantee of at least one poly(A) track presence). I was not able to get the alleles from Ensembl's biomart as it is not providing the alleles of variants from COSMIC database at the time of writing (and those were crucial for the analysis!). From COSMIC genes cDNA sequences only headers were used in order to map the names of transcripts between Ensembl and COSMIC databases (worth noting, the distribution of data between transcripts from COSMIC were later in the analysis treated as uninformative, due to purported lack of proper curation i.e. there are serious premises multiple transcript-specific entries were nonetheless assigned to canonical transcript when incorporated to the COSMIC database).


#### General workflow

##### Discovery of Poly(A) related variants

SNP variants are fetched (hsapiens_snp) from Ensembl's biomart from genes such that poly(A) is guaranteed to occur (PATACSDB, human).

The polyadenylate track [poly(A)] was defined (the same as in PATACSDB) as 12A-1 (sequence of twelve subsequent adenines with one mismatch allowed). The code responsible for poly(A) detection is located in `poly_a` module.

##### Filling up gaps (getting more data about variants)
Variants retrieved from biomart and populated with additional data from VCF files are sorted by genes. Later in surroundings of each variant the search for poly(A) is performed - both in original and mutated sequence.


#### Future development ideas

It would be beneficial to split the final association table to separate tables for different variant consequences (synonymous, stop gained, coding sequence and missense) and even to include less restrictive filters on consequences (i.e. take all variants) but then have multiple tables showing associations with respect to consequences of variants. Some types of variant's consequences might be used as control - for example we are expecting to find no correlation in UTR related variants.

##  Summary of set up

### Dependencies

To use biomart API, the parser uses [fork of pip-distributed biomart](https://github.com/krassowski/biomart) package (which is 2 features ahead of the original repository, pull request has been already sent).
To download it, type following from snp_aaa directory:
```
git clone https://github.com/krassowski/biomart
```

Other dependencies (vcf and tabix parser) might be installed easily with use of pip:

```
python2 -m pip install -r requirements.txt
```
or
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

All files from Ensembl should be placed in `ensembl` directory. Those files could be downloaded using `download_ensembl.sh` bash script (note that it requires edits when using different genome assembly GRCh versions).

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

```bash
cd ncbi
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20161122.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20161122.vcf.gz.tbi
mv All_20161122.vcf.gz 00-All.vcf.gz
mv All_20161122.vcf.gz.tbi 00-All.vcf.gz.tbi
```

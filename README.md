# Analysis of SNPs affecting poly(A) tracks and their impact on protein expression in Human

The repository contains multiple Python modules, useful when analyzing data from Ensembl, Cosmic, NCBI and other services while searching for information about effect of particular SNPs.
The pipeline is written in Python 2.7.
The major part of the analysis revolves around data from [PATACSDB](https://peerj.com/articles/cs-45/) and expands the research published in Science Advances: ["Translational control by lysine-encoding A-rich sequences"](http://advances.sciencemag.org/content/1/6/e1500154).

### Analyses
All analyses run on genome assembly GRCh37 and use Ensembl release 88 unless stated otherwise.

Output of all analyses will be written to `.txt` files in `reports` directory.

The primary goal of analysis was to investigate effects of variants which elongate or shorten poly(A) motifs on:
* expression levels,
* copy number variations,

and other quantitative protein-level determinants.


#### Poly(A)

Aim: Select only poly(A) related variants (such that make the track longer, shorter or does not change its length).

Full genome:
    Data sources:
        - Raw Ensembl MySQL import data files (genes and variants)
```bash
./snp_parser.py --variants ensembl --report poly_aaa
```

Only genes known to have poly_aaa:
    Data sources:
        - Ensembl's biomart (variants from dbSNP, Cosmic and others),
        - PATACSDB (gene names)
```bash
./snp_parser.py --variants biomart --report poly_aaa
```


#### Expression vs poly(A)


Aim: Compare length / lengthening of poly(A) tracks by a variant with variant's effect sizes as defined in GTEx.

Data source: GTEx

Takeaway: not enough poly(A) related variants in GTEx as for 2017 to perform this analysis.

##### Using local copy of GTEx

Load GTEx expression database (use it just once, generates GTEx key-value database from raw GTEx data):

```bash
./snp_parser.py gtex --reload
```

Perform analysis:

```bash
./snp_parser.py --report poly_aaa_vs_expression
```

##### Using Ensembl API:

```bash
./snp_parser.py --report gtex_over_api
```

Results from offline version and from API may vary,
with the API version returning more hits
(which - strangely - are not present even on online GTEx page).

#### SPIDEX vs poly(A)

Aim: Compare length / lengthening of poly(A) tracks by a variant with its predicted impact on splicing inclusion index.

Data source: SPIDEX

Check effect of mutations changing length of poly_aaa in spidex database
```bash
./snp_parser.py --variants ensembl --report poly_aaa_vs_spidex
```

Full genome:
```bash
./snp_parser.py --variants ensembl --report poly_aaa_vs_spidex
```

Only genes known to have poly_aaa (this might be very slow):
```bash
./snp_parser.py --variants biomart --report poly_aaa_vs_spidex
```


#### GTEx vs SPIDEX

Aim: Verification of an assumption about extended predictive capabilities of SPIDEX database.

```bash
./snp_parser.py -no_variants --report gtex_on_spidex
```


##### Motifs of mutations having the same effect determined in GTEx and predicted in SPIDEX:

Aim: Find sequence motifs (if such exist) of mutations which are determined/predicted to change
expression in the same way bey both: SPIDEX and GTEx

Use DREME (control=shuffled input sequences):
```bash
./snp_parser.py -n --report gtex_on_spidex_motifs_dreme
```

Use DREME (control=whole sequences of transcripts of variants where mutation occurred):
```bash
./snp_parser.py -n --report gtex_on_spidex_motifs_dreme_with_control
```

Only prepare files for online MEME discriminative analysis:
```bash
./snp_parser.py -n --report gtex_on_spidex_motifs_meme_online
```


#### CNV vs poly(A) [outdated]

Aim: Compare aggregations of mutations causing lengthening/shortening of poly(A) tracks and net effect of such groups with copy number expression of relevant region as reported by Cosmic.

Data sources: Cosmic v79

Genome: GRCh38 (Ensembl 87)

Output: 'Poly_A_and_expression_table.txt' file in format:

```
# Poly A and expression table
#Gene	CNV+	CNV-	AAA+	AAA-
ENSG00000130066    4    1    2    1
```

Where:

* CNV+/- represents count of samples where the increase or decrease of expression was called in COSMIC database (compare the example data with corresponding COSMIC entry: [SAT1 ENSG00000130066](https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=SAT1#cnv_t), `CN type` column)
* AAA+/- shows number of variants that create/elongate/shorten to poly(A) track.

From COSMIC genes cDNA sequences only headers were used in order to map the names of transcripts between Ensembl and COSMIC databases (worth noting, the distribution of data between transcripts from COSMIC were later in the analysis treated as uninformative, due to purported lack of proper curation i.e. there are serious premises multiple transcript-specific entries were nonetheless assigned to canonical transcript when incorporated to the COSMIC database).

Following data were used:
* [gene expression level 3 data](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/CosmicCompleteGeneExpression.tsv.gz) - from the TCGA portal
* [all copy number abberations](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/CosmicCompleteCNA.tsv.gz) - used to map gene expression from samples to gene transcripts
* [all COSMIC genes cDNA sequences](http://cancer.sanger.ac.uk/cosmic/files?data=/files/grch38/cosmic/v79/All_COSMIC_Genes.fasta.gz) - names mappings


##### ZCRB1:c.411G>A [to be finished]

Aim: As above, but taking into account well tested mutation with known poly(A) track shortening effect.

Data source: Genomic coordinates from Cosmic

Genome: GRCh37 (Ensembl 75)

#### Miscellaneous

To get stats on ensembl variants, use:

```bash
./snp_parser.py ensembl --transcript_stats
```


## Data sources in depth

### Poly(A) tracks

The polyadenylate track [poly(A)] is defined (the same as in PATACSDB) as
12A-1 (sequence of twelve subsequent adenines with one mismatch allowed).

Search for poly(A) is performed locally:
for each variant a surrounding sequence is fetched (base on its CDS start and CDS end coordinates),
then in that sequence poly(A) presence is checked. Later original sequence is replaced with mutated one
and the search is repeated. If either search ends with a hit, the variant is classified as poly(A) related.

The code responsible for poly(A) detection is located in `poly_a` module.

### Sequences

Transcript CDS sequences were downloaded from Ensembl as fasta files.

### Variants

Following variants properties are determined by a source-specific 'variants_getter' script:
* identifiers,
* location (genomic and in cDNA),
* reference allele,
* lists of affected transcripts

Then all variants are parsed to retrieve alternative alleles (based on proper VCF file),
find out surrounding sequence (based on Ensembl CDS transcript database) and check reference correctness.

During parsing non-poly(A) related variants are rejected if requested, to speed up computations.

#### Alternative alleles

For some variants Ensembl does not provide alternative alleles correctly;
for example COSMIC mutations have alelle string like: "A/COSMIC_MUTATION",
whereas "fully" described mutations have strings like: "A/C".

In order to retrieve lacking alternative alleles and make sure
that all others are correct VCF files from three organizations:
Ensembl, Cosmic and NCBI
were used.

Some of the VCF files cover more than one source of mutations i.e.:
* NCBI VCF has both dbSNP and ClinVar variants,
* Ensembl covers ESP and HGMD-PUBLIC

#### Raw Ensembl

* All variants with consequences (so_term): coding_sequence_variant and all of its subterms.

Variants will be filtered at loading to reduce the tremendous quantity;
An option `--early_selection` defines criteria for filtering:
* `all_potential_poly_aaa` will load all variants which might affect or lay nearby poly(A) sequence
* `spidex_poly_aaa` will load only such poly(A) related variants which are single 1-1 substitutions, like variants analysed in SPIDEX database.
```bash
./snp_parser.py --variants ensembl ensembl --early_selection all_potential_poly_aaa
```


Consequences can be adjusted using `--consequences` option:
```bash
./snp_parser.py --variants ensembl ensembl --consequences stop_lost,stop_gained
```

#### Ensembl Biomart

* All variants with consequences (so_term): coding_sequence_variant.
* Only SNPs from genes present in [PATACSDB](http://sysbio.ibb.waw.pl/patacsdb) (Homo sapiens dataset),
so presence of at least one poly(A) track was guaranteed.


#### Future development ideas

It would be beneficial to split the final association table to separate tables for different variant consequences (synonymous, stop gained, coding sequence and missense) and even to include less restrictive filters on consequences (i.e. take all variants) but then have multiple tables showing associations with respect to consequences of variants. Some types of variant's consequences might be used as control - for example we are expecting to find no correlation in UTR related variants.

##  Installation & setup

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

`samtools` or at least `htslib` is required:
```
conda install -c bioconda samtools=1.4.1
```

Installing numba will speed up computations a lot:

```
conda install numba
```

Pigz provides parallel decompression capabilities:
```
conda install pigz
```

### Databases download

You will be securely asked for you Cosmic credentials to fetch necessary files from their FTP server.
```
./download_databases.sh
```
To adjust downloaded files, you need to manually modify relevant scripts, eg. to analyse only somatic mutations one should change `Homo_sapiens.vcf.gz` to `Homo_sapiens_somatic.vcf.gz`
in the download script.


## About

The code in the repository as available before 2016-06-30 was written during 60-hours internship in Institute of Biochemistry and Biophysics Polish Academy of Sciences, under supervision and guidance of Paweł Szczęsny.

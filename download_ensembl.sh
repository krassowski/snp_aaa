#!/usr/bin/env bash
#mkdir -p esp
#cd esp
#wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
#tar -zxvf ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz

ensembl_version=88
assembly=37

cd ensembl
if [ -d "av$ensembl_version" ]
then
    echo "Ensembl data from $ensembl_version release are already downloaded"
else
    mkdir -p "v$ensembl_version"
    cd "v$ensembl_version"

    if (( $ensembl_version < 76 ))
    then
        assembly="$assembly.$ensembl_version"
    fi

    wget ftp://ftp.ensembl.org/pub/grch$assembly/release-$ensembl_version/fasta/homo_sapiens/cds/Homo_sapiens.GRCh$assembly.cds.all.fa.gz
    wget ftp://ftp.ensembl.org/pub/grch$assembly/release-$ensembl_version/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh$assembly.cdna.all.fa.gz
    # download chromosomes, without pathes

    if (( $assembly > 37 ))
    then
        wget ftp://ftp.ensembl.org/pub/release-$ensembl_version/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz
        wget ftp://ftp.ensembl.org/pub/release-$ensembl_version/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz.tbi
    else
        wget ftp://ftp.ensembl.org/pub/grch$assembly/release-$ensembl_version/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz
        wget ftp://ftp.ensembl.org/pub/grch$assembly/release-$ensembl_version/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz.tbi
        #if (( $ensembl_version < 76 ))
        #then
        #    echo "Indexing vcf.gz file ..."
        #    ../../vcf_to_tabix.sh Homo_sapiens.vcf.gz
        #    mv Homo_sapiens.vcf.gz Homo_sapiens.vcf.gz.unsorted
        #    mv Homo_sapiens.vcf.gz.bgz Homo_sapiens.vcf.gz
        #    mv Homo_sapiens.vcf.gz.bgz.tbi Homo_sapiens.vcf.gz.tbi
        #fi
    fi


    #wget ftp://ftp.ensembl.org/pub/grch$assembly/release-$ensembl_version/fasta/homo_sapiens/dna/Homo_sapiens.GRCh$assembly.dna.chromosome.X.fa.gz
    #wget ftp://ftp.ensembl.org/pub/grch$assembly/release-$ensembl_version/fasta/homo_sapiens/dna/Homo_sapiens.GRCh$assembly.dna.chromosome.Y.fa.gz
    wget "ftp://ftp.ensembl.org/pub/grch$assembly/release-$ensembl_version/fasta/homo_sapiens/dna/Homo_sapiens.GRCh$assembly.dna.chromosome.*.fa.gz" -R "Homo_sapiens.GRCh$assembly.dna.chromosome.H*.fa.gz"
fi

shopt -s nullglob
for i in *.fa.gz;
do
	echo "Recompressing $i";
	gunzip $i -c | bgzip > $i.bgz;
	mv $i.bgz $i;
done
cd ..

#!/usr/bin/env bash
wget http://meme-suite.org/meme-software/4.11.3/meme_4.11.3_1.tar.gz
tar -xvzf meme_4.11.3_1.tar.gz
cd meme_4.11.3
./configure --prefix=$HOME/meme --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt
make
make test
make install
export PATH=$HOME/meme/bin:$PATH 


git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git checkout release/88
perl INSTALL.pl --ASSEMBLY GRCh37 --VERSION 88

ln -s ~/.vep/homo_sapiens/88_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

cd
cd .vep
gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz -c | bgzip > Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz.bgz
mv Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz.bgz Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

sudo apt-get install python3-dev graphviz libgraphviz-dev pkg-config
# assuming Ubuntu
pip3 install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/"
pip3 install -r requirements.txt

wget 'https://sourceforge.net/projects/song/files/Sequence%20Ontology/so_2_5_3/so_2_5_3.obo/download' -O so_2_5_3.obo

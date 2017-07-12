#!/usr/bin/env bash
grch=37
version=81
mkdir -p cosmic
cd cosmic
mkdir -p "v$version"
cd "v$version"
read -p "Please, enter your Cosmic account email: " email
read -p "Password: " -s pass
export SSHPASS="$pass"
sshpass -e sftp -oBatchMode=no -b - "$email"@sftp-cancer.sanger.ac.uk <<CMDS
get /files/grch$grch/cosmic/v$version/VCF/CosmicCodingMuts.vcf.gz
bye
CMDS

../../vcf_to_tabix.sh CosmicCodingMuts.vcf.gz



#!/bin/bash

# define colors
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'

# error to terminate if failed
error_cheker(){
    lastexit=$1
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR in exit value of last process${NC}"
        exit
    fi
}

# create env from yml file 
echo -e "${YEL}Creating wes_gatk ENV${NC}"
conda env create -f wes_gatk.yml
error_cheker $?

# create dir for all dbs and sources
mkdir -p $HOME/annDB/{annovar_source,bcftools}

# Nirvana
echo -e "${YEL}Cloning Nirvana from github${NC}"
cd $HOME/annDB
git clone https://github.com/Illumina/Nirvana.git
error_cheker $?

cd Nirvana
error_cheker $?

echo -e "${YEL}Installing Nirvana${NC}"
dotnet build -c Release
error_cheker $?

echo -e "${YEL}Downloading Nirvana DB${NC}"
mkdir -p ~/annDB/Nirvana/DB && cd ~/annDB/Nirvana
dotnet bin/Release/net*/Downloader.dll --ga GRCh38 -o DB
error_cheker $?

# Annovar
cd ~/annDB/annovar_source
echo -e "${YEL}Downloading annovar${NC}"
echo -e "${GRE}Annovar was registerred by a.ahmad@nu.edu.eg ${NC}"
wget -c http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz 
error_cheker $?
echo -e "${YEL}Installing annovar and dbs${NC}"
tar xvfz annovar.latest.tar.gz
cd annovar
for db in refGene avsnp150 gme gnomad211_exome cosmic70 dbnsfp31a_interpro 1000g2015aug clinvar_20221231 
do
    perl annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $db humandb/
    #error_cheker $?
done


# genome GFF for bcftools consequence 
echo -e "${YEL}Downloading Genome GFF for bcftools annotation${NC}"
cd ~/annDB/bcftools/hg38/
wget -c https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz

# downloading ref genome 
### TO DO: ( needs conversion from 2bit to fa ) ###
mkdir -p ~/annDB/ref/fa
cd ~/annDB/ref/fa
echo -e "${YEL}Downloading refrecne genome from : http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit${NC}"
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit
error_cheker $?


# downloading VEP db 
echo -e "${YEL}Downloading veb cache${NC}"
cd ~
vep_install -a cf -s homo_sapiens -y GRCh38 --CONVERT
error_cheker $?

echo -e "${GRE}seems that everything went well :D${NC}"

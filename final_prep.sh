#!/bin/bash

YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'

getdir=${PWD}

error_cheker(){
    lastexit=$1
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR in exit value of last process${NC}"
        exit
    fi
}


### you can download test samples from : https://zenodo.org/record/3243160/files/ ###

#### script ####
# mkdir samples
# cd samples
# for sample in father mother proband
# do
#     wget https://zenodo.org/record/3243160/files/${sample}_R1.fq.gz
#     wget https://zenodo.org/record/3243160/files/${sample}_R2.fq.gz
# done

## test genome :

# wget https://zenodo.org/record/3243160/files/hg19_chr8.fa.gz

# these samples are a trio subsampled to chr 8 only 



echo -e "${YEL}Creating wes_gatk ENV${NC}"
conda env create -f wes_gatk.yml
error_cheker $?


mkdir -p $HOME/annDB/{annovar_source,bcftools}

echo -e "${YEL}Cloning Nirvana from github${NC}"
cd $HOME/annDB
git clone https://github.com/Illumina/Nirvana.git
error_cheker $?

cd Nirvana
error_cheker $?

echo -e "${YEL}Reactivating ENV${NC}"
mamba create -n dotnet -y -c conda-forge dotnet=6.0.408
conda activate dotnet

echo -e "${YEL}Installing Nirvana${NC}"
dotnet build -c Release
error_cheker $?

echo -e "${YEL}Downloading Nirvana DB${NC}"
mkdir -p ~/annDB/Nirvana/DB && cd ~/annDB/Nirvana
dotnet bin/Release/net*/Downloader.dll --ga GRCh38 -o DB
error_cheker $?



cd ~/annDB/annovar_source
echo -e "${YEL}Downloading annovar${NC}"
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



echo -e "${YEL}Downloading Genome GFF for bcftools annotation${NC}"
cd ~/annDB/bcftools/hg38/
wget -c https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz



echo -e "${YEL}Downloading Nirvana DB${NC}"
mkdir -p ~/annDB/Nirvana/DB && cd ~/annDB/Nirvana
dotnet bin/Release/net*/Downloader.dll --ga GRCh38 -o DB
error_cheker $?


mkdir -p ~/annDB/ref/fa
cd ~/annDB/ref/fa
echo -e "${YEL}Downloading refrecne genome from : http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit${NC}"
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit
error_cheker $?

echo -e "${YEL}Downloading veb cache${NC}"
cd ~
mamba create -n vep -c bioconda ensembl-vep 
conda activate vep
mamba install -c conda-forge perl-compress-raw-zlib=2.202
vep_install -a cf -s homo_sapiens -y GRCh38 --CONVERT
error_cheker $?

echo -e "${GRE}seems that everything went well :D${NC}"

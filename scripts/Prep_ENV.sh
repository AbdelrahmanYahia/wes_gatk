#!/bin/bash

# define colors
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'
annovar_link="" ## to obtain the link, you need to register at "https://www.openbioinformatics.org/annovar/annovar_download_form.php"

script_dir="$(dirname "$(readlink -f "$0")")"
script_dir="${script_dir%'/scripts'}"

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
mamba env create -f ${script_dir}/workflow/env/wes_gatk.yml
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
mamba create -n dotnet -y -c conda-forge dotnet=6.0.408
conda activate dotnet
dotnet build -c Release
#dotnetPath=$(which dotnet)
error_cheker $?

echo -e "${YEL}Downloading Nirvana DB${NC}"
mkdir -p ~/annDB/Nirvana/DB && cd ~/annDB/Nirvana
dotnet bin/Release/net*/Downloader.dll --ga GRCh38 -o DB
error_cheker $?

# Annovar
cd ~/annDB/annovar_source
echo -e "${YEL}Downloading annovar${NC}"
wget -c "${annovar_link}" 
error_cheker $?
echo -e "${YEL}Installing annovar and dbs${NC}"
tar xvfz annovar.latest.tar.gz
# cd annovar
# for db in refGene avsnp150 gme gnomad211_exome cosmic70 dbnsfp31a_interpro 1000g2015aug clinvar_20221231 
# do
#     perl annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $db humandb/
#     #error_cheker $?
# done

# ## NOTE: Please modify this file to have the proper path, then uncomment and run
chmod 777 ${script_dir}/scripts/download_annovar_db.sh 
bash ${script_dir}/scripts/download_annovar_db.sh  ~/annDB/annovar_source/annovar

## not used!
# # genome GFF for bcftools consequence 
# echo -e "${YEL}Downloading Genome GFF for bcftools annotation${NC}"
# cd ~/annDB/bcftools/hg38/
# wget -c https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz

# # downloading ref genome 
# mkdir -p ~/annDB/ref/fa
# cd ~/annDB/ref/fa
# echo -e "${YEL}Downloading refrecne genome from : http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit${NC}"
# wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit
# error_cheker $?
# twoBitToFa hg38.analysisSet.2bit hg38.analysisSet.fa

# Download the GATK Resource bundle: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811
mkdir -p ~/annDB/broad_hg38 && cd ~/annDB/broad_hg38
source $script_dir/scripts/gatk_download_data.sh
# wget https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz


## not used!
# # downloading VEP db 
# echo -e "${YEL}Downloading veb cache${NC}"
# cd ~
# mamba create -n vep -c bioconda ensembl-vep 
# conda activate vep
# mamba install -c conda-forge perl-compress-raw-zlib=2.202
# vep_install -a cf -s homo_sapiens -y GRCh38 --CONVERT
# error_cheker $?

# echo -e "${GRE}seems that everything went well :D${NC}"

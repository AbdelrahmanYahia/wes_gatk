#!/bin/bash
# isntall GUAP source code in home

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

echo -e "${YEL}Creating wes_gatk ENV${NC}"
# Create a conda environment
conda create -n wes_gatk -y
error_cheker $?

echo -e "${YEL}Activating ENV${NC}"
eval "$(conda shell.bash hook)"
conda activate wes_gatk

echo -e "${YEL}Installing Packages${NC}"
conda install -n wes_gatk -y -c bioconda bwa fastqc samtools vcftools gatk4 picard qualimap bcftools multiqc bedtools ensembl-vep 
error_cheker $?

conda install -n wes_gatk -y -c conda-forge pip
error_cheker $?

conda install -n wes_gatk -y -c anaconda pandas psutil pyyaml argparse tqdm

# Install dotnet
conda install -n wes_gatk -y -c conda-forge dotnet
error_cheker $?


echo -e "${YEL}Cloning Nirvana from github${NC}"

# Clone and build Nirvana
cd ~
git clone https://github.com/Illumina/Nirvana.git
error_cheker $?

cd Nirvana
error_cheker $?

echo -e "${YEL}Reactivating ENV${NC}"
conda deactivate 
conda activate wes_gatk


echo -e "${YEL}Installing Nirvana${NC}"
dotnet build -c Release
error_cheker $?


echo -e "${YEL}Downloading Nirvana DB${NC}"
# Download the database
cd ~
mkdir -p Nirvana/DB
cd Nirvana
dotnet bin/Release/net*/Downloader.dll --ga GRCh38 -o DB
error_cheker $?

echo -e "${GRE}seems that everything went well :D${NC}"

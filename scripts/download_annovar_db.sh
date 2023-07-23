#!/bin/bash

# color codes for command line 
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
CYN='\x1b[1;36m'
NC='\e[0m'

annovar_path=$1

try_download(){
    db_Name=$1
    echo -e "${YEL}Trying to download: ${NC}${db_Name} ..." 
    ${annovar_path}/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $db_Name ${annovar_path}/humandb/
    lastexit=$?
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR IN:${NC} Downloading ${db_Name} from ANNOVAR's server, trying UCSC..."
        ${annovar_path}/annotate_variation.pl -downdb -buildver hg38 $db_Name ${annovar_path}/humandb/
        lastexit=$?
        if [[ $(( lastexit )) -ne 0 ]];then
            echo -e "${RED}DOWNLOADING ${db_Name} FAILED!${NC}"
        else
            echo -e "${GRE}DOWNLOADING ${db_Name} SUCCESS!${NC}"
        fi
    else
        echo -e "${GRE}DOWNLOADING ${db_Name} SUCCESS!${NC}"
    fi
}

dbs_list_annovar=( 
    refGene 
    ensGene 
    esp6500siv2_all 
    avsnp150 
    dbnsfp42c 
    1000g2015aug
)

## source: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/

dbs_list_UCSC=( 
    phastConsElements100way 
    cytoBand
    wgRna genomicSuperDups dgvMerged 
    gwasCatalog 
    snp151 
    keggPathway keggMapDesc 
    clinvarMain clinvarCnv
    encRegTfbsClustered encTfChipPkENCFF998KDQ
)


dbs_wgENCODE=(
    wgEncodeTreatment
    wgEncodeRegDnaseClustered
    wgEncodeGencodePdbV43
    wgEncodeGencodePseudoGeneV43
    wgEncodeRegDnaseUwHacHotspot
)

other_UCSC_dbs=(
    orfeomeMrna
    lincRNAsTranscripts
)


for i in ${dbs_list_annovar[@]} ; do 
    echo -e "${YEL}Checking: ${NC}${i} ..."
    [ -n "$(find ${annovar_path}/humandb/ -maxdepth 1 -name "hg38_${i}*")" ] || try_download $i
done

for i in ${dbs_list_UCSC[@]} ; do 
    echo -e "${YEL}Checking: ${NC}${i} ..."
    [ -n "$(find ${annovar_path}/humandb/ -maxdepth 1 -name "hg38_${i}*")" ] || try_download $i
done


### you can turn those of if needed

for i in ${dbs_wgENCODE[@]} ; do 
    echo -e "${YEL}Checking: ${NC}${i} ..."
    [ -n "$(find ${annovar_path}/humandb/ -maxdepth 1 -name "hg38_${i}*")" ] || try_download $i
done

for i in ${other_UCSC_dbs[@]} ; do 
    echo -e "${YEL}Checking: ${NC}${i} ..."
    [ -n "$(find ${annovar_path}/humandb/ -maxdepth 1 -name "hg38_${i}*")" ] || try_download $i
done

 





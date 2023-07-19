#!/bin/bash

# color codes for command line 
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
CYN='\x1b[1;36m'
NC='\e[0m'

try_download(){
    db_Name=$1
    echo -e "${YEL}Trying to download: ${NC}${db_Name} ..." 
    ~/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $db_Name ~/annovar/humandb/
    lastexit=$?
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR IN:${NC} Downloading ${db_Name} from ANNOVAR's server, trying UCSC..."
        ~/annovar/annotate_variation.pl -downdb -buildver hg38 $db_Name ~/annovar/humandb/
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

dbs_list_annovar=( refGeneWithVer knownGene ensGene 
esp6500siv2_all avsnp150 dbnsfp42c 1000g2015aug)

dbs_list_UCSC=( phastConsElements46way cytoBand tfbsConsSites
wgRna targetScanS genomicSuperDups dgvMerged gwasCatalog 
wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75 
wgEncodeRegTfbsClustered snp153 
wgEncodeUwDnaseSeqHotspotsRep2Gm12878 exac10 )

for i in ${dbs_list_annovar[@]} ; do 
    echo -e "${YEL}Checking: ${NC}${i} ..."
    [ -n "$(find ~/annovar/humandb/ -maxdepth 1 -name "hg38_${i}*")" ] || try_download $i
done

for i in ${dbs_list_UCSC[@]} ; do 
    echo -e "${YEL}Checking: ${NC}${i} ..."
    [ -n "$(find ~/annovar/humandb/ -maxdepth 1 -name "hg38_${i}*")" ] || try_download $i
done
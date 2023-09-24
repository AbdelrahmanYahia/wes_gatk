#!/bin/bash
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'
# usage message
__usage="
Usage: $(basename $0) [OPTIONS]

extract_genes.sh -v <inputvcf> -i <gene_info> -o <outputdir> -g <genes> -t INT 

Options:
  -v <str>                              Input VCF file path
  -i <str>                              Input gene info file path
  -o <str>                              Output directory path
  -g <str>                              Comma-delimited list of genes to extract
  -t <int>                              Number of threads              default = 2
  -h                                    Help message ( This message )
"

usage() {
    echo -e "$__usage" 1>&2; exit 1
}

THREADS=2

while getopts hi:o:g:t:v: OPTION; do
    case "${OPTION}" in
        v) INPUT=${OPTARG};; 
        i) gene_info=${OPTARG};; 
        o) OUTPUT=${OPTARG};;
        t) THREADS=${OPTARG};;
        h) usage;;
        g) GENES=${OPTARG};;
        \?) usage;;
    esac
done

if [ -z "${INPUT}" ] || [ -z "${OUTPUT}" ] || [ -z "${GENES}" ] || [ -z "${gene_info}" ]; then
    echo -e "${RED}Missing required argument!${NC}"
    usage
fi

mkdir -p "$OUTPUT"

error_cheker(){
    lastexit=$1
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR in exit value of last process${NC}"
        exit
    fi
}

extract_genes() {
    id=$1
    out=$2

    interval=$(grep -i "Name=$id;" $gene_info | awk -F"\t" '{print "chr"$1":"$4-1000"-"$5+1000}')
    error_cheker $?

    echo -e "${YEL}Extracting ${id} [$interval] from file...${NC}"

    bcftools view -Ov "$INPUT" -r "$interval" --threads $THREADS --output-file "$out"
    error_cheker $?

}


IFS=',' read -ra GENE_ARRAY <<< "$GENES"
for gene in "${GENE_ARRAY[@]}"; do
    output_vcf="$OUTPUT/${gene}.vcf"
    extract_genes "$gene" "$output_vcf"
done

output_concat="$OUTPUT/all_genes.vcf"
bcftools concat -o "$output_concat" "$OUTPUT"/*.vcf
error_cheker $?

echo -e "${GRE}Genes extracted and concatenated successfully.${NC}"
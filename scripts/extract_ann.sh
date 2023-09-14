## Extract annotation of a target gene
mkdir -p ens_ann && cd ens_ann
wget https://ftp.ensembl.org/pub/release-91/gff3/homo_sapiens/Homo_sapiens.GRCh38.91.gff3.gz
zcat Homo_sapiens.GRCh38.91.gff3.gz | grep "ID=gene" > gene.info
#webget https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
conda activate dotnet
snv="/mnt/gs21/scratch/mansourt/wes_gatk/output/05_Annotation/Nirvana/snvs/Annotation.json.gz"
Jasix="$HOME/annDB/Nirvana/bin/Release/net6.0/Jasix.dll"
export_Nirvana="/mnt/gs21/scratch/mansourt/wes_gatk/scripts/export_Nirvana_csv.py"
outdir="/mnt/gs21/scratch/mansourt/wes_gatk/output/05_Annotation/Nirvana/snvs/extracted_genes" && mkdir -p $outdir

id=ATG101
interval=$(grep -i "Name=$id;" gene.info | awk -F"\t" '{print "chr"$1":"$4-1000"-"$5+1000}')
echo $interval
dotnet "$Jasix" -i "$snv" -q "$interval" -t | tr -d '\n' | sed 's/{  \"chromosome/\n{\"chromosome/g' | sed 's/ //g' | sed "s/]}]}]}/]}]}\n]}/" | gzip > $id.json.gz
python3 "$export_Nirvana" -i $id.json.gz -o "$outdir"/ann.tsv --extract-genes "$id","test" --extracted-genes-outdir "$outdir"


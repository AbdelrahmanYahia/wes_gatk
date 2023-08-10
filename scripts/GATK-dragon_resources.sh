destination=$1
echo "$destination"
if ( -z "$destination"  ) then 
    cmd=""
    echo "no dir" 
else
    cmd="--directory-prefix $destination"
    mkdir -p "$destination"
    echo "found dir"
fi





#wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/GRCh38_gencode.v27.refFlat.txt"
wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.UD"
wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.mu"
#wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.V"
wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.bed"
#wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.exome_calling_regions.v1.README.sh"
#wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.exome_calling_regions.v1.UD"
#wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.exome_calling_regions.v1.bed"
#wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.exome_calling_regions.v1.mu"
#wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/README"
wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt"

# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G_phase3_v4_20130502.sites.hg38.vcf"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G_phase3_v4_20130502.sites.hg38.vcf.idx"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/GRCh38.primary_assembly.genome.fa"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dict"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.amb"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.ann"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.bwt"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.pac"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.sa"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.ref_cache.tar.gz"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.tile_db_header.vcf"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.vid"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/WholeGenomeShotgunContam.vcf"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/WholeGenomeShotgunContam.vcf.idx"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/autosomes-1kg-minusNA12878-ALL.vcf"
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/TileDB/Homo_sapiens_assembly38.tile_db_header.vcf" 
# wget $cmd -c "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/TileDB/Homo_sapiens_assembly38.vid" 

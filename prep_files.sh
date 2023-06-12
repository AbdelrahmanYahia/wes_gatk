python3 ./wes.py WES \
    -i /home/marc/Desktop/data/samples/complex_wes/ \
    -o /home/marc/Desktop/output_new_wes_gatk \
    --reference-fasta /home/marc/Desktop/data/refs/fa/hg38/hg38.fa \
    --bed-file /home/marc/Desktop/data/refs/BED/Twist_Exome_Core_Covered_Targets_hg38.bed \
    --gff-file /home/marc/Desktop/data/refs/gtf/Homo_sapiens.GRCh38.109.gff3.gz \
    --nirvana-path /home/marc/Nirvana \
    --annovar-path /home/marc/annovar \
    --known-variants /home/marc/Desktop/data/refs/vcf/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz \
    --reference-index /home/marc/Desktop/data/refs/indexes/bwa/hg38/hg38 \
    --generate-confs-only \
    --threads 12
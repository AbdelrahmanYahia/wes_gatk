
import os
import pandas as pd

# samples = ["father", "mother", "proband"]
samples = ["proband", "mother", "father"]
aligners = ["bwa"]
variant_callers = ["GATK"]

RS = "R"
EXT = "fastq"

samples_dir = "/home/marc/WES_GATK/test/samples"
out_dir = "/home/marc/WES_GATK/test/snake_out2"
ref_bwa = "/home/marc/Desktop/data/refs/indexes/bwa/hg38/hg38"
ref_bwa_path = "/home/marc/Desktop/data/refs/indexes/bwa"
ref_fasta = "~/Desktop/data/refs/fa/hg38/hg38.fa"
ref_fasta_path = "/home/marc/Desktop/data/refs/fa/hg38"
bed_file = "/home/marc/Desktop/data/refs/BED/Twist_Exome_Core_Covered_Targets_hg38.bed"
known_variants = "/home/marc/Desktop/data/refs/vcf/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
nirvana_path = "/home/marc/Nirvana"
Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll"
Nirvana_supplementray = f"{nirvana_path}/Data/SupplementaryAnnotation/GRCh38/"
Nirvana_ref = f"{nirvana_path}/Data/References/Homo_sapiens.GRCh38.Nirvana.dat"
Nirvana_cache = f"{nirvana_path}/Data/Cache/GRCh38/Both"
gff = "/home/marc/Desktop/data/refs/gtf/Homo_sapiens.GRCh38.109.gff3.gz"
env_path = "/home/marc/WES_GATK/test_workflow/envs"


workdir: out_dir
include: "utils.smk"
include: "alignment.smk"
include: "bam_processing.smk"
include: "QC.smk"
include: "variant_calling.smk"
include: "variant_processing.smk"
include: "annotation.smk"


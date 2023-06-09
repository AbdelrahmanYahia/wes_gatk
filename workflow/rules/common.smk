import os
import pandas as pd

RS = "R"
EXT = "fastq"

samples_dir = "/home/marc/WES_GATK/test/samples"

out_dir = config["output"]


ref_bwa = config["reference_index"]

ref_bwa_path = confi["reference_output_path"]
ref_prefix = confi["reference_output_prefix"]

ref_fasta = config["refence_fasta"]
ref_fasta_path = os.path.dirname(ref_fasta)

bed_file = config["bed_file"]

known_variants = config["known_variants"]


nirvana_path = config["nirvana_path"]
annovar_dir = config["annovar_path"]

Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll"
Nirvana_supplementray = f"{nirvana_path}/DB/SupplementaryAnnotation/GRCh38/"
Nirvana_ref = f"{nirvana_path}/DB/References/Homo_sapiens.GRCh38.Nirvana.dat"
Nirvana_cache = f"{nirvana_path}/DB/Cache/GRCh38/Both"


gff = config["gff_file"]


workdir: out_dir
include: "utils.smk"
include: "alignment.smk"
include: "bam_processing.smk"
include: "QC.smk"
include: "variant_calling.smk"
include: "variant_processing.smk"
include: "annotation.smk"


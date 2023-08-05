import os
import pandas as pd

# get options
PATH = config["path"]
EXTT = config["ext"]
RS = config["R"]
if config["decompress"]:
    EXT = EXTT.replace(".gz","")
else:
    EXT = EXTT

working_dir = config["working_dir"]
sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file)

files_R1s = set(SampleTable.iloc[:, 1-1])
files_R2s = set(SampleTable.iloc[:, 8-1])
samples = set(SampleTable.iloc[:, 3-1]) # sample full name
units = set(SampleTable.iloc[:, 2-1])
samples_IDs = set(SampleTable.iloc[:, 4-1])
# print(SampleTable)
units = (
    pd.read_csv('samples.tsv', sep="\t", dtype={"sample_id": str, "unit": str})
    .set_index(["sample_id", "unit"], drop=False)
    .sort_index()
)

ALL_THREADS = config["threads"]
MEM = config["total_mem"]
GUAP_FOLDER = config["GUAP_DIR"]
R = [1, 2]
source = PATH

working_dir = config["working_dir"]
source_dir = config["GUAP_DIR"]
common_rules = config["common_rules"]

samples_dir = config["input"]

out_dir = config["output"]

ref_bwa = config["reference_index"]

ref_bwa_path = config["reference_output_path"]
ref_prefix = config["reference_output_prefix"]

ref_fasta = config["reference_fasta"]
ref_fasta_path = os.path.dirname(ref_fasta)

bed_file = config["bed_file"]

unique_values = set(line.split('\t')[0] for line in open(bed_file))

known_variants_snps = config["known_variants_snps"]
known_variants_indels = config["known_variants_indels"]
known_variants_indels2 = config["known_variants_indels2"]

nirvana_path = config["nirvana_path"]
annovar_dir = config["annovar_path"]

Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll"
Nirvana_supplementray = f"{nirvana_path}/DB/SupplementaryAnnotation/GRCh38/"
Nirvana_ref = f"{nirvana_path}/DB/References/Homo_sapiens.GRCh38.Nirvana.dat"
Nirvana_cache = f"{nirvana_path}/DB/Cache/GRCh38/Both"

gff = config["gff_file"]

call_cnv = config["call_CNV"]

include: "01_utils.smk"
include: "02_samples_preprocessing.smk"
include: "03_alignment.smk"
include: "04_bam_processing.smk"
include: "05_variant_calling.smk"
include: "06_variant_processing.smk"
include: "07_annotation.smk"
include: "08_CNV.smk"
include: "09_CNV_postprocessing.smk"  
include: "10_QC.smk"
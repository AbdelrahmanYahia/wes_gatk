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

TAIL = f"_{config['tail']}"
sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file)

files_R1s = list(SampleTable.iloc[:, 0])
files_R2s = list(SampleTable.iloc[:, 14])
samples = list(SampleTable.iloc[:, 1]) # sample full name
library_index = list(SampleTable.iloc[:, 3])
samples_IDs = list(SampleTable.iloc[:, 2])
samples_names = list(SampleTable.iloc[:, 1])



ALL_THREADS = config["threads"]
MEM = config["total_mem"]
GUAP_FOLDER = config["GUAP_DIR"]
R = [1, 2]
source = PATH

working_dir = config["working_dir"]
source_dir = config["GUAP_DIR"]
common_rules = config["common_rules"]


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


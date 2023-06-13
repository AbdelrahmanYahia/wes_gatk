## THIS WORKFLOW IS STILL UNDER DEVELOPMENT AND IS NOT READY YET!

## Introduction

**Note: This workflow is still under development.**

WES GATK is a flexible and user-friendly whole exome sequencing workflow based on [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). It is designed for processing Illumina WES short reads data and features automatic sample table generation, Snakemake configuration file, and simplified workflow execution.

## Workflow

![Workflow](workflow.png)

## Quick Start

### Installation

Clone the repository:
```
git clone https://github.com/AbdelrahmanYahia/wes_gatk.git
```


If you prefer to download the dependencies manually, you can find them [here].

### Step-by-Step Installation Guide

#### Register ANNOVAR

To obtain the link, you need to register at [Annovar website](https://www.openbioinformatics.org/annovar/annovar_download_form.php).

Make sure you have Conda and Mamba installed. If not, follow these steps:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh # Follow the prompts
```


Install Mamba:
```
conda install mamba -n base -c conda-forge

```

Restart your terminal.

#### Prepare the Environment

Change the permissions for the `Prep_ENV.sh` file:
```
chmod 777 wes_gatk/scripts/Prep_ENV.sh
```


Create the environment, install the tools, and download the annotations database:
```
./wes_gatk/scripts/Prep_ENV.sh ANNOVAR_LINK
```


This process may take some time.

You can also download all the required reference files using `wes_gatk/scripts/gatk_download_data.sh`:
```
chmod 777 wes_gatk/scripts/gatk_download_data.sh
bash wes_gatk/scripts/gatk_download_data.sh DOWNLOAD_DIR
```



### Running the Analysis

To start the analysis, activate the `wes_gatk` environment and run the `wes.py` file:
```
conda activate wes_gatk
python3 wes.py WES --help
```

You can also use the Python file to generate the sample table and the config file, and then run Snakemake independently. Modify the parameters according to your needs:
```
conda activate wes_gatk
python3 wes.py WES \
  --input PTH/to/samples \
  --output PTH/to/outdir \
  --reference-fasta broad_hg38/Homo_sapiens_assembly38.fasta \
  --bed-file exome_bed/S07604715_Padded.bed \
  --gff-file broad_hg38/Homo_sapiens.GRCh38.109.gff3.gz \
  --nirvana-path ~/Nirvana \
  --annovar-path ~/annovar_source/annovar \
  --known-variants broad_hg38/1000G_omni2.5.hg38.vcf.gz \
  --reference-index broad_hg38/Homo_sapiens_assembly38.fasta \
  --generate-confs-only
```

To run a Snakemake dry-run:
```
conda activate wes_gatk
snakemake \
  --snakefile wes_gatk/workflow/Snakefile \
  -c PTH/to/outdir/config.yml \
  -n -j THREADS
```
### Advanced Parameters

For advanced usage, you can refer to the following command-line options:
usage: Basic Run Usage example:
    guap WES -i indir -o outdir --bed-file file --reference-fasta fasta.fasta --reference-index indexpath 
        
```
options:
  -h, --help            show this help message and exit

basic config:
  -i in path, --input in path
                        Input directory path
  -o out path, --output out path
                        Output directory path

Workflow configure:
  --threads N           Number of total threads to use [default = all]
  --reference-fasta path/to/file.fa
                        path to reference fasta file
  --bed-file path       bed file path
  --gff-file path       gff file path
  --nirvana-path path   Path for Nirvana
  --annovar-path path   Path for annovar
  --generate-confs-only
                        Generate sample table and config file only

Aligner configuration:
  --threads-index N     Number of threads to use during indexing ref [default
                        = 4]
  --threads-align N     Number of threads to use during sample alignment
                        [default = 4]
  --aligner-extra-args '-args'
                        Extra arguments for aligner, use it with no spaces and
                        add = ( --aligner-extra-args='-arg1 -arg2' )
  --reference-index path/to/ref
                        path to reference index
  --reference-output-path path/to/ref
                        path to reference index
  --reference-output-prefix path/to/ref
                        path to reference index
  --index-fasta         Index fasta file

Variant caller configuration:
  --known-variants path
                        path to reference fasta file

Snakemake Options:
  --dry-run             performs snakemake dry run
  --export-dag          performs snakemake dry run and exports DAG
  --smk-extra-args ='-args'
                        A string value of extra args for snakemake(must be
                        used with = with no spaces (--smk-extra-args='-arg1
                        -arg2'))
  --parse-snakemake-output
                        prints progress bar instead of snakemake regular
                        output

Annotation configuration:
  --annovar-protocol str
                        Annovar Protocol defaults: refGene,avsnp150,clinvar_20
                        221231,cosmic70,dbnsfp31a_interpro,EAS.sites.2015_08,E
                        UR.sites.2015_08,gme,gnomad211_exome,SAS.sites.2015_08
  --annovar-operation str
                        Annovar Protocol defaults: g,f,f,f,f,f,f,f,f,f

Other:
  --continue            continue analysis when re-run
  --overwrite           overwrite output dir if exsits
  -n str, --name str    Name of files [ default = guap_run[date time] ]
  --verbose             verbose
  --quit                print many output
  --print-last-run      Prints last run on screen

```

import pandas as pd
import argparse
import json
import gzip
import sys
import os

parser = argparse.ArgumentParser(description='Exports Nirvana output to CSV')
parser.add_argument(
    '--input', '-i',
    help='path to Nirvana gz json', 
    required=True
    )

parser.add_argument(
    '--output', '-o', 
    help='output dir', 
    required=True
    )  

parser.add_argument(
    '--jasix-path', '-j', 
    help='full path to jasix.dll', 
    default="$HOME/annDB/Nirvana/bin/Release/net6.0/Jasix.dll"
    )  

parser.add_argument(
    '--gene-info', '-g', 
    help='full path to gene.info file',
    required=True
    ) 

parser.add_argument(
    '--compine-genes', 
    help='Compine all genes in one file', 
    action='store_true'
    )  

parser.add_argument(
    '--genes', 
    help='comma seprated gene list', 
    default=None,
    required='--gene-list' not in sys.argv
    )  

parser.add_argument(
    '--padding', 
    help='Padding for gene region', 
    default=1000
    )  

parser.add_argument(
    '--gene-list', 
    help='File containing gene each gene in a line', 
    default=None,
    required='--genes' not in sys.argv
    )  

args = parser.parse_args()

file = args.input
Jasix = os.path.expanduser(args.jasix_path)
gene_info = args.gene_info
padding = args.padding

directory_path = args.output
os.makedirs(directory_path, exist_ok=True)

def extract_gene_region(Ann_file, gene_info, id, padding=1000):
    interval = os.popen(f"grep -i 'Name={id};' {gene_info} | awk -F'\t' '{{print \"chr\"$1\":\"$4-{padding}\"-\"$5+{padding}}}'").read().strip()
    if interval == '':
        print("didn't find Gene in Gene list!")
        return None
    else:
        print(F"\nExporting {id} with padding [ {padding} ] ({interval})  from Annotation file...\n")
        dotnet_command = f"dotnet {Jasix} -i {file} -q {interval}"
        dotnet_output = os.popen(dotnet_command)
        dotnet_output_str = dotnet_output.read().replace("\n","")
        dotnet_output_exitcode = dotnet_output.close()
        if dotnet_output_exitcode is not None:
            print(f"ERROR with dotnet command! maybe didn't find the gene {id}")
            return None
        else:
            formatted_output = dotnet_output_str.replace("{  \"chromosome", "\n{\"chromosome")
            formatted_output = formatted_output.replace(" ", "")
            formatted_output = formatted_output.replace('{"positions":[\n{', '[{').replace("}]}]}]}","}]}]}]")

            return pd.DataFrame(json.loads(formatted_output))

def process_region(*args):
    counter = 0
    for df in args:
        if counter == 0:
            my_df = df
        else:
            my_df = pd.concat([my_df, df], axis=0, ignore_index=True)
        counter += 1
        
    return my_df

def process_dataframe(df,target_gene_symbol=None):
    df = pd.DataFrame(df)
    try:
        variants_df = df['variants'].apply(pd.Series)
        variants_df2 = variants_df[0].apply(pd.Series)
        final_df = pd.concat([df, variants_df2], axis=1).drop(['variants', 'samples'], axis=1)
        if target_gene_symbol is not None:
            final_df.to_csv(f"{directory_path}/{target_gene_symbol}.tsv",sep='\t',index=False)
        else:
            final_df.to_csv(f"{directory_path}/Multiple_genes.tsv",sep='\t',index=False)
        return final_df
    

    except:
        print("No match found for any of the provided genes!")



def print_stats(final_df,target_gene_symbol=None):
    try:
        total_rows = len(final_df)
    except:
        exit(1)
    total_variants = total_rows
    try:
        variants_with_annotations = total_rows - final_df['dbsnp'].isna().sum()
    except:
        variants_with_annotations = 0
    try:
        variants_with_regulatoryRegions = total_rows - final_df['regulatoryRegions'].isna().sum()
    except:
        variants_with_regulatoryRegions = 0
    try:
        variants_with_inlowcomplexity = total_rows - final_df['inLowComplexityRegion'].isna().sum()
    except:
        variants_with_inlowcomplexity = 0
    novel_variants = total_variants - variants_with_annotations
    percentage_novel = round(((novel_variants / total_rows) * 100),2)
    if target_gene_symbol is not None:
        print(f"#--------------------{target_gene_symbol}--------------------#")
    print(f"Total Variants: {total_variants}")
    print(f"No. Variants with dbSNP Annotations: {variants_with_annotations}")
    print(f"No. Variants affecting regualtory regions: {variants_with_regulatoryRegions}")
    print(f"No. Variants in low complexity regions: {variants_with_inlowcomplexity}")
    print(f"Percentage of Unknown Variants (dbsnp only): {percentage_novel} %\n")



def parse_genes_arg(genes_arg):
    if ',' in genes_arg:
        return genes_arg.split(',')

if args.gene_list is not None:
    print("Sorry still under development, use --genes instead!")
    exit(0)

Genes  =  parse_genes_arg(args.genes)
egenes = []

if args.compine_genes:
    for gene in Genes:
        egenes.append(extract_gene_region(file, gene_info, gene, padding))
    combined = process_region(*egenes)
    final_df = print_stats(process_dataframe(combined))
else:
    for gene in Genes:
        egene = extract_gene_region(file, gene_info, gene, padding)
        processed = process_dataframe(egene,target_gene_symbol=gene)
        final_df = print_stats(processed,target_gene_symbol=gene)
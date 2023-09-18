import pandas as pd
import argparse
import json
import gzip
import sys
import os

parser = argparse.ArgumentParser(description='Exports Nirvana output to CSV')
parser.add_argument('--input', '-i', help='path to Nirvana gz json', required=True)
parser.add_argument('--output', '-o', help='output file name path (tsv)', required=True)  
parser.add_argument('--export-genes', help='Export genes df  from Nirvana json file', action='store_true')  
parser.add_argument('--extract-genes', help='comma seprated gene list or file with gene symbols', default=None)  
parser.add_argument('--extracted-genes-outdir', help='Output dir for extracted genes', required='--extract-genes' in sys.argv)  

args = parser.parse_args()

file = args.input

header = ''
positions = []
genes = []
is_header_line = True
is_position_line = False
is_gene_line = False
gene_section_line = '],"genes":['
end_line = ']}'

with gzip.open(file, 'rt') as f:
    position_count = 0
    gene_count = 0
    for line in f:
        trim_line = line.strip()
        if is_header_line:
            ## only keep the "header" field content from the line
            header = trim_line[10:-14]
            is_header_line = False
            is_position_line = True
            continue
        if trim_line == gene_section_line:
            is_gene_line = True
            is_position_line = False
            continue
        elif trim_line == end_line:
            break
        else:
            if is_position_line:
                ## remove the trailing ',' if there is
                positions.append(trim_line.rstrip(','))
                position_count += 1
            if is_gene_line:
                ## remove the trailing ',' if there is
                genes.append(trim_line.rstrip(','))
                gene_count += 1

print ('number of positions:', position_count)
print ('number of genes:', gene_count)
print('Processing file, it might take some time...')

mylist = []
for position in positions:
    position_dict = json.loads(position)
    mylist.append(position_dict)
df = pd.DataFrame(mylist)
variants_df = df['variants'].apply(pd.Series)
variants_df2 = variants_df[0].apply(pd.Series)
final_df = pd.concat([df, variants_df2], axis=1).drop(['variants', 'samples'], axis=1)

final_df.to_csv(args.output,sep='\t',index=False)

if args.export_genes:
    print('exporting genes...')
    mylist2 = []
    for gene in genes:
        position_dict = json.loads(gene)
        mylist2.append(position_dict)
    genedf = pd.DataFrame(mylist2)
    directory = os.path.dirname(args.output)
    genedf.to_csv(f"{directory}/Genes.tsv",sep='\t',index=False)


genes_arg = args.extract_genes

if genes_arg != None:
    def parse_genes_arg(genes_arg):
        if ',' in genes_arg:
            return genes_arg.split(',')
        else:
            try:
                with open(genes_arg, 'r') as file:
                    return [line.strip() for line in file.readlines()]
            except FileNotFoundError:
                raise argparse.ArgumentTypeError(f"File '{genes_arg}' not found.")
            
    def filter_by_gene_symbol(row):
        if isinstance(row, list):
            for d in row:
                if isinstance(d, dict) and d.get('hgnc') == target_gene_symbol:
                    return True
        return False  


    def is_novel(row):
        return pd.isna(row['dbsnp'])
    directory_path = args.extracted_genes_outdir
    os.makedirs(directory_path, exist_ok=True)
    Genes  =  parse_genes_arg(genes_arg)
    for target_gene_symbol in Genes:
        filtered_df = final_df[final_df['transcripts'].apply(filter_by_gene_symbol)]
        if len(filtered_df) > 0:
            total_rows = len(filtered_df)
            total_variants = total_rows
            variants_with_annotations = total_rows - filtered_df['dbsnp'].isna().sum()
            novel_variants = total_variants - variants_with_annotations
            percentage_novel = (novel_variants / total_rows) * 100
            print()
            print(f"#--------------------{target_gene_symbol}--------------------#")
            print(f"Total Variants: {total_variants}")
            print(f"Variants with Annotations: {variants_with_annotations}")
            print(f"Variants without dbSNP Annotations: {novel_variants}")
            print(f"Percentage of Novel Variants: {percentage_novel} %")
            filtered_df.to_csv(f"{directory_path}/{target_gene_symbol}.tsv",sep='\t',index=False)
        else:
            print()
            print(f"Didn't find any variants in Gene: {target_gene_symbol}")
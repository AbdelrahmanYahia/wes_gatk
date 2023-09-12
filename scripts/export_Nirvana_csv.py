import pandas as pd
import argparse
import json
import gzip


parser = argparse.ArgumentParser(description='Exports Nirvana output to CSV')
parser.add_argument('--input', '-i', help='path to Nirvana gz json', required=True)
parser.add_argument('--output', '-o', help='output file name path', required=True)  

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

mylist = []
for position in positions:
    position_dict = json.loads(position)
    mylist.append(position_dict)
df = pd.DataFrame(mylist)
variants_df = df['variants'].apply(pd.Series)
variants_df2 = variants_df[0].apply(pd.Series)
final_df = pd.concat([df, variants_df2], axis=1).drop(['variants', 'samples'], axis=1)

final_df.to_csv(args.output,sep='\t',index=False)
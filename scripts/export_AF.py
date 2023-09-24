import matplotlib.pyplot as plt
import argparse
import os


parser = argparse.ArgumentParser(description='Exports Nirvana output to CSV')
parser.add_argument(
    '--input', '-i',
    help='path to vcf', 
    required=True
    )

parser.add_argument(
    '--output', '-o', 
    help='output file', 
    required=True
    )  

args = parser.parse_args()
infile = args.input
outfile = args.output

output = os.popen(f"bcftools stats {infile} | grep '^AF' | cut -f3").read().strip()

allele_frequencies = [float(line) for line in output.split('\n')]

plt.bar(range(len(allele_frequencies)), allele_frequencies, color='blue')
plt.xlabel('Index')
plt.ylabel('Allele Frequency')
plt.title('Allele Frequency Distribution')

plt.savefig(outfile)
print("DONE exporting file")
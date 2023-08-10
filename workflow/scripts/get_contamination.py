import csv
import sys
import argparse

# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f

parser = argparse.ArgumentParser(description='Get contamination from selfSM file')
parser.add_argument('--input', help='selfSM file', required=True)
parser.add_argument('--output', help='output file', required=True)  
args = parser.parse_args()

# colors 
red = "\033[1;31m"
green = "\033[0;32m"
yellow = "\033[1;33m"
blue = "\033[1;34m"
nc = "\033[0m"
with open(args.input) as selfSM:
    reader = csv.DictReader(selfSM, delimiter='\t')
    i = 0
    for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
        # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
        # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
        # vcf and bam.
            sys.stderr.write(f"{red}Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).{nc}")
            sys.exit(1)
        print(float(row["FREEMIX"])/0.75)
        with open(args.output, 'w') as output:
            output.write(str(float(row["FREEMIX"])/0.75))

        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
            sys.stderr.write(f"{red}Found {i} rows in .selfSM file. Was expecting exactly 1. This is an error{nc}")
            sys.exit(2)
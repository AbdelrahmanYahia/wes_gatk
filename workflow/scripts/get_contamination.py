import sys
import csv
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Creates contamintaion file')
parser.add_argument('--output', help='file', required=True)  
parser.add_argument('--threads', help='number of threads', required=True)
parser.add_argument('--prefix', help='prefix for files', required=True)
parser.add_argument('--svdprefix', help='prefix for files', required=True)
parser.add_argument('--get_con_file', help='path to get_contamination.py', required=True)
parser.add_argument('--target_bam', help='path to target bam file', required=True)
parser.add_argument('--ref', help='path to reference file', required=True)

args = parser.parse_args()

threads = args.threads
prefix = args.prefix
svdprefix = args.svdprefix
get_con_file = args.get_con_file
target_bam = args.target_bam
ref = args.ref
output = args.output

varify_outputname = prefix + ".selfSM"

red = "\033[1;31m"
green = "\033[0;32m"
yellow = "\033[1;33m"
blue = "\033[1;34m"
nc = "\033[0m"

p = subprocess.Popen("verifybamid2 --NumPC 4 --NumThread {threads} --Output {prefix} --BamFile {target_bam} --Reference {ref} --SVDPrefix {svdprefix}", 
                     shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)

stdout, stderr= p.communicate()

if p.returncode != 0:
    with open(varify_outputname, 'w') as f:
        pass

    with open(output, 'w') as f:
        f.write("0")

else:
    with open(varify_outputname) as selfSM:
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

import pandas as pd
import threading
import argparse
import json
import gzip
import sys
import os
import concurrent.futures
from queue import Queue

parser = argparse.ArgumentParser(description='Exports Nirvana output to CSV')
parser.add_argument('--input', '-i', help='path to Nirvana gz json', required=True)
parser.add_argument('--output', '-o', help='output file name path (tsv)', required=True)
parser.add_argument('--export-genes', help='Export genes df from Nirvana json file', action='store_true')
parser.add_argument('--extract-genes', help='comma-separated gene list or file with gene symbols', default=None)
parser.add_argument('--extracted-genes-outdir', help='Output dir for extracted genes', required='--extract-genes' in sys.argv)
parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads to use for processing')
parser.add_argument('--chunks', '-c', type=int, default=10000, help='Number of lines to process at the same time, [default = 10,000]')

args = parser.parse_args()

def process_nirvana_chunk(chunk, output_file_prefix):
    header = ''
    positions = []
    genes = []
    is_header_line = True
    is_position_line = False
    is_gene_line = False
    gene_section_line = '],"genes":['
    end_line = ']}'

    position_count = 0
    gene_count = 0

    for line in chunk:
        trim_line = line.strip()
        if is_header_line:
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
                positions.append(trim_line.rstrip(','))
                position_count += 1
            if is_gene_line:
                genes.append(trim_line.rstrip(','))
                gene_count += 1

    mylist = []
    for position in positions:
        position_dict = json.loads(position)
        mylist.append(position_dict)
    df = pd.DataFrame(mylist)
    variants_df = df['variants'].apply(pd.Series)
    variants_df2 = variants_df[0].apply(pd.Series)
    final_df = pd.concat([df, variants_df2], axis=1).drop(['variants', 'samples'], axis=1)

    output_file = f"{output_file_prefix}_positions.tsv"
    final_df.to_csv(output_file, sep='\t', index=False, mode='a', header=False)  # Append to the output file

    if args.export_genes:
        mylist2 = []
        for gene in genes:
            position_dict = json.loads(gene)
            mylist2.append(position_dict)
        genedf = pd.DataFrame(mylist2)
        gene_output_file = f"{output_file_prefix}_genes.tsv"
        genedf.to_csv(gene_output_file, sep='\t', index=False, mode='a', header=False)  # Append to the output file

    print(f"Processed {position_count} positions and {gene_count} genes in {output_file}")

def process_file_threaded(input_file, output_file_prefix, chunk_size, line_queue):
    print(f"Processing file: {input_file}")
    with gzip.open(input_file, 'rt') as f:
        chunk = []
        line_count = 0
        for line in f:
            line_count += 1
            line_queue.put(line)
            if line_queue.qsize() >= chunk_size:
                chunk = []
                while not line_queue.empty():
                    chunk.append(line_queue.get())
                process_nirvana_chunk(chunk, output_file_prefix)
        while not line_queue.empty():
            chunk.append(line_queue.get())
        if chunk:  # Process any remaining lines in the last chunk.
            process_nirvana_chunk(chunk, output_file_prefix)

    print(f"Finished processing {input_file} with {line_count} lines")

def main():
    input_file = args.input
    output_file_prefix = args.output
    chunk_size = args.chunks
    line_queue = Queue()

    # Create a thread pool
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        # Function to read the file and put lines into the queue
        def read_and_enqueue():
            with gzip.open(input_file, 'rt') as file:
                for line in file:
                    line_queue.put(line)

        # Start the reading thread
        reading_thread = threading.Thread(target=read_and_enqueue)
        reading_thread.start()

        # Start worker threads to process the chunks
        for _ in range(args.threads):
            executor.submit(process_file_threaded, input_file, output_file_prefix, chunk_size, line_queue)

        # Wait for all threads to finish
        executor.shutdown()
        reading_thread.join()

    print("All threads have finished processing.")

if __name__ == "__main__":
    main()
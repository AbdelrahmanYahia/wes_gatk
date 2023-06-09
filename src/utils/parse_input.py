"""This script parses the input dir and creates the sample table and config files"""

## need to add lane information recognition
## need to add validate input
import re 
import sys
import yaml 
import shutil
import pandas as pd
from collections import defaultdict
# GUAP modules import
from .globals import *
import subprocess
from tqdm import tqdm

def check_extension(df): # takes pandas df and returns string
    # checks all files have same extension from pandas df, to use in generete sample table function
    uniques = df['ext'].unique()
    if len(uniques) > 1:
        glogger.prnt_fatel(f"Your input directory has multible fastq file extensions, please check directory.")

    else:
        return uniques[0]


def check_PE(df): # takes pandas df and returns string
    """checks all files either single or paried ended from pandas df, to use in generete sample table function"""
    uniques = df['PE'].unique()
    if len(uniques) > 1:
        glogger.prnt_fatel(f"Your input directory has both Paired and single files, please check directory.")
    else:
        return uniques[0]


def check_R(df): # takes pandas df and returns string
    """checks all files have same naming patterns from pandas df, to use in generete sample table function"""
    uniques = df['read_num'].unique()
    if len(uniques) > 1:
        glogger.prnt_fatel(f"Your input directory has multible fastq file naming patterns, please check directory.")
    else:
        return uniques[0].replace("1","")


def check_pattern(df): # takes pandas df and returns a string 
    """checks all files have same naming patterns from pandas df, to use in generete sample table function"""
    uniques = df['matched_pattern'].unique()
    if len(uniques) > 1:
        glogger.prnt_fatel(f"Your input directory has multible fastq file naming patterns, please check directory.")
    else:
        return uniques[0]
    

def recogize_pattern(file_name): # takes string of fastq file name and returns dict with read info and id
    """ using re to recognize the naming pattern of samples (illumina, srr and general naming patten)"""
    # naming pattern for re 
    patterns = { # ! fix (_|\.) group for R pattern in dict config !
        "illumina": "(((.+)_(S\d+)_(L00\d))_(R1|R2|r1|r2|read1|read2)_(00\d)\.(fastq\.gz|fastq|fq\.gz|fq))",
        "SRR": "(((SRR)(\d+))(_|\.)(1|2|R1|R2|r1|r2|read1|read2)\.(fastq\.gz|fastq|fq\.gz|fq))",
        "general": "(((.+))(_|\.)(1|2|R1|R2|r1|r2|read1|read2)\.(fastq\.gz|fastq|fq\.gz|fq))"
    }

    matched_pattern = None
    ## loop on pattern to and checks whichs one matches 
    ## starting with illumina because general would match any ways
    ## breaks once successful 
    for ptrn_name, pattern in patterns.items():
        try:
            matched = re.match(pattern, file_name) 
        except:
            continue
        
        if bool(matched) :
            matched_pattern = ptrn_name
            break
            
        else:
            continue

    # Extracts sample information
    if matched_pattern == "illumina":
        file_name, sample_name, sample_id, sample_number, read_num, lane, tail, ext = matched.groups()[0], matched.groups()[1], matched.groups()[2], matched.groups()[3], matched.groups()[5], matched.groups()[4], matched.groups()[6], matched.groups()[7]

    elif matched_pattern == "SRR":
        glogger.prnt_fatel(f"{RED}{matched_pattern}{NC} is currntly not supported sample naming pattern\nonly {GRE}'Illumina'{NC} naming pattern is supported at the moment")
        file_name, sample_name, sample_id, sample_number, read_num, lane, tail, ext = matched.groups()[0], matched.groups()[1], matched.groups()[3], "", matched.groups()[5], "", "", matched.groups()[6]

    elif matched_pattern == "general":
        glogger.prnt_fatel(f"{RED}{matched_pattern}{NC} is currntly not supported sample naming pattern\nonly {GRE}'Illumina'{NC} naming pattern is supported at the moment")
        file_name, sample_name, sample_id, sample_number, read_num, lane, tail, ext = matched.groups()[0], matched.groups()[1], matched.groups()[1], "", matched.groups()[4], "", "", matched.groups()[5]

    else:
        glogger.prnt_fatel(f"{RED}Your Samples Pattern is an unfamiler pattern.{NC}\nPlease contact my Developpers and they will look into it :D")
        file_name = sample_name = sample_id = sample_number = read_num = lane = tail = ext = None

    # Returns a dictionary of sample information
    return {
        "file_name": file_name,
        "sample_name": sample_name,
        "sample_id": sample_id,
        "sample_number": sample_number,
        "read_num": read_num,
        "lane": lane,
        "tail": tail,
        "ext": ext,
        "matched_pattern": ptrn_name
    }

def qiime_table_checker(df):
    if ("sample-id" or "id" or "sampleid" or "sample id") in df.columns:
        if "#q2:types" in (df.iloc[:1]).values.tolist()[0]:
            return True
        else:
            return False
    else:
        return False

def check_metadata(args):
    # samples dataframe and output path declaration 
    samples = parse_samples(os.path.abspath(args.input))
    samples_IDs = list(samples.iloc[:, 2])
    outpath = os.path.abspath(args.output)

    # check metadata file 
    if args.metadata is None or args.create_metadata_table:
        if args.create_metadata_table:
            pass
        else:
            # prompt the user to create metadata file 
            create_meta = input(f"{RED}No metadata file supplied,{NC} do I create an empty one with sample IDs? (y/n) ")
            if create_meta == ( 'y' or 'Y'):
                pass
            else:
                # exiting 
                glogger.prnt_fatel("No metadata supplied or created!")
        # creating metadata file 
        glogger.prnt_warning(f"I will create one at '{args.output}', {YEL}\nplease re-run the analysis with: \n -m {args.output}/sample-metada.tsv after modifing the file\n{NC}(check https://docs.qiime2.org/2021.11/tutorials/metadata/)")
        # creating empty data frame to store metadata info
        header = pd.DataFrame({"1": ["#q2:types", "categorical"]}).T
        header.columns = ["sample-id", "condition"]

        new_samples = {}
        c = 2
        # get sample ids
        for sample in samples_IDs:
            new_samples[c] = [sample, "BLANK"]
            c += 1
        # store sample IDs and BLANK at new df, and export
        samples_df = pd.DataFrame(new_samples).T
        samples_df.columns = ["sample-id", "condition"]
        samples_df = samples_df.sort_values(["sample-id"])
        metadata_file = pd.concat([header, samples_df])
        # Create the directory if it doesn't exist
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        metadata_file.to_csv(outpath+"/"+"sample-metadata.tsv",sep='\t',index=False) 
        args.metadata = f"{outpath}/sample-metadata.tsv"
        glogger.prnt_fatel(f"{GRE}Metadata empty file created, {RED}Exiting...{NC}")

    else:
        if not os.path.isfile(args.metadata):
            glogger.prnt_fatel(f"{RED}{args.metadata} Doesn't exist!{NC}")
        else:
            # checking metadata file extension 
            metadataname, metadataextension = os.path.splitext(args.metadata)
            if metadataextension == (".tsv"):
                the_file = pd.read_csv(args.metadata,sep="\t")
            elif metadataextension == (".csv"):
                the_file = pd.read_csv(args.metadata,sep="\t")
            else:
                glogger.prnt_fatel(f"{RED}{args.metadata} Have strange extension (use: tsv or csv){NC}")

            rest_of_cols = the_file.columns[1:]
            # checks all samples are present in the file with proper names
            if False in samples['sample_id'].isin(the_file.T.iloc[0]).tolist():
                glogger.prnt_fatel(f"{RED} Please check {args.metadata} sample IDs!\n{YEL}Note: {NC}You can re-run your code WITHOUT suppling a metadata file and \nGUAP will ask to create an empty one for you with the samples IDs.")
     
            else:
                # check qiime2 format and modifing if not
                if qiime_table_checker(the_file):
                    if metadataextension == (".csv"):
                        the_file.to_csv(outpath+"/"+"sample-metadata.tsv",sep='\t',index=False) 
                        args.metadata = "{outpath}/sample-metadata.tsv"
                else:
                    print(f"{RED}Note:{NC}Modifing sample metadata file to be compatible with QIIME2")
                    columns_names = rest_of_cols.to_list()
                    columns_names.insert(0, 'sample-id')
                    columns_names
                    t = f'{(len(columns_names) -1 )* "categorical,"}'.split(",")[:-1]
                    t.insert(0, "#q2:types")
                    header_N = pd.DataFrame({1: t}).T
                    header_N.columns = columns_names
                    the_file.columns = columns_names
                    last_file = pd.concat([header_N, the_file])
                    last_file.to_csv(outpath+"/"+"sample-metadata.tsv",sep='\t',index=False) 
                    args.metadata = f"{outpath}/sample-metadata.tsv"

def parse_samples(inpath): # takes path return contains fastq files, returns df contains sample information
    ## takes input path
    ## gets the file names containg fastq and fq
    ## performs the recogize_pattern function to 
    ## capture sample information and stores it in 
    ## pandas df
    # input path to absolute path
    path = os.path.abspath(inpath)
    # list all files
    all_files = os.listdir(path)
    samples = defaultdict(dict)
    # takes fastq files only
    for file_name in all_files:
        if os.path.isfile(path + "/" + file_name) and ("fastq" in file_name or "fq" in file_name):
            # Captures the file path and name
            filename, file_extension = os.path.splitext(file_name)
            if "fastq" in filename or "fq" in filename:
                filename, new_ext = os.path.splitext(filename)
                file_extension = new_ext + file_extension
            # recogize_pattern function returns a dictitionary with sample names, id, and read information
            sample_info = recogize_pattern(file_name)

            # get only forward reads and replace the read number to get R2
            # appends sample information to a dict of dicts
            if "1" in sample_info["read_num"]:
                read_2 = sample_info["read_num"].replace("1","2")
                if sample_info["matched_pattern"] == "illumina":
                    read_1 = f"{sample_info['read_num']}_{sample_info['tail']}.f"
                    read_2 = f"{read_2}_{sample_info['tail']}.f"
                else:
                    read_1 = f"{sample_info['read_num']}.f"
                    read_2 = f"{read_2}.f"

                f2 = file_name.replace(read_1, read_2)
                if f2 in all_files:
                    sample_info["file2"] = f2
                    sample_info["PE"] = True
                    samples[sample_info["sample_id"]] = sample_info

                else:
                    sample_info["file2"] = ""
                    sample_info["PE"] = False

    # converts the dict to pandas df and returns the df
    m_samples = samples
    samples= pd.DataFrame(samples).T
    samples = samples.sort_values(by=['sample_id'])

    return samples

# sample table info
### file_name sample_name sample_id read_num lane tail ext matched_pattern file2 PE
### uniques = df['PE'].unique()

def parse_input_args(args): # takes args (object) returns dict of args information
    global glogger
    # set vars
    global verbose
    try:
        verbose = args.verbose
    except:
        verbose = False
    
    glogger.create_console_handler(verbose=verbose)
    glogger.prnt_info("Logger is set to INFO")
    outpath = os.path.abspath(args.output)
    path = os.path.abspath(args.input)

    # create output path if doesn't exsit 
    if not os.path.exists(outpath):
        os.mkdir(outpath)
        glogger.create_file_handler(f"{outpath}/main_log.txt")
        glogger.prnt_info(f"Created Working Dir @ {outpath}")

    else:        
        try:
            if args.overwrite:
                # Remove directory and all its contents
                shutil.rmtree(outpath)
                os.mkdir(outpath)
                glogger.create_file_handler(f"{outpath}/main_log.txt")
                glogger.prnt_warning(f"Overwrited working dir {outpath}")
            else:
                glogger.create_file_handler(f"{outpath}/main_log.txt")

        except:
            glogger.create_file_handler(f"{outpath}/main_log.txt")


    # creates a dict with input args 
    all_args = vars(args)

    if all_mem < 7:
        glogger.prnt_error("Your System doesn't have enough memory to run the analyis")

    # validate samples
    samples = parse_samples(args.input)

    # check n of threads and use all threads if not supplied 
    if args.threads is None:
        glogger.prnt_warning(f"You did not specify any number of threads, I will use all ({all_threads}). ")
        args.threads = all_threads

    elif args.threads < 4:
        glogger.prnt_error(f"GUAP can't use less than 4 threads, would you please change the number of threads?  ")
        new_threads = int(input("Number of threads: (4 and above, enter any thing else to exit) "))
        try:
            int(new_threads)
            if int(new_threads) < 4:
                glogger.prnt_fatel(f"Value is smaller than 4. exiting...")
            else:
                args.threads = new_threads
        except ValueError:
            glogger.prnt_fatel(f"exiting...")


    ext = str(check_extension(samples))
    PE = bool(check_PE(samples))
    R = str(check_R(samples))
    compressed = False
    EXT = ext
    pattern = str(check_pattern(samples))
    tail = "001" if pattern == "illumina" else ""

    # to perform gunzipping 
    if ".gz" in ext:
        compressed = True
        EXT = ext.replace(".gz","")
    
    # check if analysis run before and created sample table 
    if os.path.exists(outpath+"/"+"samples.tsv"):
        glogger.prnt_warning(f"Found an exsiting sample.tsv file in output directory, will not override.")
    else:
        samples.to_csv(outpath+"/"+"samples.tsv",sep='\t')  
    
    extra_info = {
        "path": path,
        "working_dir": outpath,
        "ext": ext,
        "tail": tail,
        "R": R,
        "naming_pattern": pattern,
        "R1_pattern": f"_{R}1_{tail}.{EXT}",
        "R2_pattern": f"_{R}2_{tail}.{EXT}",
        "compressed" : compressed,
        "total_mem": all_mem,
        "GUAP_DIR": GUAP_DIR,
        "common_rules": f"{GUAP_DIR}/workflows/common/rules/"
    }

    if "decompress" not in all_args:
        all_args.update({"decompress":False})
    all_args.update(extra_info)
    # create config file 
    with open(f'{outpath}/config.yaml', 'w') as yaml_file:
        yaml.safe_dump(all_args, yaml_file, default_flow_style=False, sort_keys=False)
    
    # Check if a command was given as an argument
    if len(sys.argv) > 1:
        # Get the first argument (excluding the script name)
        command = " ".join(sys.argv[1:])
        
        # Write the command to a file
        with open(f"{outpath}/command.txt", "w") as f:
            f.write(command)
        # Write the command to a lastcommand file
        with open(f"{GUAP_DIR}/.last_run.txt", "w") as f:
            f.write(command)
    return all_args


def process_snakemake_standard_output(snakemake_cmd, outfilelog):
    time_stamp_pattern = r'\[\w{3} \w{3} \d{2} \d{2}:\d{2}:\d{2} \d{4}\]'
    rule_name_pattern = r'rule (\w+):$'
    jobid_pattern = r'    jobid: (\d+)'
    done_status_pattern = r'(\d+) of (\d+) steps \((\d+)\%\) done'
    total_pattern = r'total\s+(\d+)\s+\d+\s+\d+'
    error_pattern = r'.*one of the commands exited with non-zero exit code.*'
    rule_error_pattern = r"Error in rule .*"
    exiting_message = r"Exiting because a job execution failed. Look above for error message"
    no_more_jobs = r"Nothing to be done (all requested files are present and up to date)."
    running_jobs = {}
    finished_jobs = {}
    failed_jobs = {}
    proc = subprocess.Popen(f"{snakemake_cmd} ", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    def check_job_id(line, rule_name):
        job_id_match = re.search(jobid_pattern, line)
        if job_id_match:
            job_number = job_id_match.group(1)
            job_number = int(job_number)
            running_jobs[job_number] = rule_name
            progress_bar.set_description(f"Performing {rule_name} job no. {job_number}")
            check_next_line(line)
        elif re.search(time_stamp_pattern, line):
            check_next_line(line)
        else:
            try:
                next_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                outfile.write(next_line + "\n")
                check_job_id(next_line,rule_name)
            except StopIteration:
                    progress_bar.close()

    def check_next_line(line):
        global rule_name, job_number
        try:
            next_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
            outfile.write(next_line + "\n")

            rule_name_match = re.search(rule_name_pattern, next_line)
            finished_match = re.search(r'Finished job (\d+)\.', next_line)
            finished_status_match = re.search(done_status_pattern, next_line)
            exiting_message_match = re.search(exiting_message,next_line)
            nothing_match = re.search(no_more_jobs,next_line)
            rule_error_pattern_match = re.search(r"Error in rule (\w+).*", next_line)
            if rule_error_pattern_match:
                # glogger.prnt_fatel(f"Error in {rule_error_pattern_match.group(1)}")
                current_rule_error = rule_error_pattern_match.group(1)
                try:
                    following_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                    outfile.write(following_line + "\n")
                    job_id_match = re.search(jobid_pattern, following_line)
                    if job_id_match:
                        job_number = job_id_match.group(1)
                        job_number = int(job_number)
                        failed_jobs[job_number] = current_rule_error
                        tqdm.write(f"{YEL}job {job_number}, of Rule ({current_rule_error}) had and {RED}Error{NC}")
                        check_next_line(following_line)
                    else:
                        progress_bar.close()
                        glogger.prnt_fatel(f"Something seems to be wrong!\nCheck log for errors as I found an error in rule: {current_rule_error}")
                except StopIteration:
                    progress_bar.close()
                    glogger.prnt_fatel(f"{RED_}Something went wrong!{NC}")
            # print(next_line)
            if re.search(time_stamp_pattern, next_line):
                check_next_line(next_line)
            elif rule_name_match:
                rule_name = rule_name_match.group(1)
                check_job_id(next_line,rule_name)
            elif finished_match:
                finished_job_number = finished_match.group(1)
                finished_jobs[rule_name] = finished_job_number
                check_next_line(next_line)
            elif finished_status_match:
                current_finished_job = finished_status_match.group(1)
                total_jobs = finished_status_match.group(2)
                percentage = finished_status_match.group(3)
                percentage = int(percentage)
                if percentage == 100:
                    progress_bar.set_description(f"{GRE}All process finished :D {NC}")
                    progress_bar.update((int(percentage)- progress_bar.n))
                else:
                    progress_bar.update((int(percentage)- progress_bar.n))
                # print(f"Progress update: {finished_job_number}")
                check_next_line(next_line)
            elif exiting_message_match:
                tqdm.write(f"{YEL}An error occurred. Stopping progress bar and exiting.{NC}")
                try:
                    following_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                    outfile.write(following_line + "\n")
                    last_line = re.search(r'^Complete log:.*', following_line)
                    if last_line:
                        progress_bar.close()
                        glogger.prnt_fatel(f"A job or more has Failed, check log file\nFailed jobs:\n{failed_jobs}")
                    else:
                        progress_bar.close()
                        glogger.prnt_fatel(f"Something seems to be wrong!\nCheck log for errors as I found an error in rule: {current_rule_error}")
                except StopIteration:
                    progress_bar.close()
                    glogger.prnt_fatel(f"{RED_}Something went wrong!{NC}")
            elif nothing_match:
                tqdm.write(f"{GRE}Nothing to be done (all requested files are present and up to date){NC}")
                progress_bar.close()

            else:
                check_next_line(line)
        except StopIteration:
            progress_bar.close()

    def check_time_stamp(line):
        global rule_name, job_number
        # Check if line matches timestamp pattern
        timestamp_match = re.search(time_stamp_pattern, line)
        if timestamp_match:
            # A new timestamp means a new job is starting or has finished
            rule_name = None
            job_number = None
            outfile.write(line + "\n")
            check_next_line(line)
        else:
            outfile.write(line + "\n")
            check_next_line(line)

    rule_name = None
    job_number = None
    with open(outfilelog, "w") as outfile:
        with tqdm(total=100,desc="Workflow running", leave=True, bar_format="{l_bar}{bar}| [ Elapsed: {elapsed} ]") as progress_bar:
            for line in iter(proc.stdout.readline, b''):
                line = line.decode('utf-8').rstrip()
                try:
                    check_time_stamp(line)
                except Exception as E:
                    glogger.prnt_fatel(f"Error in checking time stamp:\n{RED_}{E}{NC}")
  


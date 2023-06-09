from ..workflows import *

class WES(WorkflowCli):
    name = 'WES'
    help = '''guap WES -i/--input dir'''
    usage = f"""{YEL}Basic Run Usage example:{NC}
    guap WES -i indir -o outdir --bed-file file \
            --reference-fasta fasta.fasta \
            --reference-index indexpath 
        """

    def add_arguments(self, parser):

        # basic configuration
        basic_conf = parser.add_argument_group(f'{CYN}basic config{NC}')

        basic_conf.add_argument(
            '-i', '--input', 
            help="Input directory path", 
            metavar='in path', 
            type=os.path.abspath, 
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)
        ) 

        basic_conf.add_argument(
            '-o', '--output', 
            help= "Output directory path", 
            metavar='out path', 
            type=os.path.abspath, 
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)
        )  


        # workflow configure
        workflow_conf = parser.add_argument_group(f'{CYN}Workflow configure{NC}')

        workflow_conf.add_argument(
            '--threads', 
            metavar = "N",
            help= "Number of total threads to use [default = all]", 
            type=int,
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)

        )

        workflow_conf.add_argument(
            '--reference-fasta',
            metavar='path/to/file.fa',
            type=os.path.abspath,
            help="path to reference fasta file",
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)

        )

        workflow_conf.add_argument(
            '--bed-file', 
            help='bed file path', 
            metavar='path',
            type=os.path.abspath,
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)

        )

        workflow_conf.add_argument(
            '--gff-file', 
            help='gff file path', 
            metavar='path',
            type=os.path.abspath,
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)

        )

        workflow_conf.add_argument(
            '--nirvana-bath', 
            help='Path for Nirvana', 
            metavar='path',
            type=os.path.abspath,
            default="$HOME/annDB/Nirvana"

        )

        workflow_conf.add_argument(
            '--annovar-bath', 
            help='Path for annovar', 
            metavar='path',
            type=os.path.abspath,
            default="$HOME/annDB/annovar_source/annovar"
        )

        workflow_conf.add_argument(
            '--generate-confs-only', 
            help='Generate sample table and config file only', 
            action='store_true'
        )


        # qc_conf = parser.add_argument_group(f'{BLU}QC configuration{NC}')

        # qc_conf.add_argument(
        #     '--trimmomatic',
        #     dest='trimmomatic',
        #     action='store_true',
        #     help="Use trimmomatic"
        # )

        # qc_conf.add_argument(
        #     '--trim-t', 
        #     help= "Number of threads to use during trim step", 
        #     type=int ,
        #     metavar = "N",
        #     default= 4 
        # )

        # qc_conf.add_argument(
        #     "--trim-min-length", 
        #     type=int,
        #     default=30,
        #     metavar = "N",
        #     help='trimmomatic min length [default = 30]'
        # )

        # qc_conf.add_argument(
        #     "--slidingwindow-size", 
        #     type=int,
        #     default=4,
        #     metavar = "N",
        #     help='trimmomatic sliding window size [default = 4]'
        # )

        # qc_conf.add_argument(
        #     "--slidingwindow-quality", 
        #     type=int,
        #     default=10,
        #     metavar = "N",
        #     help='trimmomatic sliding window quality score [default = 10]'
        # )

        # qc_conf.add_argument(
        #     '--trimmomatic-extra-args',
        #      type=str,
        #     metavar="='-args'",
        #     help="A string value of extra args for trimmomatic (must be used with = with no spaces (--trimmomatic-extra-args='-arg1 -arg2'))",
        #     default=""
        # )


        # qc_conf.add_argument(
        #     '--skip-QC',
        #      action='store_true',
        #      help="Skipp Fastqc step"
        # )

        # qc_conf.set_defaults(trimmomatic=False)

        aligner_conf = parser.add_argument_group(f'{BLU}Aligner configuration{NC}')

        # aligner_conf.add_argument(
        #     '--aligner', 
        #     help = "Choose aligner [default = bwa]",
        #     choices=["bwa", "bowtie2"], 
        #     metavar = "bwa|bowtie2",
        #     type=str,
        #     default='bwa'
        # )

        aligner_conf.add_argument(
            '--threads-index', 
            help= "Number of threads to use during indexing ref [default = 4]", 
            type=int, 
            default= 4,
            metavar = "N",
        )

        aligner_conf.add_argument(
            '--threads-align', 
            help= "Number of threads to use during sample alignment [default = 4]", 
            type=int, 
            default= 4,
            metavar = "N",
        )

        aligner_conf.add_argument(
            '--aligner-extra-args', 
            help = "Extra arguments for aligner, use it with no spaces and add = ( --aligner-extra-args='-arg1 -arg2' ) ",
            type=str,
            metavar = "'-args'",
            default=""
        )

        aligner_conf.add_argument(
            '--reference-index',
            metavar='path/to/ref', 
            type=os.path.abspath,
            help="path to reference index"
        )

        aligner_conf.add_argument(
            '--index-fasta', 
            help='Index fasta file', 
            action='store_true', 
        )

        aligner_conf.set_defaults(index_fasta=False)

        variant_caller_conf = parser.add_argument_group(f'{BLU}Variant caller configuration{NC}')

        # variant_caller_conf.add_argument(
        #         '--variant-caller', 
        #         help = "Choose variant caller", 
        #         choices=["GATK", "mpileup", "lofreq"], 
        #         type=str, 
        #         default='GATK'
        # )

        # variant_caller_conf.add_argument(
        #         '--bqsr-mem', 
        #         help = "base quality score recalipration memory for each sample", 
        #         choices=["GATK", "mpileup", "lofreq"], 
        #         type=str, 
        #         default='GATK'
        # )

        variant_caller_conf.add_argument(
            '--known-variants',
            metavar='path', 
            type=os.path.abspath, 
            help="path to reference fasta file",
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)

        )

        # variant_caller_conf.add_argument(
        #     '--caller-extra-args', 
        #     help = "Extra arguments for caller, use it with no spaces and add = ( --caller-extra-args='-arg1 -arg2' ) ",
        #     type=str,
        #     metavar = "'-args'",
        #     default=""
        # )

        # variant_caller_conf.add_argument(
        #     '--threads-calling', 
        #     help= "Number of threads to use during variant calling [default = 4]", 
        #     type=int, 
        #     default= 4,
        #     metavar = "N",
        # )

        # Snakemake Options
        snakemake_options = parser.add_argument_group(f'{CYN}Snakemake Options{NC}')

        snakemake_options.add_argument(
            '--dry-run', 
            action='store_true', 
            help="performs snakemake dry run"
        )

        snakemake_options.add_argument(
            '--export-dag', 
            action='store_true', 
            help="performs snakemake dry run and exports DAG"
        )

        snakemake_options.add_argument(
            "--smk-extra-args", 
            metavar="='-args'",
            help="A string value of extra args for snakemake(must be used with = with no spaces (--smk-extra-args='-arg1 -arg2'))",
            default="", type=str
        )

        snakemake_options.add_argument(
            "--parse-snakemake-output", 
            action='store_true', 
            help="prints progress bar instead of snakemake regular output "
        )

        annotation_conf = parser.add_argument_group(f'{BLU}Annotation configuration{NC}')

        annotation_conf.add_argument(
            "--annovar-protocol", 
            default="refGene,avsnp150,clinvar_20221231,cosmic70,dbnsfp31a_interpro,EAS.sites.2015_08,EUR.sites.2015_08,gme,gnomad211_exome,SAS.sites.2015_08", 
            metavar = 'str',
            help=f"Annovar Protocol\n defaults:\nrefGene,avsnp150,clinvar_20221231,cosmic70,dbnsfp31a_interpro,EAS.sites.2015_08,EUR.sites.2015_08,gme,gnomad211_exome,SAS.sites.2015_08"
        )

        annotation_conf.add_argument(
            "--annovar-operation", 
            default="g,f,f,f,f,f,f,f,f,f", 
            metavar = 'str',
            help=f"Annovar Protocol\n defaults:\ng,f,f,f,f,f,f,f,f,f"
        )

        # other options
        other_conf = parser.add_argument_group(f'{CYN}Other{NC}')

        other_conf.add_argument(
            '--continue', 
            action='store_true', 
            help="continue analysis when re-run"
        )

        other_conf.add_argument(
            '--overwrite', 
            action='store_true', 
            help="overwrite output dir if exsits"
        )

        other_conf.add_argument(
            '-n',"--name", 
            default=f"guap_run[{os.environ['start_time']}]", 
            metavar = 'str',
            help=f"Name of files [ default = guap_run[date time] ]"
        )

        other_conf.add_argument(
            '--verbose', 
            action='store_true', 
            help="verbose"
        )

        other_conf.add_argument(
            '--quit', 
            dest='verbose', 
            action='store_false', 
            help="print many output"
        )

        other_conf.add_argument(
            '--print-last-run', 
            action='store_true', 
            help="Prints last run on screen",
        )


    def run(self, args):

        workflow = "WES"
        if args.print_last_run:
            with open(f"{GUAP_DIR}/.last_run.txt", 'r') as last_run:
                lines = last_run.readlines()
            last_command = lines[0]
            print(f"guap {last_command}")
            exit()
        all_args = parse_input_args(args)
        if all_args["generate_confs_only"]:
            print(f"{PRP}{runtime.elapsed()}{NC}")
        else:
            glogger.prnt_fatel(f"{RED}Sorry still under dev :(\n{YEL}Make sure to include: {GRE}--generate-confs-only{NC}")
            snakemake_cmd = smk_cmd(all_args, f"{workflow}")
            try:
                if all_args['export_dag'] is True and all_args['dry_run'] != True:
                    if all_args['continue']:
                        subprocess.run(f"snakemake --snakefile '{GUAP_DIR}/workflows/{workflow}/Snakefile' --configfile '{all_args['working_dir']}/config.yaml' -j {all_args['threads']} --progress {all_args['smk_extra_args']}", shell=True)
                    else:
                        subprocess.run(snakemake_cmd, shell=True)
                        subprocess.run(f"snakemake --snakefile '{GUAP_DIR}/workflows/{workflow}/Snakefile' --configfile '{all_args['working_dir']}/config.yaml' -j {all_args['threads']} --progress {all_args['smk_extra_args']}", shell=True)

                    print(f"{PRP}{runtime.elapsed()}{NC}")
                elif all_args['dry_run']:
                    try:
                        subprocess.run(snakemake_cmd, shell=True)
                    except Exception as E:
                        glogger.prnt_fatel(f"Error in Dry Run:\n{E}")
                    subprocess.run(f"{snakemake_cmd} -q -n --rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'", shell=True)
                    print(f"{PRP}{runtime.elapsed()}{NC}") 

                else:

                    subprocess.run(f"{snakemake_cmd} -q -n --rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'", shell=True)

                    try:
                        if args.parse_snakemake_output:
                            process_snakemake_standard_output(snakemake_cmd, f"{all_args['working_dir']}/output.log")
                        else:
                            subprocess.run(f"{snakemake_cmd}")
                    except Exception as E:
                        glogger.prnt_error(f"Error in snakemake parsing:\n{RED_}{E}{NC}\n trying normal run using subprocess.run()")
                        try:
                            subprocess.run(f"{snakemake_cmd}")
                        except Exception as e:
                            glogger.prnt_fatel(f"Error in snakemake run:\n{RED_}{e}{NC}")


                    print(f"{PRP}{runtime.elapsed()}{NC}") 
            except Exception as E:
                glogger.prnt_fatel(f"Error in snakemake run:\n{RED_}{E}{NC}")
                print(f"{PRP}{runtime.elapsed()}{NC}")  

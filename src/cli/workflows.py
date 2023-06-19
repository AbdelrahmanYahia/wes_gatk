import os
import sys
import argparse
import subprocess
from ..utils.globals import *
from ..utils.parse_input import parse_input_args, process_snakemake_standard_output

def smk_cmd(all_args,workflow):
    smk_cmd = ""

    if all_args['dry_run'] is True:
        smk_cmd += " -n"
        glogger.prnt_warning("performing dry run")
        if all_args['export_dag'] is True:
            smk_cmd += f" --rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'"
            glogger.prnt_warning("exporting dag")
        else:
            pass
    elif all_args['export_dag'] is True and all_args['dry_run'] != True:
        smk_cmd += f"--rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'"

    else:
        smk_cmd = smk_cmd + f"{all_args['smk_extra_args']}"


    snakemake_cmd = f"snakemake --snakefile '{GUAP_DIR}/workflows/{workflow}/Snakefile' --configfile '{all_args['working_dir']}/config.yaml' -j {all_args['threads']} {smk_cmd}"
    return snakemake_cmd

# Create a custom action to preserve quoted strings
class StoreQuotedString(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # If the value is surrounded by quotes, remove them
        if isinstance(values, str) and values.startswith('"') and values.endswith('"'):
            values = values[1:-1]
        setattr(namespace, self.dest, values)

# base class for subcommands (workflows)
class WorkflowCli:
    def __init__(self,subparser):
        self.parser = subparser.add_parser(self.name, help=self.help, usage=self.usage)
        self.add_arguments(self.parser)

    def add_arguments(self, parser):
        pass
    
    def run(self, args, workflow):
        # if args.print_last_run:
        #     with open(f"{GUAP_DIR}/.last_run.txt", 'r') as last_run:
        #         lines = last_run.readlines()
        #     last_command = lines[0]
        #     print(f"guap {last_command}")
        #     exit()
        all_args = parse_input_args(args)
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
                process_snakemake_standard_output(snakemake_cmd, f"{all_args['working_dir']}/output.log")
                print(f"{PRP}{runtime.elapsed()}{NC}") 
        except Exception as E:
            glogger.prnt_fatel(f"Error in snakemake run:\n{RED_}{E}{NC}")
            print(f"{PRP}{runtime.elapsed()}{NC}") 

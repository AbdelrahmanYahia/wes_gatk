import argparse
from ..utils.globals import *
from .workflow_cmd_classes.WES_cmd import WES

# main cli class
class Cli:
    def __init__(self):
        # create argument parser
        self.parser = argparse.ArgumentParser(description="GUAPtoolkit for Genomics analysis",
                                              prog="guap", 
                                              epilog = '''This tool is still under development''',
                                              usage=f"guap WES [--input <dir>] [--output <dir>] [...]")

        # setup workflows parsers
        subparsers = self.parser.add_subparsers(title=f"{NC_}Workflows{NC}", dest='workflow',
                                                description="choose workflow", 
                                                metavar="")
        WESworkflowCli = WES(subparsers)
        
        args = self.parser.parse_args()

        if args.workflow == None:
            WESworkflowCli.run(args)

        elif args.workflow == 'WES':
            WESworkflowCli.run(args)

        else:
            print("unknown arg!")


if __name__ == '__main__':
    cli = Cli()

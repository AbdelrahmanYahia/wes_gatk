#!/bin/env python

import datetime 
now = datetime.datetime.now()
start_time = now.strftime("%Y-%m-%d %H:%M:%S")
import os
os.environ['start_time'] = start_time
import sys

### TODO: should modify this for stability 
sys.setrecursionlimit(10000000)
from src.cli import cli
from src.utils.globals import *

if __name__ == '__main__':
    cli = cli.Cli()
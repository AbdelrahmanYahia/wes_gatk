import os
import psutil
import multiprocessing
from .runtime import RunTime

# cols
RED = "\033[;31;1m"
RED_ = "\033[;31;4m"
GRE = "\033[;32;1m"
YEL = "\033[;33;1m"
BLU = "\033[;34;1m"
PRP = "\033[;35;1m"
CYN = "\033[;36;1m"
BLD = "\033[;37;1m"
GRY = "\033[;30;1m"
NC = "\033[;39;0m"
NC_ = "\033[;39;4m"

start_time = os.environ['start_time']
runtime = RunTime(start_time)

# performance
all_threads = multiprocessing.cpu_count()
all_mem = int( (psutil.virtual_memory().total ) / 1000000000 )

# other vars
GUAP_DIR = (os.path.abspath(__file__)).replace("src/utils/globals.py","")
global_vars = {}

from .Logger import GLogger

glogger = GLogger()

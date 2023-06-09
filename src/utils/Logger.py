import logging
from .globals import *


class CustomFormatter(logging.Formatter):
    ### code from:
    ## https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output

    format = "[%(asctime)s] %(levelname)s: %(message)s"
    datefmtstr = "%Y-%m-%d %H:%M:%S"

    FORMATS = {
        logging.DEBUG: format.replace('%(levelname)s', GRY + '%(levelname)s' + NC),
        logging.INFO: format.replace('%(levelname)s', GRY + '%(levelname)s' + NC),
        logging.WARNING: format.replace('%(levelname)s', YEL + '%(levelname)s' + NC),
        logging.ERROR: format.replace('%(levelname)s', RED + '%(levelname)s' + NC),
        logging.CRITICAL: format.replace('%(levelname)s', RED_ + '%(levelname)s' + NC),
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt=self.datefmtstr)
        return formatter.format(record)
    

class GLogger:
    def __init__(self):
        # create logger attr and set level
        self.logger = logging.getLogger('GUAP_logger')
        self.formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    def create_console_handler(self, verbose=False):
        # create handler
        self.ch = logging.StreamHandler()

        # set level according to user verbose option
        if verbose:
            self.ch.setLevel(logging.DEBUG)
        else:
            self.ch.setLevel(logging.WARNING)

        self.ch.setFormatter(CustomFormatter())

        # add handler
        self.logger.addHandler(self.ch)

    def create_file_handler(self, file):
        # create handler
        self.fh = logging.FileHandler(file)
        self.fh.setLevel(logging.DEBUG)
        self.fh.setFormatter(self.formatter)

        self.logger.addHandler(self.fh)


    def prnt_info(self, str):
        self.logger.info(f"{str}")


    def prnt_warning(self, str):
        self.logger.warning(f"{str}")


    def prnt_error(self, str):
        self.logger.error(f"{str}")


    def prnt_fatel(self, str):
        self.logger.fatal(f"{str}")
        print(f"{PRP}{runtime.elapsed()}{NC}")
        exit(1)


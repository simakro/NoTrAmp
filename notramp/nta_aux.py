from os import path
from collections import Counter
import logging
import logging.config


log_file_path = path.join(path.dirname(__file__),  "resources", "logging.conf")
logging.config.fileConfig(log_file_path , disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def fastq_autoscan(read_file):
    """Scanning readfile to determine filetype"""
    logger.debug("Scanning readfile to determine filetype")
    with open(read_file, "r", encoding="utf-8") as rfl:
        line_ct = 0
        fastq = False
        while line_ct < 100:
            for line in rfl:
                line_ct += 1
                if line.startswith("@"):
                    fastq = True
                    break
    return fastq


class BedColumnError(Exception):
    """raised when number of columns in bed file is not matching"""

    def __init__(self):
        self.message = "Not all lines in primer bed have same number of" \
                    " columns. File may be corrupt, please check integrity" \
                    " of this file. Columns may not be empty. Take care that" \
                    " headers or comments begin with any of the following:" \
                    " 'track', 'browser', '#'."
        super().__init__(self.message)


def bed_scan(bed_file):
    """Scan bed-file to get col-num"""
    logger.debug("Scanning bed-file")
    col_num = 6
    headers = []
    with open(bed_file, "r", encoding="utf-8") as bfl:
        line_ct = 0
        line_lens = []
        for line in bfl:
            line_ct += 1
            lst = line.strip()
            if lst.startswith("#"):
                headers.append("#")
            elif lst.startswith("browser"):
                headers.append("browser")
            elif lst.startswith("track"):
                headers.append("track")
            else:
                if len(lst) > 0:
                    lsp = lst.split("\t")
                    line_lens.append(len(lsp))
        len_set = set(line_lens)
        if len(len_set) == 1:
            col_num = list(len_set)[0]
        else:
            raise BedColumnError
        if len(headers) == 0:
            headers = False
    return col_num, headers


def bed_chk_header(line, headers):
    head = False
    for string in headers:
        if line.startswith(string):
            head = True
    return head

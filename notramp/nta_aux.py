from os import path
import logging
import logging.config
import psutil


logger = logging.getLogger(__name__)

class BedColumnError(Exception):
    """raised when number of columns in bed file is not matching"""

    def __init__(self):
        self.message = "Not all lines in primer bed have same number of" \
                    " columns. File may be corrupt, please check integrity" \
                    " of this file. Columns may not be empty. Take care that" \
                    " headers or comments begin with any of the following:" \
                    " 'track', 'browser', '#'."
        super().__init__(self.message)


# class BedColumnError(Exception):
#     """raised when number of columns in bed file is not matching"""

#     def __init__(self):
#         self.message = "Not all lines in primer bed have same number of" \
#                     " columns. File may be corrupt, please check integrity" \
#                     " of this file. Columns may not be empty. Take care that" \
#                     " headers or comments begin with any of the following:" \
#                     " 'track', 'browser', '#'."
#         super().__init__(self.message)


def fastq_autoscan(read_file):
    """Scanning readfile to determine filetype"""
    logger.info("Scanning readfile to determine filetype")
    rf_size = path.getsize(read_file)
    if rf_size == 0:
        logger.warning("Read file appears to be completely empty!")
    elif rf_size < 100:
        logger.warning("Read file appears to be almost empty!")
    else:
        logger.info("Read file size is %s", rf_size)
    with open(read_file, "r", encoding="utf-8") as rfl:
        line_ct = 0
        fastq = False
        # while line_ct < 100:
        for line in rfl:
            line_ct += 1
            if line.startswith("@"):
                fastq = True
                break
            if line_ct > 100:
                break
    return fastq


def chk_mem_fit(kw):
    """Check for available memory"""
    logger.info("Checking for available memory")
    fastq = fastq_autoscan(kw["reads"])
    rf_size = path.getsize(kw["reads"])
    if fastq:
        load_size = rf_size if kw["fq_out"] else rf_size/2
    else:
        load_size = rf_size
    avail_mem = psutil.virtual_memory()[1]

    if load_size > 0:
        if avail_mem / load_size > 4:
            return True
        else:
            return False
    else:
        return True


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
    """Check if line conforms to header type in headers list"""
    head = False
    for string in headers:
        if line.startswith(string):
            head = True
    return head

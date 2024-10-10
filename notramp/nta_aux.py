from os import path
import logging
import logging.config
import traceback as trace

try:
    import psutil
except ModuleNotFoundError:
    print(
        "WARNING: Module psutil could not be be found. Running NoTrAmp without "
        "psutil can increase runtime significantly!"
        )


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


class PrimerSchemeError(Exception):
    """raised Naming scheme does not match primer name in bed file"""

    def __init__(self, obj_type):
        self.message = f"Error during {obj_type}-objects generation. Naming scheme does not" \
                f" match primer names in bed file. Please ensure that you are us" \
                f"ing the right combination of naming scheme and primer bed-file" \
                f" and that your primers consistently follow this scheme."
        super().__init__(self.message)


def exception_handler(exception_type, exception, tb):
    logger.error(f"{exception_type.__name__}, {exception}")
    tb_str = "".join(trace.format_list(trace.extract_tb(tb)))
    logger.info(tb_str)


def fastq_autoscan(read_file):
    """Scanning readfile to determine filetype"""
    logger.info("Scanning readfile to determine filetype")
    rf_size = path.getsize(read_file)
    if rf_size == 0:
        logger.warning(f"File {read_file} appears to be completely empty!")
    elif rf_size < 100:
        logger.warning(f"File {read_file} appears to be almost empty!")
    else:
        logger.info("Read file size is %s bytes", rf_size)
    with open(read_file, "r", encoding="utf-8") as rfl:
        line_ct = 0
        fastq = False
        for line in rfl:
            line_ct += 1
            if line.startswith("@"):
                fastq = True
                break
            if line_ct > 100:
                break
    return fastq


# def fastq_gen():
#     while True:
#         yield "header"
#         yield "seq"
#         yield "plus"
#         yield "qstring"


# def fastq_chk(fastq_file):
#     """Check integrity of nucleotide fastq file"""
#     logger.info("Checking integrity of nucleotide fastq")
#     ln_ct = 0
#     ln_type = fastq_gen()
#     bases = ["A", "T", "C", "G", "N"]
#     errors = 0
#     good_reads = 0

#     with open(fastq_file, "r", encoding="utf-8") as fq:
#         curr_read = {}
#         for line in fq:
#             ln_ct += 1
#             supposed_ln = next(ln_type)
#             if supposed_ln == "qstring":
#                 curr_read[supposed_ln] = line
#                 tests = [
#                     curr_read["header"].startswith("@"),
#                     curr_read["seq"][0].upper() in bases,
#                     curr_read["plus"] == "+\n",
#                     len(curr_read["qstring"]) == len(curr_read["seq"]),
#                 ]
#                 if all(tests):
#                     good_reads += 1
#                 else:
#                     print(f"Format error in fastq in lines {ln_ct-3}-{ln_ct}")
#                     errors += 1 
#                 curr_read = {}
#             else:
#                 curr_read[supposed_ln] = line
#     logger.info(f"Found {errors} format errors in {fastq_file}")
#     logger.info(f"Found {good_reads} good reads in {fastq_file}")
#     if errors > 0:
#         pass


def chk_fqblock_integrity(title, seq, qual):
    """Check integrity of all entries in one fastq read-block"""
    bases = "ATGCN"
    if all([title, seq, qual]):
        tests = [
                    title.startswith("@"),
                    seq[0].upper() in bases,
                    # all([c.upper() in bases for c in set(seq)]),
                    len(qual) == len(seq),
                    ]
        return all([tests])
    elif any([title, seq, qual]):
        return False
    else:
        return None


def fastq_iterator(handle):
    """Fast fastq iterator with integrity check"""
    while True:
        title = handle.readline().strip()
        seq = handle.readline().strip()
        handle.readline()  # Skip '+' line
        qual = handle.readline().strip()
        chk = chk_fqblock_integrity(title, seq, qual)
        if chk:
            yield title, seq, qual
        elif chk == False:
            logger.warning("Encountered irregularity in fastq file")
        else: # if chk == None:
            logger.info("Reached end of file")
            break


def chk_mem_fit(kw):
    """Check for available memory"""
    logger.info("Checking for available memory")
    fastq = kw["fastq_in"]
    rf_size = path.getsize(kw["reads"])
    if fastq:
        load_size = rf_size if kw["fq_out"] else rf_size/2
    else:
        load_size = rf_size
    avail_mem = psutil.virtual_memory()[1]

    if load_size > 0:
        if avail_mem - load_size > avail_mem*0.2:
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

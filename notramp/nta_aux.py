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
        self.message = "Not all lines in primer bed have same number of colum" \
            "ns. File may be corrupt, please check integrity of this file. Co" \
            "lumns may not be empty. Take care that headers or comments begin" \
            " with any of the following: 'track', 'browser', '#'."
        super().__init__(self.message)


class PrimerSchemeError(Exception):
    """raised when naming scheme does not match primer name in bed file"""

    def __init__(self, obj_type):
        self.message = f"Error during {obj_type}-objects generation. Naming " \
            f"scheme does not match primer names in bed file. Please ensure " \
            f"that you are using the right combination of naming scheme and " \
            f"primer bed-file and that your primers consistently follow this" \
            f" scheme."
        super().__init__(self.message)


class LoadReadsError(Exception):
    """raised upon failure of loading reads into memory"""

    def __init__(self):
        self.message = "An error occured while trying to load reads. Aborting" \
            "NoTrAmp. This can be caused by corrupt input files. Please check" \
            " log for details."
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


def analyze_read_file(file_path, fastq=False):
    """Parse entire fasta/q file and report all irregularities and errors"""
    logger.info("Analyzing readfile to identify potential errors.")
    
    def test_header(line, kwargs):
        correct_ind = line.startswith("@" if kwargs["fastq"] else ">")
        return correct_ind
    
    def test_seq(line, kwargs):
        bases = "ATGCN"
        header_inds = [">", "@"]
        found = set(line)
        checked = [char in bases for char in found]
        if all(checked):
            return True
        else:
            unexpect_chars = [char for char in found if char not in bases]
            if line[0] not in header_inds:
                logger.warning(
                    f"Found unexpected chars {unexpect_chars} in seq-line "
                    f"{kwargs['l_ct']}"
                    )
            return False
    
    def test_plus(line, kwargs):
        kwargs = kwargs
        return line == "+"
    
    def test_qual(line, kwargs):
        return len(line) == kwargs["len_seq"]

    test_funcs = {
        "header": test_header,
        "seq": test_seq,
        "+": test_plus,
        "qual": test_qual
    }    
    err_lines = {
        "empty": [],
        "incorrect": []
    }
    expect_fq = {
        "header": "seq",
        "seq": "+", 
        "+": "qual",
        "qual": "header"
        }
    expect_fa = {
        "header": "seq",
        "seq": "header"
        }

    l_ct = 0
    exp_curr = "header"
    len_last_seq = 0
    line_dct = expect_fq if fastq else expect_fa
    with open(file_path, "r", encoding="utf-8") as f:
        for l in f:
            l = l.strip()
            l_ct += 1
            if len(l) == 0:
                err_lines["empty"].append(l_ct)
            if exp_curr == "seq":
                len_last_seq = len(l)
            kwargs = {"len_seq": len_last_seq, "fastq": fastq, "l_ct": l_ct}
            chk = test_funcs[exp_curr](l, kwargs)
            if not chk:
                err_lines["incorrect"].append(l_ct)
            exp_curr = line_dct[exp_curr]
    logger.warning(f"Detected {len(err_lines['incorrect'])} irregular lines in the input sequence file.")
    logger.warning(f"Detected {len(err_lines['empty'])} empty lines in the input sequence file.")
    logger.info(f"Analyzed file has {l_ct} lines.") 


def chk_fablock_integrity(title, seq):
    """Check integrity of all entries in one fasta read-block"""
    bases = "ATGCN"
    if all([title, seq]):
        tests = [
                    title.startswith(">"),
                    seq[0].upper() in bases,
                    ]
        return all([tests])
    elif any([title, seq]):
        return False, False
    else:
        return None, None


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
        elif chk is False:
            logger.warning("Encountered irregularity in fastq file")
            yield chk
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

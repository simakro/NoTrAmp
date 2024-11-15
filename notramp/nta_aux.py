from os import path, getcwd
import sys
import logging
import logging.config
import traceback as trace
import subprocess as sp
from time import strftime, localtime

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
            "lumns may not be empty. Take care that headers or comments (if p" \
            "resent) begin with any of the following: 'track', 'browser', '#'."
        super().__init__(self.message)


class TargetReferenceError(Exception):
    """raised if more than one sequence entry is found in the reference fasta"""

    def __init__(self):
        self.message = "More than one sequence entry was found in the referen" \
            "ce fasta. Multiple reference sequences are currently not support" \
            "ed in NoTrAmp and could lead to unexpected results. Amplicon til" \
            "ing schemes typically aim at a single longer contiguous sequence" \
            ". Please run tiling panels for independent target sequences sepa" \
            "retly. If there is interest to combine multiple targets that can" \
            " not be represented in a single fasta entry, please contact the " \
            "developer/maintainer to request this feature."
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


class AmpliconGenerationError(Exception):
    """raised if amplicon generation fails or yields unexpected results"""

    def __init__(self, expect_amps, gen_amps):
        self.message = f"Amplicon generation yielded unexpected results. Expe" \
            f"cted {expect_amps} amps. Amplicons generated: {gen_amps}. Pleas" \
            f"e check if you are using the right Primer(.bed) and Primer-sche" \
            f"me(.json) files. Aborting NoTrAmp."
        super().__init__(self.message)


class SeqReadFileAnalyzer():
    """Class for analyis of and error detection in fasta/q files"""

    def __init__(self, file_path, fastq=False):
        self.file_path = file_path
        self.fastq = fastq
        self.test_funcs = {
                            "header": self.test_header,
                            "seq": self.test_seq,
                            "+": self.test_plus,
                            "qual": self.test_qual
                        }
        self.err_lines = {
                            "empty": [],
                            "header": [],
                            "seq": [],
                            "+": [],
                            "qual": []
                        }
        self.expect_fq = {
                            "header": "seq",
                            "seq": "+",
                            "+": "qual",
                            "qual": "header"
                            }
        self.expect_fa = {
                            "header": "seq",
                            "seq": "header"
                            }
        self.l_ct = 0
        self.last_line = ""
        self.len_last_seq = 0
        self.read_names = []
        self.messages = {"infos": [], "warnings": [], "errors": []}
        self.analyze_read_file()

    def test_header(self, line):
        correct_ind = line.startswith("@" if self.fastq else ">")
        if correct_ind:
            self.read_names.append(line.split(" ")[0])
        return correct_ind

    def test_seq(self, line):
        bases = "ATGCN"
        header_inds = [">", "@"]
        found = set(line)
        checked = [char in bases for char in found]
        if all(checked):
            return True
        else:
            unexpect_chars = [char for char in found if char not in bases]
            if line[0] not in header_inds:
                logger.info(
                    f"Found unexpected chars {unexpect_chars} in seq-line "
                    f"{self.l_ct}"
                    )
            return False

    def test_plus(self, line):
        return line == "+"

    def test_qual(self, line):
        chk = len(line) == self.len_last_seq
        if not chk:
            if self.last_line == "+":
                logger.info(
                    f"Quality string in line {self.l_ct} does not have the same"
                    f" length as previous seq-line."
                    )
        return chk

    def prepare_report(self):
        unique_names = set(self.read_names)
        if not len(unique_names) == len(self.read_names):
            duplicate_ct = len(self.read_names) - len(unique_names)
            self.messages["errors"].append(
                f"File contains at least {duplicate_ct} duplicate read-name(s)"
            )
        if not self.fastq:
            self.err_lines.pop("qual")
            self.err_lines.pop("+")
        for l_type in self.err_lines:
            level = "warnings" if len(l_type) > 0 else "infos"
            self.messages[level].append(
                f"Found {len(self.err_lines[l_type])} erroneous {l_type} lines "
                f"in the input sequence file."
            )
        total_err = sum([len(self.err_lines[k]) for k in self.err_lines])
        self.messages["warnings"].append(
                f"Detected a total of {total_err} errors in the input file."
            )
        self.messages["infos"].append(f"Analyzed file has {self.l_ct} lines.")

    def analyze_read_file(self):
        logger.info("Analyzing read file!")
        exp_curr = "header"
        line_dct = self.expect_fq if self.fastq else self.expect_fa
        with open(self.file_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                self.l_ct += 1
                if len(line) == 0:
                    self.err_lines["empty"].append(self.l_ct)
                if exp_curr == "seq":
                    self.len_last_seq = len(line)
                chk = self.test_funcs[exp_curr](line)
                if not chk:
                    self.err_lines[exp_curr].append(self.l_ct)
                exp_curr = line_dct[exp_curr]
                self.last_line = line
        self.prepare_report()

    def report_results(self):
        reporter = {
            "infos": logger.info,
            "warnings": logger.warning,
            "errors": logger.error,
        }
        for key in self.messages:
            for msg in self.messages[key]:
                reporter[key](msg)


# class NonPyDepChecker():
#     """Class for self testing of NoTrAmp"""
#     def __init__(self):
#         self.deps = ["minimap2"]
#         self.cmds = {"minimap2": "--version"}
#         self.vers = {}
    
    # def test_dep(self, dep):
    #     try:
    #         dep_ver = sp.run([self.dep, self.cmds[self.dep]], stdout=sp.PIPE)
    #         dep_ver = dep_ver.stdout.replace(
    #                                         b',', b';'
    #                                     ).replace(
    #                                             b'\n', b''
    #                                         ).replace(
    #                                                 b'\r', b''
    #                                             ).decode("utf-8")
    #         logger.info(f"{dep} installed and responsive")
    #         logger.info(dep_ver)
    #         # if dep=="minimap2":
    #         #     chk_minimap_version(dep_ver, args)
    #     except sp.CalledProcessError as ce:
    #         logger.error(f"{dep} check CalledProcessError: {ce}")
    #         logger.error(
    #             f"The {dep} executable is installed but appears not"
    #             f" to be working properly. Check logs for more info."
    #             )
    #         logger.error(dep_unavail_msg)
    #         # args.__dict__[dep.split(".")[0]] = False
    #     except FileNotFoundError as fe:
    #         logger.error(f"Dependency check {dep}: FileNotFoundError: {fe}")
    #         pf = sys.platform



class SelfTest():
    """Class for self testing of NoTrAmp"""

    def __init__(self):
        self.pkg_path = path.split(path.abspath(__file__))[0]
        self.python = sys.executable
        self.py_ok = self._chk_pyver()
        self.main = path.join(self.pkg_path, "notramp_main.py")
        self.mode = "-a"
        self.ref = path.join(
            self.pkg_path, "test_self",
            "MN908947.3_SARS-CoV2_Wuhan-reference.fasta"
            )
        self.reads = path.join(
            self.pkg_path, "test_self",
            "ARTIC_v5.3.2_Amplicons_Wuhan1_reference.fasta"
            )
        self.prim_scheme = path.join(
            self.pkg_path, "test_self",
            "artic_nCoV_scheme_v5.3.2.json"
            )
        self.prim_bed = path.join(
            self.pkg_path, "test_self",
            "ARTIC_SARS-CoV-2_v5.3.2_400.primer.bed"
            )
        self.outdir = path.join(
            getcwd(),
            f"{strftime('%Y%m%d_%H-%M-%S', localtime())}_notramp-selftest"
            )
        self.run_selftest()
        self.log = path.join(self.outdir, "notramp.log")
        self.print_log()

    def run_selftest(self):
        if self.py_ok:
            try:
                stdo_log = sp.run(
                        [
                            self.python,
                            self.main,
                            self.mode,
                            "--figures",
                            "-p", self.prim_bed,
                            "-g", self.ref,
                            "-r", self.reads,
                            "-n", self.prim_scheme,
                            "-o", self.outdir
                        ],
                        check=True,
                        capture_output=True
                        ).stdout.decode()
                return stdo_log
            except sp.CalledProcessError as err:
                log_sp_error(err, "SelfTest NoTrAmp subprocess failed.")
        else:
            logger.error(
                "The python version you are using is older than the minimum re"
                "quirement of 3.6. Please use a newer version or set up a ded"
                "icated environment with a more recent version of python."
                )

    def print_log(self):
        with open(self.log, "r") as log:
            for line in log:
                print(line.strip())

    def _chk_pyver(self):
        version_str = sys.version
        version = version_str.split(" ")[0].split(".")
        major, minor = version[0], version[1]
        checks = [False, False]
        if int(major) >= 3:
            checks[0] = True
        if int(minor) >= 6:
            checks[1] = True
        if all(checks):
            return version_str
        else:
            return False

    def _get_main(self):
        return path.join(self.pkg_path, "notramp_main.py")


def exception_handler(exception_type, exception, tb):
    logger.error(f"{exception_type.__name__}, {exception}")
    tb_str = "".join(trace.format_list(trace.extract_tb(tb)))
    logger.info(tb_str)


def log_sp_error(error, message):
    """custom logging of subprocess errors"""
    logger.error(message)
    print("Check notramp.log for error details")
    logger.error(error.stdout.decode('utf-8'))
    logger.error("Exiting notramp")
    sys.exit()


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


def chk_fablock_integrity(title, seq):
    """Check integrity of all entries in one fasta read-block"""
    bases = "ATGCN"
    if all([title, seq]):
        tests = [
                    title.startswith(">"),
                    seq[0].upper() in bases
                    ]
        return tests
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
        handle.readline()  # Skip over plus-line
        qual = handle.readline().strip()
        chk = chk_fqblock_integrity(title, seq, qual)
        if chk:
            yield title, seq, qual
        elif chk is False:
            logger.warning("Encountered irregularity in fastq file")
            yield chk
        else:  # if chk == None:
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

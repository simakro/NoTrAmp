# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.


import argparse
import logging
import logging.config
import subprocess
from time import perf_counter
import random
import json
import sys
import os
from os import path
from collections import defaultdict

if __name__ == "__main__":
    import nta_aux as aux
    import amp_cov as amp_cov
    import map_trim as map_trim
    from defaults import DEFAULTS
    from version import __version__
else:
    import notramp.nta_aux as aux
    import notramp.amp_cov as amp_cov
    import notramp.map_trim as map_trim
    from notramp.defaults import DEFAULTS
    from notramp.version import __version__


logger = logging.getLogger(__name__)
sys.excepthook = aux.exception_handler


def get_arguments():
    """get arguments from the command line"""

    usg_msg = "\nnotramp (-a | -c | -t) -p PRIMERS -r READS -g REFERENCE"\
        " [optional arguments]"

    parser = argparse.ArgumentParser(
        prog="notramp", usage=usg_msg,
        description="NoTrAmp is a Tool for read-depth normalization and"
        " trimming of amplicon reads generated with long read technologies"
        " (ONT/PacBio). It can be used in amplicon-tiling approaches to cap"
        " coverage of each amplicon and to trim amplicons to their appropriate"
        " length removing barcodes, adpaters and primers (if desired) in a"
        " single clipping step.",
        add_help=False
        )
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument(
        "-p", "--primers", required=True,
        help="Path to primer bed-file (primer-names must adhere to a consistent"
        " naming scheme see readme)"
        )
    required_args.add_argument(
        "-r", "--reads", required=True,
        help="Path to sequencing reads fasta"
        )
    required_args.add_argument(
        "-g", "--reference", required=True,
        help="Path to reference (genome)"
        )

    read_input = required_args.add_mutually_exclusive_group(required=True)
    read_input.add_argument(
        "-a", "--all",
        help="Perform read depth normalization by coverage-capping/downsampling"
        " first, then clip the normalized reads. (mut.excl. with -c, -t)",
        default=False,  action='store_true'
        )
    read_input.add_argument(
        "-c", "--cov",
        help="Perform only read-depth normalization/downsampling."
        " (mut.excl. with -a, -t)",
        default=False, action='store_true'
        )
    read_input.add_argument(
        "-t", "--trim",
        help="Perform only trimming to amplicon length (excluding primers"
        " by default; to include primers set --incl_prim flag)."
        " (mut.excl. with -a, -c)",
        default=False, action='store_true'
        )

    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument(
        "-h", "--help",
        action="help",
        help="Print help message and exit"
        )
    optional_args.add_argument(
        "-o", dest='out_dir',
        default=False,
        help="Optionally specify a directory for saving of outfiles. If this "
        "argument is not given, out-files will be saved in the directory where"
        " the input reads are located. [default=False]"
        )
    optional_args.add_argument(
        "-m", dest='max_cov',
        default=DEFAULTS["max_cov"], type=int,
        help=f"Provide threshold for maximum read-depth per amplicon as integer"
        f" value. [default={DEFAULTS['max_cov']}]"
        )
    optional_args.add_argument(
        "--incl_prim", dest='incl_prim',
        default=False, action='store_true',
        help="Set this flag if you want to include the primer sequences in the "
        "trimmed reads. By default primers are removed together with all "
        "overhanging sequences like barcodes and adapters."
        )
    optional_args.add_argument(
        "-s", dest='seq_tec',
        default=DEFAULTS['seq_tec'],
        help=f"Specify long-read sequencing technology (ont/pb)."
        f" [default={DEFAULTS['seq_tec']}]"
        )
    optional_args.add_argument(
        "-n", dest='name_scheme',
        default=DEFAULTS["name_scheme"],
        help=f"Provide path to json-file containing a naming scheme which is "
        f"consistently used for all primers.[default={DEFAULTS['name_scheme']}]"
        )
    optional_args.add_argument(
        "--set_min_len", dest='set_min_len',
        default=False, type=bool,
        help="Set a minimum required length for alignments of reads to "
        "amplicon. If this is not set the min_len will be 0.5*average_amp_len."
        " If amplicon sizes are relatively homogenous this parameter is not"
        " required [default=False]"
        )
    optional_args.add_argument(
        "--set_max_len", dest='set_max_len',
        default=False, type=bool,
        help="Set a maximum allowed length for alignments of reads to"
        " amplicon. If this is not set the max_len will be 1.5*average_amp_len."
        " If amplicon sizes are relatively homogenous this parameter is not "
        " required [default=False]"
        )
    optional_args.add_argument(
        "--set_margins", dest='margins',
        default=DEFAULTS["margins"], type=int,
        help=f"Set length of tolerance margins for sorting of mappings to "
        f" amplicons. [default={DEFAULTS['margins']}]"
        )
    optional_args.add_argument(
        "--figures", dest='figures',
        default=False, action="store", const=20, nargs='?',
        help="Set to generate figures of input and output read_counts. Availabl"
        "e for --all and --cov modes. You can optionally provide a value to dra"
        "w a red helper line in the output read plot, showing a threshold, e.g."
        " min. required reads. [default=False; default_threshold=20]"
        )
    optional_args.add_argument(
        "--fastq", dest='fq_out',
        default=False, action="store_true",
        help="Set this flag to request output in fastq format. By default"
        " output is in fasta format. Has no effect if input file is fasta."
        )
    optional_args.add_argument(
        "--split", dest='split_out',
        default=False, action="store_true",
        help="Set this flag to request output of capped, untrimmed reads split"
        " to amplicon specific files (can be a lot)."
        )
    optional_args.add_argument(
        "-v", "--version", action="version",
        version=f"notramp {__version__}",
        help="Print version and exit"
        )

    args = parser.parse_args()
    return args


class Amp:
    """class for generation of amplicon objects, able to extract, store, \
    handle and process data of attached Primer, Mapping and Read objects"""

    def __init__(self, amp_no, primers):
        self.name = self._set_name(amp_no)
        self.primers = primers
        self.start = min([primer.start for primer in primers])
        self.end = max([primer.end for primer in primers])
        self.mappings = {}
        self.reads_dct = {}
        self.selected = []

    def _set_name(self, amp_no):
        try:
            return int(amp_no)
        except ValueError as ve:
            tb = ve.__traceback__
            raise aux.PrimerSchemeError("Amplicon").with_traceback(tb)

    def fwp_boundary(self):
        """get inner boundary of forward primer"""
        return max(
            [p for p in self.primers if p.pos == p.scheme["fw_indicator"]],
            key=lambda x: x.end
        ).end

    def revp_boundary(self):
        """get inner boundary of reverse primer"""
        return min(
            [p for p in self.primers if p.pos == p.scheme["rev_indicator"]],
            key=lambda x: x.start
        ).start

    def get_max_len(self):
        """get maximum length of amplicon"""
        return self.end - self.start

    def random_sample(self, cutoff):
        """select a random subsample from all reads"""
        if len(self.reads_dct) > cutoff:
            read_names = list(self.reads_dct.keys())
            self.selected = random.sample(read_names, k=cutoff)
        else:
            self.selected = list(self.reads_dct.keys())

    def trim_clip_all(self, incl_prim):
        """trim all attached reads"""
        clip_ct_left = 0
        clip_ct_right = 0
        for rname, read in self.reads_dct.items():
            try:
                mapping = self.mappings[rname]
                if incl_prim:
                    left, right = read.trim_to_amp(
                        self.start, self.end, mapping
                    )
                else:
                    left, right = read.clip_primers(
                        self.fwp_boundary(), self.revp_boundary(), mapping
                    )
                clip_ct_left += left
                clip_ct_right += right
            except KeyError as err:
                logger.exception(err)
        return clip_ct_left, clip_ct_right


class Read:
    """class for generation of read objects, and clipping/trimming"""

    def __init__(self, header, seq, fastq=False):
        self.header = header
        self.fastq = fastq
        self.name = self.header.split("@")[1] if self.fastq else self.header.split(">")[1]
        self.seq = seq
        self.plus = False
        self.qstr = False

    @staticmethod
    def non_neg(num):
        """return 0 if result < 0"""
        return num if num >= 0 else 0

    def trim_to_amp(self, amp_start, amp_end, mapping):
        """trim to amplicon boundaries including primers"""
        qlen, samestrand = mapping.qlen, mapping.samestrand
        qstart, qend = mapping.qstart, mapping.qend
        tstart, tend = mapping.tstart, mapping.tend

        if samestrand:
            clip_left = 0
            ldiff = abs(amp_start - tstart)
            if tstart >= amp_start:
                clip_left = Read.non_neg(qstart - ldiff)
            else:
                clip_left = qstart + ldiff

            clip_right = qlen
            rdiff = abs(tend - amp_end)
            if tend <= amp_end:
                clip_right = qend + rdiff
            else:
                clip_right = Read.non_neg(qend - rdiff)
        else:
            clip_left = 0
            ldiff = abs(tend - amp_end)
            if tend <= amp_end:
                clip_left = qstart - ldiff
            else:
                clip_left = qstart + ldiff

            clip_right = qlen
            rdiff = abs(amp_start - tstart)
            if tstart >= amp_start:
                clip_right = qend + rdiff
            else:
                clip_right = qend - rdiff

        self.seq = self.seq[clip_left:clip_right]
        if self.qstr:
            self.qstr = self.qstr[clip_left:clip_right]
        return clip_left, qlen - clip_right

    def clip_primer_plus_strand(
            self, tstart, qstart, qlen, tend, qend, fwp_boundary, revp_boundary
            ):
        clip_left = 0
        if tstart <= fwp_boundary:
            add_clip = fwp_boundary - tstart
            clip_left = qstart + add_clip
        else:
            ldiff = tstart - fwp_boundary
            if qstart >= ldiff:
                clip_left = qstart - ldiff

        clip_right = qlen
        if tend >= revp_boundary:
            sub_clip = tend - revp_boundary
            clip_right = qend - sub_clip
        else:
            rdiff = revp_boundary - tend
            if qlen - qend >= rdiff:
                clip_right = qend + rdiff
        return clip_left, clip_right

    def clip_primer_minus_strand(
            self, tstart, qstart, qlen, tend, qend, fwp_boundary, revp_boundary
            ):
        clip_right = qlen
        if tstart <= fwp_boundary:
            add_clip = fwp_boundary - tstart
            clip_right = qend - add_clip
        else:
            rdiff = tstart - fwp_boundary
            if qlen - qend >= rdiff:
                clip_right = qend + rdiff

        clip_left = 0
        if tend >= revp_boundary:
            sub_clip = tend - revp_boundary
            clip_left = qstart + sub_clip
        else:
            ldiff = revp_boundary - tend
            if qstart >= ldiff:
                clip_left = qstart - ldiff
        return clip_left, clip_right

    def clip_primers(self, fwp_boundary, revp_boundary, mapping):
        """trim to amplicon boundaries excluding primers"""
        qlen, samestrand = mapping.qlen, mapping.samestrand
        qstart, qend = mapping.qstart, mapping.qend
        tstart, tend = mapping.tstart, mapping.tend

        if samestrand:
            clip_left, clip_right = self.clip_primer_plus_strand(
                tstart, qstart, qlen, tend, qend, fwp_boundary, revp_boundary
                )
        else:
            clip_left, clip_right = self.clip_primer_minus_strand(
                tstart, qstart, qlen, tend, qend, fwp_boundary, revp_boundary
                )

        self.seq = self.seq[clip_left:clip_right]
        if self.qstr:
            self.qstr = self.qstr[clip_left:clip_right]
        return clip_left, qlen - clip_right


class Mapping:
    """class for storage and handling of mapping information"""

    def __init__(
        self,
        qname,
        qlen,
        qstart,
        qend,
        samestrand,
        tname,
        tlen,
        tstart,
        tend,
        matches,
        total_bp,
        qual,
        kwattr,
    ):
        self.qname = qname
        self.qlen = int(qlen)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.samestrand = self.eval_strand(samestrand)
        self.tname = tname
        self.tlen = int(tlen)
        self.tstart = int(tstart)
        self.tend = int(tend)
        self.matches = int(matches)
        self.total_bp = int(total_bp)
        self.qual = int(qual)
        self.kwattr = kwattr
        self.gen_kw_attr()

    @staticmethod
    def eval_strand(strand_info):
        """evaluate strand info"""
        return True if strand_info == "+" else False

    def gen_kw_attr(self):
        """generate class attributes from key-worded entries in mapping output"""
        kwattr_dict = {kw.split(":")[0]: kw.split(":")[-1] for kw in self.kwattr}
        for key in kwattr_dict:
            self.__dict__[key] = kwattr_dict[key]


class Primer:
    """class for storage and handling of primer information"""

    def __init__(self, split_bed_line, naming_scheme):
        self.start = None
        self.end = None
        self.name = None
        self.add_info = None
        self.type = "primary"
        self.scheme = naming_scheme
        self.unexpected_vals = defaultdict(list)
        self.get_col_infos(split_bed_line)
        self.get_name_infos()

    def get_col_infos(self, split_bed_line):
        """generate class attributes from bed columns"""
        opt_cols = {
            "score": str,
            "strand": str,
            "thickStart": int,
            "thickEnd": int,
            "itemRgb": str,
            "blockCount": int,
            "blockSizes": str,
            "blockStarts": str,
        }
        chrom, start, end, name = split_bed_line[:4]
        try:
            self.start = int(start)
            self.end = int(end)
        except ValueError as e:
            logger.error(
                "Faulty primer bed file. Both start and end must be of type int"
                )
            logger.error(e)
        self.name = name
        add_info = {"chrom": chrom}

        add_cols = split_bed_line[4:]
        zip_col = zip(list(opt_cols.keys())[:len(add_cols)], add_cols)
        for key, val in zip_col:
            try:
                add_info[key] = opt_cols[key](val)
            except ValueError:
                self.unexpected_vals[(key, opt_cols[key])].append(val)

    def get_name_infos(self):
        """extract information from primer name"""
        try:
            lsp = self.name.split(self.scheme["sep"])
            self.__dict__["amp_no"] = lsp[self.scheme["amp_num"]]
            self.__dict__["pos"] = lsp[self.scheme["position"]]
            if "alt" in self.scheme:
                if len(lsp) == self.scheme["max_len"]:
                    self.type = "alt"
                    self.alt_name = lsp[self.scheme["alt"]]
        except Exception as e:
            tb = e.__traceback__
            raise aux.PrimerSchemeError("Primer").with_traceback(tb)


def log_sp_error(error, message):
    """custom logging of subprocess errors"""
    logger.error(message)
    print("Check notramp.log for error details")
    logger.error(error.stdout.decode('utf-8'))
    logger.error("Exiting notramp")
    sys.exit()


def create_primer_objs(primer_bed, name_scheme):
    """generate primer objects"""
    logger.info("generating primer objects")
    with open(name_scheme, "r", encoding="utf-8") as scheme:
        naming_scheme = json.loads(scheme.read())

    col_num, headers = aux.bed_scan(primer_bed)
    with open(primer_bed, "r", encoding="utf-8") as bed:
        primers = []
        for line in bed:
            if len(line.strip()) == 0:
                pass
            else:
                lst = line.strip().split("\t")
                if len(lst) == col_num:
                    prim = Primer(lst, naming_scheme)
                    primers.append(prim)
                else:
                    if headers:
                        head = aux.bed_chk_header(line, headers)
                        if head:
                            pass
                        else:
                            raise aux.BedColumnError
                    else:
                        raise aux.BedColumnError
    unexpected_bed_vals = defaultdict(list)
    for p in primers:
        for k, v in p.unexpected_vals.items():
            unexpected_bed_vals[k].extend(v)
    for k, l in dict(unexpected_bed_vals).items():
        logger.info(
            f"Expected value of type {k[1]} for optional column {k[0]}. Got"
            f" the following values instead in {len(l)} cases:\n\t\t\t\t{l}"
        )
        logger.info(
            "Values in optional columns do not affect the functionality of"
            " NoTrAmp. As long as the core values (start, end, name) are c"
            "orrect there is nothing to worry about. Frequently primer seq"
            "uences that are included in primer bed files are the reason."
            )
    return sorted(primers, key=lambda x: x.end)


def generate_amps(primers):
    """generate amplicon objects"""
    logger.info("generating amplicon objects")
    amp_nums = set([primer.amp_no for primer in primers])
    amps = []
    amp_lens = []
    for num in amp_nums:
        amp_obj = Amp(num, [primer for primer in primers if primer.amp_no == num])
        amps.append(amp_obj)
        amp_lens.append(amp_obj.get_max_len())
    av_amp_len = sum(amp_lens) / len(amp_lens)
    return sorted(amps, key=lambda x: x.name), av_amp_len


def load_fasta(rfa, reads):
    for line in rfa:
        if line.startswith(">"):
            header = line.strip().split(" ")[0]
            seq = next(rfa)
            chk = aux.chk_fablock_integrity(header, seq)
            if all(chk):
                read = Read(header, seq.strip())
            else:
                # ignores if all in block are empty (None), as can happen at eof
                if chk[0] is False:
                    raise aux.LoadReadsError()
            reads[read.name] = read
    return reads        


def load_fastq(kw, rfa, reads):
    """Load reads from fastq file"""
    logger.info("Loading reads from fastq file")
    read_ct = 0
    for block in aux.fastq_iterator(rfa):
        if len(block) == 3:
            title, seq, qual = block
            read = Read(title.split(" ")[0], seq.strip(), fastq=True)
            read_ct += 1
            if kw["fq_out"]:
                read.plus = "+\n"
                read.qstr = qual
            reads[read.name] = read
        else:
            # ignores if all in block are empty (None), as can happen at eof
            if block is False:
                raise aux.LoadReadsError()    
    return reads


def load_reads(kw, step="amp_cov"):
    """load reads into memory"""
    logger.info("loading reads")
    read_dct = {}
    fastq = kw["fastq_in"]
    if step == "map_trim":
        fastq = kw["fq_out"]
    try:
        with open(kw["reads"], "r", encoding="utf-8") as rfa:
            if fastq:
                reads = load_fastq(kw, rfa, read_dct)
            else:
                reads = load_fasta(rfa, read_dct)
    except Exception as e:
        aux.analyze_read_file(kw["reads"], fastq=fastq)
        tb = e.__traceback__
        raise aux.LoadReadsError().with_traceback(tb)
    return reads


def name_out_paf(reads, reference, mod_name):
    """construct name for output reads-file"""
    logger.debug("naming mapping file")
    read_dir, reads_file = path.split(reads)
    reads_name = ".".join(reads_file.split(".")[:-1])
    ref_name = ".".join(path.split(reference)[1].split(".")[:-1])
    paf_name = f"{reads_name}_mapto_{ref_name}.{mod_name}.paf"
    return path.join(read_dir, paf_name)


def name_out_reads(reads, suffix, outdir, kw):
    """construct name for output reads-file"""
    logger.debug("naming output reads-file")
    read_dir, reads_file = path.split(reads)
    rf_spl = reads_file.split(".")
    reads_name = ".".join(rf_spl[:-1])
    fext = "fastq" if kw["fq_out"] else "fasta"
    out_name = f"{reads_name}.{suffix}.{fext}"
    if outdir:
        out_path = path.join(outdir, out_name)
    else:
        out_path = path.join(read_dir, out_name)
    return out_path


def map_reads(reads, reference, out_paf, seq_tech="ont"):
    """mapping reads"""
    logger.info("mapping reads")
    try:
        seq_tech = "map-" + seq_tech
        subprocess.run(
            ["minimap2", "-x", seq_tech, reference, reads, "-o", out_paf, "--secondary=no"],
            check=True,
            capture_output=True
        )
    except subprocess.CalledProcessError as err:
        log_sp_error(err, "Mapping of reads to reference failed")
    return out_paf


def gen_mapping_objs(mm2_paf):
    """generate mapping objects"""
    logger.info("creating mapping objects")
    with open(mm2_paf, "r", encoding="utf-8") as paf:
        map_dct = defaultdict(list)
        for line in paf:
            if len(line.strip()) > 0:
                lisp = line.strip().split("\t")
                mapping = Mapping(*lisp[:12], lisp[12:])
                map_dct[mapping.qname].append(mapping)
    return map_dct


def filter_read_mappings(mappings, min_len, max_len):
    """filter reads mappings"""
    logger.info("filtering read mappings")
    mappings = [m for m in mappings if min_len < m.qlen < max_len]
    mappings = [m for m in mappings if m.qual == 60]
    return mappings


def create_filt_mappings(mm2_paf, av_amp_len, set_min_len, set_max_len):
    """create filtered mapping objects"""
    logger.info("creating filtered mapping objects")
    min_len = av_amp_len*0.5 if not set_min_len else set_min_len
    max_len = av_amp_len*1.5 if not set_max_len else set_max_len
    map_dct = gen_mapping_objs(mm2_paf)
    mult_maps = {n: ml for n, ml in map_dct.items() if len(ml) > 1}
    mappings = [m for k, l in map_dct.items() for m in l]
    mappings = filter_read_mappings(mappings, min_len, max_len)
    mono_mappings = [m for m in mappings if m.qname not in mult_maps]
    dual_mappings = {k: v for k, v in mult_maps.items() if len(v) == 2}
    incl_max = [
        max(dual_mappings[mname], key=lambda x: x.matches)
        for mname in dual_mappings
    ]
    incl_max = filter_read_mappings(incl_max, min_len, max_len)
    mono_mappings.extend(incl_max)
    mappings = mono_mappings
    return sorted(mappings, key=lambda x: (x.tend, x.tstart))


def run_amp_cov_cap(kw):
    """Cap coverage per amplicon and return subsampled reads"""
    logger.info("Start capping of read-depths per amplicon")
    primers = create_primer_objs(kw["primers"], kw["name_scheme"])
    out_paf = name_out_paf(kw["reads"], kw["reference"], "cap")
    mm2_paf = map_reads(kw["reads"], kw["reference"], out_paf, kw["seq_tec"])
    amps, av_amp_len = generate_amps(primers)
    mappings = create_filt_mappings(
        mm2_paf, av_amp_len, kw["set_min_len"], kw["set_max_len"]
        )
    binned = amp_cov.bin_mappings_ac(
        amps, mappings, kw["max_cov"], kw["margins"], kw["figures"], kw
        )
    cap_out = name_out_reads(kw["reads"], "cap", kw["out_dir"], kw)
    try:
        mem_fit = aux.chk_mem_fit(kw)
    except Exception:
        mem_fit = False
        logger.warning(
            "Memory check failed, likely because module psutil could not be fou"
            "nd. NoTrAmp will still work, but for safety reasons all writing an"
            "d loading of reads will be performed from file instead from memory"
            ". This can increase runtime significantly."
            )

    if mem_fit:
        loaded_reads = load_reads(kw)
        amp_cov.write_capped_from_loaded(binned, loaded_reads, cap_out, kw)
        if kw["split_out"]:
            amp_cov.write_to_split_files(binned, loaded_reads, cap_out, kw)
        del loaded_reads
    else:
        amp_cov.write_capped_from_file(binned, kw["reads"], cap_out, kw)
        if kw["split_out"]:
            logger.warning(
                "Writing to split files could not be performed, because reads w"
                "ere not loaded into memory. This can occur if the in-file is t"
                "oo large for your memory or when NoTrAmp is run without psutil"
                "."
                )
    os.remove(out_paf)
    return cap_out, primers, amps, av_amp_len


def run_map_trim(kw):
    """trimming/clipping of reads"""
    logger.info("Start trimming/clipping of reads")
    if kw["all"]:
        primers = kw["primer_objects"]
    else:
        primers = create_primer_objs(kw["primers"], kw["name_scheme"])
    amps, av_amp_len = generate_amps(primers)
    out_paf = name_out_paf(kw["reads"], kw["reference"], "trim")
    mm2_paf = map_reads(kw["reads"], kw["reference"], out_paf, kw["seq_tec"])
    mappings = create_filt_mappings(
        mm2_paf, av_amp_len, kw["set_min_len"], kw["set_max_len"]
        )
    amps_bin_maps = map_trim.bin_mappings_mt(amps, mappings, kw["margins"])
    loaded_reads = load_reads(kw, step="map_trim")
    amps_bin_reads = map_trim.load_amps_with_reads(amps_bin_maps, loaded_reads)
    clipped_out = name_out_reads(kw["reads"], "clip", kw["out_dir"], kw)
    map_trim.clip_and_write_out(
        amps_bin_reads, clipped_out, kw["incl_prim"], kw["fq_out"]
        )
    os.remove(out_paf)
    return clipped_out


def get_main_logger(outdir):
    log_conf_path = path.join(
        path.dirname(__file__),
        "resources",
        "logging.conf"
        )
    logging.config.fileConfig(
        log_conf_path,
        disable_existing_loggers=False,
        defaults={'logfilename': path.join(outdir, "notramp.log")}
    )
    logger = logging.getLogger(__name__)
    return logger


def run_notramp():
    """NoTramp main function"""
    prstart = perf_counter()
    kwargs = vars(get_arguments())

    if not kwargs["out_dir"]:
        kwargs["out_dir"] = path.split(kwargs["reads"])[0]
    if not path.exists(kwargs["out_dir"]):
        os.makedirs(kwargs["out_dir"], mode=777, exist_ok=False)
        os.chmod(kwargs["out_dir"], 0o777)

    logger = get_main_logger(kwargs["out_dir"])
    logger.info("notramp version %s", __version__)
    logger.info("notramp started with: %s", kwargs)
    pkg_dir = path.split(path.abspath(__file__))[0]
    if kwargs["name_scheme"] == DEFAULTS["name_scheme"]:
        kwargs["name_scheme"] = path.join(
                                            pkg_dir,
                                            "resources",
                                            kwargs["name_scheme"] + ".json"
                                            )

    kwargs["fastq_in"] = aux.fastq_autoscan(kwargs["reads"])
    # kwargs["fq_in"] = True if kwargs["fastq_in"] else False
    if kwargs["fq_out"]:
        if not kwargs["fastq_in"]:
            kwargs["fq_out"] = False
            logger.warning("No fastq input. Outfile-type was changed to fasta")

    if kwargs["cov"]:
        run_amp_cov_cap(kwargs)
    elif kwargs["trim"]:
        run_map_trim(kwargs)
    elif kwargs["all"]:
        capped_reads, prims, amps, amp_lens = run_amp_cov_cap(kwargs)
        kwargs["reads"] = capped_reads
        kwargs["primer_objects"] = prims
        run_map_trim(kwargs)
    else:
        print("Missing argument for notramp operation mode")

    prend = perf_counter()
    logger.info(
            "finished notramp after %s seconds of runtime",
            str(prend-prstart)
        )


if __name__ == "__main__":
    run_notramp()

import argparse
import logging
import subprocess
from time import perf_counter
import random
import json
import sys
import os
from collections import defaultdict

import map_trim
import amp_cov


def get_arguments():
    parser = argparse.ArgumentParser(description=
        "NoTrAmp is a Tool for the preprocessing of long NGS reads (ONT/PacBio) generated with amplicon-tiling approaches. "
        " It can normalize read-depth by capping coverage of each amplicon to a provided threshold and trim amplicons to their "
        " appropriate length by removing primers, barcodes and adpaters in a single clipping step.", add_help=True)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-p", "--primers", required=True, 
        help="Path to primer bed-file (primer-names must adhere to a consistent naming scheme see readme)")
    required_args.add_argument("-r", "--reads", required=True, 
        help="Path to sequencing reads fasta")
    required_args.add_argument("-g", "--reference", required=True, 
        help="Path to reference (genome)")
    
    read_input = required_args.add_mutually_exclusive_group(required=True)
    read_input.add_argument("-a", "--all", 
        help="Perform read depth normalization by coverage-capping/downsampling first, then clipp the normalized reads.", 
        default=False,  action='store_true')
    read_input.add_argument("-c", "--cov", 
        help="Perform only read-depth normalization/downsampling.", 
        default=False, action='store_true')
    read_input.add_argument("-t", "--trim", 
        help="Perform only trimming to amplicon length (excluding primers).", 
        default=False, action='store_true')


    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-o", dest='out_dir', 
        default=False, 
        help="Optionally specify a directory for saving of outfiles. If this argument is not given, "
        "out-files will be saved in the directory where the input reads are located. [default=False]")
    optional_args.add_argument("-m", dest='max_cov', 
        default=200, type=int,
        help="Provide threshold for maximum read-depth per amplicon as integer value. [default=200]")
    optional_args.add_argument("-s", dest='seq_tec', 
        default="ont", 
        help="Specify long-read sequencing technology (ont/pb). [default='ont']")
    optional_args.add_argument("-n", dest='name_scheme', 
        default="artic_nCoV_scheme", 
        help="Provide path to json-file containing a naming scheme which is consistently used for all primers. [default='artic_nCoV_scheme']")
    optional_args.add_argument("--set_min_len", dest='set_min_len', 
        default=False, type=int,
        help="Set a minimum required length for alignments of reads to amplicon. "
        "If this is not set the min_len will be 0.5 * average_amp_len. "
        "If amplicon sizes are relatively homogenous this parameter is not required [default=False]")
    optional_args.add_argument("--set_max_len", dest='set_max_len', 
        default=False, type=int,
        help="Set a maximum required length for alignments of reads to amplicon. "
        "If this is not set the max_len will be 1.5 * average_amp_len. "
        "If amplicon sizes are relatively homogenous this parameter is not required [default=False]")


    args = parser.parse_args()

    return args


class Read:
    def __init__(self, header, seq, fastq=False):
        self.header = header
        self.fastq = fastq
        self.name = self.header.split("@")[1] if self.fastq else self.header.split(">")[1]
        self.seq = seq


    def clip_primers(self, fwp_boundary, revp_boundary, mapping):
        qlen, samestrand = mapping.qlen, mapping.samestrand
        qstart, qend = mapping.qstart, mapping.qend
        tstart, tend = mapping.tstart, mapping.tend

        if samestrand:
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

        else:
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

        self.seq = self.seq[clip_left:clip_right]

        return clip_left, qlen - clip_right

class Mapping:
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

    def gen_kw_attr(self):
        kwattr_dict = {kw.split(":")[0]: kw.split(":")[-1] for kw in self.kwattr}
        for key in kwattr_dict:
            self.__dict__[key] = kwattr_dict[key]

    def eval_strand(self, strand_info):
        if strand_info == "+":
            return True
        else:
            return False

class Primer:
    def __init__(self, contig, start, end, name, score, strand, naming_scheme):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand
        self.type = "primary"
        self.scheme = naming_scheme
        self.amp_no, self.pos = self.get_name_infos()
        

    def get_name_infos(self):
        # ls = self.name.split("_")
        ls = self.name.split(self.scheme["sep"])
        if len(ls) == self.scheme["max_len"]:
            self.type = "alt"
            self.alt_name = ls[self.scheme["alt"]]
        return ls[self.scheme["amp_num"]], ls[self.scheme["position"]]

class Amp:
    def __init__(self, amp_no, primers):
        self.name = int(amp_no)
        self.primers = primers
        self.start = min([primer.start for primer in primers])
        self.end = max([primer.end for primer in primers])
        self.max_len = self.end - self.start
        self.fwp_boundary = max(
            [prim for prim in primers if prim.pos == prim.scheme["fw_indicator"]], 
            key=lambda x: x.end
        ).end
        self.revp_boundary = min(
            [prim for prim in primers if prim.pos == prim.scheme["rev_indicator"]], 
            key=lambda x: x.start
        ).start
        self.mappings = dict()
        self.read_names = list()
        self.reads = list()
        self.reads_dct = dict()
    
    
    def random_sample(self, cutoff):
        if len(self.read_names) > cutoff:
            self.selected = random.choices(self.read_names, k=cutoff)
        else:
            self.selected = self.read_names

    def primer_clipping_all(self):
        clip_ct_left = 0
        clip_ct_right = 0
        for read in self.reads_dct:
            try:
                mapping = self.mappings[read]
                left, right = self.reads_dct[read].clip_primers(
                    self.fwp_boundary, self.revp_boundary, mapping
                )
                clip_ct_left += left
                clip_ct_right += right
            except KeyError as e:
                print(f"KeyError in primer_clipping_all: {e}")
        return clip_ct_left, clip_ct_right


def log_sp_error(error, message):
    logging.error(message)
    print("Check notramp.log for error details")
    logging.error(error.stdout.decode('utf-8'))
    logging.info("Exiting integration locator")
    sys.exit()


def fastq_autoscan(read_file):
    with open(read_file, "r") as rf:
        line_ct = 0
        fastq = False
        while line_ct < 100:
            for line in rf:
                line_ct += 1
                if line.startswith("@"):
                    fastq = True
                    break
    return fastq

        
def create_primer_objs(primer_bed, name_scheme="artic_nCov_scheme"):
    with open(name_scheme, "r") as scheme:
        pd = json.loads(scheme.read())

    with open(primer_bed, "r") as bed:
        primers = list()
        for line in bed:
            prim = Primer(*line.strip().split("\t"), pd)
            primers.append(prim)
    return sorted(primers, key=lambda x: x.end)


def generate_amps(primers):
    amp_nums = set([primer.amp_no for primer in primers])
    amps = list()
    amp_lens = list()
    for num in amp_nums:
        ao = Amp(num, [primer for primer in primers if primer.amp_no == num])
        amps.append(ao)
        amp_lens.append(ao.max_len)
    av_amp_len = sum(amp_lens) / len(amp_lens)
    return sorted(amps, key=lambda x: x.name), av_amp_len


def load_reads(read_file):
    reads = dict()
    fastq = fastq_autoscan(read_file)
    header_ind = "@" if fastq else ">"
    with open(read_file, "r") as rfa:
        for line in rfa:
            if line.startswith(header_ind):
                header = line.strip().split(" ")[0]
                seq = next(rfa)
                read = Read(header, seq.strip(), fastq=fastq)
                reads[read.name] = read
    return reads


def name_out_paf(reads, reference, mod_name):
    read_dir, reads_file = os.path.split(reads)
    reads_name = ".".join(reads_file.split(".")[:-1])
    ref_name =  ".".join(os.path.split(reference)[1].split(".")[:-1])
    paf_name = f"{reads_name}_mapto_{ref_name}.{mod_name}.paf"
    return os.path.join(read_dir, paf_name)


def map_reads(reads, reference, out_paf, seq_tech="ont"):
    try:
        seq_tech = "map-" + seq_tech
        subprocess.run(["minimap2", "-x", seq_tech, reference, reads, "-o", out_paf, "--secondary=no"], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        log_sp_error(e, "Mapping of reads to reference failed")
    return out_paf


def filter_read_mappings(mappings, min_len, max_len):
    mappings = [m for m in mappings if min_len < m.qlen < max_len]
    mappings = [m for m in mappings if m.qual == 60]
    return mappings


def create_read_mappings(mm2_paf, av_amp_len, set_min_len, set_max_len):
    min_len = av_amp_len*0.5 if not set_min_len else set_min_len
    max_len = av_amp_len*1.5 if not set_max_len else set_max_len
    with open(mm2_paf, "r") as paf:
        map_dct = defaultdict(list)
        for line in paf:
            if len(line.strip()) > 0:
                ls = line.strip().split("\t")
                mapping = Mapping(*ls[:12], ls[12:])
                map_dct[mapping.qname].append(mapping)
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


if __name__ == "__main__":
    logger = logging.getLogger(name='main')
    logging.basicConfig(level=logging.DEBUG)

    start = perf_counter()

    args = get_arguments()
    print(args)
    primer = args.primers
    reads = args.reads
    ref = args.reference
    max_cov = args.max_cov
    out_dir = args.out_dir
    seq_tec = args.seq_tec
    pkg_dir = os.path.split(os.path.abspath(__file__))[0]
    name_scheme = args.name_scheme
    if name_scheme == 'artic_nCoV_scheme':
        name_scheme = os.path.join(pkg_dir, "resources", name_scheme + ".json")
    # set_min_len = int(args.set_min_len) if args.set_min_len else args.set_min_len
    # set_max_len = int(args.set_max_len) if args.set_max_len else args.set_max_len
    set_min_len = args.set_min_len
    set_max_len = args.set_max_len

    if args.cov:
        amp_cov.run_amp_cov_cap(primer, name_scheme, reads, ref, max_cov, set_min_len, set_max_len, seq_tech=seq_tec)
    elif args.trim:
        map_trim.run_map_trim(primer, name_scheme, reads, ref, set_min_len, set_max_len, seq_tech=seq_tec)
    elif args.all:
        capped_reads = amp_cov.run_amp_cov_cap(primer, name_scheme, reads, ref, max_cov, set_min_len, set_max_len, seq_tech=seq_tec)
        map_trim.run_map_trim(primer, name_scheme, capped_reads, ref, set_min_len, set_max_len, seq_tech=seq_tec)
    else:
        print("Missing argument for notramp operation mode")
    
    end = perf_counter()
    runtime = end-start
    print("runtime", runtime)






import logging
import logging.config
import subprocess
from time import perf_counter
import random
import json
import sys
import os
from collections import defaultdict
from classes import *


logging.config.fileConfig('logging.conf', disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def log_sp_error(error, message):
    """custom logging of subprocess errors"""
    logging.error(message)
    print("Check notramp.log for error details")
    logging.error(error.stdout.decode('utf-8'))
    logging.error("Exiting notramp")
    sys.exit()


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


def create_primer_objs(primer_bed, name_scheme):
    """generate primer objects"""
    logger.info("generating primer objects")
    with open(name_scheme, "r", encoding="utf-8") as scheme:
        psd = json.loads(scheme.read())

    with open(primer_bed, "r", encoding="utf-8") as bed:
        primers = []
        for line in bed:
            contig, start, end, name, score, strand = line.strip().split("\t")
            add_info = {"contig": contig, "strand": strand, "score": score}
            prim = Primer(start, end, name, add_info, psd)
            primers.append(prim)
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


def load_reads(read_file):
    """load reads into memory"""
    logger.info("loading reads")
    reads = {}
    fastq = fastq_autoscan(read_file)
    header_ind = "@" if fastq else ">"
    with open(read_file, "r", encoding="utf-8") as rfa:
        for line in rfa:
            if line.startswith(header_ind):
                header = line.strip().split(" ")[0]
                seq = next(rfa)
                read = Read(header, seq.strip(), fastq=fastq)
                reads[read.name] = read
    return reads


def name_out_paf(reads, reference, mod_name):
    """construct name for output reads-file"""
    logger.debug("naming mapping file")
    read_dir, reads_file = os.path.split(reads)
    reads_name = ".".join(reads_file.split(".")[:-1])
    ref_name = ".".join(os.path.split(reference)[1].split(".")[:-1])
    paf_name = f"{reads_name}_mapto_{ref_name}.{mod_name}.paf"
    return os.path.join(read_dir, paf_name)


def name_out_reads(reads, suffix, outdir):
    """construct name for output reads-file"""
    logger.debug("naming output reads-file")
    read_dir, reads_file = os.path.split(reads)
    rf_spl = reads_file.split(".")
    reads_name = ".".join(rf_spl[:-1])
    out_name = f"{reads_name}.{suffix}.fasta"
    if outdir:
        out_path = os.path.join(outdir, out_name)
    else:
        out_path = os.path.join(read_dir, out_name)
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
    with open(mm2_paf, "r", encoding="utf-8") as paf:
        map_dct = defaultdict(list)
        for line in paf:
            if len(line.strip()) > 0:
                lisp = line.strip().split("\t")
                mapping = Mapping(lisp[0], lisp[1:12], lisp[12:])
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
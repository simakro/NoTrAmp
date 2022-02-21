# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.

from collections import defaultdict
# import subprocess
import logging
from notramp import *
import os



class Read:
    def __init__(self, header, seq):
        self.header = header
        self.name = self.header.split(">")[1]
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


# class Mapping:
#     def __init__(
#         self,
#         qname,
#         qlen,
#         qstart,
#         qend,
#         samestrand,
#         tname,
#         tlen,
#         tstart,
#         tend,
#         matches,
#         total_bp,
#         qual,
#         kwattr,
#     ):
#         self.qname = qname
#         self.qlen = int(qlen)
#         self.qstart = int(qstart)
#         self.qend = int(qend)
#         self.samestrand = self.eval_strand(samestrand)
#         self.tname = tname
#         self.tlen = int(tlen)
#         self.tstart = int(tstart)
#         self.tend = int(tend)
#         self.matches = int(matches)
#         self.total_bp = int(total_bp)
#         self.qual = int(qual)
#         self.kwattr = kwattr
#         self.gen_kw_attr()

#     def gen_kw_attr(self):
#         kwattr_dict = {kw.split(":")[0]: kw.split(":")[-1] for kw in self.kwattr}
#         for key in kwattr_dict:
#             self.__dict__[key] = kwattr_dict[key]

#     def eval_strand(self, strand_info):
#         if strand_info == "+":
#             return True
#         else:
#             return False


# class Primer:
#     def __init__(self, contig, start, end, name, score, strand):
#         self.contig = contig
#         self.start = int(start)
#         self.end = int(end)
#         self.name = name
#         self.score = int(score)
#         self.strand = strand
#         self.type = "primary"
#         self.amp_no, self.pos = self.get_name_infos()

#     def get_name_infos(self):
#         ls = self.name.split("_")
#         if len(ls) == 4:
#             self.type = "alt"
#             self.alt_name = ls[3]
#         return ls[1], ls[2]


# class Amp:
#     def __init__(self, amp_no, primers):
#         self.name = int(amp_no)
#         self.primers = primers
#         self.start = min([primer.start for primer in primers])
#         self.end = max([primer.end for primer in primers])
#         self.max_len = self.end - self.start
#         self.fwp_boundary = max(
#             [prim for prim in primers if prim.pos == "LEFT"], key=lambda x: x.end
#         ).end
#         self.revp_boundary = min(
#             [prim for prim in primers if prim.pos == "RIGHT"], key=lambda x: x.start
#         ).start
#         self.mappings = dict()
#         self.reads = dict()

#     def primer_clipping_all(self):
#         clip_ct_left = 0
#         clip_ct_right = 0
#         for read in self.reads:
#             try:
#                 mapping = self.mappings[read]
#                 left, right = self.reads[read].clip_primers(
#                     self.fwp_boundary, self.revp_boundary, mapping
#                 )
#                 clip_ct_left += left
#                 clip_ct_right += right
#             except KeyError as e:
#                 print(f"KeyError in primer_clipping_all: {e}")
#         return clip_ct_left, clip_ct_right




# def create_primer_objs(primer_bed):
#     with open(primer_bed, "r") as bed:
#         primers = list()
#         for line in bed:
#             prim = Primer(*line.strip().split("\t"))
#             primers.append(prim)
#     return sorted(primers, key=lambda x: x.end)


# def generate_amps(primers):
#     amp_nums = set([primer.amp_no for primer in primers])
#     amps = list()
#     for num in amp_nums:
#         ao = Amp(num, [primer for primer in primers if primer.amp_no == num])
#         amps.append(ao)
#     return sorted(amps, key=lambda x: x.name)


# def map_reads(reads, reference, out_paf, seq_tech="map-ont"):
#     try:
#         subprocess.run(["minimap2", "-x", seq_tech, reference, reads, "-o", out_paf, "--secondary=no"], check=True, capture_output=True)
#     except subprocess.CalledProcessError as e:
#         log_sp_error(e, "Mapping of reads to reference failed")
#     return out_paf


# def filter_read_mappings(mappings):
#     mappings = [m for m in mappings if 300 < m.qlen < 600]
#     mappings = [m for m in mappings if m.qual == 60]
#     return mappings


# def create_read_mappings(mm2_paf):
#     with open(mm2_paf, "r") as paf:
#         map_dct = defaultdict(list)
#         for line in paf:
#             if len(line.strip()) > 0:
#                 ls = line.strip().split("\t")
#                 mapping = Mapping(*ls[:12], ls[12:])
#                 map_dct[mapping.qname].append(mapping)
#         mult_maps = {n: ml for n, ml in map_dct.items() if len(ml) > 1}
#         mappings = [m for k, l in map_dct.items() for m in l]
#         mappings = filter_read_mappings(mappings)
#         mono_mappings = [m for m in mappings if m.qname not in mult_maps]
#         dual_mappings = {k: v for k, v in mult_maps.items() if len(v) == 2}
#         incl_max = [
#             max(dual_mappings[mname], key=lambda x: x.matches)
#             for mname in dual_mappings
#         ]
#         incl_max = filter_read_mappings(incl_max)
#         mono_mappings.extend(incl_max)
#         mappings = mono_mappings
#     return sorted(mappings, key=lambda x: (x.tend, x.tstart))


def bin_mappings(amp_bins, mappings):
    binned = list()
    na = list()
    while len(amp_bins) > 0:
        if len(mappings) > 0:
            if mappings[0].tend <= amp_bins[0].end + 5:
                if mappings[0].tstart >= amp_bins[0].start - 5:
                    amp_bins[0].mappings[mappings[0].qname] = mappings[0]
                    mappings.pop(0)
                else:
                    na.append(mappings[0].qname)
                    mappings.pop(0)
            else:
                binned.append(amp_bins[0])
                amp_bins.pop(0)
        else:
            break
    return binned


def load_reads(read_fasta, amp_bins):
    reads = dict()
    with open(read_fasta, "r") as rfa:
        for line in rfa:
            if line.startswith(">"):
                # print(line)
                header = line.strip().split(" ")[0]
                seq = next(rfa)
                read = Read(header, seq.strip())
                reads[read.name] = read

    for amp in amp_bins:
        amp.reads_dct = {k: v for k, v in reads.items() if k in amp.mappings}


    return amp_bins


def clip_and_write_out(amp_bins, clipped_out):
    with open(clipped_out, "w") as out:
        clip_ct = {"left": 0, "right": 0}
        for amp in amp_bins:
            left, right = amp.primer_clipping_all()
            clip_ct["left"] += left
            clip_ct["right"] += right
            for read in amp.reads_dct:
                out.write(amp.reads_dct[read].header + "\n")
                out.write(amp.reads_dct[read].seq + "\n")
    print(
        f"{clip_ct['left']} bases were clipped from the fw/left-primer side of reads and "
        f"{clip_ct['right']} bases were clipped from the rev/right-primer side of reads"
    )


def name_out_paf(reads, reference):
    read_dir, reads_file = os.path.split(reads)
    reads_name = reads_file.split(".")[:-1]
    ref_name = os.path.split(reference)[1].split(".")[:-1]
    paf_name = f"{reads_name}_mapto_{ref_name}.paf"
    return os.path.join(read_dir, paf_name)


def name_clipped(reads):
    read_dir, reads_file = os.path.split(reads)
    rf_spl = reads_file.split(".")
    reads_name, ext = rf_spl[:-1], rf_spl[-1]
    clipped_name = f"{reads_name}.clip.{ext}"
    return os.path.join(read_dir, clipped_name)


def run_map_trim(primer_bed, reads, reference, seq_tech="map-ont"):
    primers = create_primer_objs(primer_bed)
    amps = generate_amps(primers)
    out_paf = name_out_paf(reads, reference)
    mm2_paf = map_reads(reads, reference, out_paf, seq_tech=seq_tech)
    mappings = create_read_mappings(mm2_paf)
    amps_bin_maps = bin_mappings(amps, mappings)
    amps_bin_reads = load_reads(reads, amps_bin_maps)
    clipped_out = name_clipped(reads)
    clip_and_write_out(amps_bin_reads, clipped_out)
    os.remove(out_paf)
    return clipped_out
    

if __name__ == "__main__":
    import sys

    primer_bed = sys.argv[1]
    reads = sys.argv[2]
    reference = sys.argv[3]

    run_map_trim(primer_bed, reads, reference, seq_tech="map-ont")





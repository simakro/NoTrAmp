# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.

import logging
import notramp  as nta
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


def name_clipped(reads):
    read_dir, reads_file = os.path.split(reads)
    rf_spl = reads_file.split(".")
    reads_name, ext = ".".join(rf_spl[:-1]), rf_spl[-1]
    clipped_name = f"{reads_name}.clip.{ext}"
    return os.path.join(read_dir, clipped_name)


def run_map_trim(primer_bed, name_scheme, reads, reference, seq_tech="map-ont"):
    primers = nta.create_primer_objs(primer_bed, name_scheme=name_scheme)
    amps = nta.generate_amps(primers)
    out_paf = nta.name_out_paf(reads, reference, "trim")
    mm2_paf = nta.map_reads(reads, reference, out_paf, seq_tech=seq_tech)
    mappings = nta.create_read_mappings(mm2_paf)
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





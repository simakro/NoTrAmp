# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.

from collections import defaultdict
import notramp as nta
import os


def bin_mappings(amp_bins, mappings):
    binned = list()
    na = list()
    while len(amp_bins) > 0:
        if len(mappings) > 0:
            if mappings[0].tend <= amp_bins[0].end + 5:
                if mappings[0].tstart >= amp_bins[0].start - 5:
                    amp_bins[0].read_names.append(mappings[0].qname)
                    mappings.pop(0)
                else:
                    na.append(mappings[0].qname)
                    mappings.pop(0)
            else:
                binned.append(amp_bins[0])
                amp_bins.pop(0)
        else:
            break

    for bin in binned:
        bin.random_sample(200)
    #     print(bin.name, len(bin.read_names), "selected:", len(bin.selected))
    # print("na", len(na))

    return binned


def write_capped_reads(binned, reads, fa_out):
    all_picks = ["@" + name for amp in binned for name in amp.selected]
    with open(reads, "r") as fq, open(fa_out, "w") as fa:
        for line in fq:
            if line.startswith("@"):
                readname = line.split(" ")[0]
                if readname in all_picks:
                    readname = readname.replace("@", ">")
                    fa.write(readname + "\n")
                    fa.write(next(fq))


def name_capped(reads):
    read_dir, reads_file = os.path.split(reads)
    rf_spl = reads_file.split(".")
    reads_name = ".".join(rf_spl[:-1]) # ext = rf_spl[-1]
    capped_name = f"{reads_name}.cap.fasta"
    return os.path.join(read_dir, capped_name)


def run_amp_cov_cap(primer_bed, name_scheme, reads, reference, seq_tech="map-ont"):
    primers = nta.create_primer_objs(primer_bed, name_scheme=name_scheme)
    out_paf = nta.name_out_paf(reads, reference, "cap")
    mm2_paf = nta.map_reads(reads, reference, out_paf, seq_tech=seq_tech)
    amps = nta.generate_amps(primers)
    mappings = nta.create_read_mappings(mm2_paf)
    binned = bin_mappings(amps, mappings)
    fa_out = name_capped(reads)
    write_capped_reads(binned, reads, fa_out)
    os.remove(out_paf)
    return fa_out


if __name__ == "__main__":
    import sys

    primer_bed = sys.argv[1]
    reads = sys.argv[2]
    reference = sys.argv[3]

    run_amp_cov_cap(primer_bed, reads, reference, seq_tech="map-ont")


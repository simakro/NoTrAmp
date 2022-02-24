# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.

from collections import defaultdict
import notramp as nta
import os
import psutil
import logging

logger = logging.getLogger(name=__name__)
# logger.setLevel(logging.INFO)
# # logging.basicConfig(filename='notramp.log', level=logging.DEBUG, format='%(asctime)s:%(levelname)s:%(message)s:%(lineno)s')
# # logger = nta.create_logger()
logger.info("From amp_cov")
logger.warning("From amp_cov")

def bin_mappings(amp_bins, mappings, max_cov):
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
            binned.append(amp_bins[0])
            break
    

    num_selected = 0
    for bin in binned:
        bin.random_sample(max_cov)
        print(bin.name, len(bin.read_names), "selected:", len(bin.selected))
        num_selected += len(bin.selected)
    print("na", len(na))
    print("Total reads selected:", num_selected)

    return binned


def write_capped_from_file(binned, reads, fa_out):
    fastq = nta.fastq_autoscan(reads)
    hi = "@" if fastq else ">"
    all_picks = [hi + name for amp in binned for name in amp.selected]
    with open(reads, "r") as fi, open(fa_out, "w") as fa:
        for line in fi:
            if line.startswith(hi):
                readname = line.split(" ")[0]
                if readname in all_picks:
                    if fastq:
                        readname = readname.replace("@", ">")
                    fa.write(readname + "\n")
                    fa.write(next(fi))


def write_capped_from_loaded(binned, loaded_reads, fa_out):
    all_picks = [name for amp in binned for name in amp.selected]
    with open(fa_out, "w") as fa:
        for name in all_picks:
            try:
                header = ">" + name + "\n"
                seq = loaded_reads[name].seq + "\n"
                fa.write(header)
                fa.write(seq)
            except KeyError as e:
                print(f"Error: read {name} was not found in loaded reads")


def chk_mem_fit(read_path):
    fastq = nta.fastq_autoscan(read_path)
    rf_size = os.path.getsize(read_path)
    load_size = rf_size if not fastq else rf_size/2
    avail_mem = psutil.virtual_memory()[1]
    if avail_mem / load_size > 4:
        return True
    else:
        return False


def run_amp_cov_cap(**kw):
    primers = nta.create_primer_objs(kw["primers"], kw["name_scheme"])
    out_paf = nta.name_out_paf(kw["reads"], kw["reference"], "cap")
    mm2_paf = nta.map_reads(kw["reads"], kw["reference"], out_paf, kw["seq_tec"])
    amps, av_amp_len = nta.generate_amps(primers)
    mappings = nta.create_read_mappings(mm2_paf, av_amp_len, kw["set_min_len"], kw["set_max_len"])
    binned = bin_mappings(amps, mappings, kw["max_cov"])
    fa_out = nta.name_out_reads(kw["reads"], "cap", kw["out_dir"])
    mem_fit = chk_mem_fit(kw["reads"])
    if mem_fit:
        loaded_reads = nta.load_reads(kw["reads"])
        write_capped_from_loaded(binned, loaded_reads, fa_out)
        del(loaded_reads)
    else:
        write_capped_from_file(binned, kw["reads"], fa_out)
    os.remove(out_paf)
    return fa_out

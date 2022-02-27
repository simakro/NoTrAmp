# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.


from os import path
import logging
import logging.config
import psutil

if __name__ == "amp_cov":
    from nta_aux import fastq_autoscan
else:
    from notramp.nta_aux import fastq_autoscan


log_file_path = path.join(path.dirname(__file__), "resources", "logging.conf")
logging.config.fileConfig(log_file_path, disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def bin_mappings(amp_bins, mappings, max_cov, margins):
    """sort mappings to amplicons"""
    logger.info("sorting mappings to amplicons")
    binned = []
    not_av = []
    while len(amp_bins) > 0:
        if len(mappings) > 0:
            if mappings[0].tend <= amp_bins[0].end + margins:
                if mappings[0].tstart >= amp_bins[0].start - margins:
                    amp_bins[0].reads_dct[mappings[0].qname] = mappings[0].qname
                    mappings.pop(0)
                else:
                    not_av.append(mappings[0].qname)
                    mappings.pop(0)
            else:
                binned.append(amp_bins[0])
                amp_bins.pop(0)
        else:
            binned.append(amp_bins[0])
            break

    num_selected = 0
    for amp_bin in binned:
        amp_bin.random_sample(max_cov)
        amn = amp_bin.name
        rct = str(len(amp_bin.reads_dct))
        slr = str(len(amp_bin.selected))
        logger.info("Amp_%s %s selected: %s", amn, rct, slr)
        num_selected += len(amp_bin.selected)
    return binned


def write_capped_from_file(binned, reads, fa_out):
    """write subsample to outfile using reads-file as source"""
    logger.info("writing subsample to outfile using reads-file as source")
    fastq = fastq_autoscan(reads)
    hin = "@" if fastq else ">"
    all_picks = [hin + name for amp in binned for name in amp.selected]
    with open(reads, "r", encoding="utf-8") as rfi, open(fa_out, "w", encoding="utf-8") as fao:
        for line in rfi:
            if line.startswith(hin):
                readname = line.split(" ")[0]
                if readname in all_picks:
                    if fastq:
                        readname = readname.replace("@", ">")
                    fao.write(readname + "\n")
                    fao.write(next(rfi))


def write_capped_from_loaded(binned, loaded_reads, fa_out):
    """write subsample to outfile using loaded reads as source"""
    logger.info("writing subsample to outfile using loaded reads as source")
    all_picks = [name for amp in binned for name in amp.selected]
    with open(fa_out, "w", encoding="utf-8") as fao:
        for name in all_picks:
            try:
                header = ">" + name + "\n"
                seq = loaded_reads[name].seq + "\n"
                fao.write(header)
                fao.write(seq)
            except KeyError:
                logging.exception("Error: read %s was not found in loaded reads", name)


def chk_mem_fit(read_path):
    """Check for available memory"""
    logger.info("Checking for available memory")
    fastq = fastq_autoscan(read_path)
    rf_size = path.getsize(read_path)
    load_size = rf_size if not fastq else rf_size/2
    avail_mem = psutil.virtual_memory()[1]
    if avail_mem / load_size > 4:
        return True
    else:
        return False

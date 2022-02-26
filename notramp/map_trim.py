# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.

"""module for trimming of reads to amplicon size"""

import os
import logging
import logging.config
import notramp as nta


logging.config.fileConfig("logging.conf", disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def bin_mappings(amp_bins, mappings, margins):
    """sort mappings to amplicons"""
    logger.info("sorting mappings to amplicons")
    binned = []
    not_av = []
    while len(amp_bins) > 0:
        if len(mappings) > 0:
            if mappings[0].tend <= amp_bins[0].end + margins:
                if mappings[0].tstart >= amp_bins[0].start - margins:
                    amp_bins[0].mappings[mappings[0].qname] = mappings[0]
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
    return binned


def load_amps_with_reads(amp_bins, loaded_reads):
    """load amplicon objects with reads"""
    logger.info("loading amplicon objects with reads")
    for amp in amp_bins:
        amp.reads_dct = {k: v for k, v in loaded_reads.items() if k in amp.mappings}
    return amp_bins


def clip_and_write_out(amp_bins, clipped_out, incl_prim):
    """Trim reads to amplicon boundaries (including/excluding primers)"""
    incl_excl = "excluding" if not incl_prim else "including"
    logger.info("Trimming reads to amplicons boundaries %s primers", incl_excl)
    with open(clipped_out, "w", encoding="utf-8") as out:
        clip_ct = {"left": 0, "right": 0}
        clipped_out = 0
        for amp in amp_bins:
            left, right = amp.trim_clip_all(incl_prim)
            clip_ct["left"] += left
            clip_ct["right"] += right
            for read in amp.reads_dct:
                out.write(amp.reads_dct[read].header + "\n")
                out.write(amp.reads_dct[read].seq + "\n")
                clipped_out += 1
    logger.info("%s reads were processed for clipping", clipped_out)
    logger.info("%d bases were clipped from the fw/left-primer side of reads.", clip_ct["left"])
    logger.info("%d bases were clipped from the rev/right-primer side of reads", clip_ct["right"])


def run_map_trim(**kw):
    """trimming/clipping of reads"""
    logger.info("Start trimming/clipping of reads")
    primers = nta.create_primer_objs(kw["primers"], kw["name_scheme"])
    amps, av_amp_len = nta.generate_amps(primers)
    out_paf = nta.name_out_paf(kw["reads"], kw["reference"], "trim")
    mm2_paf = nta.map_reads(kw["reads"], kw["reference"], out_paf, kw["seq_tec"])
    mappings = nta.create_filt_mappings(mm2_paf, av_amp_len, kw["set_min_len"], kw["set_max_len"])
    amps_bin_maps = bin_mappings(amps, mappings, kw["margins"])
    loaded_reads = nta.load_reads(kw["reads"])
    amps_bin_reads = load_amps_with_reads(amps_bin_maps, loaded_reads)
    clipped_out = nta.name_out_reads(kw["reads"], "clip", kw["out_dir"])
    clip_and_write_out(amps_bin_reads, clipped_out, kw["incl_prim"])
    os.remove(out_paf)
    return clipped_out
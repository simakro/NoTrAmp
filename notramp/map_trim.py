# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.


import logging


logger = logging.getLogger(__name__)


def bin_mappings_mt(amp_bins, mappings, margins):
    """sort mappings to amplicons (map_trim)"""
    logger.info("sorting mappings to amplicons (map_trim)")
    binned_ct = 0
    binned = []
    not_av = []
    if len(amp_bins) == 0:
        logger.warning(
            "No mappings of any reads were created! No trimming will occur."
            )
    while len(amp_bins) > 0:
        if len(mappings) > 0:
            if mappings[0].tend <= amp_bins[0].end + margins:
                if mappings[0].tstart >= amp_bins[0].start - margins:
                    amp_bins[0].mappings[mappings[0].qname] = mappings[0]
                    binned_ct += 1
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
    logger.info(
        f"{binned_ct} reads were sorted to bins. "
        f"{len(not_av)} could not be sorted to an amplicon"
        )
    if len(mappings) > 0:
        perc_binned = (binned_ct/len(mappings))*100
        if perc_binned < 75:
            logger.warning(
                "A high proportion (>25%) of reads could not be sorted to ampli"
                "cons. This can occur if you don't use the right amplicon tilin"
                "g scheme (primer bed-file) or can be an indication of sequenci"
                "ng issues."
                )
    return binned


def load_amps_with_reads(amp_bins, loaded_reads):
    """load amplicon objects with reads"""
    logger.info("loading amplicon objects with reads")
    for amp in amp_bins:
        amp.reads_dct = {k: v for k, v in loaded_reads.items() if k in amp.mappings}
    return amp_bins


def clip_and_write_out(amp_bins, clipped_out, incl_prim, output_fq):
    """Trim reads to amplicon boundaries (including/excluding primers)"""
    incl_excl = "excluding" if not incl_prim else "including"
    logger.info("Trimming reads to amplicon boundaries %s primers", incl_excl)
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
                if output_fq:
                    out.write("+\n")
                    out.write(amp.reads_dct[read].qstr + "\n")
                clipped_out += 1
    logger.info("%s reads were processed for clipping", clipped_out)
    logger.info(
        "%d bases were clipped from the fw/left-primer side of reads.",
        clip_ct["left"]
        )
    logger.info(
        "%d bases were clipped from the rev/right-primer side of reads",
        clip_ct["right"]
        )

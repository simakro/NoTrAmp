# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.


import logging


logger = logging.getLogger(__name__)

def bin_mappings_ac(amp_bins, mappings, max_cov, margins, figures, kw):
    """sort mappings to amplicons (amp_cov)"""
    logger.info("sorting mappings to amplicons (amp_cov)")
    binned_ct = 0
    binned = []
    not_av = []
    if len(mappings) == 0:
        logger.warning(
            "No mappings of any reads were created!"
            " No capping or trimming will occur."
            )
    logger.info("mappings of %s reads available", str(len(mappings)))

    for m in mappings:
        if m.tend <= amp_bins[0].end + margins:
            if m.tstart >= amp_bins[0].start - margins:
                amp_bins[0].reads_dct[m.qname] = m.qname
                binned_ct += 1
            else:
                not_av.append(m.qname)
        else:
            binned.append(amp_bins[0])
            amp_bins.pop(0)
            # to avoid systematic loss of reads when switching to next amp_bin
            # check if the current mapping fits into the next (or later) bin
            switch_binned = False
            for bin in amp_bins:
                if m.tend <= bin.end + margins:
                    if m.tstart >= bin.start - margins:
                        bin.reads_dct[m.qname] = m.qname
                        binned_ct += 1
                        switch_binned = True
                        break
            if not switch_binned:
                not_av.append(m.qname)
            if len(amp_bins) == 0:
                # find a way to add potential remaining mapping to not_av
                break

    if len(amp_bins) > 0:
        binned.extend(amp_bins)
    
    # not_av.extend([m.qname for m in mappings])
    mappings = []
    logger.info(
        f"{binned_ct} reads were sorted to bins. "
        f"{len(not_av)} could not be sorted to an amplicon."
        )

    num_selected = 0
    avail_reads = {}
    select_reads = {}
    for amp_bin in binned:
        amp_bin.random_sample(max_cov)
        amn = amp_bin.name
        rct = len(amp_bin.reads_dct)
        slr = len(amp_bin.selected)
        avail_reads[amn] = rct
        select_reads[amn] = slr
        logger.info(
            "Amp_%s reads available %s selected: %s", amn, str(rct), str(slr)
            )
        num_selected += len(amp_bin.selected)
    logger.info(
        "A total of %s reads was selected.", str(sum(select_reads.values()))
        )
    if figures:
        try:
            import plot_amp_reads
            plot_amp_reads. gen_overview_fig(avail_reads, select_reads, kw)
        except:
            logger.warning(
                "Figures could not be created. The functionality of NoTrAmp wil"
                "l not be affected otherwise."
                )
    
    logger.debug("%s reads could not be binned to an Amp", str(len(not_av)))
    return binned


def write_read_to_file(rname, read_fobj, out_fobj, infile_fastq, outfile_type):
    """Write reads from in_file to out_file"""
    if infile_fastq:
        if outfile_type == "fasta":
            rname = rname.replace("@", ">")
            out_fobj.write(rname + "\n")
            out_fobj.write(next(read_fobj))
        else:
            try:
                out_fobj.write(rname + "\n")
                out_fobj.write(next(read_fobj))
                out_fobj.write(next(read_fobj))
                out_fobj.write(next(read_fobj))
            except StopIteration:
                logger.exception("Unexpected eof during writing of fastq")
    else:
        out_fobj.write(rname + "\n")
        out_fobj.write(next(read_fobj))


def write_capped_from_file(binned, reads, fa_out, kw):
    """write subsample to outfile using reads-file as source"""
    logger.info("writing subsample to outfile using reads-file as source")
    fastq = kw["fq_in"]
    hin = "@" if fastq else ">"
    all_picks = [hin + name for amp in binned for name in amp.selected]
    with open(reads, "r", encoding="utf-8") as rfi, open(fa_out, "w", encoding="utf-8") as ofo:
        for line in rfi:
            if line.startswith(hin):
                readname = line.split(" ")[0]
                if readname in all_picks:
                    write_read_to_file(readname, rfi, ofo, fastq, kw["fq_out"])


def write_capped_from_loaded(binned, loaded_reads, out_file, kw):
    """write subsample to outfile using loaded reads as source"""
    logger.info("writing subsample to outfile using loaded reads as source")
    all_picks = [name for amp in binned for name in amp.selected]
    with open(out_file, "w", encoding="utf-8") as outf:
        for name in all_picks:
            try:
                if not kw["fq_out"]:
                    header = ">" + name + "\n"
                    seq = loaded_reads[name].seq + "\n"
                    outf.write(header)
                    outf.write(seq)
                else:
                    fq_read = loaded_reads[name]
                    header = "@" + name + "\n"
                    outf.write(header)
                    outf.write(fq_read.seq + "\n")
                    outf.write("+\n")
                    outf.write(fq_read.qstr + "\n")
            except KeyError:
                logging.exception(
                    "Error: read %s was not found in loaded reads", name
                    )

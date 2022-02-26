# NoTrAmp
Normalization and Trimming of long-read (ONT, PB) amplicon sequencing data


NoTrAmp is a Tool for read-depth normalization and trimming of amplicon reads generated with long read technologies (ONT/PacBio).
It can be used in amplicon-tiling approaches to cap coverage of each amplicon and to trim amplicons to their
appropriate length removing barcodes, adpaters and primers (if desired) in a single clipping step.

## Table of Contents

- [Commmand line Options](#options)
- [Dependencies](#depend)

## <a name="options"></a>Commmand line Options
```sh
usage: notramp.py [-h] -p PRIMERS -r READS -g REFERENCE (-a | -c | -t) [-o OUT_DIR] [-m MAX_COV] [--incl_prim] [-s SEQ_TEC] [-n NAME_SCHEME] [--set_min_len SET_MIN_LEN] [--set_max_len SET_MAX_LEN]
                  [--set_margins MARGINS]

NoTrAmp is a Tool for read-depth normalization and trimming of amplicon reads generated with long read technologies (ONT/PacBio). It can be used in amplicon-tiling approaches to cap coverage of each
amplicon and to trim amplicons to their appropriate length removing barcodes, adpaters and primers (if desired) in a single clipping step.

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -p PRIMERS, --primers PRIMERS
                        Path to primer bed-file (primer-names must adhere to a consistent naming scheme see readme)
  -r READS, --reads READS
                        Path to sequencing reads fasta
  -g REFERENCE, --reference REFERENCE
                        Path to reference (genome)
  -a, --all             Perform read depth normalization by coverage-capping/downsampling first, then clip the normalized reads. (mut.excl. with -c, -t)
  -c, --cov             Perform only read-depth normalization/downsampling. (mut.excl. with -a, -t)
  -t, --trim            Perform only trimming to amplicon length (excluding primers by default; to include primers set --incl_prim flag). (mut.excl. with -a, -c)

Optional arguments:
  -o OUT_DIR            Optionally specify a directory for saving of outfiles. If this argument is not given, out-files will be saved in the directory where the input reads are located. [default=False]
  -m MAX_COV            Provide threshold for maximum read-depth per amplicon as integer value. [default=200]
  --incl_prim           Set to False if you want to include the primer sequences in the trimmed reads. By default primers are removed together with all overhanging sequences. [default=False]
  -s SEQ_TEC            Specify long-read sequencing technology (ont/pb). [default='ont']
  -n NAME_SCHEME        Provide path to json-file containing a naming scheme which is consistently used for all primers. [default='artic_nCoV_scheme']
  --set_min_len SET_MIN_LEN
                        Set a minimum required length for alignments of reads to amplicon. If this is not set the min_len will be 0.5*average_amp_len. If amplicon sizes are relatively homogenous this
                        parameter is not required [default=False]
  --set_max_len SET_MAX_LEN
                        Set a maximum required length for alignments of reads to amplicon. If this is not set the max_len will be 1.5*average_amp_len. If amplicon sizes are relatively homogenous this
                        parameter is not required [default=False]
  --set_margins MARGINS
                        Set length of tolerance margins for sorting of mappings to amplicons. [default=5]

## <a name="depend"></a>Requirements/Dependencies
- Python 3.x
- minimap2
- psutil
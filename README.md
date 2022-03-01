# NoTrAmp
Normalization and Trimming of long-read (ONT, PB) amplicon sequencing data


NoTrAmp is a Tool for read-depth normalization and trimming of amplicon reads generated with long read technologies (ONT/PacBio).
It can be used in amplicon-tiling approaches to cap coverage of each amplicon and to trim amplicons to their
appropriate length removing barcodes, adpaters and primers (if desired) in a single clipping step.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Dependencies](#depend)


## <a name="install"></a>Installation
install with pip:
```sh
pip install notramp
```
install with conda:
```sh
conda create -n notramp
conda activate notramp
conda install -c simakro notramp
```

or

```
conda create -n notramp -c simakro notramp
conda activate notramp
```


## <a name="usage"></a>Usage

USAGE:

install notramp package and run:
```sh
notramp (-a | -c | -t) -p PRIMERS -r  READS -g REFERENCE [optional arguments]
```

or download source from github and run from package dir:
```sh
notramp_main.py (-a | -c | -t) -p PRIMERS -r READS -g  REFERENCE [optional arguments]
```

All arguments in detail:
```sh
NoTrAmp is a Tool for read-depth normalization and trimming of amplicon reads generated with long 
read technologies (ONT/PacBio). It can be used in amplicon-tiling approaches to cap coverage of
each amplicon and to trim amplicons to their appropriate length removing barcodes, adpaters and 
primers (if desired) in a single clipping step.

Required arguments:
  -p PRIMERS, --primers PRIMERS
                        Path to primer bed-file (primer-names must adhere to a consistent naming
                         scheme see readme)
  -r READS, --reads READS
                        Path to sequencing reads fasta
  -g REFERENCE, --reference REFERENCE
                        Path to reference (genome)
  -a, --all             Perform read depth normalization by coverage-capping/downsampling first, 
                        then clip the normalized reads. (mut.excl. with -c, -t)
  -c, --cov             Perform only read-depth normalization/downsampling. (mut.excl. with -a, -t)
  -t, --trim            Perform only trimming to amplicon length (excluding primers by default; to 
                        include primers set --incl_prim flag). (mut.excl. with -a, -c)

Optional arguments:
  -h, --help            Print help message and exit
  -o OUT_DIR            Optionally specify a directory for saving of outfiles. If this argument is 
                        not given, out-files will be saved in the directory where the input reads 
                        are located. [default=False]
  -m MAX_COV            Provide threshold for maximum read-depth per amplicon as integer value. 
                        [default=200]
  --incl_prim           Set this flag if you want to include the primer sequences in the trimmed 
                        reads. By default primers are removed together with all overhanging 
                        sequences like barcodes and adapters.
  -s SEQ_TEC            Specify long-read sequencing technology (ont/pb). [default='ont']
  -n NAME_SCHEME        Provide path to json-file containing a naming scheme which is consistently 
                        used for all primers. [default='artic_nCoV_scheme']
  --set_min_len SET_MIN_LEN
                        Set a minimum required length for alignments of reads to amplicon. If this 
                        is not set the min_len will be 0.5*average_amp_len. If amplicon sizes are 
                        relatively homogenous this parameter is not required [default=False]
  --set_max_len SET_MAX_LEN
                        Set a maximum required length for alignments of reads to amplicon. If this 
                        is not set the max_len will be 1.5*average_amp_len. If amplicon sizes are 
                        relatively homogenous this parameter is not required [default=False]
  --set_margins MARGINS
                        Set length of tolerance margins for sorting of mappings to amplicons. 
                        [default=5]
  -v, --version         Print version and exit
  ```

## <a name="depend"></a>Requirements/Dependencies
- Python 3.x
- minimap2
- psutil
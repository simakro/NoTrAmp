![PyPI Version](https://img.shields.io/pypi/v/notramp)
![PyPI monthly Downloads](https://img.shields.io/pypi/dm/notramp.svg?label=PyPI)
[![PyPI total Downloads pepy](https://static.pepy.tech/badge/notramp)](https://pepy.tech/project/notramp)
![bioconda version](https://img.shields.io/conda/v/bioconda/notramp.svg?label=Bioconda)
![bioconda total Downloads biolabel](https://img.shields.io/conda/dn/bioconda/notramp.svg?label=Bioconda)
[![Py versions](https://img.shields.io/pypi/pyversions/notramp.svg)](https://pypi.org/project/notramp)
![GitHub actions standard workflow status](https://img.shields.io/github/actions/workflow/status/simakro/notramp/python-app_test_NoTramp.yml?label=tests)


# NoTrAmp
Normalization and Trimming of amplicon sequencing data

## Table of Contents
- [Installation](#installation)
- [Introduction](#intro)
- [Usage](#usage)
- [Output](#output)
- [Primer naming schemes](#namescheme)
- [Dependencies](#depend)

## <a name="intro"></a>Introduction

NoTrAmp is a Tool for super fast trimming and read-depth normalization of amplicon reads.
It is designed to be used in amplicon-tiling panels (or similar multiplexed amplicon sequencing approaches) to cap coverage of each amplicon and to trim amplicons to their
appropriate length removing barcodes, adpaters and primers (if desired) in a single clipping step.  

Amplicon-tiling schemes are employed to target and amplify specific sequences and enable coverage of longer regions of DNA with small, contiguous segments using overlapping amplicons. 
This approach is particularly useful for detection of mutations, characterization of genetic variation and allows generation of high quality assemblies from low input, fragmented DNA. 
It is frequently utilized for the sequencing of viral genomes and has seen extensive use during the SARS-CoV2 pandemic or during Ebola outbreaks [Citations, links to ARTIC], but is also very useful for exploration of specific genomic loci at high resolution in bacteria or eukaryotes.  

Amplicon-tiling protocols include amplification of the target sequences in separate multiplex PCRs build on (typically) two complementary primer pools.
The performance of individual amplicons in these multiplex PCRs can be vastly different, resulting in large variations of read counts for different regions of the target sequence.
The necessity to accumulate enough reads at weak amplicons usually results in amassing orders of magnitude more reads than required at the more efficient amplicons.
This net overproduction increases the data load and can significantly slow down downstream processes. Additionaly, adapters and barcodes that are attached to DNA fragments during sequencing library preparation, as well as the PCR primers, which could otherwise conceal mutations/variations, need to be removed for downstream processing sequencing.  
NoTrAmp addresses these issues by limiting the read depth at each amplicon to a set count and performs extremely fast one-step trimming, by removing primers, barcodes and adapters in the same clipping operation.

NoTrAmp is suitable for use with both long (e.g. ONT/PacBio) and short reads (e.g Illumina). However, when using reads that are significantly shorter than amplicon sizes, you should adjust the minimum required alignment length using the --set_min_len argument (see below).

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
usage: 
notramp (-a | -c | -t) -p PRIMERS -r READS -g REFERENCE [optional arguments]

NoTrAmp is a Tool for read-depth normalization and trimming of reads in amplicon-tiling approaches. It trims amplicons to their appropriate length removing barcodes, adpaters and
primers (if desired) in a single clipping step and can be used to cap coverage of all amplicons at a chosen value.

Required arguments:
  -p PRIMERS, --primers PRIMERS
                        Path to primer bed-file (primer-names must adhere to a consistent naming scheme see readme)
  -r READS, --reads READS
                        Path to sequencing reads fasta
  -g REFERENCE, --reference REFERENCE
                        Path to reference (genome) fasta file. Must contain only one target sequence. Multiple target sequences are not currently supported.
  -a, --all             Perform read depth normalization by coverage-capping/downsampling first, then clip the normalized reads. (mut.excl. with -c, -t)
  -c, --cov             Perform only read-depth normalization/downsampling. (mut.excl. with -a, -t)
  -t, --trim            Perform only trimming to amplicon length (excluding primers by default; to include primers set --incl_prim flag). (mut.excl. with -a, -c)

Optional arguments:
  -h, --help            Print help message and exit
  -o OUT_DIR            Optionally specify a directory for saving of outfiles. If this argument is not given, out-files will be saved in the directory where the input reads are located.
                        [default=False]
  -m MAX_COV            Provide threshold for maximum read-depth per amplicon as integer value. [default=200]
  --incl_prim           Set this flag if you want to include the primer sequences in the trimmed reads. By default primers are removed together with all overhanging sequences like
                        barcodes and adapters.
  -s SEQ_TEC            Specify long-read sequencing technology (ont/pb). [default=ont]
  -n NAME_SCHEME        Provide path to json-file containing a naming scheme which is consistently used for all primers.[default=artic_nCoV_scheme_v5.3.2]
  --set_min_len SET_MIN_LEN
                        Set a minimum required length for alignments of reads to amplicons. If this is not set the min_len will be 0.8*shortest_amp_len. When using reads that are
                        shorter than amplicon sizes use this argument to adjust. For long reads this option is usually not required.
  --set_max_len SET_MAX_LEN
                        Set a maximum allowed length for alignments of reads to amplicon. If this is not set the max_len will be 1.2*longest_amp_len. The default setting normally
                        doesn't need to be changed.
  --set_margins MARGINS
                        Set length of tolerance margins for sorting of mappings to amplicons. [default=5]
  --figures [FIGURES]   Set to generate figures of input and output read_counts. Available for --all and --cov modes. You can optionally provide a value to draw a red helper line in the
                        output read plot, showing a threshold, e.g. min. required reads. [default=False; default_threshold=20]
  --fastq               Set this flag to request output in fastq format. By default output is in fasta format. Has no effect if input file is fasta.
  --split               Set this flag to request output of capped, untrimmed reads split to amplicon specific files (can be a lot).
  --selftest            Run a selftest of NoTrAmp using included test-data. Overrides all other arguments and parameters. Useful for checking how NoTrAmp runs in your environment.
  -v, --version         Print version and exit
  ```


## <a name="output"></a>Output
NoTrAmp by default generates a separate read file as output for capping and trimming. Capped untrimmed reads are contained in a file ending on ".cap.fasta". 
Clipped reads are stored in a file ending on ".clip.fasta". If both capping and trimming were selected, trimmed versions of the capped reads, are written to  "YourFileName.cap.clip.fasta". If quality information is required downstream 
in your workflow, you can request output in fastq format, by setting the --fastq flag. It is recommended that quality control and 
filtering of data is performed before running NoTramp.  
Additionaly a log-file ("notramp.log") is generated, that also contains detailed information about processed and selected reads, read coverage/amplicon and trimmed bases. 
A visual representation (see below) of input and output reads can also be requested by setting the --figures flag.
<p align="center">
  <img src="https://github.com/simakro/NoTrAmp/blob/main/notramp/resources/notramp_amplicon_coverage_large_test_data_thr30x.png" width="600" height="auto" align="center"/>
</p>
Upper plot: Input reads (before)  
Lower plot: Capped output reads (after)

## <a name="namescheme"></a>Primer naming schemes
NoTramp requires primers in multiplex amplicon tiling panels to follow a consistent scheme.
A primer name must consist of a number of fields delineated by a defined separator.
Two fields are mandatory: amplicon-number and primer-position(FW/REV).
NoTramp offers a high degree of flexibility for naming of primers by allowing the user to specify their own naming scheme.
This is provided as json file, which defines the information contained in each field of the primer name.
The convention specified in this naming scheme must be followed consistently throughout the entire panel.

Minimal primer naming scheme containing all required key/value pairs (see also "notramp/resources/minimal_scheme.json"):
```sh
{
    "sep": "_", 
    "min_len": 2, 
    "max_len": 2, 
    "amp_num": 0, 
    "position": 1, 
    "fw_indicator": "fw",
    "rev_indicator": "rev"
}
```
|"Keyword"|"Description"|
|:---|:---|
|"sep" | Separator used to delineate fields in primer name|
|"min_len" | minimal number of fields in primer names [must be an int]|
|"max_len" | maximum number of fields in primer names (if no alternative primer are in the panel, the same as min_len) [must be an int]|
|"amp_num" | 0 based index of the field containing the amplicon-number [must be an int]|
|"position" | 0 based index of the field containing information on primer position in the amplicon (e.g. left/right or fw/rev) [must be an int]|
|"fw_indicator" | indicator used to identify directionality of the primer; can be anything as long as consistent; typical indicators: "fw", "FW", "left", "start", "+"|
|"rev_indicator" | indicator used to identify directionality of the primer; can be anything as long as consistent; typical indicators: "rev", "REV", "right", "end", "-"|

Examples for primers named after a minimal scheme:
```sh
  1_fw 
  1_rev
  2_fw
  2_rev
  3_fw
  3_rev
  ...
```

Use of such naming schemes enforces primers to consistently have names with the same number of fields, which always carry the same kind of information.
However, there can be one exception from the same number of fields rule.

Sometimes alternative primer pairs for the "same" amplicon are used to boost underperforming amplicons.
These are typically primer sequences that are shifted just a couple of bases to the left or right of the original/primary primer pair,
in the hope of increasing the yield for this region. The products generated under involvement of such alternative primers
cover the same region of the target sequence, with only slight deviations on the rims/fringes/edges/corners, in the parts overlapping
with neighbouring amplicons.
Alternative primers are required to carry the same amplicon-number as primary ones, but must have some indicator to be distinguished from those.
NoTramp accounts for this irregularity by allowing for alternative primers to be labeled with an "alt" field, which must be supplied as a postfix.
Therefore in panels containing alternative primers, max_len must min_len + 1.
If your panel includes alternative primers, you have to add the "alt" keyword to your naming scheme, with its corresponding index as value, which must be the last possible field (max_len-1).
The alt indicator you use in the name can be anything you want (e.g. "v2", "2", "b", "alt" etc.), but the keyword in the json dict must be "alt".

In addition to the optional "alt" field, any number of custom fields can be added to the naming scheme, if your primer contains additionals field.
These could typically be a common name or pool indicator. These fields can be added for completeness, but don't have an effect on the inner workings of NoTramp.

Generic primer scheme (see also "notramp/resources/generic_scheme.json"):
```sh
{
    "sep": "_", 
    "min_len": 3, 
    "max_len": 4, 
    "root_name": 0,
    "amp_num": 1, 
    "position": 2, 
    "alt": 3, 
    "fw_indicator": "FW",
    "rev_indicator": "REV"
}
```
|"Keyword"|"Description"|
|:---|:---|
|"alt" |0 based index of the field containing the alternative primer indicator|   
|"root_name"|custom field used for a common name |

Examples for primers complying to the generic scheme above:
```sh
Target-Gene_1_FW
Target-Gene_1_REV
Target-Gene_2_FW
Target-Gene_2_REV
Target-Gene_2_FW_v2
Target-Gene_2_REV_v2
Target-Gene_3_FW
Target-Gene_3_REV
...
```

<!-- use v5 scheme here instead show scheme and examples (excerpt from artic-ncov2019 v5 primer panel)-->
If no custom scheme is supplied, the current default is the Artic SARS-CoV2 v5.3.2 primer scheme:
Generic primer scheme (see also "notramp/resources/artic_nCoV_scheme_v5.json"):
```sh
{
    "sep": "_", 
    "min_len": 4, 
    "max_len": 5, 
    "root_name": 0,
    "amp_scale": 1, 
    "amp_num": 2, 
    "position": 3,
    "iteration": 4,
    "alt": 5, 
    "fw_indicator": "LEFT", 
    "rev_indicator": "RIGHT"
}
```


## <a name="depend"></a>Requirements/Dependencies
required:   
- Python >= 3.7
- minimap2  

recommended:   
- psutil  

optional:   
- matplotlib (for figures if desired)
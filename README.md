# NoTrAmp
Normalization and Trimming of long-read (ONT, PB) amplicon sequencing data


NoTrAmp is a Tool for read-depth normalization and trimming of amplicon reads generated with long read technologies (ONT/PacBio).
It can be used in amplicon-tiling approaches to cap coverage of each amplicon and to trim amplicons to their
appropriate length removing barcodes, adpaters and primers (if desired) in a single clipping step.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Primer naming schemes](#namescheme)
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
|---|---|
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
"alt" = 0 based index of the field containing the alternative primer indicator   
"root_name" = custom field used for a common name   

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
- Python 3.x
- minimap2
- psutil
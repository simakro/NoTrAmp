import argparse
import logging
import subprocess
import random
import sys
import os



    # mm2_paf = sys.argv[1]
    # primer_bed = sys.argv[2]
    # reads = sys.argv[3]
    # clipped_out = reads + "_primer_clipped"

def get_arguments():
    parser = argparse.ArgumentParser(description="This is a Tool to detect the number and location of integrations of a"
                                                 " transgenic construct or other inserted DNA sequence in a genome from"
                                                 " long-read whole-genome sequencing data.", add_help=True)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-c", "--construct", required=True, help="Path to fasta file containing the sequence"
                               " of the integrating construct")
    required_args.add_argument("-g", "--genome", required=True, help="Path to reference genome")
    required_args.add_argument("-o", "--outdir", required=True, help="Path for output directory")
    
    read_input = required_args.add_mutually_exclusive_group(required=True)
    read_input.add_argument("-f", "--reads_fasta", help="Path to a single fasta/q file containing "
                            "the sequencing reads")
    read_input.add_argument("-d", "--reads_dir", help="Path to a single fasta/q file containing the "
                            "sequencing reads")

    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("--min_read_supp", dest='min_supp', default=3, type=ranged_typechk(1,1000,int),
                        help="Set threshold for minimum read support for integration validation")
    optional_args.add_argument("--spec", dest='species', default=None, 
                        help="Select species (e.g. \"human\") to enable output of species specific ideogram and annotations")
    optional_args.add_argument("--ava_clust", dest='cluster_reads', default=False, action='store_true', 
                        help="Cluster reads based on all-vs-all alignment. Comparison of read-clusters to integration-sites "
                        "gained by  alignment of reads to the reference genome may enable identification of potential "
                        "duplicates that could arise from integration into highly repetetive genomic loci.")
    optional_args.add_argument("--verbosity", dest='verbosity', default=1, type=ranged_typechk(-1,100,int), 
                        help="Adjust level of output to the terminal")
    optional_args.add_argument("--file_check", dest='check_reads', default="quick",
                        help="Adjust thoroughness of file-check; quick/extensive [quick]")
    optional_args.add_argument("--include_Nx", dest='Nx', default=98, type=ranged_typechk(1,100,float),
                        help="Choose cutoff for contigs to be included in figures in graphical report. "
                        "Contigs comprising at least x percent of bases in RefGenome [98] will be included."
                        "This parameter can be helpful to avoid graphs becoming overcrowded when working with "
                        "reference genomes that contain large amounts of alternative loci or a high number of "
                        "miniscule unplaced contigs. Don't worry, this applies only to the figures 1-3 in the "
                        "pdf-report. Thus, no sequences will be ignored during integration search, i.e. there is no "
                        "risk of missing out on integration sites. Species specific ideogram is not affected.")
    optional_args.add_argument("--bait", dest='bait', default=False,
                        help="This argument can be supplied in addition to the construct/-c argument. "
                        "It is required if the integration construct shares partial similarity with the reference sequence. "
                        "Intloc will inform you if the use of a bait sequence is necessary. "
                        "Provide the path to a fasta file containing only part of the integrating construct as bait for "
                        "identification of integration evidence reads."
                        "The bait sequence must not share similarity with reference. [False]")
    optional_args.add_argument("--cleanup", dest='cleanup_aux', default=True,
                        help="Cleanup auxiliary files [True]. Relevant for debugging only.")
    optional_args.add_argument("--override", dest='override', default=False, action='store_true', help=argparse.SUPPRESS)
    optional_args.add_argument("--gui", dest='gui', default=False, action='store_true', help=argparse.SUPPRESS)
    optional_args.add_argument("--testing", dest='testing', default=False, action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args()

    return args


def log_sp_error(error, message):
    logging.error(message)
    print("Check intlog.log for error details")
    logging.error(error.stdout.decode('utf-8'))
    logging.info("Exiting integration locator")
    sys.exit()



class Mapping:
    def __init__(
        self,
        qname,
        qlen,
        qstart,
        qend,
        samestrand,
        tname,
        tlen,
        tstart,
        tend,
        matches,
        total_bp,
        qual,
        kwattr,
    ):
        self.qname = qname
        self.qlen = int(qlen)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.samestrand = self.eval_strand(samestrand)
        self.tname = tname
        self.tlen = int(tlen)
        self.tstart = int(tstart)
        self.tend = int(tend)
        self.matches = int(matches)
        self.total_bp = int(total_bp)
        self.qual = int(qual)
        self.kwattr = kwattr
        self.gen_kw_attr()

    def gen_kw_attr(self):
        kwattr_dict = {kw.split(":")[0]: kw.split(":")[-1] for kw in self.kwattr}
        for key in kwattr_dict:
            self.__dict__[key] = kwattr_dict[key]

    def eval_strand(self, strand_info):
        if strand_info == "+":
            return True
        else:
            return False

class Primer:
    def __init__(self, contig, start, end, name, score, strand):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand
        self.type = "primary"
        self.amp_no, self.pos = self.get_name_infos()

    def get_name_infos(self):
        ls = self.name.split("_")
        if len(ls) == 4:
            self.type = "alt"
            self.alt_name = ls[3]
        return ls[1], ls[2]

class Amp:
    def __init__(self, amp_no, primers):
        self.name = int(amp_no)
        self.primers = primers
        self.start = min([primer.start for primer in primers])
        self.end = max([primer.end for primer in primers])
        self.max_len = self.end - self.start
        self.fwp_boundary = max(
            [prim for prim in primers if prim.pos == "LEFT"], key=lambda x: x.end
        ).end
        self.revp_boundary = min(
            [prim for prim in primers if prim.pos == "RIGHT"], key=lambda x: x.start
        ).start
        self.mappings = dict()
        self.reads_dct = dict()
    
    
    def random_sample(self, cutoff):
        if len(self.read_names) > cutoff:
            self.selected = random.choices(self.read_names, k=cutoff)
        else:
            self.selected = self.read_names

    def primer_clipping_all(self):
        clip_ct_left = 0
        clip_ct_right = 0
        for read in self.reads_dct:
            try:
                mapping = self.mappings[read]
                left, right = self.reads_dct[read].clip_primers(
                    self.fwp_boundary, self.revp_boundary, mapping
                )
                clip_ct_left += left
                clip_ct_right += right
            except KeyError as e:
                print(f"KeyError in primer_clipping_all: {e}")
        return clip_ct_left, clip_ct_right


def create_primer_objs(primer_bed):
    with open(primer_bed, "r") as bed:
        primers = list()
        for line in bed:
            prim = Primer(*line.strip().split("\t"))
            primers.append(prim)
    return sorted(primers, key=lambda x: x.end)


def generate_amps(primers):
    amp_nums = set([primer.amp_no for primer in primers])
    amps = list()
    for num in amp_nums:
        ao = Amp(num, [primer for primer in primers if primer.amp_no == num])
        amps.append(ao)
    return sorted(amps, key=lambda x: x.name)


def name_out_paf(reads, reference, mod_name):
    read_dir, reads_file = os.path.split(reads)
    reads_name = reads_file.split(".")[:-1]
    ref_name = os.path.split(reference)[1].split(".")[:-1]
    paf_name = f"{reads_name}_mapto_{ref_name}.{mod_name}.paf"
    return os.path.join(read_dir, paf_name)


def map_reads(reads, reference, out_paf, seq_tech="map-ont"):
    try:
        subprocess.run(["minimap2", "-x", seq_tech, reference, reads, "-o", out_paf, "--secondary=no"], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        log_sp_error(e, "Mapping of reads to reference failed")
    return out_paf


def filter_read_mappings(mappings):
    mappings = [m for m in mappings if 300 < m.qlen < 600]
    mappings = [m for m in mappings if m.qual == 60]
    return mappings


def create_read_mappings(mm2_paf):
    with open(mm2_paf, "r") as paf:
        map_dct = defaultdict(list)
        for line in paf:
            if len(line.strip()) > 0:
                ls = line.strip().split("\t")
                mapping = Mapping(*ls[:12], ls[12:])
                map_dct[mapping.qname].append(mapping)
        mult_maps = {n: ml for n, ml in map_dct.items() if len(ml) > 1}
        mappings = [m for k, l in map_dct.items() for m in l]
        mappings = filter_read_mappings(mappings)
        mono_mappings = [m for m in mappings if m.qname not in mult_maps]
        dual_mappings = {k: v for k, v in mult_maps.items() if len(v) == 2}
        incl_max = [
            max(dual_mappings[mname], key=lambda x: x.matches)
            for mname in dual_mappings
        ]
        incl_max = filter_read_mappings(incl_max)
        mono_mappings.extend(incl_max)
        mappings = mono_mappings
    return sorted(mappings, key=lambda x: (x.tend, x.tstart))




if __name__ == "__main__":

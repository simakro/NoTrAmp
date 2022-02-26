import logging
import logging.config
# from time import perf_counter
import random
# from collections import defaultdict

class Amp:
    """class for generation of amplicon objects, able to extract, store, \
    handle and process data of attached Primer, Mapping and Read objects"""

    def __init__(self, amp_no, primers):
        self.name = int(amp_no)
        self.primers = primers
        self.start = min([primer.start for primer in primers])
        self.end = max([primer.end for primer in primers])
        self.mappings = {}
        self.reads_dct = {}
        self.selected = []

    def fwp_boundary(self):
        """get inner boundary of forward primer"""
        return max(
            [prim for prim in self.primers if prim.pos == prim.scheme["fw_indicator"]],
            key=lambda x: x.end
        ).end

    def revp_boundary(self):
        """get inner boundary of reverse primer"""
        return min(
            [prim for prim in self.primers if prim.pos == prim.scheme["rev_indicator"]],
            key=lambda x: x.start
        ).start

    def get_max_len(self):
        """get maximum length of amplicon"""
        return self.end - self.start

    def random_sample(self, cutoff):
        """select a random subsample from all reads"""
        if len(self.reads_dct) > cutoff:
            read_names = list(self.reads_dct.keys())
            self.selected = random.choices(read_names, k=cutoff)
        else:
            self.selected = list(self.reads_dct.keys())

    def trim_clip_all(self, incl_prim):
        """trim all attached reads"""
        clip_ct_left = 0
        clip_ct_right = 0
        for rname, read in self.reads_dct.items():
            try:
                mapping = self.mappings[rname]
                if incl_prim:
                    left, right = read.trim_to_amp(
                        self.start, self.end, mapping
                    )
                else:
                    left, right = read.clip_primers(
                        self.fwp_boundary(), self.revp_boundary(), mapping
                    )
                clip_ct_left += left
                clip_ct_right += right
            except KeyError as err:
                logging.exception(err)
        return clip_ct_left, clip_ct_right


class Read:
    """class for generation of read objects, and clipping/trimming"""

    def __init__(self, header, seq, fastq=False):
        self.header = header
        self.fastq = fastq
        self.name = self.header.split("@")[1] if self.fastq else self.header.split(">")[1]
        self.seq = seq

    @staticmethod
    def non_neg(num):
        """return 0 if result < 0"""
        return num if num >= 0 else 0

    def trim_to_amp(self, amp_start, amp_end, mapping):
        """trim to amplicon boundaries including primers"""
        qlen, samestrand = mapping.qlen, mapping.samestrand
        qstart, qend = mapping.qstart, mapping.qend
        tstart, tend = mapping.tstart, mapping.tend

        if samestrand:
            clip_left = 0
            ldiff = abs(amp_start - tstart)
            if tstart >= amp_start:
                clip_left = Read.non_neg(qstart - ldiff)
            else:
                clip_left = qstart + ldiff

            clip_right = qlen
            rdiff = abs(tend - amp_end)
            if tend <= amp_end:
                clip_right = qend + rdiff
            else:
                clip_right = Read.non_neg(qend - rdiff)
        else:
            clip_left = 0
            ldiff = abs(tend - amp_end)
            if tend <= amp_end:
                clip_left = qstart - ldiff
            else:
                clip_left = qstart + ldiff

            clip_right = qlen
            rdiff = abs(amp_start - tstart)
            if tstart >= amp_start:
                clip_right = qend + rdiff
            else:
                clip_right = qend - rdiff

        self.seq = self.seq[clip_left:clip_right]
        return clip_left, qlen - clip_right

    def clip_primers(self, fwp_boundary, revp_boundary, mapping):
        """trim to amplicon boundaries excluding primers"""
        qlen, samestrand = mapping.qlen, mapping.samestrand
        qstart, qend = mapping.qstart, mapping.qend
        tstart, tend = mapping.tstart, mapping.tend

        if samestrand:
            clip_left = 0
            if tstart <= fwp_boundary:
                add_clip = fwp_boundary - tstart
                clip_left = qstart + add_clip
            else:
                ldiff = tstart - fwp_boundary
                if qstart >= ldiff:
                    clip_left = qstart - ldiff

            clip_right = qlen
            if tend >= revp_boundary:
                sub_clip = tend - revp_boundary
                clip_right = qend - sub_clip
            else:
                rdiff = revp_boundary - tend
                if qlen - qend >= rdiff:
                    clip_right = qend + rdiff

        else:
            clip_right = qlen
            if tstart <= fwp_boundary:
                add_clip = fwp_boundary - tstart
                clip_right = qend - add_clip
            else:
                rdiff = tstart - fwp_boundary
                if qlen - qend >= rdiff:
                    clip_right = qend + rdiff

            clip_left = 0
            if tend >= revp_boundary:
                sub_clip = tend - revp_boundary
                clip_left = qstart + sub_clip
            else:
                ldiff = revp_boundary - tend
                if qstart >= ldiff:
                    clip_left = qstart - ldiff

        self.seq = self.seq[clip_left:clip_right]
        return clip_left, qlen - clip_right


class Mapping:
    """class for storage and handling of mapping information"""

    def __init__(self, qname, posattr, kwattr):
        self.qname = qname
        self.posattr = posattr
        self.kwattr = kwattr
        self.gen_pos_attr()
        self.gen_kw_attr()

    @staticmethod
    def eval_strand(strand_info):
        """evaluate strand info"""
        return True if strand_info == "+" else False

    def gen_pos_attr(self):
        """generate class attributes from positional info in mapping output"""
        posattr_dict = {
            "qlen": int,
            "qstart": int,
            "qend": int,
            "samestrand": Mapping.eval_strand,
            "tname": str,
            "tlen": int,
            "tstart": int,
            "tend": int,
            "matches": int,
            "total_bp": int,
            "qual": int
        }
        zip_attr = zip(posattr_dict.keys(), self.posattr)
        for key, val in zip_attr:
            self.__dict__[key] = posattr_dict[key](val)

    def gen_kw_attr(self):
        """generate class attributes from key-worded entries in mapping output"""
        kwattr_dict = {kw.split(":")[0]: kw.split(":")[-1] for kw in self.kwattr}
        for key in kwattr_dict:
            self.__dict__[key] = kwattr_dict[key]


class Primer:
    """class for storage and handling of primer information"""

    def __init__(self, start, end, name, add_info, naming_scheme):
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.add_info = add_info
        self.type = "primary"
        self.scheme = naming_scheme
        self.get_name_infos()

    def get_name_infos(self):
        """extract information from primer name"""
        lsp = self.name.split(self.scheme["sep"])
        if len(lsp) == self.scheme["max_len"]:
            self.type = "alt"
            self.alt_name = lsp[self.scheme["alt"]]
        self.__dict__["amp_no"] = lsp[self.scheme["amp_num"]]
        self.__dict__["pos"] = lsp[self.scheme["position"]]

    def get_len(self):
        """get primer length"""
        return self.end - self.start + 1
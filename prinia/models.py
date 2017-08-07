from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip)

__author__ = 'ahbbollen'

from datetime import datetime


class Variant(object):
    """
    This class describes a variant returned by LOVD
    """

    miracle_id = ''
    datum = ''
    refseq_build = ''
    chromosome = ''
    genomic_id_ncbi = ''
    position_g_start = ''
    position_g_end = ''
    variant_on_genome = ''
    is_present_father = ''
    is_present_mother = ''
    gene_id = ''
    gene_name = ''
    transcript_id_ncbi = ''
    variant_on_transcript_dna = ''
    variant_on_transcript_rna = ''
    variant_on_transcript_protein = ''
    allele = ''
    variant_on_genome_origin = ''
    in_gene_panel = ''
    confirm_in_lab = True

    @classmethod
    def from_dict(cls, dictionary, comments):
        results = []
        for n in range(comments["n_line"]):
            x = Variant()

            x.miracle_id = comments['id_miracle']
            if "active_gene_panel" in comments:
                x.datum = comments['active_gene_panel']
            else:
                x.datum = str(datetime.utcnow())

            for k in dictionary.keys():
                try:
                    setattr(x, k, dictionary[k][n])
                except:
                    pass

            x.variant_on_genome = dictionary["VariantOnGenome/DNA"][n]
            x.variant_on_transcript_dna = dictionary["VariantOnTranscript/DNA"][n]
            x.variant_on_transcript_rna = dictionary["VariantOnTranscript/RNA"][n]
            x.variant_on_transcript_protein = dictionary["VariantOnTranscript/Protein"][n]
            x.variant_on_genome_origin = dictionary["VariantOnGenome/Genetic_origin"][n]

            if x.confirm_in_lab == "0":
                x.confirm_in_lab = False
            elif x.confirm_in_lab == "1":
                x.confirm_in_lab = True

            results.append(x)

        return results


class Primer(object):
    """
    This class represents a Primer we produce
    """

    chromosome = ''
    position = ''
    left = ''
    right = ''
    left_pos = ''
    right_pos = ''
    left_len = ''
    right_len = ''
    left_gc = ''
    right_gc = ''
    fragment_sequence = ''
    contains_SNP = False
    contains_freq_SNP = False
    snp_freq = 0.0
    blathits_left = []
    blathits_right = []
    left_name = ''
    right_name = ''

    def __init__(self, **kwargs):
        for kw, kv in kwargs.items():
            setattr(self, kw, kv)

    @classmethod
    def from_p3(cls, forward, reverse, fragment=None, chromosome=None, position=None):
        # line looks like this:
        #    0 CAGCACTGCTTGAGGGGAA                0 19  0 57.89 59.926  0.00  0.00 36.40  1.074
        f_contents = [x for x in forward.strip().split(" ") if x]
        r_contents = [x for x in reverse.strip().split(" ") if x]
        primer = cls()
        primer.chromosome = chromosome
        primer.position = position
        primer.left = f_contents[1]
        primer.right = r_contents[1]
        primer.left_pos = f_contents[2]
        primer.right_pos = r_contents[2]
        primer.left_len = f_contents[3]
        primer.right_len = r_contents[3]
        primer.left_gc = f_contents[5]
        primer.right_gc = r_contents[5]
        primer.fragment_sequence = fragment
        return primer


class Region(object):

    def __init__(self, *args, **kwargs):
        self.chr = kwargs['chromosome']
        self.acc_nr = kwargs['acc_nr']

        self.start = kwargs['start']
        self.stop = kwargs['stop']

        self.padding_left = kwargs['padding_left']
        self.padding_right = kwargs['padding_right']
        if "ref" in kwargs:
            self.ref = kwargs['ref']
        else:
            self.ref = 'NA'

        if "other" in kwargs:
            self.other_information = kwargs['other']
        else:
            self.other_information = "NA"

    def __len__(self):
        return int(self.stop) - int(self.start)

    def size(self, padded=False):
        if padded:
            return int(self.stop_w_padding) - int(self.start_w_padding)
        else:
            return self.__len__()

    @property
    def start_w_padding(self):
        return int(self.start) - int(self.padding_left)

    @property
    def stop_w_padding(self):
        return int(self.stop) + int(self.padding_right)

    @classmethod
    def from_dict(cls, dictionary):
        for k in dictionary.keys():
            setattr(cls, k, dictionary[k])

    @classmethod
    def from_variant(cls, variant, padding_l=0, padding_r=0):
        return cls(chromosome=variant.chromosome, acc_nr=variant.genomic_id_ncbi,
                   start=variant.position_g_start, stop=variant.position_g_end,
                   padding_left=padding_l,
                   padding_right=padding_r)

    @classmethod
    def from_bed(cls, line, reference, padding_l=0, padding_r=0):
        contents = line.strip().split("\t")

        if len(contents) > 3:
            other = ".".join(contents[3:])
        else:
            other = "NA"

        return cls(chromosome=contents[0], acc_nr='NA',
                   start=int(contents[1]), stop=int(contents[2]), padding_left=padding_l,
                   padding_right=padding_r, ref=reference, other=other)

    @classmethod
    def cut(cls, region, max_size, padded=False):
        """
        Cut object to maximum size
        Will cut on RIGHT side
        :param region: Region to cut
        :param max_size: maximum size wished
        :param padded: whether to use padded size
        :return: (cut, remainder)
        """

        if region.size(padded=padded) <= max_size:
            return None

        if not padded:
            remainder_start = region.start + max_size + 1
            cut_stop = region.start + max_size
        else:
            remainder_start = region.start_w_padding + max_size + 1
            cut_stop = region.start_w_padding + max_size
            if cut_stop < region.stop:
                raise ValueError("Stop would be before target")

        remainder = cls(chromosome=region.chr, acc_nr=region.acc_nr,
                        start=remainder_start,
                        stop=region.stop,
                        padding_left=region.padding_left,
                        padding_right=region.padding_right,
                        ref=region.ref,
                        other=region.other_information)
        cut = cls(chromosome=region.chr, acc_nr=region.acc_nr,
                  start=region.start,
                  stop=cut_stop,
                  padding_left=region.padding_left,
                  padding_right=region.padding_right,
                  ref=region.ref,
                  other=region.other_information)

        return cut, remainder

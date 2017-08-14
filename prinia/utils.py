from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip)
__author__ = 'ahbbollen'

from datetime import datetime
import hashlib
import sys

import vcf

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import NucleotideAlphabet


def _is_vcf_version_at_least_0_6_8():
    """
    The behaviour of vcfReader.fetch changed significantly from version 0.6.8 onwards
    :return: boolean
    """
    major, minor, patch = vcf.VERSION.split(".")
    if int(major) == 0 and int(minor) == 6 and int(patch) >= 8:
        return True
    elif int(major) == 0 and int(minor) > 6:
        return True
    elif int(major) > 0:
        return True
    return False


NEW_VCF = _is_vcf_version_at_least_0_6_8()

def calc_gc(sequence):
    """
    Calculates GC percentage from DNA sequence
    Does not consider IUPAC ambiguity codes
    :param sequence:
    :return: float
    """
    if len(sequence) == 0:
        raise ValueError("Sequence must have minimum length 1")
    gc = 0
    for char in sequence:
        if char.upper() == 'G' or char == 'C':
            gc += 1
    return float(gc)/len(sequence) * 100.0


class NoPrimersException(Exception):
    pass


def datehash():
    a = datetime.utcnow()
    st = "{:f}".format((a-datetime(1970, 1, 1)).total_seconds())
    if sys.version_info[0] == 2:
        hash = hashlib.md5(st).hexdigest()
    elif sys.version_info[0] == 3:
        hash = hashlib.md5(st.encode("utf-8")).hexdigest()
    else:
        raise NotImplementedError
    return hash[:5]


def primer_to_seq_record(primer):
    """
    Create seqrecord from a primer
    :param primer: primer
    :return: 2-tuple of (SeqRecord forward, SeqRecord reverse)
    """
    id_str = datehash()
    forward_id = "{0}/1".format(id_str)
    reverse_id = "{0}/2".format(id_str)
    forward_seq = Seq(primer.left, alphabet=NucleotideAlphabet)
    reverse_seq = Seq(primer.right, alphabet=NucleotideAlphabet)
    forward_q = [40 for _ in range(len(primer.left))]
    reverse_q = [40 for _ in range(len(primer.right))]
    forward_r = SeqRecord(forward_seq, id=id_str, description=forward_id)
    forward_r.letter_annotations['phred_quality'] = forward_q
    reverse_r = SeqRecord(reverse_seq, id=id_str, description=reverse_id)
    reverse_r.letter_annotations['phred_quality'] = reverse_q
    return forward_r, reverse_r


def generate_fastq_from_primers(primers, forward_path, reverse_path):
    """
    Generate paired-end fastq files from list of primers
    :param primers:
    :return: 2-tuple of (path_R1, path_R2)
    """
    forward_handle = open(forward_path, "w")
    reverse_handle = open(reverse_path, "w")

    seqs = [primer_to_seq_record(x) for x in primers]
    forwards = [s[0] for s in seqs]
    reverses = [s[1] for s in seqs]

    SeqIO.write(forwards, forward_handle, "fastq")
    SeqIO.write(reverses, reverse_handle, "fastq")

    forward_handle.close()
    reverse_handle.close()

    return forward_path, reverse_path




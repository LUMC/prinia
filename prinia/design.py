from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip)
__author__ = 'ahbbollen'

from tempfile import TemporaryFile, NamedTemporaryFile
from subprocess import check_call
import os
import re
import warnings
from time import sleep

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

from fastools.fastools import get_reference
from pyfaidx import Fasta

from .models import Primer, Region, BlatLine
from .utils import NoPrimersException, calc_gc, NEW_VCF

PRIMER3_SCRIPT = os.path.join(os.path.join(os.path.dirname(__file__) ,"static"), 'getprimers.sh')


def get_sequence_fasta(region, reference=None, padding=True):
    ref = Fasta(reference)
    if "chr" not in list(ref.keys())[0] and "chr" in region.chr:
        chrom = region.chr.split("chr")[1]
    elif "chr" not in region.chr and "chr" in list(ref.keys())[0]:
        chrom = "chr" + region.chr
    else:
        chrom = region.chr

    if not padding:
        return ref[chrom][region.start:region.stop].seq
    else:
        return ref[chrom][region.start_w_padding:region.stop_w_padding].seq


def run_primer3(sequence, region, padding=True,
                primer3_script=PRIMER3_SCRIPT, product_size="200-450",
                n_primers=4, prim3_exe=None):
    if padding:
        target_start = region.padding_left
        target_len = len(sequence) - region.padding_left - region.padding_right
    else:
        target_start = 1
        target_len = len(sequence)

    target = ",".join(map(str, [target_start, target_len]))
    opt_size = str(sum(map(int, product_size.split("-"))) / 2)

    args = [primer3_script, product_size, target, sequence, target, opt_size, prim3_exe]
    retval = check_call(args)
    if retval != 0:
        raise ValueError("Primer3 crashed")

    # now read example.for and example.rev

    with open("example.for", "rb") as forward, open("example.rev", "rb") as reverse:
        forwards = _sanitize_p3(forward)
        reverses = _sanitize_p3(reverse)

    forwards, reverses = _get_shortest(forwards, reverses, n_primers)

    if len(forwards) == 0 or len(reverses) == 0:
        raise NoPrimersException("No acceptable primers could be found. Try increasing the padding")

    primers = [Primer.from_p3(x, y, sequence, region.chr, region.start) for x, y in zip(forwards, reverses)]
    return primers


def blat_primers(primers, blat_exe=None, ref=None):
    tmp = NamedTemporaryFile()
    records = []
    for primer in primers:
        records += [SeqRecord(Seq(primer.left, generic_dna), primer.left)]
        records += [SeqRecord(Seq(primer.right, generic_dna), primer.right)]

    SeqIO.write(records, tmp, "fasta")

    tmp.seek(0)
    out = NamedTemporaryFile()
    args = [blat_exe, '-stepSize=5', '-repMatch=2253', '-minScore=0', '-minIdentity=0', ref, tmp.name, out.name]
    retval = check_call(args)
    if retval != 0:
        raise ValueError("Blat crashed")

    blat_out = out.readlines()
    out.close()
    tmp.close()
    return blat_out


def find_best(primers, blat_out, accept_snp=False, field=None, max_freq=None, dbsnp=None):
    # should match minimally 90%
    # minimum amount of hits
    # correct strand
    sanitized = _sanitize_blat(blat_out)
    primers = _primers_with_blatlines(primers, sanitized)

    # check for existence of correct strand
    n_prims = []
    for primer in primers:
        any_left = any([x.strand == "+" for x in primer.blathits_left])
        any_right = any([x.strand == "-" for x in primer.blathits_right])
        if any_left and any_right:
            n_prims.append(primer)
    primers = n_prims

    # check for match percentage
    n_prims = []
    for primer in primers:
        any_left = any([float(x.Q_size)/float(x.match) >= 0.9 for x in primer.blathits_left])
        any_right = any([float(x.Q_size)/float(x.match) >= 0.9 for x in primer.blathits_right])
        if any_left and any_right:
            n_prims.append(primer)
    primers = n_prims

    blat_lens = [len(x.blathits_left) + len(x.blathits_right) for x in primers]
    # get primers with least hits
    primers = [x for x in primers if (len(x.blathits_right) + len(x.blathits_left)) == min(blat_lens)]

    # get correct chromosome and location
    for primer in primers:
        if "chr" in primer.chromosome and "chr" not in primer.blathits_left[0].T_name:
            chrom = primer.chromosome.split("chr")[1]
        elif "chr" not in primer.chromosome and "chr" in primer.blathits_left[0].T_name:
            chrom = "chr" + primer.chromosome
        else:
            chrom = primer.chromosome
        the_lines_left = [x for x in primer.blathits_left if x.T_name == chrom]
        the_lines_right = [x for x in primer.blathits_right if x.T_name == chrom]
        if len(the_lines_left) == 0 or len(the_lines_right) == 0:
            raise NoPrimersException("No primers could be detected")
        primer.left_pos = int(the_lines_left[0].T_start)
        primer.right_pos = int(the_lines_right[0].T_start)

    # find snps
    if not accept_snp and max_freq:
        primers = [find_snps(x, dbsnp, field) for x in primers]
        primers = [x for x in primers if x.snp_freq <= max_freq]

    if len(primers) == 0:
        raise NoPrimersException("No suitable primers could be detected")
    return primers[0]


def _freq_in_query(query, field):
    """
    Get amount of records and frequencies for a query (= reader iterator)
    :param query: VCF reader iterator
    :param field: the field to fetch
    :return:
    """
    n = 0
    freqs = []
    for record in query:
        n += 1
        try:
            freqs += list(map(float, record.INFO[field]))
        except ValueError:
            freqs += [0]
        except KeyError:
            freqs += [0]
    return n, freqs


def find_snps(primer, db_snp=None, field="AF"):
    try:
        import vcf
    except:
        return None

    reader = vcf.Reader(filename=db_snp)

    left_end = int(primer.left_pos) + len(primer.left)
    right_end = int(primer.right_pos) + len(primer.right)

    contigs = list(reader.contigs.keys())
    if "chr" in contigs[0] and not "chr" in primer.chromosome:
        chrom = "chr" + primer.chromosome
    elif "chr" in primer.chromosome and not "chr" in contigs[0]:
        chrom = primer.chromosome.split("chr")[1]
    else:
        chrom = primer.chromosome

    if NEW_VCF:
        query = reader.fetch(chrom, int(primer.left_pos), left_end)
    else:
        query = reader.fetch(chrom, int(primer.left_pos) + 1, left_end)

    left_n, left_freqs = _freq_in_query(query, field)

    # can't have two pysam queries open at the same time, so have to do this this way
    if NEW_VCF:
        right_query = reader.fetch(chrom, int(primer.right_pos), right_end)
    else:
        right_query = reader.fetch(chrom, int(primer.right_pos) + 1, right_end)

    right_n, right_freqs = _freq_in_query(right_query, field)
    n = left_n + right_n
    freqs = right_freqs + left_freqs

    if n > 0:
        primer.contains_SNP = True
        primer.contains_freq_SNP = True
        if len(freqs) > 0:
            primer.snp_freq = max(freqs)
    return primer


def _sanitize_p3(handle):
    # sanitize p3 output
    # first line should always be removed:
    _ = next(handle)
    return [x for x in handle if "#" not in x.decode()]


def _sanitize_blat(blat_out):
    # sanitize blat output
    # lines should start with a number
    regex = re.compile("^\d+\t")
    return [BlatLine(*x.split("\t")) for x in blat_out if regex.match(x.decode())]


def _get_shortest(x, y, max_v):
    small = min(len(x), len(y), max_v)
    return x[:small], y[:small]


def _primers_with_blatlines(primers, blat_lines):
    for primer in primers:
        correct_lines_left = [x for x in blat_lines if x.Q_name == primer.left]
        correct_lines_right = [x for x in blat_lines if x.Q_name == primer.right]
        primer.blathits_left = [x for x in correct_lines_left]
        primer.blathits_right = [x for x in correct_lines_right]

    return primers


def _minimize_hits(blat_lines):
    groups = {}
    for line in blat_lines:
        q_name = int(line.strip().split("\t")[9])
        try:
            groups[q_name] += [line]
        except KeyError:
            groups[q_name] = [line]

    mini = min(map(len, groups.values()))
    minimized = []
    for x in groups:
        if len(groups[x]) == mini:
            minimized += groups[x]
    return minimized


def chop_region(region, size):
    """
    Chop a region in multiple regions with max size `size`
    Child regions inherit the padding of their parents
    :param region: the region to be chopped
    :param size: integer
    :return: List[Region]
    """
    if len(region) <= size:
        return [region]
    first = Region(chromosome=region.chr, start=region.start, stop=region.start+size,
                   acc_nr=region.acc_nr, padding_left=region.padding_left, padding_right=region.padding_right,
                   other=region.other_information)
    regions = [first]
    while sum([len(x) for x in regions]) < len(region):
        last = regions[-1]
        if last.stop + size >= region.stop:
            stop = region.stop
        else:
            stop = last.stop + size

        nex = Region(chromosome=region.chr, start=last.stop, stop=stop, acc_nr=region.acc_nr,
                     padding_left=region.padding_left, padding_right=region.padding_right,
                     other=region.other_information)
        regions.append(nex)
    return regions


def get_primer_from_region(region, reference, product_size, n_prims,
                           blat_exe, primer3_exe, dbsnp=None, field=None,
                           max_freq=None, strict=False, min_margin=10):

    min_length, max_length = list(map(int, product_size.split("-")))
    regions = chop_region(region, min_length)
    if any([x.size(True) > max_length for x in regions]):
        warnings.warn("Current padding results in larger regions than allowed")

    primers = []
    return_regions = []

    # TODO change such that blat is only run ONCE

    for reg in regions:
        sequence = get_sequence_fasta(reg, reference=reference)
        raw_primers = run_primer3(sequence, reg, padding=True,
                                  product_size=product_size,
                                  n_primers=n_prims, prim3_exe=primer3_exe)
        blat_out = blat_primers(raw_primers, blat_exe=blat_exe, ref=reference)

        primer = find_best(raw_primers, blat_out, dbsnp=dbsnp, field=field, max_freq=max_freq)

        # filter out primers too close to the variant
        if not (int(primer.left_pos) + int(primer.left_len) + min_margin) < int(reg.start):
            continue
        if not (int(primer.right_pos) - min_margin) > int(reg.stop):
            continue

        n_region = Region(start=int(primer.left_pos), stop=int(primer.right_pos)+len(primer.right),
                          chromosome=primer.chromosome, padding_left=0, padding_right=0,
                          acc_nr="NA", other="NA")
        primer.fragment_sequence = get_sequence_fasta(n_region, reference=reference, padding=False)
        if strict and len(primer.fragment_sequence) > max_length:
            pass
        else:
            primers.append(primer)
            return_regions.append(reg)
    if len(primers) == 0 or len(return_regions) == 0:
        raise NoPrimersException
    return return_regions, primers


def create_left_prim(primer, reference):
    primer_r_pos = primer.position + int(primer.right_pos)
    n_region = Region(start=primer_r_pos - len(primer.right), stop=primer_r_pos,
                      chromosome=primer.chromosome, padding_left=0, padding_right=0,
                      acc_nr="NA", other="NA")
    next_left_seq = get_sequence_fasta(n_region, reference=reference, padding=False)
    next_left_pos = primer_r_pos - len(primer.right)
    next_left_gc = calc_gc(next_left_seq)
    next_left = Primer()
    next_left.left = next_left_seq
    next_left.left_gc = next_left_gc
    next_left.left_pos = next_left_pos
    next_left.left_len = len(next_left_seq)
    next_left.left_name = '.'.join([primer.chromosome, str(primer.position)]) + "_left"
    next_left.chromosome = primer.chromosome
    next_left.position = primer_r_pos

    return next_left


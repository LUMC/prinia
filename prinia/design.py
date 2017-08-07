from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip)
__author__ = 'ahbbollen'

from tempfile import NamedTemporaryFile
from subprocess import check_call
import os
import warnings

from pyfaidx import Fasta
from pysam import AlignmentFile

from .models import Primer, Region
from .utils import NoPrimersException, calc_gc, NEW_VCF, generate_fastq_from_primers

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


def aln_primers(primers, bwa_exe=None, samtools_exe=None, ref=None, output_bam=None):
    """
    Align primers with BWA. 
    This only works with BWA-ALN due to the short read length. 
    :param primers: List of primers
    :param bwa_exe: Path to BWA
    :param samtools_exe: Path to samtools
    :param ref: Path to reference fasta
    :param output_bam: Path to final bam file
    :return: instance of pysam.AlignmentFile 
    """
    fq1 = NamedTemporaryFile()
    fq2 = NamedTemporaryFile()
    sai1 = NamedTemporaryFile()
    sai2 = NamedTemporaryFile()

    generate_fastq_from_primers(primers, fq1.name, fq2.name)
    aln_args1 = [bwa_exe, 'aln', ref, fq1.name]
    aln_args2 = [bwa_exe, 'aln', ref, fq2.name]

    r = check_call(aln_args1, stdout=sai1)
    if r != 0:
        raise ValueError("bwa aln crashed with error code {0}".format(r))
    r = check_call(aln_args2, stdout=sai2)
    if r != 0:
        raise ValueError("bwa aln crashed with error code {0}".format(r))

    sam = NamedTemporaryFile()
    sam_args = [bwa_exe, 'sampe', ref, sai1.name, sai2.name, fq1.name, fq2.name]
    r = check_call(sam_args, stdout=sam)
    if r != 0:
        raise ValueError("bwa sampe crashed with error code {0}".format(r))

    bam = NamedTemporaryFile()
    bam_args = [samtools_exe, 'view', '-Shb', sam.name]
    r = check_call(bam_args, stdout=bam)
    if r != 0:
        raise ValueError("samtools view crashed with error code {0}".format(r))

    final_args = [samtools_exe, 'sort', "-f", bam.name, output_bam]
    r = check_call(final_args)
    if r != 0:
        raise ValueError("samtools sort crashed with error code {0}".format(r))

    index_args = [samtools_exe, "index", output_bam]
    r = check_call(index_args)
    if r != 0:
        raise ValueError("samtools index crashed with error code {0}".format(r))

    # cleanup
    for i in [fq1, fq2, sai1, sai2, sam, bam]:
        i.close()

    return AlignmentFile(output_bam, "rb")


def _get_match_fraction(aligned_segment):
    """Get the fraction of a read matching the reference"""
    matching_bases = aligned_segment.cigartuples[0][1]
    return float(matching_bases)/aligned_segment.query_length


def _read_on_same_chrom(aligned_segment, variant=None, region=None):
    """Return true if read is on same chromosome as variant or region"""
    read_chrom = aligned_segment.reference_name.split("chr")[-1]
    if variant is not None:
        variant_chrom = variant.chromosome.split("chr")[-1]
    elif region is not None:
        variant_chrom = region.chr.split("chr")[-1]
    else:
        raise ValueError
    return read_chrom == variant_chrom


def _has_alternative_alignments(aligned_segment):
    """Return true if read has alternative alignments"""
    tags = aligned_segment.get_tags()
    return "XA" in [x[0] for x in tags]


def create_primer_from_pair(read1, read2, position=0):
    """Create primer from read pair"""
    return Primer(
        chromosome=read1.reference_name,
        position=position,
        left=read1.query_sequence,
        right=read2.query_sequence,
        left_pos=read1.reference_start,
        right_pos=read2.reference_start,
        left_len=read1.query_length,
        right_len=read2.query_length,
        left_gc=calc_gc(read1.query_sequence),
        right_gc=calc_gc(read2.query_sequence)
    )


def create_primers_bwa(bam_handle, variant=None, region=None):
    """
    Find correct pairs in a bam file containing primers for a variant or region    
    This requires to hold the bam file in memory. Use only with small bam files
    :param bam_handle: pysam.AlignmentFile instance
    :param variant: variant
    :param region: region 
    :return: generator of primers
    """
    if variant is not None and region is not None:
        raise ValueError("Either variant _or_ region must be set")
    elif variant is None and region is None:
        raise ValueError("Either variant _or_ region must be set")

    pairs = {}  # bucket to hold pairs
    for read in bam_handle:
        if read.is_read1:
            if read.query_name in pairs:
                pairs[read.query_name].update({"r1": read})
            else:
                pairs[read.query_name] = {"r1": read}
        elif read.is_read2:
            if read.query_name in pairs:
                pairs[read.query_name].update({"r2": read})
            else:
                pairs[read.query_name] = {"r2": read}

    for pair in pairs.values():
        # must match at least 90%
        # must be on same chromosome as variant
        # may not have alternative alignments
        # pair must be complete
        if "r1" not in pair or "r2" not in pair:
            continue
        if _get_match_fraction(pair['r1']) < 0.9:
            continue
        if _get_match_fraction(pair['r2']) < 0.9:
            continue

        if variant is not None:
            if not _read_on_same_chrom(pair['r1'], variant=variant):
                continue
            if not _read_on_same_chrom(pair['r2'], variant=variant):
                continue
        else:
            if not _read_on_same_chrom(pair['r1'], region=region):
                continue
            if not _read_on_same_chrom(pair['r2'], region=region):
                continue

        if _has_alternative_alignments(pair['r1']):
            continue
        if _has_alternative_alignments(pair['r2']):
            continue

        if variant is not None:
            yield create_primer_from_pair(pair['r1'], pair['r2'], variant.position_g_start)
        else:
            position = int(region.start) + int((float((int(region.stop) - int(region.start)))/2))
            yield create_primer_from_pair(pair['r1'], pair['r2'], position)


def find_best_bwa(bam, variant=None, region=None, accept_snp=False, field=None, max_freq=None, dbsnp=None):
    primers = []
    for primer in create_primers_bwa(bam, variant=variant, region=region):
        if not accept_snp and max_freq is not None:
            prim = find_snps(primer, dbsnp, field)
            if prim.snp_freq <= max_freq:
                primers.append(prim)
        else:
            primers.append(primer)
    if len(primers) == 0:
        raise NoPrimersException("No suitable primers could be detected")
    else:
        return primers


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


def _get_shortest(x, y, max_v):
    small = min(len(x), len(y), max_v)
    return x[:small], y[:small]


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
                           bwa_exe, samtools_exe, primer3_exe,
                           output_bam=None, dbsnp=None, field=None,
                           max_freq=None, strict=False, min_margin=10):

    min_length, max_length = list(map(int, product_size.split("-")))
    regions = chop_region(region, min_length)
    if any([x.size(True) > max_length for x in regions]):
        warnings.warn("Current padding results in larger regions than allowed")

    primers = []
    return_regions = []

    for reg in regions:
        sequence = get_sequence_fasta(reg, reference=reference)
        raw_primers = run_primer3(sequence, reg, padding=True,
                                  product_size=product_size,
                                  n_primers=n_prims, prim3_exe=primer3_exe)
        bam = aln_primers(raw_primers, bwa_exe=bwa_exe, samtools_exe=samtools_exe,
                          ref=reference, output_bam=output_bam)
        prims = find_best_bwa(bam, region=reg, dbsnp=dbsnp, field=field, max_freq=max_freq)
        tmp_reg = []
        tmp_prim = []

        # filter out primers too close to the variant
        for primer in prims:
            if not (int(primer.left_pos) + int(primer.left_len) + min_margin) < int(reg.start):
                continue
            if not (int(primer.right_pos) - min_margin) > int(reg.stop):
                continue

            n_region = Region(start=int(primer.left_pos), stop=int(primer.right_pos)+len(primer.right),
                              chromosome=primer.chromosome, padding_left=0, padding_right=0,
                              acc_nr="NA", other="NA")
            primer.fragment_sequence = get_sequence_fasta(n_region, reference=reference, padding=False)
            if strict and len(primer.fragment_sequence) > max_length:
                continue
            else:
                tmp_prim.append(primer)
                tmp_reg.append(reg)
        if len(tmp_prim) > 0 and len(tmp_prim) > 0:
            primers.append(tmp_prim[0])
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

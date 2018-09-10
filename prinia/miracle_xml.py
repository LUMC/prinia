"""
Parsing of Miracle XML files
"""
__author__ = 'ahbbollen'

import os
from io import StringIO
from datetime import datetime
from lxml import etree

from .models import Variant, Primer
from .utils import calc_gc, datehash

DEFAULT_XSD = os.path.join(os.path.join(os.path.dirname(__file__), "static"),
                           'variant.xsd')
REGION_XSD = os.path.join(os.path.join(os.path.dirname(__file__), "static"),
                          'region.xsd')

try:
    from itertools import zip_longest as longzip
except ImportError:
    from itertools import izip_longest as longzip


def miracle_to_primer_and_var(xml, xsd=DEFAULT_XSD):
    with open(xsd, "rb") as xsd_handle:
        schema_handle = etree.XMLSchema(etree.parse(xsd_handle))

    with open(xml, "rb") as xml_handle:
        doc = etree.parse(xml_handle)
        schema_handle.assertValid(doc)
        root = doc.getroot()

    uitslagen = root.find('UITSLAGEN')
    mid = uitslagen.find('MIRACLE_ID').text
    paneldatum = uitslagen.find('ANALYSE').find('PANELDATUM').text
    uitslags = uitslagen.find('ANALYSE').findall('UITSLAG')

    variants = [Variant() for _ in uitslags]
    primers = [Primer() for _ in uitslags]

    for item, var, prim in zip(uitslags, variants, primers):
        chr = item.find('VARIANT').find('GEN').find('CHROMOSOOM').text
        gen_code = item.find('VARIANT').find('GEN').find('CODE').text
        gen_name = item.find('VARIANT').find('GEN').find('NAAM').text

        build = item.find('VARIANT').find('BUILD').text
        genome_ref = item.find('VARIANT').find('GENOOM_REF_SEQ').text
        alt_genome = item.find('VARIANT').find('WIJZIGING_GENOOM').text
        tr_ref = item.find('VARIANT').find('TR_REF_SEQ').text
        alt_gen = item.find('VARIANT').find('WIJZIGING_GEN').text
        alt_rna = item.find('VARIANT').find('WIJZIGING_RNA').text
        alt_prot = item.find('VARIANT').find('WIJZIGING_EIWIT').text

        prims = item.find('VARIANT').find("PRIMERS")

        primer_f_code = prims.find("PRIMER_F").find("CODE").text
        primer_f_seq = prims.find("PRIMER_F").find("SEQUENTIE").text

        primer_r_code = prims.find("PRIMER_R").find("CODE").text
        primer_r_seq = prims.find("PRIMER_R").find("SEQUENTIE").text

        var.miracle_id = mid
        var.datum = paneldatum
        var.refseq_build = build
        var.chromosome = chr
        var.genomic_id_ncbi = genome_ref
        var.variant_on_genome = alt_genome

        # is present codes are those used in MAGPIE
        if item.find('PATERNAAL').text == 'N':
            var.is_present_father = '1'
        elif item.find('PATERNAAL').text == 'J':
            var.is_present_father = '6'
        else:
            var.is_present_father = 'NA'

        if item.find('MATERNAAL').text == 'N':
            var.is_present_mother = '1'
        elif item.find('MATERNAAL').text == 'J':
            var.is_present_mother = '6'
        else:
            var.is_present_mother = '1'

        var.gene_id = gen_code
        var.gene_name = gen_name
        var.transcript_id_ncbi = tr_ref
        var.variant_on_transcript_dna = alt_gen
        var.variant_on_transcript_rna = alt_rna
        var.variant_on_trascript_protein = alt_prot

        if item.find("UITSLAG_CODE").text == 'HET':
            var.allele = 'Heterozygous'
        elif item.find('UITSLAG_CODE').text == 'HOM':
            var.allele = 'Homozygous'
        # FIXME: should include cases of maternally hemizygous and
        # paternally hemizygous
        else:
            var.allele = 'NA'

        if item.find("DE_NOVO").text == 'J':
            var.variant_on_genome_origin = 'De novo'
        else:
            var.variant_on_genome_origin = 'NA'

        prim.chromosome = chr
        prim.left = primer_f_seq
        prim.right = primer_r_seq
        prim.left_gc = calc_gc(primer_f_seq)
        prim.right_gc = calc_gc(primer_r_seq)
        prim.left_name = primer_f_code
        prim.right_name = primer_r_code

    return variants, primers


def common_xml(miracle_id, paneldatum, panel=True, old=True):
    """
    Creates common XML structure
    :param miracle_id: miracle id
    :param paneldatum: panel date
    :param panel: boolean whether PANEL or EXOMEWIDE
    :return: results object in XML tree
    """
    document = etree.Element("DOCUMENT")

    type = etree.SubElement(document, "TYPE")
    type.text = "PRIM_MIR"

    version = etree.SubElement(document, "VERSIE")

    if old:
        version.text = "2.0"
        analysis = etree.SubElement(document, "ANALYSE")
    else:
        version.text = "2.1"
        analysis = etree.SubElement(document, "ANALYSEGROEP")

    analysis_id = etree.SubElement(analysis, "ID")
    analysis_id.text = "1234"  # TODO: this is currently unknown

    date = etree.SubElement(analysis, "DATUM")
    now = datetime.utcnow()
    date.text = now.strftime("%Y-%m-%d %H:%M:%S")

    method = etree.SubElement(analysis, "METHODE")
    if panel:
        method.text = "PANEL"
    else:
        method.text = "EXOMEWIDE"

    paneldate = etree.SubElement(analysis, "PANELDATUM")
    paneldate.text = paneldatum

    results = etree.SubElement(analysis, "RESULTS")

    m_id = etree.SubElement(results, "MIRACLE_ID")
    m_id.text = miracle_id

    return document, results


def primer_to_xml(root, prim, var):

    primer_f = etree.SubElement(root, "PRIMER_F")
    f_code = etree.SubElement(primer_f, "CODE")
    f_code.text = var.gene_id + "_F." + datehash()
    f_seq = etree.SubElement(primer_f, "SEQUENTIE")
    f_seq.text = prim.left
    f_loc = etree.SubElement(primer_f, "COORDINATE")
    f_loc.text = str(prim.left_pos)

    primer_r = etree.SubElement(root, "PRIMER_R")
    r_code = etree.SubElement(primer_r, "CODE")
    r_code.text = var.gene_id + "_R." + datehash()
    r_seq = etree.SubElement(primer_r, "SEQUENTIE")
    r_seq.text = prim.right
    r_loc = etree.SubElement(primer_r, "COORDINATE")
    r_loc.text = str(prim.right_pos)

    return root


def vars_and_primers_to_xml(variants, primers, xml_path=None, xsd=DEFAULT_XSD,
                            old=True):
    panel = int(variants[0].in_gene_panel) == 1

    if old:
        xsd = os.path.join(os.path.join(os.path.dirname(__file__), 'static'),
                           'variant.xsd')
    else:
        xsd = os.path.join(os.path.join(os.path.dirname(__file__), 'static'),
                           'variant_new.xsd')

    document, results = common_xml(variants[0].miracle_id, variants[0].datum,
                                   panel, old=old)
    i = 1

    for var, prim in longzip(variants, primers, fillvalue=None):
        uitslag = etree.SubElement(results, "UITSLAG")
        variant = etree.SubElement(uitslag, "VARIANT")

        chrom = etree.SubElement(variant, "CHROMOSOOM")
        chrom.text = var.chromosome
        coord_from = etree.SubElement(variant, "COORDINATE_FROM")
        coord_from.text = var.position_g_start
        coord_to = etree.SubElement(variant, "COORDINATE_TO")
        coord_to.text = var.position_g_end

        gene = etree.SubElement(variant, "GEN")

        chr = etree.SubElement(gene, "CHROMOSOOM")
        if var.chromosome.startswith("chr"):
            chr.text = var.chromosome.split("chr")[-1]
        else:
            chr.text = var.chromosome

        code = etree.SubElement(gene, "CODE")
        code.text = var.gene_id

        name = etree.SubElement(gene, "NAAM")
        name.text = var.gene_name

        loc = etree.SubElement(gene, "LOCATIE")  # noqa
        build = etree.SubElement(variant, "BUILD")
        build.text = var.refseq_build

        genome_ref = etree.SubElement(variant, "GENOOM_REF_SEQ")
        genome_ref.text = var.genomic_id_ncbi

        variation = etree.SubElement(variant, "WIJZIGING_GENOOM")
        variation.text = var.variant_on_genome

        transcript = etree.SubElement(variant, "TR_REF_SEQ")
        transcript.text = var.transcript_id_ncbi

        variation_gene = etree.SubElement(variant, "WIJZIGING_GEN")
        variation_gene.text = var.variant_on_transcript_dna

        variation_rna = etree.SubElement(variant, "WIJZIGING_RNA")
        variation_rna.text = var.variant_on_transcript_rna

        variation_prot = etree.SubElement(variant, "WIJZIGING_EIWIT")
        variation_prot.text = var.variant_on_transcript_protein

        prims = etree.SubElement(variant, "PRIMERS")

        if prim is not None:
            frag_len = etree.SubElement(prims, "FRAGMENT_LENGTH")
            # fragment is the region between the end of the forward primer,
            #  and the start of the reverse primer

            frag_len.text = str(int(prim.right_pos) - (int(prim.left_pos) +
                                                       len(prim.left)))

            gc_perc = etree.SubElement(prims, "GC_PERC")

            try:
                gc_perc.text = str(int(calc_gc(prim.fragment_sequence)))
            except ValueError:
                gc_perc.text = '0'
            _ = primer_to_xml(prims, prim, var)  # noqa
        else:
            comment = etree.SubElement(uitslag, "OPMERKING")
            comment.text = "NO PRIMERS FOUND"

        denovo = etree.SubElement(uitslag, "DE_NOVO")
        if var.variant_on_genome_origin == "De novo":
            denovo.text = 'J'
        else:
            denovo.text = 'N'

        father = etree.SubElement(uitslag, "PATERNAAL")
        if var.is_present_father == '6' or var.is_present_father == '5':
            father.text = "J"
        else:
            father.text = "N"

        mother = etree.SubElement(uitslag, "MATERNAAL")
        if var.is_present_mother == '6' or var.is_present_mother == '5':
            mother.text = "J"
        else:
            mother.text = "N"

        allele = etree.SubElement(uitslag, "UITSLAG_CODE")
        if var.allele.upper() == 'HETEROZYGOUS':
            allele.text = 'HET'
        elif var.allele.upper() == 'HOMOZYGOUS':
            allele.text = 'HOM'
        elif var.allele.upper() == 'HEMIZYGOUS':
            allele.text = 'HEM'
        elif var.allele.upper() == "MATERNAL":
            allele.text = "HET"
        elif var.allele.upper() == "PATERNAL":
            allele.text = "HET"
        else:
            raise ValueError("Unknown zygosity {0}".format(var.allele.upper()))

        confirmation = etree.SubElement(uitslag, "BEVESTIGEN")
        print(var.confirm_in_lab)
        if var.confirm_in_lab:
            confirmation.text = "J"
        else:
            confirmation.text = "N"

        i += 1

    xml = etree.tostring(document, pretty_print=True, xml_declaration=True)
    xml_plain = etree.tostring(document)
    with open(xsd, "rb") as xsd_handle:
        xml_schema = etree.XMLSchema(etree.parse(xsd_handle))
        xml_schema.assertValid(etree.parse(StringIO(xml_plain.decode())))
    return xml


def regions_and_primers_to_xml(variants, primers, xml_path=None,
                               xsd=REGION_XSD, sample='', build="GRCh37",
                               old=True):
    if old:
        xsd = os.path.join(os.path.join(os.path.dirname(__file__), 'static'),
                           'region.xsd')
    else:
        xsd = os.path.join(os.path.join(os.path.dirname(__file__), 'static'),
                           'region_new.xsd')
    document, results = common_xml(sample, '', old=old)
    i = 1

    for region, prim in zip(variants, primers):
        gap = etree.SubElement(results, "GAP")

        g_build = etree.SubElement(gap, "BUILD")
        g_build.text = build

        chrom = etree.SubElement(gap, "CHROMOSOOM")

        if region.chr.startswith("chr"):
            chrom.text = region.chr.split("chr")[-1]
        else:
            chrom.text = region.chr

        gene = etree.SubElement(gap, "GEN")
        code = etree.SubElement(gene, "CODE")
        if region.other_information != "NA":
            try:
                code.text = region.other_information.split(
                    ","
                )[0].split("|")[0]
            except IndexError:
                code.text = ""
        else:
            code.text = ""
        name = etree.SubElement(gene, "NAAM")
        name.text = code.text

        coord_start = etree.SubElement(gap, "COORDINATE_FROM")
        coord_start.text = str(region.start)
        coord_end = etree.SubElement(gap, "COORDINATE_TO")
        coord_end.text = str(region.stop)

        primers = etree.SubElement(gap, "PRIMERS")
        frag_length = etree.SubElement(primers, "FRAGMENT_LENGTH")
        frag_length.text = str(len(prim.fragment_sequence))
        gc_perc = etree.SubElement(primers, "GC_PERC")
        gc_perc.text = str(calc_gc(prim.fragment_sequence))

        primer_f = etree.SubElement(primers, "PRIMER_F")
        f_code = etree.SubElement(primer_f, "CODE")
        if "NA" in region.other_information:
            f_code.text = prim.left_name + "." + datehash()
        else:
            f_code.text = (region.other_information.split("|")[0] + "_F." +
                           datehash())
        f_seq = etree.SubElement(primer_f, "SEQUENTIE")
        f_seq.text = prim.left
        f_coord = etree.SubElement(primer_f, "COORDINATE")
        f_coord.text = str(prim.left_pos)

        primer_r = etree.SubElement(primers, "PRIMER_R")
        r_code = etree.SubElement(primer_r, "CODE")
        if "NA" in region.other_information:
            r_code.text = prim.right_name + "." + datehash()
        else:
            r_code.text = (region.other_information.split("|")[0] + "_R." +
                           datehash())
        r_seq = etree.SubElement(primer_r, "SEQUENTIE")
        r_seq.text = prim.right
        r_coord = etree.SubElement(primer_r, "COORDINATE")
        r_coord.text = str(prim.right_pos)

        i += 1

    xml = etree.tostring(document, pretty_print=True, xml_declaration=True)
    xml_plain = etree.tostring(document)
    with open(xsd, "rb") as xsd_handle:
        xml_schema = etree.XMLSchema(etree.parse(xsd_handle))
        xml_schema.assertValid(etree.parse(StringIO(xml_plain.decode())))
    return xml

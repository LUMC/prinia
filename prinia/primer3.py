"""
prinia.primer3
~~~~~~~~~~~~~~

:copyright: (c) 2017-2018 Sander Bollen
:copyright: (c) 2017-2018 Leiden University Medical Center
:license: MIT
"""
from pathlib import Path
from tempfile import NamedTemporaryFile
from subprocess import check_call
from typing import Optional
import os
import json

import re

from jsonschema import validate

from .models import Primer


SETTINGS_SCHEMA = (Path(__file__).parent / Path("static") /
                   Path("primer3_settings_schema.json"))


def parse_settings(settings_file: Optional[Path] = None) -> dict:
    """
    Parse settings file to settings dictionary.

    Parameters not specified in settings_file will be taken from default values
    in schema. If settings_file is None, all values will be the default values
    in schema
    :param settings_file: Optional path to settings file
    :return: dict with settings
    :raises: ValidationError if settings_file does not conform to schema
    """
    with SETTINGS_SCHEMA.open("r") as schema_handle:
        schema_dict = json.load(schema_handle)

    defaults = {k: v['default'] for k, v in schema_dict['properties'].items()}

    if settings_file is None:
        return defaults

    with settings_file.open("r") as settings_handle:
        settings_dict = json.load(settings_handle)

    validate(settings_dict, schema_dict)  # raises ValidationError if failure

    default_keys = defaults.keys() - settings_dict.keys()
    default_dict = {k: defaults[k] for k in default_keys}
    generated_dict = {**settings_dict, **default_dict}

    # if range is set, but optimum size is not, calculate opt size
    if ('primer_product_opt_size' in default_dict and
            "primer_product_size_range" in settings_dict):
        product_range = settings_dict['primer_product_size_range']
        min_size, max_size = product_range.split("-")
        opt_size = int((int(min_size) + (int(max_size) - int(min_size)) // 2))
        generated_dict['primer_product_opt_size'] = opt_size

    return generated_dict


class Primer3(object):

    def __init__(self, primer3_exe, template, target, excluded_region,
                 opt_prim_length=25, opt_gc_perc=50, min_melting_t=58,
                 max_melting_t=62, min_product_size=200,
                 max_product_size=600):
        self.primer3_exe = primer3_exe
        self.template = template
        self.target = target
        self.excluded_region = excluded_region
        self.opt_prim_length = opt_prim_length
        self.opt_gc_perc = opt_gc_perc
        self.min_melting_t = min_melting_t
        self.max_melting_t = max_melting_t
        self.min_product_size = min_product_size
        self.max_product_size = max_product_size

    @property
    def range(self):
        return "{0}-{1}".format(self.min_product_size, self.max_product_size)

    @property
    def opt_size(self):
        return int(
            self.min_product_size +
            (float((self.max_product_size - self.min_product_size)) / 2)
        )

    def create_config(self, handle):
        """Create config for primer3."""

        cfg_str = "SEQUENCE_ID=example\n" \
                  "SEQUENCE_TEMPLATE={seq}\n" \
                  "SEQUENCE_TARGET={tar}\n" \
                  "SEQUENCE_EXCLUDED_REGION={exc}\n" \
                  "PRIMER_TASK=pick_detection_primers\n" \
                  "PRIMER_PICK_LEFT_PRIMER=1\n" \
                  "PRIMER_PICK_INTERNAL_OLIGO=0\n" \
                  "PRIMER_PICK_RIGHT_PRIMER=1\n" \
                  "PRIMER_MIN_GC=20.0\n" \
                  "PRIMER_INTERNAL_MIN_GC=20.0\n" \
                  "PRIMER_OPT_GC_PERCENT={gc}\n" \
                  "PRIMER_MAX_GC=80.0\n" \
                  "PRIMER_INTERNAL_MAX_GC=80.0\n" \
                  "PRIMER_WT_GC_PERCENT_LT=0.0\n" \
                  "PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0\n" \
                  "PRIMER_GC_CLAMP=0\n" \
                  "PRIMER_MAX_END_GC=5\n" \
                  "PRIMER_OPT_SIZE={size}\n" \
                  "PRIMER_MIN_SIZE={isize}\n" \
                  "PRIMER_MAX_SIZE={asize}\n" \
                  "PRIMER_MAX_NS_ACCEPTED=0\n" \
                  "PRIMER_PRODUCT_SIZE_RANGE={range}\n" \
                  "PRIMER_PRODUCT_OPT_SIZE={osize}\n" \
                  "PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.1\n" \
                  "PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.1\n" \
                  "P3_FILE_FLAG=1\n" \
                  "SEQUENCE_INTERNAL_EXCLUDED_REGION=37,21\n" \
                  "PRIMER_EXPLAIN_FLAG=1\n" \
                  "PRIMER_MIN_TM={it}\n" \
                  "PRIMER_MAX_TM={at}\n" \
                  "PRIMER_NUM_RETURN=200\n" \
                  "=".format(seq=self.template, tar=self.target,
                             exc=self.excluded_region,
                             gc=self.opt_gc_perc,
                             size=self.opt_prim_length,
                             isize=self.opt_prim_length-5,
                             asize=self.opt_prim_length+5,
                             range=self.range,
                             osize=self.opt_size,
                             it=self.min_melting_t,
                             at=self.max_melting_t
                             )

        handle.write(cfg_str)

    def run(self):
        cfg = NamedTemporaryFile(delete=False, mode="w")
        out = NamedTemporaryFile()

        self.create_config(cfg)
        cfg.close()
        args = [self.primer3_exe, "-output", out.name, cfg.name]

        _ = check_call(args)  # noqa

        retval = []
        with open(out.name) as handle:
            for l in handle:
                retval.append(l.strip())

        out.close()
        os.remove(cfg.name)

        return retval


def _parse_single_pair_id(id, lines):
    """Parse single pair id"""
    left_seq_key = "PRIMER_LEFT_{0}_SEQUENCE=".format(id)
    right_seq_key = "PRIMER_RIGHT_{0}_SEQUENCE=".format(id)
    left_pos_key = "PRIMER_LEFT_{0}=".format(id)
    right_pos_key = "PRIMER_RIGHT_{0}=".format(id)
    left_gc_key = "PRIMER_LEFT_{0}_GC_PERCENT".format(id)
    right_gc_key = "PRIMER_RIGHT_{0}_GC_PERCENT".format(id)

    left_seq_line = [x for x in lines if x.startswith(left_seq_key)]
    right_seq_line = [x for x in lines if x.startswith(right_seq_key)]
    left_pos_line = [x for x in lines if x.startswith(left_pos_key)]
    right_pos_line = [x for x in lines if x.startswith(right_pos_key)]
    left_gc_line = [x for x in lines if x.startswith(left_gc_key)]
    right_gc_line = [x for x in lines if x.startswith(right_gc_key)]

    d = {}
    for name, val in zip(
            ["left", "right", "left_gc", "right_gc"],
            [left_seq_line, right_seq_line, left_gc_line, right_gc_line]
    ):
        if len(val) != 1:
            raise ValueError("Value for {n} "
                             "must have exactly one match".format(n=name))
        d[name] = val[0].split("=")[-1]

    for name, val in zip(["left_pos", "right_pos"],
                         [left_pos_line, right_pos_line]):
        if len(val) != 1:
            raise ValueError("Vale for {n} "
                             "must have exactly one match".format(n=name))
        d[name] = val[0].split("=")[-1].split(",")[0]

    return Primer(**d)


def parse_primer3_output(lines):
    """Parse primer3 output to Primers"""
    seq_re = re.compile('^PRIMER_LEFT_(\d+)_SEQUENCE=.+$')

    matches = [seq_re.match(x) for x in lines]
    pair_ids = [int(x.group(1)) for x in matches if x is not None]

    for id in pair_ids:
        yield _parse_single_pair_id(id, lines)

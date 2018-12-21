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

from jinja2 import Template
from jsonschema import validate

from .models import Primer


SETTINGS_SCHEMA = (Path(__file__).parent / Path("static") /
                   Path("primer3_settings_schema.json"))
TEMPLATE = Path(__file__).parent / Path("static") / Path("primer3_conf.j2")


def parse_settings(settings_file: Optional[Path] = None,
                   settings_schema: Path = SETTINGS_SCHEMA) -> dict:
    """
    Parse settings file to settings dictionary.

    Parameters not specified in settings_file will be taken from default values
    in schema. If settings_file is None, all values will be the default values
    in schema
    :param settings_file: Optional path to settings file
    :param settings_schema: Path to schema json for settings. Default is
    supplied schema in static/primer3_settings_schema.json
    :return: dict with settings
    :raises: ValidationError if settings_file does not conform to schema
    """
    with settings_schema.open("r") as schema_handle:
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


def create_primer3_config(settings_dict: dict,
                          sequence_template: str,
                          sequence_target: str,
                          excluded_region: str,
                          template_path: Path = TEMPLATE,
                          thermodynamic_params: Optional[Path] = None)-> str:
    """Create primer3 configuration from jinja2 template and settings_dict"""
    with template_path.open("r") as thandle:
        template = Template(thandle.read(), trim_blocks=True)

    # primer3 wants paths to end in slashes
    if thermodynamic_params is not None:
        therm_params = str(thermodynamic_params)
        if not therm_params.endswith("/"):
            therm_params += "/"
    else:
        therm_params = None

    return template.render(
        {
            "seq": sequence_template,
            "target": sequence_target,
            "excluded_region": excluded_region,
            "settings": settings_dict,
            "thermodynamic_params": therm_params
        }
    )


class Primer3(object):

    def __init__(self, primer3_exe: str, template: str,
                 target: str, excluded_region: str,
                 settings_dict: dict,
                 thermodynamic_params: Optional[Path] = None):
        self.primer3_exe = primer3_exe
        self.template = template
        self.target = target
        self.excluded_region = excluded_region
        self.settings_dict = settings_dict
        self.thermodynamic_params = thermodynamic_params

    def run(self):
        cfg = NamedTemporaryFile(delete=False, mode="w")
        out = NamedTemporaryFile()

        settings = create_primer3_config(
            self.settings_dict, self.template,
            self.target, self.excluded_region,
            thermodynamic_params=self.thermodynamic_params
        )

        cfg.write(settings)
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
    seq_re = re.compile(r'^PRIMER_LEFT_(\d+)_SEQUENCE=.+$')

    matches = [seq_re.match(x) for x in lines]
    pair_ids = [int(x.group(1)) for x in matches if x is not None]

    for id in pair_ids:
        yield _parse_single_pair_id(id, lines)

from builtins import (map, open, str, zip)

__author__ = 'ahbbollen'

import re
from collections import OrderedDict

from .models import Variant


def var_from_lovd(path):
    d = OrderedDict()
    comments = {}

    seen_header = False
    data_lines = 0
    regex = re.compile('^(("\{\{[\w\/]+\}\}"|\{\{[\w\/]+\}\})\t)+')

    with open(path, "rb") as handle:
        for line in handle:
            line = line.strip().decode()

            if line.startswith("#") and "=" in line:
                # we have a header/comment line
                contents = list(map(str.strip, line.split("#")[1].split("=")))
                comments[contents[0]] = contents[1]

            elif not seen_header and regex.match(line):
                seen_header = True
                contents = line.split("\t")
                for header in contents:
                    d[header.strip('""').strip("{{").strip("}}")] = []

            elif not seen_header and not regex.match(line):
                raise ValueError("Malformed file")

            elif seen_header:
                data_lines += 1
                contents = line.split("\t")
                if len(contents) != len(d.keys()):
                    raise ValueError("Number of data fields does not "
                                     "match number of headers.")
                for header, field in zip(d.keys(), contents):
                    if field.strip('""'):
                        d[header].append(field.strip('""'))
                    else:
                        d[header].append(u'NA')

            else:
                raise ValueError("Something odd happened")

    comments["n_line"] = data_lines

    return Variant.from_dict(d, comments)

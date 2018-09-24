"""
test_primer3.py
~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from pathlib import Path
from pytest import fixture

from prinia.primer3 import parse_settings


data_dir = Path(__file__).parent / Path("data")

@fixture
def valid_settings():
    return data_dir / Path("valid_settings.json")


@fixture
def valid_partial_settings():
    return data_dir / Path("valid_partial_settings.json")


def test_none_settings():
    parsed = parse_settings()
    assert parsed == {
        "primer_min_gc": 20,
        "primer_internal_min_gc": 20,
        "primer_opt_gc_percent": 50,
        "primer_max_gc": 80,
        "primer_internal_max_gc": 80,
        "primer_wt_gc_percent_lt": 0,
        "primer_internal_wt_gc_percent_lt": 0,
        "primer_wt_gc_percent_gt": 0,
        "primer_internal_wt_gc_percent_gt": 0,
        "primer_gc_clamp": 0,
        "primer_max_end_gc": 5,
        "primer_opt_size": 25,
        "primer_min_size": 20,
        "primer_max_size": 30,
        "primer_max_ns_accepted": 0,
        "primer_product_size_range": "200-450",
        "primer_product_opt_size": 325,
        "primer_pair_wt_product_size_gt": 0.1,
        "primer_pair_wt_product_size_lt": 0.1,
        "primer_min_tm": 58,
        "primer_max_tm": 62,
        "primer_num_return": 200
    }


def test_complete_settings(valid_settings):
    parsed = parse_settings(valid_settings)

    # TODO: work out


def test_partial_settings(valid_partial_settings):
    parsed = parse_settings(valid_partial_settings)

    # TODO: work out

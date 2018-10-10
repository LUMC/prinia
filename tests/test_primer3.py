"""
test_primer3.py
~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from pathlib import Path
from shutil import which
import pytest
from prinia.primer3 import parse_settings, create_primer3_config, Primer3


data_dir = Path(__file__).parent / Path("data")


@pytest.fixture
def valid_settings():
    return data_dir / Path("valid_settings.json")


@pytest.fixture
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
        "primer_num_return": 200,
        "primer_opt_tm": 60,
        "primer_pair_max_diff_tm": 100,
        "primer_max_hairpin_th": 47
    }


def test_complete_settings(valid_settings):
    parsed = parse_settings(valid_settings)
    assert parsed == {
        "primer_min_gc": 30,
        "primer_internal_min_gc": 30,
        "primer_opt_gc_percent": 51,
        "primer_max_gc": 81,
        "primer_internal_max_gc": 81,
        "primer_wt_gc_percent_lt": 0.12,
        "primer_internal_wt_gc_percent_lt": 0.12,
        "primer_wt_gc_percent_gt": 0.12,
        "primer_internal_wt_gc_percent_gt": 0.12,
        "primer_gc_clamp": 1,
        "primer_max_end_gc": 0,
        "primer_opt_size": 21,
        "primer_min_size": 19,
        "primer_max_size": 29,
        "primer_max_ns_accepted": 1,
        "primer_product_size_range": "300-500",
        "primer_product_opt_size": 300,
        "primer_pair_wt_product_size_gt": 0.12,
        "primer_pair_wt_product_size_lt": 0.12,
        "primer_min_tm": 55,
        "primer_max_tm": 75,
        "primer_num_return": 100,
        "primer_opt_tm": 60,
        "primer_pair_max_diff_tm": 100,
        "primer_max_hairpin_th": 47
    }


def test_partial_settings(valid_partial_settings):
    parsed = parse_settings(valid_partial_settings)

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
        "primer_gc_clamp": 1,
        "primer_max_end_gc": 0,
        "primer_opt_size": 21,
        "primer_min_size": 19,
        "primer_max_size": 29,
        "primer_max_ns_accepted": 1,
        "primer_product_size_range": "300-500",
        "primer_product_opt_size": 400,
        "primer_pair_wt_product_size_gt": 0.12,
        "primer_pair_wt_product_size_lt": 0.12,
        "primer_min_tm": 55,
        "primer_max_tm": 75,
        "primer_num_return": 100,
        "primer_opt_tm": 60,
        "primer_pair_max_diff_tm": 100,
        "primer_max_hairpin_th": 47
    }


configuration_data = [
    (None, data_dir / Path("primer3_conf_default.txt")),
    (
        data_dir / Path("valid_settings.json"),
        data_dir / Path("primer3_conf_full_settings.txt")
    ),
    (
        data_dir / Path("valid_partial_settings.json"),
        data_dir / Path("primer3_conf_partial_settings.txt")
    )
]


@pytest.mark.parametrize("settings_json, conf", configuration_data)
def test_primer3_configuration(settings_json, conf):
    parsed_settings = parse_settings(settings_json)
    generated_conf = create_primer3_config(parsed_settings, "ACTG",
                                           "50-60", "50-60")
    with conf.open("r") as chandle:
        expected_conf = chandle.read().strip()

    assert generated_conf == expected_conf


@pytest.fixture
def rCRS():
    rCRS_path = data_dir / Path("rCRS.fa")
    with rCRS_path.open("r") as handle:
        lines = handle.readlines()
    return "".join(list(map(str.strip, lines[1:])))


@pytest.mark.parametrize("settings_json, ignored", configuration_data)
def test_primer3_run_non_failing(settings_json, ignored, rCRS):
    settings_dict = parse_settings(settings_json)
    # assumes primer3_core is on the PATH
    primer3_exe = which("primer3_core")
    p = Primer3(primer3_exe, rCRS, "3200,300", "3300,50", settings_dict)
    p.run()

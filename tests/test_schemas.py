"""
test_schemas.py
~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from pathlib import Path
import json

from jsonschema.validators import Draft4Validator

from pytest import fixture
import pytest


@fixture
def settings_schema():
    return (Path(__file__).parent.parent / Path("prinia") / Path("static")
            / Path("primer3_settings_schema.json"))


data_dir = Path(__file__).parent / Path("data")

setting_jsons_data = [
    (data_dir / Path("valid_settings.json"), True),
    (data_dir / Path("valid_partial_settings.json"), True),
    (data_dir / Path("valid_partial_settings2.json"), True),
    (data_dir / Path("settings_wrong_types.json"), False),
    (data_dir / Path("settings_invalid_values.json"), False)
]


def test_settings_schema_is_loadable_json(settings_schema):
    with settings_schema.open("r") as handle:
        json.load(handle)


def test_settings_schema_is_valid_schema(settings_schema):
    with settings_schema.open("r") as handle:
        s = json.load(handle)
    Draft4Validator.check_schema(s)


@pytest.mark.parametrize("settings_json, should_be_valid", setting_jsons_data)
def test_settings_json_valid(settings_schema, settings_json, should_be_valid):
    with settings_schema.open("r") as shandle:
        schema = json.load(shandle)

    schema_instance = Draft4Validator(schema)

    with settings_json.open("r") as jhandle:
        j = json.load(jhandle)
        assert schema_instance.is_valid(j) == should_be_valid

"""
test_utils.py
~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
import pytest

from prinia.utils import is_at_least_version_samtools

samtools_version_data = [
    ("1.9-33-g2d34e15", (1, 3), True),
    ("1.2-33-g2d34e15", (1, 3), False),
    ("1.9", (1, 3), True),
    ("1.2", (1, 3), False),
    ("samtools 1.9-33-g2d34e15", (1, 3), True),
    ("samtools 1.9", (1, 3), True),
    ("samtools 1.2", (1, 3), False)
]


@pytest.mark.parametrize("version_str, version_tupl, expected",
                         samtools_version_data)
def test_is_at_least_version_samtools(version_str, version_tupl, expected):
    val = is_at_least_version_samtools(version_str, version_tupl)
    assert val == expected

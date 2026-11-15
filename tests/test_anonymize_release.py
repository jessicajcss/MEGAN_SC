"""tools/anonymize_release.py: only the value of "ee_project" changes in a provenance file."""
import os, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tools'))
from anonymize_release import anonymize_bytes


def test_only_the_project_value_changes():
    raw = b'{\n "year": 2019,\n "engine": "/content/drive/MyDrive/megan-sc",\n "ee_project": "some-project-01",\n "started": "x"\n}'
    out = anonymize_bytes(raw)
    assert out == raw.replace(b'"some-project-01"', b'"withheld"')


def test_one_field_is_required():
    with pytest.raises(ValueError):
        anonymize_bytes(b'{"year": 2019}')

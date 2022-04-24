import io
from unittest import TestCase

import pytest
import requests
from Bio.PDB import MMCIFIO, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from mmcifix import dict_to_file, fix_dict
from mmcifix.fixers import FixAuthSeqId, FixLabelSeqId, FixDatabaseId
from mmcifix.util import find_changes


def is_biopython_parseable(d: dict) -> bool:
    struct_id = d["_entry.id"][0]
    output = io.StringIO()
    dict_to_file(d, output)
    output.seek(0)
    try:
        MMCIFParser().get_structure(struct_id, output)
        return True
    except ValueError as e:
        raise e
    return False


@pytest.fixture(scope="module")
def alphafill_P04406() -> dict:
    res = requests.get("https://alphafill.eu/v1/aff/P04406", allow_redirects=True)
    return MMCIF2Dict(io.StringIO(res.text))


def test_fix_label_seq_id(alphafill_P04406: dict):
    fixed = FixLabelSeqId().run(alphafill_P04406)

    # Check that we have fixed the unknown ids
    assert all([it not in [".", "?"] for it in fixed['_atom_site.label_seq_id']])


def test_fix_auth_seq_id(alphafill_P04406):
    fixed = FixAuthSeqId().run(alphafill_P04406)

    # Check that we have fixed the unknown ids
    assert all([it not in [".", "?"] for it in fixed['_atom_site.auth_seq_id']])

def test_fix_database_id(alphafill_P04406):
    assert alphafill_P04406["_database_2.database_id"] == ["AF"]

    fixed = FixDatabaseId().run(alphafill_P04406)

    # Check that we have fixed the unknown ids
    assert fixed["_database_2.database_id"] != ["AF"]
    assert isinstance(fixed["_database_2.database_id"], list)

def test_fix_alphafill(alphafill_P04406):
    fixed = fix_dict(alphafill_P04406, fixers=[
        "label_seq_id",
        "auth_seq_id",
    ])
    # Check that we have fixed the unknown ids
    assert all([it not in [".", "?"] for it in fixed['_atom_site.label_seq_id']])
    assert all([it not in [".", "?"] for it in fixed['_atom_site.auth_seq_id']])

    # Check that biopython can now parse it
    assert is_biopython_parseable(fixed)


def test_find_changes():
    assert list(find_changes([1, 2, 3])) == [1, 2]
    assert list(find_changes([1, 1, 2, 2, 3, 3])) == [2, 4]

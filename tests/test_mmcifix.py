import io
from unittest import TestCase

import pytest
import requests
from Bio.PDB import MMCIFIO, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path

from mmcifix import dict_to_file, fix_dict
from mmcifix.fixers import FixAuthSeqId, FixLabelSeqId, FixDatabaseId, FixAsymIdForPdb
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
    with (Path(__file__).parent / "AF-P04406-F1-model_v1.cif").open() as fp:
        return MMCIF2Dict(fp)


@pytest.fixture(scope="module")
def alphafill_P27037() -> dict:
    with (Path(__file__).parent / "AF-P27037-F1-model_v1.cif").open() as fp:
        return MMCIF2Dict(fp)


@pytest.fixture(scope="module")
def alphafill_P27037_one_ligand() -> dict:
    with (Path(__file__).parent / "P27037_one_ligand.cif").open() as fp:
        return MMCIF2Dict(fp)


def test_fix_label_seq_id(alphafill_P04406: dict):
    fixed = FixLabelSeqId().run(alphafill_P04406)

    # Check that we have fixed the unknown ids
    assert all([it not in [".", "?"] for it in fixed["_atom_site.label_seq_id"]])


def test_fix_auth_seq_id(alphafill_P04406):
    fixed = FixAuthSeqId().run(alphafill_P04406)

    # Check that we have fixed the unknown ids
    assert all([it not in [".", "?"] for it in fixed["_atom_site.auth_seq_id"]])


def test_fix_database_id(alphafill_P04406):
    assert alphafill_P04406["_database_2.database_id"] == ["AF"]

    fixed = FixDatabaseId().run(alphafill_P04406)

    # Check that we have fixed the unknown ids
    assert fixed["_database_2.database_id"] != ["AF"]
    assert isinstance(fixed["_database_2.database_id"], list)


def test_fix_asym_id(alphafill_P27037_one_ligand: dict):
    cif = alphafill_P27037_one_ligand
    assert "AD" in cif["_atom_site.auth_asym_id"]
    assert len(set(cif["_atom_site.auth_asym_id"])) == 2

    fixed = FixAsymIdForPdb().run(cif)

    # AD should no longer be present
    assert "AD" not in fixed["_atom_site.auth_asym_id"]
    assert "AD" not in fixed["_atom_site.label_asym_id"]

    assert len(set(fixed["_atom_site.auth_asym_id"])) == 2
    new_id = next(iter(set(fixed["_atom_site.auth_asym_id"]) - {"A"}))
    # The replacement ID should only be 1 character long
    assert len(new_id) == 1
    # The new ID should be deterministic
    assert new_id == "B"


def test_fix_alphafill(alphafill_P04406):
    fixed = fix_dict(
        alphafill_P04406,
        fixers=[
            "label_seq_id",
            "auth_seq_id",
        ],
    )
    # Check that we have fixed the unknown ids
    assert all([it not in [".", "?"] for it in fixed["_atom_site.label_seq_id"]])
    assert all([it not in [".", "?"] for it in fixed["_atom_site.auth_seq_id"]])

    # Check that biopython can now parse it
    assert is_biopython_parseable(fixed)

def test_fix_alt_id(alphafill_P04406):
    fixed = fix_dict(
        alphafill_P04406,
        fixers=[
            "alt_id",
        ],
    )
    assert all([it == "." for it in fixed["_atom_site.label_alt_id"]])


def test_find_changes():
    assert list(find_changes([1, 2, 3])) == [1, 2]
    assert list(find_changes([1, 1, 2, 2, 3, 3])) == [2, 4]

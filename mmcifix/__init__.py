import io
import sys
import typing
from math import inf
from typing import Union, Iterable, IO

from Bio.PDB import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from more_itertools import windowed
import click


def try_parse_int(s: str, default: Union[int, float]) -> Union[int, float]:
    """
    Tries to parse a string as an integer. If it fails, returns the default value
    """
    try:
        as_int = int(s)
        return as_int
    except ValueError:
        return default


def find_changes(l: list) -> Iterable[int]:
    """
    Returns a list of positions where the list has changed values
    """
    end = object()
    for pos, (a, b) in enumerate(windowed(l, 2, fillvalue=end)):
        if a != b:
            yield pos
        if b == end:
            break


def fix_cif_list(l: list, segments: Iterable[int]) -> list[str]:
    """
    Replaces any ? or . character in a cif file list with an unused integer, which is
    incremented whenever the index goes past a segment (generally when the chain ID changes)
    """
    as_integers = [try_parse_int(s, -inf) for s in l]
    replacement: int = max(as_integers)
    new = []
    for i, value in enumerate(l):

        if i in segments:
            # Increment the replacement value whenever we find a new ligand
            replacement = replacement + 1

        if value in [".", "?"]:
            # Replace empty values with an unused integer
            new.append(str(replacement))
        else:
            new.append(value)

    return new


def fix_auth_seq_id(struct: dict) -> dict:
    result = struct.copy()
    ligand_changes = list(find_changes(struct['_atom_site.label_asym_id']))
    result["_atom_site.auth_seq_id"] = fix_cif_list(struct["_atom_site.auth_seq_id"], ligand_changes)
    return result


def fix_label_seq_id(struct: dict) -> dict:
    result = struct.copy()
    ligand_changes = list(find_changes(struct['_atom_site.label_asym_id']))
    result["_atom_site.label_seq_id"] = fix_cif_list(struct["_atom_site.label_seq_id"], ligand_changes)
    return result


FIXERS = {
    "label_seq_id": fix_label_seq_id,
    "auth_seq_id": fix_auth_seq_id,
}


def fix_dict(structure: dict, fixers: list) -> dict:
    """
    Fixes an mmCIF file which has already been converted to dict
    """
    parsed = structure

    for fixer in fixers:
        parsed = FIXERS[fixer](parsed)

    return parsed


def fix_file(in_file: IO[str], out_file: IO[str], fixers: list[str]):
    """
    Fix an mmCIF file, and write the fixed structure to out_file
    """
    parsed = MMCIF2Dict(in_file)
    fixed = fix_dict(parsed, fixers)
    dict_to_file(fixed, out_file)

def dict_to_file(d: dict, file: IO[str]):
    mmcifio = MMCIFIO()
    mmcifio.set_dict(d)
    mmcifio.save(file)

@click.command("mmcifix")
@click.option("--fixer", type=click.Choice(FIXERS.keys()), multiple=True)
@click.argument("structure", type=click.File())
def main(structure: IO[str], fixer: list[str]):
    fix_file(in_file=structure, out_file=sys.stdout, fixers=fixer)

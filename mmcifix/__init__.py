import sys
from typing import IO

import click
from Bio.PDB import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from mmcifix.fixers import FIXERS


def fix_dict(structure: dict, fixers: list) -> dict:
    """
    Fixes an mmCIF file which has already been converted to dict
    """
    parsed = structure

    for fixer in fixers:
        parsed = FIXERS[fixer]().run(parsed)

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

@click.command("mmcifix", help="Fixes one or more issues with an mmCIF file")
@click.option("--fixer", type=click.Choice(FIXERS.keys()), multiple=True, help="Apply a fixer")
@click.argument("structure", type=click.File())
def main(structure: IO[str], fixer: list[str]):
    fix_file(in_file=structure, out_file=sys.stdout, fixers=fixer)

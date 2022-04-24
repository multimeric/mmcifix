from math import inf
from typing import Union, Iterable

from more_itertools import windowed


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
            yield pos + 1
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


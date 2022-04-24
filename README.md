# mmcifix

A collection of "fixer" functions for mmCIF files that can fix common issues with files in that format

## CLI

```
Usage: mmcifix [OPTIONS] STRUCTURE

  Fixes one or more issues with an mmCIF file

Options:
  --fixer [auth_seq_id|label_seq_id|database_id]
                                  Apply a fixer
  --help                          Show this message and exit.
```

## Python API

```python
from mmcifix import fix_file

with open("some_input.cif") as in_fp, open("some_output.cif") as out_fp:
    fix_file(in_file=in_fp, out_file=out_fp, fixers=["auth_seq_id"])
```

## Changes

### 0.2.1

* Fix bug in `database_id` fixer

### 0.2.0

* Convert all fixers into classes
* Add `database_id` fixer
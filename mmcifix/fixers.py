import string

from mmcifix.util import find_changes, fix_cif_list
from abc import abstractmethod, ABC
from inflection import underscore

FIXERS = {}

class FixerBase(ABC):
    @classmethod
    def name(cls) -> str:
        """
        The name of this fixer
        """
        return underscore(cls.__name__).lstrip("fix_")

    def description(self) -> str:
        """
        The description of this fixer
        """
        return ""

    def needed(self, struct: dict) -> bool:
        """
        Returns true if the given fixer is necessary for the given cif
        """
        return False

    @abstractmethod
    def run(self, struct: dict) -> dict:
        """
        Applies the fixer to a cif and returns the transformed version
        """

def fixer(fixer: FixerBase):
    """
    Decorator for registering fixers
    """
    FIXERS[fixer.name()] = fixer
    return fixer

@fixer
class FixAuthSeqId(FixerBase):
    def run(self, struct: dict) -> dict:
        result = struct.copy()
        ligand_changes = list(find_changes(struct['_atom_site.label_asym_id']))
        result["_atom_site.auth_seq_id"] = fix_cif_list(struct["_atom_site.auth_seq_id"], ligand_changes)
        return result

@fixer
class FixLabelSeqId(FixerBase):
    def run(self, struct: dict) -> dict:
        result = struct.copy()
        ligand_changes = list(find_changes(struct['_atom_site.label_asym_id']))
        result["_atom_site.label_seq_id"] = fix_cif_list(struct["_atom_site.label_seq_id"], ligand_changes)
        return result

@fixer
class FixDatabaseId(FixerBase):
    def run(self, struct: dict) -> dict:
        result = struct.copy()
        result["_database_2.database_id"] = ["PDB"]
        return result

@fixer
class FixAsymIdForPdb(FixerBase):
    """
    Renames the chain IDs to single-letter IDs for PDB compatability
    """
    def run(self, struct: dict) -> dict:
        used_chains = set(struct["_atom_site.auth_asym_id"]) | set(struct["_atom_site.label_asym_id"])
        valid_ids = set(string.ascii_uppercase)
        invalid_ids = used_chains - valid_ids

        if len(invalid_ids) == 0:
            # If all our IDs are already single letters, nothing needs to be done
            return struct

        available_ids = valid_ids - used_chains

        # We convert each invalid ID to an invalid ID via a dictionary
        mapping = dict(zip(sorted(invalid_ids), sorted(available_ids)))

        # Translate each invalid ID using the mapping
        output = struct.copy()
        output["_atom_site.auth_asym_id"] = [mapping[it] if it in mapping else it for it in output["_atom_site.auth_asym_id"]]
        output["_atom_site.label_asym_id"] = [mapping[it] if it in mapping else it for it in output["_atom_site.auth_asym_id"]]
        return output

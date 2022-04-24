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

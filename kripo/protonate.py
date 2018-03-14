import logging
import subprocess
from typing import Tuple, List

import pybel
from atomium.files.pdb import Pdb
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.structures import Model, Atom
from rdkit.Chem import MolToMolBlock, Mol, MolFromMolBlock, SanitizeMol

logger = logging.getLogger()


def protonate_protein(pdb_block: str, timeout: int=300, flags: List[str]=('-OH', '-HIS', '-NOHETh',)) -> Tuple[str, str]:
    """Passes a pdb block through reduce program to hydrogenate it.

    See http://kinemage.biochem.duke.edu/software/reduce.php .

    Expects `reduce` program to be in the PATH.

    Args:
        pdb_block: PDB block to hydrogenate
        timeout: Number of seconds to wait for completion
        flags: Flags to pass to reduce

    Returns:
        The hydrogenated pdb block and the standard error output of reduce
    """
    args = ['reduce'] + list(flags) + ['-']

    proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        out, err = proc.communicate(pdb_block.encode(), timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, err = proc.communicate()

    return out.decode(), err.decode()


def protonate_molecule(mol_in: Mol, ph=7.4) -> Mol:
    molblock_in = MolToMolBlock(mol_in)
    babel_mol = pybel.readstring('mol', molblock_in)
    babel_mol.OBMol.AddHydrogens(False, True, ph)
    molblock_out = babel_mol.write('mol')
    mol = MolFromMolBlock(molblock_out, removeHs=False, sanitize=False)
    try:
        SanitizeMol(mol)
    except ValueError:
        # Try again, but without ph correction
        babel_mol = pybel.readstring('mol', molblock_in)
        babel_mol.OBMol.AddHydrogens(False, False)
        molblock_out = babel_mol.write('mol')
        mol = MolFromMolBlock(molblock_out, removeHs=False, sanitize=False)
        SanitizeMol(mol)
    return mol


def fill_serial_numbers(pdb: Pdb):
    """The reduce program can add hydrogens to ligands, those hydrogens will have no atom serial numbers
    RDKit will give parse error on a PDB block with atoms without an atom serial number.

    This method adds serial ids and bonds those hydrogens to their heavy atom based on it's name

    Args:
        pdb: The pdb to fill serial ids in

    """
    model = pdb.model()
    max_serial_number = max([a.atom_id() for a in model.atoms() if a.atom_id()])
    for mol in model.molecules(generic=True):
        for a in mol.atoms(element='H'):
            if a.atom_id() != 0:
                continue
            max_serial_number += 1
            # a.atom_id() is not a setter, so set it using private prop
            a._id = max_serial_number
            hgrp = a.name()[1]
            bonded = False
            for heavy in mol.atoms(exclude='H'):
                oname = heavy.name()
                if len(oname) > 1 and oname[1] == hgrp:
                    logger.info('Binding {0}:{1} with {2}:{3}'.format(
                        a.atom_id(), a.name(), heavy.atom_id(), heavy.name()
                    ))
                    heavy.bond(a)
                    bonded = True
            if not bonded:
                logger.warning('Unable to bind {0}:{1} to heavy atom'.format(
                    a.atom_id(), a.name()
                ))


def protonate_pdb(pdb: Pdb) -> Pdb:
    """Hydrogenate a first model of PDB

    1. Passes generic molecules to protonate_ligand,  ligands which already contain hydrogens are skipped
    2. Take H in protonated ligand block and add to existing pdb
    3. Pass pdb to protonate_protein

    Args:
        pdb: The pdb to hydrogenate

    Returns:
        A pdb with hydrogens added
    """
    # Protonate whole pdb
    unprotonated_block = pdb.to_file_string()
    protonated_block = protonate_protein(unprotonated_block)
    protonated_pdb = pdb_dict_to_pdb(pdb_string_to_pdb_dict(protonated_block[0]))
    fill_serial_numbers(protonated_pdb)
    return protonated_pdb

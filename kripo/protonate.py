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


def protonate_ligand(pdb_block: str, ph=7.4) -> str:
    """Passes a pdb block through open babel to hydrogenate it.

    Uses the cli equivalent of `babel in.pdb out.pdb -h -p -ipdb  -opdb`.

    Args:
        pdb_block: PDB block to hydrogenate
        ph: Add hydrogens appropriate for this pH

    Returns:
        The hydrogenated pdb block
    """
    if pdb_block == '':
        return ''
    mol = pybel.readstring('pdb', pdb_block)
    mol.OBMol.AddHydrogens(False, True, ph)
    return mol.write('pdb')


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


def protonate(pdb: Pdb, aa=True, het=True) -> Pdb:
    """Hydrogenate a first model of PDB

    1. Passes generic molecules to protonate_ligand,  ligands which already contain hydrogens are skipped
    2. Take H in protonated ligand block and add to existing pdb
    3. Pass pdb to protonate_protein

    Args:
        pdb: The pdb to hydrogenate
        aa: If true then protonate protein
        het: If true then protonate ligands

    Returns:
        A pdb with hydrogens added
    """
    if het:
        model = pdb.model()

        # Protonate ligands
        ligands_model = Model()
        [ligands_model.add_molecule(l) for l in model.molecules(generic=True) if len(l.atoms(element='H')) == 0]
        ligands_block = ligands_model.to_file_string('pdb')
        protonated_ligands_block = protonate_ligand(ligands_block)
        protonated_ligands_pdb = pdb_dict_to_pdb(pdb_string_to_pdb_dict(protonated_ligands_block))
        protonated_ligands_model = protonated_ligands_pdb.model()

        # Add hydrogens of protonated_ligands_pdb to model
        add_hydrogens_from_ligand2model(model, protonated_ligands_model)

    if aa:
        # Protonate whole pdb
        unprotonated_block = pdb.to_file_string()
        protonated_block = protonate_protein(unprotonated_block)
        protonated_pdb = pdb_dict_to_pdb(pdb_string_to_pdb_dict(protonated_block[0]))
        fill_serial_numbers(protonated_pdb)
        return protonated_pdb
    else:
        return pdb


def add_hydrogens_from_ligand2model(model: Model, protonated_ligands_model: Model):
    """The ligand is protonated by openbabel, sadly openbabel renumbers the atom serial number.

    This method creates hydrogens in the model based on the ones in the protonated ligand model

    Args:
        model: Model of a PDB with unprotonated ligands
        protonated_ligands_model: Model with protonated ligands

    """
    max_serial_number = max([a.atom_id() for a in model.atoms() if a.atom_id()])
    h_atom_id = max_serial_number + 1
    for l in protonated_ligands_model.molecules(generic=True):
        for a in l.atoms():
            for ba in a.bonded_atoms():
                if ba.element() == 'H':
                    newh = Atom(element='H',
                                x=ba.location()[0],
                                y=ba.location()[1],
                                z=ba.location()[2],
                                atom_id=h_atom_id,
                                name=ba.name())
                    mol = model.molecule(name=l.name())
                    heavy = mol.atom(name=a.name())
                    model.add_atom(newh)
                    mol.add_atom(newh)
                    heavy.bond(newh)
                    # newh.bond(heavy)
                    h_atom_id += 1

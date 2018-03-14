import logging
from typing import List, Dict

import atomium
from atomium.structures import Model
from atomium.structures.molecules import Molecule
from atomium.files.pdb import Pdb
from rdkit.Chem import Mol

from .ligand import Ligand
from .protonate import protonate_pdb, protonate_molecule

logger = logging.getLogger(__name__)

"""Unwanted heteros, like solvents, metals, sugars, etc."""
UNWANTED_HETEROS = {
    # REMOVE SOLVENTS, METALS AND ALT LOCATIONS
    'MG', 'K', 'MN', 'NA', 'ZN', 'MG', 'CA', 'CD', 'HG', 'FE', 'CL', 'CU', '4MO', 'AU', 'IOD', 'SM', 'AF3', 'MOO',
    'HOH', 'ACE', 'DTT', 'SO4', 'PO3', 'PO4', 'MG8', 'BU3', 'GOL', 'DMS', 'FMT', 'LPA', 'F3S', 'ACE', 'NH2', 'CO',
    # HEM as seed
    'EDO',
    '1FH', '2FH', 'CCH', 'CLN', 'COH', 'DDH', 'DEU', 'DHE', 'FDD', 'FDE', 'FEC', 'FMI', 'HAS', 'HE6',
    'HEA', 'HEB', 'HEG', 'HEM', 'HEO', 'HEV', 'HIF', 'HKL', 'HME', 'HNI', 'MNH', 'MNR', 'MP1', 'VEA',
    'VER', 'ZEM', 'ZNH', 'HMG', '6HE', '7HE', 'HCO', 'HDD', 'HDM', 'HE5', 'HEC', 'HES', 'HFM', 'PC3',
    # S-F cluster
    'B51', 'SF4', 'WCC', 'NFS', 'CFM',
    'SUC', 'GAL', 'SUC', 'NAG', '1G6', '5AX', 'A2G', 'AB6', 'AHO', 'AMU', 'AMV', 'AS5', 'BGN', 'BM3',
    'CBS', 'CTO', 'DR2', 'EAG', 'FSM', 'HSQ', 'LXB', 'LXZ', 'M5G', 'MA8', 'MAG', 'MGC', 'MMR', 'MUB',
    'NA1', 'NAA', 'NAG', 'NDG', 'NGA', 'NGZ', 'NLC', 'TNR', 'TNY', 'WZ5', 'GAL', 'BGC', '1GL', '289',
    '293', '2DG', '2GS', '3MG', 'AAL', 'ABE', 'ACX', 'AHR', 'ALL', 'AMG', 'ARA', 'ARW', 'AXR', 'B2G',
    'B4G', 'B7G', 'B8D', 'BCD', 'BDF', 'BDR', 'BGC', 'BGL', 'BHE', 'BHG', 'BMA', 'BNG', 'BXP', 'BXX',
    'BXY', 'CBI', 'CDR', 'CE5', 'CE6', 'CE8', 'CEX', 'CEY', 'CTR', 'CTT', 'DA8', 'DDA', 'DDL',
    'DEG', 'DEL', 'DFR', 'DLF', 'DLG', 'DMU', 'DOM', 'DR5', 'DRI', 'DXI', 'EBG', 'EPG', 'FCA',
    'FCB', 'FRU', 'FUB', 'FUC', 'FUL', 'G4D', 'G6D', 'GAL', 'GL0', 'GLA', 'GLC', 'GLD', 'GMH',
    'GU8', 'GU9', 'GUP', 'GXL', 'GYP', 'GZL', 'HSG', 'HSH', 'HSY', 'JZR', 'KHO', 'LAK', 'LAT',
    'LB2', 'LBT', 'LDY', 'LFR', 'LMT', 'LMU', 'LRH', 'LXC', 'M13', 'M3M', 'M5S', 'MAB', 'MAL',
    'MAN', 'MBG', 'MDA', 'MDM', 'MFB', 'MFU', 'MGL', 'MLR', 'MMA', 'MRP', 'MT7', 'MTT', 'MXY',
    'MXZ', 'NGR', 'OPM', 'PSV', 'QKH', 'RAE', 'RAM', 'RAO', 'RGG', 'RIB', 'RIP', 'RM4', 'SHD',
    # sugar
    'SUC', 'SWE', 'TRE', 'TYV', 'UMQ', 'XLF', 'XLM', 'XYP', 'XYS', 'XYZ',
}

# PROTEIN PREPARATION CONSTANTS
"""Maximum size for ligand structures"""
LIGAND_MAX_MASS = 800  #
"""Minimum size for ligand structures"""
LIGAND_MIN_MASS = 50  #
"""Maximum allowed distance of ligand to protein"""
MAX_CONTACT_DISTANCE = 2.5


class NoLigands(ValueError):
    pass


class PdbDumpError(TypeError):
    pass


def ligands(pdb: Pdb, ligand_expo: Dict[str, Mol]) -> List[Ligand]:
    """Ligands of a pdb

    Args:
        pdb: The pdb
        ligand_expo: Dictionary with molecules of ligand expo

    Raises:
        NoLigands: When PDB has no ligands

    Returns:
        List of ligands, ordered by name
    """
    model = pdb.model()
    ligs = {}
    for amol in model.molecules(generic=True):
        amol_id = amol.molecule_id()
        lig_id = pdb.code().lower() + '_' + amol.name() + '_1_' + amol_id[0] + '_' + amol_id[1:]
        try:
            lig = ligand_expo[lig_id]
            plig = protonate_molecule(lig)
            ligs[lig_id] = Ligand(amol, plig)
        except KeyError:
            logger.warning('Unable to find {0} in ligand expo db, skipping'.format(lig_id))
            pass

    if not ligs:
        raise NoLigands()
    return sorted(ligs.values(), key=lambda l: l.name())


def pdb_from_file(path, hydrogenate=True, clean=True) -> Pdb:
    """Construct Pdb object from a PDB file

    Args:
        path: Path to PDB file
        hydrogenate: Whether to add hydrogens
        clean: Whether to remove unwanted molecules

    Returns:
        Protein databank entry
    """
    pdb = atomium.pdb_from_file(path)
    return pdb_from_atomium_pdb(pdb, hydrogenate, clean)


def pdb_from_atomium_pdb(pdb: Pdb, hydrogenate=True, clean=True) -> Pdb:
    """Construct a PDB entry from a Atomium pdb

    Args:
        pdb: Atomium PDB entry
        hydrogenate: Whether to add hydrogens
        clean: Whether to remove unwanted molecules

    Returns:
        A PDB entry which can optionally be hydrogenated and have it unwanted molecules removed.
    """
    if clean:
        remove_unwanted_molecules(pdb)
    if hydrogenate:
        try:
            pdb = protonate_pdb(pdb)
        except TypeError as e:
            raise PdbDumpError(pdb) from e
    if clean:
        remove_non_contacting_molecules(pdb)
    return pdb


def remove_unwanted_molecules(pdb: Pdb):
    """Remove unwanted molecules from model

    Cleans pdb by removing molecules which:
        * Have name in UNWANTED_HETEROS list
        * Is out side LIGAND_MAX_MASS..LIGAND_MIN_MASS mass range
        * Have name already processed (aka removes duplicates)

    Removing is done in-place.

    Args:
        pdb: Atomium PDB entry containing possible unwanted molecules

    """
    unique_names = set()
    model = pdb.model()
    for mol in sorted(model.molecules(generic=True), key=lambda m: m.molecule_id()):
        is_unwanted = mol.name() in UNWANTED_HETEROS
        in_mass_range = LIGAND_MIN_MASS < mol.mass() < LIGAND_MAX_MASS
        seen_before = mol.name() in unique_names
        if is_unwanted or not in_mass_range or seen_before:
            model.remove_molecule(mol)
        else:
            unique_names.add(mol.name())


def remove_non_contacting_molecules(pdb: Pdb):
    """Remove unwanted molecules from model

    Cleans pdb by removing molecules which:
        * Is more then MAX_CONTACT_DISTANCE away from protein

    Removing is done in-place.

    Args:
        pdb: Atomium PDB entry containing possible unwanted molecules

    """
    model = pdb.model()
    for mol in sorted(model.molecules(generic=True), key=lambda m: m.molecule_id()):
        in_contact_with_protein = ligand_contacts_protein(mol, model)
        if not in_contact_with_protein:
            try:
                model.remove_molecule(mol)
            except KeyError:
                # in 1efr was unable to delete atom with key 22969, ignore error
                pass


def ligand_contacts_protein(ligand: Molecule, model: Model) -> bool:
    """Determines if molecule is in contact with protein

    Args:
        ligand: The molecule
        model: Model containing protein chains

    Returns:
        True if in contact and false when not in contact
    """
    dist = 999.0
    for latom in ligand.atoms():
        for chain in model.chains():
            for patom in chain.atoms():
                ndist = latom.distance_to(patom)
                if ndist < dist:
                    dist = ndist
    return dist < MAX_CONTACT_DISTANCE

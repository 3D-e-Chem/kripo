from typing import List

import atomium
from atomium.structures.molecules import Molecule

from .ligand import Ligand
from .protonate import protonate

"""Unwanted heteros, like solvents, metals, sugars, etc."""
UNWANTED_HETEROS = {
    # REMOVE SOLVENTS, METALS AND ALT LOCATIONS
    'MG', 'K', 'MN', 'NA', 'ZN', 'MG', 'CA', 'CD', 'HG', 'FE', 'CL', 'CU', '4MO', 'AU', 'IOD', 'SM', 'AF3', 'MOO',
    'HOH', 'ACE', 'DTT', 'SO4', 'PO3', 'PO4', 'MG8', 'BU3', 'GOL', 'DMS', 'FMT', 'LPA', 'F3S', 'ACE', 'NH2',
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


class Pdb:
    """Protein databank entry

    Attributes:
        path: File path of entry
        pdb (atomium.files.pdb.Pdb): Atomium PDB object
        model (atomium.structures.models.Model): Atomium model of first model in pdb

    """
    def __init__(self, path, hydrogenate=True, clean=True):
        """Construct PDB from file path

        Protonates the protein and ligands.

        Cleans pdb by removing molecules which:
        * Have name in UNWANTED_HETEROS list
        * Is out side LIGAND_MAX_MASS..LIGAND_MIN_MASS mass range
        * Is more then MAX_CONTACT_DISTANCE away from protein
        * Have name already processed (aka removes duplicates)

        Args:
            path (str): Path to PDB file

        """
        self.path = path
        self.pdb = atomium.pdb_from_file(path)
        self.model = self.pdb.model()
        if hydrogenate:
            self.model = protonate(self.pdb.model())
        if clean:
            self._clean()

    def _clean(self):
        unique_names = set()
        for mol in sorted(self.model.molecules(generic=True), key=lambda m: m.molecule_id()):
            is_unwanted = mol.name() in UNWANTED_HETEROS
            in_mass_range = LIGAND_MIN_MASS < mol.mass() < LIGAND_MAX_MASS
            in_contact_with_protein = self._min_distance(mol) < MAX_CONTACT_DISTANCE
            seen_before = mol.name() in unique_names
            if is_unwanted or not in_mass_range or not in_contact_with_protein or seen_before:
                self.model.remove_molecule(mol)
            else:
                unique_names.add(mol.name())

    def _min_distance(self, ligand: Molecule):
        """Calculate the minimum distance any atom of a ligand is from any atom of the protein chains"""
        dist = 999.0
        for latom in ligand.atoms():
            for chain in self.model.chains():
                for patom in chain.atoms():
                    ndist = latom.distance_to(patom)
                    if ndist < dist:
                        dist = ndist
        return dist

    def ligands(self) -> List[Ligand]:
        """Ligands of pdb

        Returns:
            Unique list of ligands

        Raises:
            NoLigands: When PDB has no ligands, after cleaning

        """
        ligs = {mol.name(): Ligand(mol) for mol in self.model.molecules(generic=True, water=False)}
        if not ligs:
            raise NoLigands()
        return list(ligs.values())

    def pdb_block(self):
        return self.model.to_file_string('pdb')

    def code(self) -> str:
        """Code of pdb

        Returns:
            Code of pdb

        """
        return self.pdb.code()

    def scop_classification(self):
        # TODO see if needed and if so then implement
        pass


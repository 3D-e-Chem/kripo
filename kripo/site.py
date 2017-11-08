from collections import Counter

from atomium.structures.atoms import Atom
from atomium.structures.chains import Site
from rdkit.Chem import AddHs, MolFromPDBBlock, Mol
from rdkit.Chem.AllChem import EmbedMolecule, ETKDG


def protonate(unprotonated_site: Site) -> Site:
    unprotonated_pdb_block = to_pdb_block(unprotonated_site)

    unprotonated_mol = MolFromPDBBlock(unprotonated_pdb_block)

    protonated_mol = protonate_rdkit_molecule(unprotonated_mol)

    return protonate_site_from_rdkit(unprotonated_site, protonated_mol)


def to_pdb_block(site: Site):
    block = site.to_file_string('pdb')
    return block


def protonate_rdkit_molecule(unprotonated_mol: Mol) -> Mol:
    protonated_mol = AddHs(unprotonated_mol)
    EmbedMolecule(protonated_mol, ETKDG())
    # TODO check positions of non-hydrogens have not been changed
    return protonated_mol


def protonate_site_from_rdkit(unprotonated_site: Site, protonated_mol: Mol) -> Site:
    site = unprotonated_site
    # In PDB each atom must have a unique id, use max+1 for new hydrogens
    hatom_id = max([a.atom_id() for a in site.atom().model().atoms()]) + 1

    positions = protonated_mol.GetConformer().GetPositions()

    for atom in protonated_mol.GetAtoms():
        if not atom.GetPDBResidueInfo():
            # atom is an added hydrogen

            # Add hydrogen atom
            hatom_pos = positions[atom.GetIdx()]
            hatom = Atom(element=atom.GetSymbol(),
                         x=hatom_pos[0],
                         y=hatom_pos[1],
                         z=hatom_pos[2],
                         atom_id=hatom_id)
            site.add_atom(hatom)
            hatom_id += 1

            # Bond hydrogen atom to residue atom
            res_atom_id = atom.GetBonds()[0].GetBeginAtom().GetPDBResidueInfo().GetSerialNumber()
            res_atom = site.atom(atom_id=res_atom_id)
            res_atom.bond(hatom)

    return site


def chain_of_site(site: Site):
    # A site can be in different chains, take the chain most residues belong to
    most_common = Counter([r.chain().chain_id() for r in site.residues()]).most_common(1)
    return most_common[0][0]

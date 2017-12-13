# -*- coding: utf-8 -*-

"""Console script for kripo."""

import click

from kripodb.db import FingerprintsDb, FragmentsDb
from kripodb.pharmacophores import PharmacophoresDb

from .fingerprint.threepoint import BIT_INFO
from .fragment import Fragment
from .ligand import Ligand
from .pharmacophore import from_fragment, NoFeatures
from .pdb import pdb_from_file, ligands
from .site import chain_of_site


@click.group()
@click.version_option()
def main():
    """Key Representation of Interaction in POckets"""
    pass


@main.command()
@click.argument('pdbs', type=click.File('rt'))
@click.argument('fragments', type=click.Path(writable=True, dir_okay=False))
@click.argument('pharmacophores', type=click.Path(writable=True, dir_okay=False))
@click.argument('fingerprints', type=click.Path(writable=True, dir_okay=False))
def generate(pdbs, fragments, pharmacophores, fingerprints):
    """Generate fragments, pharmacophores and fingerprints for given pdb files

    * PDBS, Name of file with a PDB filename on each line

    * FRAGMENTS, Fragments database output filename

    * PHARMACOPHORES, Pharmacophores database output file name

    * FINGERPRINTS, Fingerprints database output file name
    """
    fragments_db = FragmentsDb(fragments)
    pharmacophores_db = PharmacophoresDb(pharmacophores, mode='w')
    pharmacophore_points = pharmacophores_db.points
    fingerprints_db = FingerprintsDb(fingerprints)
    fingerprints_dict = fingerprints_db.as_dict(len(BIT_INFO))

    for pdb_fn in pdbs:
        pdb_fn = pdb_fn.strip()
        generate_from_pdb(pdb_fn, fragments_db, pharmacophore_points, fingerprints_dict)

    fragments_db.close()
    pharmacophores_db.close()
    fingerprints_db.close()


def generate_from_pdb(pdb_fn, fragments_db, pharmacophore_points, fingerprints_dict):
    """Generate pharmacophore fingerprints from a pdb and store them

    Args:
        pdb_fn (str): Filename of PDB file
        fragments_db (kripodb.db.FragmentsDb): Fragments database
        pharmacophore_points (kripodb.pharmacophores.PharmacophorePointsTable): Pharmacophores database
        fingerprints_dict (kripodb.db.IntbitsetDict): Fingerprints db dictionary

    """
    click.echo('Parsing {0}'.format(pdb_fn))
    pdb = pdb_from_file(pdb_fn)
    for ligand in ligands(pdb):
        for frag_nr, fragment in enumerate(ligand.fragments(), 1):
            generate_from_fragment(pdb,
                                   ligand,
                                   fragment,
                                   frag_nr,
                                   fragments_db,
                                   pharmacophore_points,
                                   fingerprints_dict)


def generate_from_fragment(pdb, ligand, fragment, frag_nr, fragments_db, pharmacophore_points, fingerprints_dict):
    try:
        frag_id = build_frag_id(pdb, ligand, frag_nr)
        fragment.name = frag_id
        click.echo('Generating pharmacophore fingerprint for {0}'.format(frag_id))

        pharmacophore = from_fragment(fragment)
        fingerprint = pharmacophore.fingerprint()

        add_fragment2db(pdb, ligand, frag_nr, fragment, fragments_db)
        add_pharmacophore2db(pharmacophore_points, frag_id, pharmacophore)
        fingerprints_dict[frag_id] = fingerprint
    except NoFeatures:
        msg = 'Fragment {0} of ligand {1} of pdb {2} ' \
              'contains no pharmacophore features, skipping'.format(frag_nr, ligand.id(), pdb.code())
        click.echo(msg)


def build_frag_id(thepdb, ligand, frag_nr):
    pdb_code = thepdb.code().lower()
    het_code = ligand.name()
    return pdb_code + '_' + het_code + '_frag' + str(frag_nr)


def add_fragment2db(thepdb, ligand: Ligand, frag_nr, fragment: Fragment, fragments_db):
    # TODO move to kripodb, should not use sql here
    frag_id = fragment.name
    pdb_code = thepdb.code().lower()
    het_code = ligand.name()
    het_chain = ligand.chain()
    het_seq_nr = ligand.seq_nr()
    # A site can be in different chains, take the chain most residues belong to
    prot_chain = chain_of_site(fragment.site())
    hash_code = fragment.hash_code()
    atom_codes = ','.join(fragment.atom_names())
    nr_r_groups = fragment.nr_r_groups()

    fragment_sql = '''INSERT INTO fragments (
               frag_id,
               pdb_code,
               prot_chain,
               het_code,
               frag_nr,
               atom_codes,
               hash_code,
               het_chain,
               het_seq_nr,
               nr_r_groups
           ) VALUES (
               :frag_id,
               :pdb_code,
               :prot_chain,
               :het_code,
               :frag_nr,
               :atom_codes,
               :hash_code,
               :het_chain,
               :het_seq_nr,
               :nr_r_groups
           )'''
    fragment_row = {
        'frag_id': frag_id,
        'pdb_code': pdb_code,
        'prot_chain': prot_chain,
        'het_code': het_code,
        'frag_nr': frag_nr,
        'atom_codes': atom_codes,
        'hash_code': hash_code,
        'het_chain': het_chain,
        'het_seq_nr': het_seq_nr,
        'nr_r_groups': nr_r_groups,
    }
    fragments_db.cursor.execute(fragment_sql, fragment_row)
    fragments_db.commit()

    fragments_db.add_molecule(fragment.molecule)


def add_pharmacophore2db(pharmacophore_points, frag_id, pharmacophore):
    points = []
    for feature in pharmacophore.features:
        point = (feature.kind, feature.position[0], feature.position[1], feature.position[2])
        points.append(point)
    point_ids = range(len(points))
    pharmacophore_points.add_fragment(frag_id, point_ids, points)

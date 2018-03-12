# -*- coding: utf-8 -*-

"""Console script for kripo."""
import gzip
from sqlite3 import IntegrityError

import click

from kripodb.db import FingerprintsDb, FragmentsDb, FastInserter
from kripodb.pharmacophores import PharmacophoresDb
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier
from rdkit.Chem.rdmolops import SanitizeMol

from .ligandexpodb import LigandExpoDb
from .pharmacophore import Feature, Pharmacophore
from .generator import generate_from_pdb
from .fingerprint.threepoint import BIT_INFO, from_pharmacophore
from .pdb import PdbDumpError, NoLigands, UNWANTED_HETEROS


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
@click.option('--fuzzy-factor',
              type=int,
              help='Number of bins below/above actual bin to include in fingerprint',
              default=1,
              show_default=True)
@click.option('--fuzzy-shape',
              type=click.Choice(('all', 'one', 'v1')),
              help='Shape of bins around actual bin',
              default='all',
              show_default=True)
@click.option('--fragmentation/--no-fragmentation',
              help='Fragment ligands',
              default=True,
              show_default=True)
@click.option('--ligand-expob',
              help='Ligand expo database file name',
              type=click.Path(dir_okay=False),
              default='ligand-expo.db')
def generate(pdbs, fragments, pharmacophores, fingerprints, fuzzy_factor, fuzzy_shape, fragmentation, ligand_expo):
    """Generate fragments, pharmacophores and fingerprints for given pdb files

    * PDBS, Name of file with a PDB filename on each line

    * FRAGMENTS, Fragments database output filename

    * PHARMACOPHORES, Pharmacophores database output file name

    * FINGERPRINTS, Fingerprints database output file name
    """

    with FragmentsDb(fragments) as fragments_db, \
            PharmacophoresDb(pharmacophores, mode='a') as pharmacophores_db, \
            FingerprintsDb(fingerprints) as fingerprints_db, \
            LigandExpoDb(ligand_expo) as ligand_expo_db:
        pharmacophore_points = pharmacophores_db.points
        fingerprints_dict = fingerprints_db.as_dict(len(BIT_INFO))
        ligand_expo_dict = ligand_expo_db.as_dict()
        for pdb_fn in pdbs:
                pdb_fn = pdb_fn.strip()
                try:
                    generate_from_pdb(pdb_fn, fragments_db, pharmacophore_points, fingerprints_dict, ligand_expo_dict, fuzzy_factor, fuzzy_shape, fragmentation)
                except PdbDumpError:
                    msg = 'Unable to dump {0}, skipping'.format(pdb_fn)
                    click.secho(msg, bold=True)
                except NoLigands:
                    msg = 'No ligands found in {0}, skipping'.format(pdb_fn)
                    click.secho(msg, bold=True)


@main.group(name='pharmacophores')
def pharmacophores_group():
    """Commands on pharmacophores"""
    pass


@pharmacophores_group.command(name='fingerprints')
@click.argument('pharmacophores', type=click.Path(dir_okay=False))
@click.argument('fingerprints', type=click.Path(dir_okay=False, writable=True))
@click.option('--fuzzy_factor',
              type=int,
              help='Number of bins below/above actual bin to include in fingerprint',
              default=1,
              show_default=True)
@click.option('--fuzzy_shape',
              type=click.Choice(('all', 'one', 'v1')),
              help='Shape of bins around actual bin',
              default='all',
              show_default=True)
def pharmacophore2fingerprints(pharmacophores, fingerprints, fuzzy_factor, fuzzy_shape):
    """Generate fingerprints from pharmacophores

    * PHARMACOPHORES, Pharmacophores input data file (*.h5)

    * FINGERPRINTS, Fingerprints output data file

    """
    with PharmacophoresDb(pharmacophores, mode='r') as pharmacophores_db, \
            FingerprintsDb(fingerprints) as fingerprints_db:
        fingerprints_dict = fingerprints_db.as_dict(len(BIT_INFO))
        for frag_id, points in pharmacophores_db:
            features = [Feature(p[0], (p[1], p[2], p[3])) for p in points]
            pharmacophore = Pharmacophore(features)
            fingerprint = from_pharmacophore(pharmacophore, fuzzy_factor=fuzzy_factor, fuzzy_shape=fuzzy_shape)
            fingerprints_dict[frag_id] = fingerprint


@main.group(name='ligands')
def ligands_group():
    """Commands on ligands"""
    pass


@ligands_group.command('import')
@click.argument('ligandssdf', type=click.File('rb'))
@click.argument('ligandsdb', type=click.Path(dir_okay=False, writable=True))
def import_ligands(ligandsdb, ligandssdf):
    """Convert ligand expo sdf to database

    On http://ligand-expo.rcsb.org/ld-download.html download http://ligand-expo.rcsb.org/dictionaries/all-sdf.sdf.gz
    and use as LIGANDSSDF
    """
    sdf_fn = gzip.open(ligandssdf)
    gzsuppl = ForwardSDMolSupplier(sdf_fn, sanitize=False, removeHs=False)
    with LigandExpoDb(ligandsdb) as db:
        cursor = db.cursor
        with FastInserter(cursor):
            for mol in gzsuppl:
                if mol is None:
                    continue
                mol_name = mol.GetProp('_Name')
                # PDB ID _ Component ID _ Model No. _ Chain ID _ Residue No.
                cols = mol_name.split('_')
                if cols[1] in UNWANTED_HETEROS:
                    click.secho('Unwanted ligand ' + cols[1] + ', skipping', fg='yellow')
                    continue
                try:
                    SanitizeMol(mol)
                except ValueError:
                    click.secho('Unable to sanitize ' + mol_name + ', skipping', fg='yellow')
                    continue
                lig_id = cols[0] + '_' + cols[1] + '_' + cols[2] + '_' + cols[3] + '_' + cols[4]
                try:
                    cursor.execute('INSERT INTO ligands VALUES (?, ?)', (lig_id, mol,))
                except IntegrityError:
                    click.secho('Duplicate ' + mol_name + ', skipping', fg='yellow')

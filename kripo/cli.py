# -*- coding: utf-8 -*-

"""Console script for kripo."""
import gzip
from io import BytesIO
from sqlite3 import IntegrityError

import click
import requests

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
@click.option('--ligand-expo',
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


    Any ligand with zero atoms is downloaded from https://www.rcsb.org/pdb/download/downloadLigandFiles.do
    """
    sdf_fn = gzip.open(ligandssdf)
    gzsuppl = ForwardSDMolSupplier(sdf_fn, sanitize=False, removeHs=False)
    with LigandExpoDb(ligandsdb) as db:
        cursor = db.cursor
        ids2fetch = {}
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
                pdb_code = cols[0]
                het_code = cols[1]
                lig_id = pdb_code + '_' + het_code + '_' + cols[2] + '_' + cols[3] + '_' + cols[4]
                if mol.GetNumAtoms() == 0:
                    if het_code in ids2fetch:
                        ids2fetch[het_code].append(pdb_code)
                    else:
                        ids2fetch[het_code] = [pdb_code]
                    click.secho('Molecule {0} contains no atoms, fetching from alternate download'.format(lig_id), fg='yellow')
                    continue
                try:
                    SanitizeMol(mol)
                except ValueError:
                    click.secho('Unable to sanitize ' + lig_id + ', skipping', fg='yellow')
                    cursor.execute('INSERT INTO corrupt_ligands VALUES (?, ?)', (lig_id, 'sanitize_error'))
                    continue
                try:
                    cursor.execute('INSERT INTO ligands VALUES (?, ?)', (lig_id, mol,))
                except IntegrityError:
                    click.secho('Duplicate ' + lig_id + ', skipping', fg='yellow')
            click.secho('Downloading sdf of molecules with no atoms from alternate location')
            hets = []
            pdbs = []
            for het, het_pdbs in sorted(ids2fetch.items(), key=lambda d: len(d[1]), reverse=True):
                hets.append(het)
                pdbs.extend(het_pdbs)
                if 3 * len(hets) + 4 * len(pdbs) > 200:
                    fetch_chunk(cursor, hets, pdbs)
                    hets = []
                    pdbs = []
            if hets:
                fetch_chunk(cursor, hets, pdbs)


def fetch_chunk(cursor, hets, pdbs):
    for mol in fetch_ligand_sdf(hets, pdbs):
        mol_name = mol.GetProp('_Name')
        # '2ZFS_12U_A_501'
        cols = mol_name.split('_')
        lig_id = cols[0].lower() + '_' + cols[1] + '_1_' + cols[2] + '_' + cols[3]
        if mol.GetNumAtoms() == 0:
            click.secho('Molecule {0} contains no atoms, skipping'.format(lig_id), fg='yellow')
            cursor.execute('INSERT INTO corrupt_ligands VALUES (?, ?)', (lig_id, 'zero_atoms'))
            continue
        try:
            SanitizeMol(mol)
        except ValueError:
            click.secho('Unable to sanitize ' + lig_id + ', skipping', fg='yellow')
            try:
                cursor.execute('INSERT INTO corrupt_ligands VALUES (?, ?)', (lig_id, 'sanitize_error'))
            except IntegrityError:
                continue
        try:
            cursor.execute('INSERT INTO ligands VALUES (?, ?)', (lig_id, mol,))
        except IntegrityError:
            click.secho('Duplicate ' + lig_id + ', skipping', fg='yellow')


def fetch_ligand_sdf(hets, pdbs):
    params = {
        'ligandIdList': ','.join(hets),
        'structIdList':  ','.join(pdbs),
        'instanceType': 'all',
        'excludeUnobserved': 'false',
        'includeHydrogens': 'true',
    }
    url = 'https://www.rcsb.org/pdb/download/downloadLigandFiles.do'
    response = requests.get(url, params=params)
    sdf_fn = BytesIO(response.content)
    gzsuppl = ForwardSDMolSupplier(sdf_fn, sanitize=False, removeHs=False)
    for mol in gzsuppl:
        if mol is None:
            continue
        yield mol

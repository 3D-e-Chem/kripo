# -*- coding: utf-8 -*-

"""Console script for kripo."""

import click

from kripodb.db import FingerprintsDb, FragmentsDb
from kripodb.pharmacophores import PharmacophoresDb

from .generator import generate_from_pdb
from .fingerprint.threepoint import BIT_INFO
from .pdb import PdbDumpError, NoLigands


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

    with FragmentsDb(fragments) as fragments_db, \
            PharmacophoresDb(pharmacophores, mode='a') as pharmacophores_db, \
            FingerprintsDb(fingerprints) as fingerprints_db:
        pharmacophore_points = pharmacophores_db.points
        fingerprints_dict = fingerprints_db.as_dict(len(BIT_INFO))
        for pdb_fn in pdbs:
                pdb_fn = pdb_fn.strip()
                try:
                    generate_from_pdb(pdb_fn, fragments_db, pharmacophore_points, fingerprints_dict)
                except PdbDumpError:
                    msg = 'Unable to dump {0}, skipping'.format(pdb_fn)
                    click.secho(msg, bold=True)
                except NoLigands:
                    msg = 'No ligands found in {0}, skipping'.format(pdb_fn)
                    click.secho(msg, bold=True)



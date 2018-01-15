# -*- coding: utf-8 -*-

"""Console script for kripo."""

import click

from kripodb.db import FingerprintsDb, FragmentsDb
from kripodb.pharmacophores import PharmacophoresDb

from .pharmacophore import Feature, Pharmacophore
from .generator import generate_from_pdb
from .fingerprint.threepoint import BIT_INFO, from_pharmacophore
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
@click.option('--fuzzy_factor',
              type=int,
              help='Number of bins below/above actual bin to include in fingerprint',
              default=1)
def generate(pdbs, fragments, pharmacophores, fingerprints, fuzzy_factor):
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
                    generate_from_pdb(pdb_fn, fragments_db, pharmacophore_points, fingerprints_dict, fuzzy_factor)
                except PdbDumpError:
                    msg = 'Unable to dump {0}, skipping'.format(pdb_fn)
                    click.secho(msg, bold=True)
                except NoLigands:
                    msg = 'No ligands found in {0}, skipping'.format(pdb_fn)
                    click.secho(msg, bold=True)


@main.group(name='pharmacophores')
def pharmacophores_group():
    pass


@pharmacophores_group.command(name='fingerprints')
@click.argument('pharmacophores', type=click.Path(dir_okay=False))
@click.argument('fingerprints', type=click.Path(dir_okay=False, writable=True))
@click.option('--fuzzy_factor',
              type=int,
              help='Number of bins below/above actual bin to include in fingerprint',
              default=1)
def pharmacophore2fingerprints(pharmacophores, fingerprints, fuzzy_factor):
    """Generate fingerprints from pharmacophores

    * PHARMACOPHORES, Pharmacophores input data file

    * FINGERPRINTS, Fingerprints output data file

    """
    with PharmacophoresDb(pharmacophores, mode='r') as pharmacophores_db, \
            FingerprintsDb(fingerprints) as fingerprints_db:
        fingerprints_dict = fingerprints_db.as_dict(len(BIT_INFO))
        for frag_id, points in pharmacophores_db:
            features = [Feature(p[0], (p[1], p[2], p[3])) for p in points]
            pharmacophore = Pharmacophore(features)
            fingerprint = from_pharmacophore(pharmacophore, fuzzy_factor=fuzzy_factor)
            fingerprints_dict[frag_id] = fingerprint

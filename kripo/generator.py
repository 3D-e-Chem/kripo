from sqlite3 import IntegrityError

import click

from kripo.fragment import Fragment
from kripo.ligand import RdkitParseError, AtomiumParseError, Ligand
from kripo.pdb import pdb_from_file, ligands
from kripo.pharmacophore import from_fragment, NoFeatures
from kripo.site import chain_of_site


def generate_from_pdb(pdb_fn, fragments_db, pharmacophore_points, fingerprints_dict, fuzzy_factor, fragmentation=True):
    """Generate pharmacophore fingerprints from a pdb and store them

    Args:
        pdb_fn (str): Filename of PDB file
        fragments_db (kripodb.db.FragmentsDb): Fragments database
        pharmacophore_points (kripodb.pharmacophores.PharmacophorePointsTable): Pharmacophores database
        fingerprints_dict (kripodb.db.IntbitsetDict): Fingerprints db dictionary
        fuzzy_factor (int): The fuzzy factor
        fragmentation (bool): When true the ligand is fragmented and the whole ligand
            and it's fragments produce pharmacophores. When false only the whole ligand produced a single pharmacophore.

    """
    click.echo('Parsing {0}'.format(pdb_fn))
    pdb = pdb_from_file(pdb_fn)
    for ligand in ligands(pdb):
        click.echo('Ligand {0}'.format(ligand.name()))
        if is_ligand_stored(fragments_db, pdb.code(), ligand.name()):
            msg = 'Ligand {0} of pdb {1} already present, skipping'.format(ligand.name(), pdb.code())
            click.secho(msg, bold=True)
            continue
        try:
            if fragmentation:
                fragments = ligand.fragments()
            else:
                fragments = [ligand.as_fragment()]
            for frag_nr, fragment in enumerate(fragments, 1):
                click.echo('Fragment {0}'.format(frag_nr))
                generate_from_fragment(pdb,
                                       ligand,
                                       fragment,
                                       frag_nr,
                                       fragments_db,
                                       pharmacophore_points,
                                       fingerprints_dict,
                                       fuzzy_factor)
        except RdkitParseError:
            msg = 'Ligand {0} of {1} was not parseable by RDKit, skipping'.format(ligand.name(), pdb.code())
            click.secho(msg, bold=True)
        except AtomiumParseError:
            msg = 'Ligand {0} of {1} was not parseable by atomium, skipping'.format(ligand.name(), pdb.code())
            click.secho(msg, bold=True)


def generate_from_fragment(pdb,
                           ligand,
                           fragment,
                           frag_nr,
                           fragments_db,
                           pharmacophore_points,
                           fingerprints_dict,
                           fuzzy_factor):
    try:
        frag_id = build_frag_id(pdb, ligand, frag_nr)
        fragment.name = frag_id
        click.echo('Generating pharmacophore fingerprint for {0}'.format(frag_id))

        pharmacophore = from_fragment(fragment)
        fingerprint = pharmacophore.fingerprint(fuzzy_factor)

        click.echo('Saving fragment/pharmacophore/fingerprint')
        add_fragment2db(pdb, ligand, frag_nr, fragment, fragments_db)
        add_pharmacophore2db(pharmacophore_points, frag_id, pharmacophore)
        fingerprints_dict[frag_id] = fingerprint
    except NoFeatures:
        msg = 'Fragment {0} of ligand {1} of pdb {2} ' \
              'contains no pharmacophore features, skipping'.format(frag_nr, ligand.id(), pdb.code())
        click.secho(msg, bold=True)
    except IntegrityError:
        msg = 'Fragment {0} of ligand {1} of pdb {2} already present, skipping'.format(frag_nr,
                                                                                       ligand.name(),
                                                                                       pdb.code())
        click.secho(msg, bold=True)


def build_frag_id(thepdb, ligand, frag_nr):
    pdb_code = thepdb.code().lower()
    het_code = ligand.name()
    return pdb_code + '_' + het_code + '_frag' + str(frag_nr)


def is_ligand_stored(fragments_db, pdb_code, het_code):
    sql = 'SELECT 1 FROM fragments WHERE pdb_code=? AND het_code=?'
    fragments_db.cursor.execute(sql, (pdb_code.lower(), het_code))
    res = fragments_db.cursor.fetchone()
    return res is not None


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

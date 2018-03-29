=====
Kripo
=====

.. image:: https://travis-ci.org/3D-e-Chem/kripo.svg?branch=master
    :target: https://travis-ci.org/3D-e-Chem/kripo

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1209378.svg
   :target: https://doi.org/10.5281/zenodo.1209378

Key Representation of Interaction in POckets, see http://dx.doi.org/10.1186/1758-2946-6-S1-O26 for more information.

Command line tool to generate Kripo fingerprints from Protein Data Bank files.
The fingerprints can be used to can be used to find proteins which can bind a fragment of a ligand.

The kripo algorithm in a nutshell is:

1. Foreach PDB
2. Foreach ligand
3. Fragment ligand into fragments
4. Determine which residues of protein are in contact with fragment
5. Translate contact residues to pharmacophoric features
6. Convert pharmacophore to three point fingerprint

It uses https://github.com/3D-e-Chem/kripodb to store output in database files.

Usage
-----

The ligands in PDB files have no bond information so we use the ligand expo (http://ligand-expo.rcsb.org) for this.
The ligand expo sdf file must be downloaded and indexed using:

.. code-block:: bash

    wget http://ligand-expo.rcsb.org/dictionaries/all-sdf.sdf.gz
    kripo ligands import all-sdf.sdf.gz ligand-expo.db

To generate pharmacophore fingerprints from 2 pdb files use:

.. code-block:: bash

    echo tests/fixtures/3HEG.pdb > pdb.list
    echo tests/fixtures/5IS0.pdb >> pdb.list
    kripo generate --ligand-expo ligand-expo.db pdb.list frags.db phars.h5 fingerprints.db

This will generate Kripo fragments/pharmacophores/fingerprints for the PDB files listed in the `pdb.list` file.

To host the Kripo web service the following steps must be done:

.. code-block:: bash

    # Fetch Uniprot, Enzyme for each PDB entry
    kripodb fragments pdb frags.db
    # Calculate similarity scores between fingerprints pairs
    kripodb fingerprints similarities --fragmentsdbfn frags.db fingerprints.db fingerprints.db similarities.h5
    # Freeze similarity pairs into a matrix which can be quickly queried for most similar
    kripodb similarities freeze similarities.h5  similarities.frozen.h5
    # Startup the web service
    kripodb serve similarities.frozen.h5 frags.db phars.h5

Other commands are described at http://kripodb.readthedocs.io/en/latest/cli.html .

Install
-------

To install run

.. code-block:: bash

    conda env create -n kripo -f environment.yml
    conda activate kripo
    python setup.py install

Conda is used to install rdkit, reduce and openbabel.

You will have `kripo` and `kripodb` commands available in your PATH.

For development of kripo code use `pip install -r requirements_dev.txt` instead of `python setup.py install`.

Reference
---------

KRIPO – a structure-based pharmacophores approach explains polypharmacological effects;
Tina Ritschel, Tom JJ Schirris, and Frans GM Russel; J Cheminform. 2014; 6(Suppl 1): O26;
Published online 2014 Mar 11; http://dx.doi.org/10.1186/1758-2946-6-S1-O26

For protonation of protein Reduce http://kinemage.biochem.duke.edu/software/reduce.php is used:
Word, et. al. (1999) Asparagine and Glutamine: Using Hydrogen Atom
Contacts in the Choice of Side-chain Amide Orientation, J. Mol. Biol. 285, 1733-1747.

For protonation of ligand OpenBabel https://doi.org/10.1186/1758-2946-3-33 is used.

Legal
-----

Kripo is copyright 2018 Vrije Unversiteit Amsterdam and Netherlands eScience Center. It is free, open source software. You can redistribute it and/or modify it under the terms of either:

- The GNU General Public License as published by the Free Software Foundation, either version 2 or, (at your option), any later version, or
- The Apache License 2.0

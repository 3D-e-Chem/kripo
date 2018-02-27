=====
Kripo
=====


.. image:: https://img.shields.io/pypi/v/kripo.svg
        :target: https://pypi.python.org/pypi/kripo

.. image:: https://img.shields.io/travis/3D-e-Chem/kripo.svg
        :target: https://travis-ci.org/3D-e-Chem/kripo

.. image:: https://readthedocs.org/projects/kripo/badge/?version=latest
        :target: https://kripo.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/3D-e-Chem/kripo/shield.svg
     :target: https://pyup.io/repos/github/3D-e-Chem/kripo/
     :alt: Updates


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

To generate pharmacophore fingerprints from 2 pdb files use:

.. code-block:: bash

    echo tests/fixtures/3HEG.pdb > pdb.list
    echo tests/fixtures/5IS0.pdb >> pdb.list
    kripo generate pdb.list frags.db phars.h5 fingerprints.db

This will generate Kripo fragments/pharmacophores/fingerprints for the PDB files listed in the `pdb.list` file.

To host the Kripo web service the following steps must be done:

.. code-block:: bash

    kripodb fragments pdb frags.db
    kripodb fingerprints similarities --fragmentsdbfn frags.db fingerprints.db fingerprints.db similarities.h5
    kripodb serve similarities.h5 frags.db phars.h5

Install
-------

.. code-block:: bash

    conda env create -f environment.yml
    python setup.py install

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

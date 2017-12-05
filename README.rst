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


Key Representation of Interaction in POckets


* Free software: Apache Software License 2.0
* Documentation: https://kripo.readthedocs.io.


Usage
-----

To generate pharmacophore fingerprints from 2 pdb files use:

.. code-block:: bash

    echo tests/fixtures/3HEG.pdb > pdb.list
    echo tests/fixtures/5IS0.pdb >> pdb.list
    kripo generate pdb.list frags.db phars.h5 fingerprints.db


This will generate Kripo fragments/pharmacophores/fingerprints for the PDB files listed in the `pdb.list` file.

Install
-------

.. code-block:: bash

    conda env create -f environment.yml
    python setup.py install

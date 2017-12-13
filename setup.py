#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0',
    'atomium',
    'pybel',  # conda install -c openbabel openbabel
    'kripodb',
    'intbitset',
]

setup_requirements = [
    'pytest-runner',
]

test_requirements = [
    'pytest',
    'numpy',
]

setup(
    name='kripo',
    version='0.1.0',
    description="Key Representation of Interaction in POckets",
    long_description=readme + '\n\n' + history,
    author="Stefan Verhoeven",
    author_email='s.verhoeven@esciencecenter.nl',
    url='https://github.com/3D-e-Chem/kripo',
    packages=find_packages(include=['kripo']),
    entry_points={
        'console_scripts': [
            'kripo=kripo.cli:main'
        ]
    },
    package_data={'kripo': [
        'data/combined.smirks',
        'data/PHARMACKEY_3PFP_25BINS_6FEATURES_new.txt.bz2',
    ]},
    include_package_data=True,
    install_requires=requirements,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='kripo',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)

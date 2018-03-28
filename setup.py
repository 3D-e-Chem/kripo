#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

exec(open('kripo/version.py').read())

requirements = [
    'Click>=6.0',
    'atomium',
    'pybel',
    'kripodb',
    'intbitset',
    'requests',
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
    version=__version__,
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
    license="Apache Software License 2.0 or GPLv2+",
    zip_safe=False,
    keywords='kripo',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)

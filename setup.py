#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = [
    "reframed>=1.2",
    "pandas>=0.20.0",
    "requests>=2.18"
]

included_files = {
    'carveme': [
        'config.cfg',
        'data/input/bigg_models.csv',
        'data/input/biomass_db.tsv',
        'data/input/manually_curated.csv',
        'data/input/media_db.tsv',
        'data/input/metabolomics_park2016.csv',
        'data/input/unbalanced_metabolites.csv',
        'data/input/bigg_proteins.faa',
        'data/input/bigg_proteins.dmnd',
        'data/input/equilibrator_compounds.tsv.gz',
        'data/input/refseq_release_201.tsv.gz',
        'data/generated/bigg_gibbs.csv',
        'data/generated/bigg_gprs.csv.gz',
        'data/generated/model_specific_data.csv.gz',
        'data/generated/universe_draft.xml.gz',
        'data/generated/universe_bacteria.xml.gz',
        'data/generated/universe_grampos.xml.gz',
        'data/generated/universe_gramneg.xml.gz',
        'data/generated/universe_archaea.xml.gz',
    ]
}


setup(
    name='carveme',
    version='1.4.0',
    description="CarveMe: automated metabolic model reconstruction",
    long_description=readme,
    author="Daniel Machado",
    author_email='cdanielmachado@gmail.com',
    url='https://github.com/cdanielmachado/carveme',
    entry_points={
        'console_scripts': [
            'build_universe=carveme.cli.build_universe:main',
            'carve=carveme.cli.carve:main',
            'gapfill=carveme.cli.gapfill:main',
            'merge_community=carveme.cli.merge_community:main',
            'benchmark=carveme.cli.benchmark:main',
        ],
    },
    #   package_dir={'':'src'},
    packages=find_packages(),
    include_package_data=True,
    package_data=included_files,
    install_requires=requirements,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='carveme',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console', 
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: Apache Software License',
    ],
    setup_requires=['setuptools_scm']
    #    test_suite='tests',
    #    tests_require=test_requirements,
)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = [
    "framed==0.5.1",
    "pandas>=0.20.0",
    "requests>=2.18"
]

script_list = [
    "scripts/carveme_init",
    "scripts/build_universe",
    "scripts/carve",
    "scripts/gapfill",
    "scripts/merge_community"
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
        'data/generated/bigg_gibbs.csv'
    ]
}


setup(
    name='carveme',
    version='1.2.1',
    description="CarveMe: automated metabolic model reconstruction",
    long_description=readme,
    author="Daniel Machado, Sergej Andrejev",
    author_email='cdanielmachado@gmail.com',
    url='https://github.com/cdanielmachado/carveme',
    scripts=script_list,
#    package_dir={'':'src'},
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
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: Apache Software License',
    ],
    setup_requires=['setuptools_scm']
#    test_suite='tests',
#    tests_require=test_requirements,
)

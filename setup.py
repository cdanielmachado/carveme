#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os
from configparser import ConfigParser
from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = [
    "reframed>=1.5.4",
    "pandas>=1.0",
]

included_files = {
    'carveme': [
        'config.cfg',
        'data/input/bigg_models.tsv',
        'data/input/biomass_db.tsv',
        'data/input/manually_curated.tsv',
        'data/input/media_db.tsv',
#        'data/input/metabolomics_park2016.csv', deleted 5cbc611af5aa265c39882f7a88bf357f3261b170
        'data/input/unbalanced_metabolites.csv',
        'data/generated/bigg_proteins.faa',
        'data/input/mnx_compounds.tsv',
        'data/input/refseq_release_201.tsv.gz',
        'data/generated/gene_annotations.tsv.gz',
#        'data/generated/bigg_gibbs.csv', # deleted c897f41d7d03c27ca12ecd9ee97337355338c378
        'data/generated/bigg_gprs.csv.gz',
        'data/generated/model_specific_data.csv.gz',
        'data/generated/bigg_universe.xml.gz',
        'data/generated/universe_bacteria.xml.gz',
        'data/generated/universe_grampos.xml.gz',
        'data/generated/universe_gramneg.xml.gz',
        'data/generated/universe_archaea.xml.gz',
        'data/benchmark/media_db.tsv',
        'data/benchmark/biolog/bsub/biolog_carbon.tsv',
        'data/benchmark/biolog/bsub/biolog_phosphorus.tsv',
        'data/benchmark/biolog/bsub/biolog_nitrogen.tsv',
        'data/benchmark/biolog/bsub/biolog_sulfur.tsv',
        'data/benchmark/biolog/ecol/biolog_carbon.tsv',
        'data/benchmark/biolog/ecol/biolog_phosphorus.tsv',
        'data/benchmark/biolog/ecol/biolog_nitrogen.tsv',
        'data/benchmark/biolog/ecol/biolog_sulfur.tsv',
        'data/benchmark/biolog/paer/biolog_carbon.tsv',
        'data/benchmark/biolog/rsol/biolog_carbon.tsv',
        'data/benchmark/biolog/rsol/biolog_phosphorus.tsv',
        'data/benchmark/biolog/rsol/biolog_nitrogen.tsv',
        'data/benchmark/biolog/rsol/biolog_sulfur.tsv',
        'data/benchmark/biolog/sone/biolog_carbon.tsv',
        'data/benchmark/biolog/sone/biolog_nitrogen.tsv',
        'data/benchmark/essentiality/bsub.tsv',
        'data/benchmark/essentiality/ecol.tsv',
        'data/benchmark/essentiality/mgen.tsv',
        'data/benchmark/essentiality/paer.tsv',
        'data/benchmark/essentiality/sone.tsv',
        'data/benchmark/fasta/Bsubtilis_168.faa',
        'data/benchmark/fasta/Paeruginosa_PAO1.faa',
        'data/benchmark/fasta/Ecoli_K12_MG1655.faa',
        'data/benchmark/fasta/M_genitalium_G37.faa',
        'data/benchmark/fasta/Rsolanacearum_GMI1000.faa',
        'data/benchmark/fasta/Soneidensis_MR1.faa',
        'data/benchmark/models/bsub.xml',
        'data/benchmark/models/rsol.xml',
        'data/benchmark/models/ecol.xml',
        'data/benchmark/models/paer.xml',
        'data/benchmark/models/sone.xml',
        'data/benchmark/models/mgen.xml',
        'data/benchmark/results/biolog.tsv',
        'data/benchmark/results/essentiality.tsv',
    ]
}
missing_files = []
for path in included_files["carveme"]:
    fullpath = os.path.join("carveme", path)
    if not os.path.exists(fullpath):
        missing_files.append(fullpath)
if missing_files:
    print("files required for install are not found:\n")
    print("\n".join(missing_files))
    raise ValueError("missing files; exiting")

config = ConfigParser()
project_dir = "carveme"
config.read(os.path.join(project_dir, 'config.cfg'))
config_files = []
for chunk in ["input", "generated"]:
    for k,v in config[chunk].items():
        vpath = os.path.join(project_dir, v)
        if k in ["folder", "diamond_db"]: continue
        if not os.path.exists(vpath):
            raise ValueError(f'file {vpath} not found')
        elif v not in included_files["carveme"]:
            raise ValueError(f'config file {vpath} not included in setup.py')
        else:
            config_files.append(v)



setup(
    name='carveme',
    version='1.6.4',
    description="CarveMe: automated metabolic model reconstruction",
    long_description=readme,
    author="Daniel Machado",
    author_email='cdanielmachado@gmail.com',
    url='https://github.com/cdanielmachado/carveme',
    entry_points={
        'console_scripts': [
            'build_universe=carveme.cli.build_universe:main',
            'curate_universe=carveme.cli.curate_universe:main',
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
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'License :: OSI Approved :: Apache Software License',
    ],
    setup_requires=['setuptools_scm']
    #    test_suite='tests',
    #    tests_require=test_requirements,
)

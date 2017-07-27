#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = [
    "framed",
    "pandas"
]

script_list = [
    "scripts/carveme_init.py",
    "scripts/build_universe.py",
    "scripts/carveme.py",
    "scripts/gapfill.py"
]

#setup_requirements = []
#test_requirements = []

setup(
    name='carveme',
    version='0.1.0',
    description="CarveMe: automated metabolic model reconstruction",
    long_description=readme,
    author="Daniel Machado, Sergej Andrejev",
    author_email='cdanielmachado@gmail.com',
    url='https://github.com/cdanielmachado/carveme',
    scripts=script_list,
    packages=find_packages(include=['src']),
    include_package_data=True,
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
#    test_suite='tests',
#    tests_require=test_requirements,
#    setup_requires=setup_requirements,
)

#!/usr/bin/env python

from carveme import config, project_dir
import argparse
import os
import pandas as pd

from carveme.universe.curation import curate_universe
from carveme.reconstruction.utils import load_biomass_db
from reframed import load_cbmodel


def curate(inputfile=None, outputfile=None, taxa=None, biomass=None, biomass_db_path=None, normalize_biomass=False):

    if inputfile:
        universe_draft = inputfile
        model_specific_data = os.path.splitext(inputfile)[0] + '.csv'
    else:
        universe_draft = project_dir + config.get('generated', 'bigg_universe')
        model_specific_data = project_dir + config.get('generated', 'model_specific_data')

    if not biomass:
        biomass = 'gramneg' if taxa == 'cyanobacteria' else taxa

    if not outputfile:
        folder = project_dir + config.get('generated', 'folder')
        outputfile = folder + f"universe_{taxa}.xml.gz"

    bigg_models = project_dir + config.get('input', 'bigg_models')
    bigg_models = pd.read_csv(bigg_models, sep='\t')

    manual_curation = project_dir + config.get('input', 'manually_curated')
    manually_curated = pd.read_csv(manual_curation, index_col=0, sep='\t')

    unbalanced = project_dir + config.get('input', 'unbalanced_metabolites')
    unbalanced = pd.read_csv(unbalanced, header=None)
    unbalanced = unbalanced[0].tolist()

#    metabolomics = project_dir + config.get('input', 'metabolomics')
#    metabolomics_data = pd.read_csv(metabolomics, index_col=1)

    try:
        model = load_cbmodel(universe_draft, flavor=config.get('sbml', 'default_flavor'))
        model_specific_data = pd.read_csv(model_specific_data)
    except IOError:
        raise IOError('Universe draft model not found. Please run --build first to download BiGG data.')

    if biomass_db_path is None:
        biomass_db_path = project_dir + config.get('input', 'biomass_library')

    biomass_db = load_biomass_db(biomass_db_path, normalize_weight=normalize_biomass, model=model)

    if biomass not in biomass_db:
        valid_ids = ','.join(biomass_db.keys())
        raise RuntimeError('Biomass identifier not in database. Currently in database: ' + valid_ids)

    biomass_eq = biomass_db[biomass]

    curate_universe(model,
                    outputfile=outputfile,
                    taxa=taxa,
                    biomass_eq=biomass_eq,
                    model_specific_data=model_specific_data,
                    bigg_models=bigg_models,
#                    metabolomics_data=metabolomics_data,
                    manually_curated=manually_curated,
                    unbalanced_metabolites=unbalanced,
                    )


def main():
    parser = argparse.ArgumentParser(description="Generate universal model to use with CarveMe")

    parser.add_argument('-i', '--input', dest='input', help="Input file")
    parser.add_argument('-o', '--output', dest='output', help="Output file")

    taxa = parser.add_mutually_exclusive_group(required=False)

    taxa.add_argument('--cyanobacteria', action='store_true',
                      help='Generate a template for cyanobacteria (includes thylakoid compartment).')

    taxa.add_argument('--archaea', action='store_true',
                      help='Generate a template for archaea (includes methanogenic reactions).')

    taxa.add_argument('--grampos', action='store_true',
                      help='Generate a template for gram-positive bacteria.')

    taxa.add_argument('--gramneg', action='store_true',
                      help='Generate a template for gram-negative bacteria.')

    parser.add_argument('--biomass', help="Advanced options: biomass equation identifier")

    parser.add_argument('--biomass-db', dest='biomass_db', help="Advanced options: biomass database file")

    parser.add_argument('--normalize-biomass', dest='normalize_biomass', action='store_true',
                        help="Advanced options: Normalize biomass dry weight to 1 gram")

    args = parser.parse_args()

    if args.cyanobacteria:
        taxa = 'cyanobacteria'
    elif args.archaea:
        taxa = 'archaea'
    elif args.grampos:
        taxa = 'grampos'
    elif args.gramneg:
        taxa = 'gramneg'
    else:
        taxa = 'bacteria'

    curate(inputfile=args.input, outputfile=args.output, taxa=taxa, biomass=args.biomass,
               biomass_db_path=args.biomass_db, normalize_biomass=args.normalize_biomass)


if __name__ == '__main__':
    main()


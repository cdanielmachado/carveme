#!/usr/bin/env python

from carveme import config, project_dir
import argparse
import os
import pandas as pd

from carveme.universe.download import download_universal_model, download_model_specific_data
from carveme.universe.curation import curate_universe
from carveme.universe.annotate import annotate_with_eQuilibrator
from carveme.reconstruction.utils import load_biomass_db
from reframed import load_cbmodel, save_cbmodel


def maincall(mode, inputfile=None, outputfile=None, biomass=None, biomass_db_path=None,
             normalize_biomass=False, taxa=None):

    if mode == 'build':

        if outputfile:
            universe_draft = outputfile
            model_specific_data = os.path.splitext(outputfile)[0] + '.csv'
            bigg_gprs = os.path.splitext(outputfile)[0] + '_gprs.csv'
            fasta_file = os.path.splitext(outputfile)[0] + '.faa'
        else:
            universe_draft = project_dir + config.get('generated', 'bigg_universe')
            model_specific_data = project_dir + config.get('generated', 'model_specific_data')
            bigg_gprs = project_dir + config.get('generated', 'bigg_gprs')
            fasta_file = project_dir + config.get('generated', 'fasta_file')

        download_universal_model(universe_draft)
        download_model_specific_data(model_specific_data, bigg_gprs, fasta_file)

    elif mode == 'annotate':

        if not inputfile:
            inputfile = project_dir + config.get('generated', 'bigg_universe')

        if not outputfile:
            outputfile = project_dir + config.get('generated', 'bigg_annotated')

        model = load_cbmodel(inputfile, flavor="bigg")

        annotate_with_eQuilibrator(model)

        save_cbmodel(model, outputfile)

    elif mode == 'curate':

        universe_draft = project_dir + config.get('generated', 'universe_draft')
        model_specific_data = project_dir + config.get('generated', 'model_specific_data')

        if not biomass:
            if taxa == 'archaea':
                biomass = 'archaea'
            else:
                biomass = config.get('universe', 'default_biomass')

        if outputfile:
            universe_final = outputfile
        else:
            tag = taxa if taxa != 'bacteria' else biomass
            universe_final = "{}{}universe_{}.xml.gz".format(project_dir, config.get('generated', 'folder'), tag)

        bigg_models = project_dir + config.get('input', 'bigg_models')
        bigg_models = pd.read_csv(bigg_models, index_col=0)

        manual_curation = project_dir + config.get('input', 'manually_curated')
        manually_curated = pd.read_csv(manual_curation, index_col=0)

        unbalanced = project_dir + config.get('input', 'unbalanced_metabolites')
        unbalanced = pd.read_csv(unbalanced, header=None)
        unbalanced = unbalanced[0].tolist()

        try:
            model = load_cbmodel(universe_draft, flavor=config.get('sbml', 'default_flavor'))
            model_specific_data = pd.read_csv(model_specific_data)
        except IOError:
            raise IOError('Universe draft not found. Please run --draft first to download BiGG data.')

        if biomass_db_path is None:
            biomass_db_path = project_dir + config.get('input', 'biomass_library')

        biomass_db = load_biomass_db(biomass_db_path, normalize_weight=normalize_biomass, model=model)

        if biomass not in biomass_db:
            valid_ids = ','.join(biomass_db.keys())
            raise RuntimeError('Biomass identifier not in database. Currently in database: ' + valid_ids)

        biomass_eq = biomass_db[biomass]
        metabolomics = project_dir + config.get('input', 'metabolomics')
        metabolomics_data = pd.read_csv(metabolomics, index_col=1)

        curate_universe(model,
                        taxa=taxa,
                        outputfile=universe_final,
                        model_specific_data=model_specific_data,
                        bigg_models=bigg_models,
                        thermodynamics_data=thermodynamics_data,
                        metabolomics_data=metabolomics_data,
                        manually_curated=manually_curated,
                        unbalanced_metabolites=unbalanced,
                        biomass_eq=biomass_eq)

    else:
        print('Unrecognized option:', mode)


def main():
    parser = argparse.ArgumentParser(description="Generate universal model to use with CarveMe")

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--build', action='store_true',
                      help='Download all data from BiGG and save SBML, FASTA, and additional files.')

    mode.add_argument('--annotate', action='store_true',
                      help='Annotate universal model with formula, charge, and Gibbs free energy.')

    mode.add_argument('--curate', action='store_true',
                      help='Run semi-automated curation of universal model.')

    parser.add_argument('-i', '--input', dest='input', help="Input file")
    parser.add_argument('-o', '--output', dest='output', help="Output file")

    taxa = parser.add_mutually_exclusive_group(required=False)

    taxa.add_argument('--cyanobacteria', action='store_true',
                      help='Generate a template for cyanobacteria (includes thylakoid compartment).')

    taxa.add_argument('--archaea', action='store_true',
                      help='Generate a template for archaea (includes methanogenic reactions).')

    parser.add_argument('--biomass', help="Advanced options: biomass equation identifier")

    parser.add_argument('--biomass-db', dest='biomass_db', help="Advanced options: biomass database file")

    parser.add_argument('--normalize-biomass', dest='normalize_biomass', action='store_true',
                        help="Advanced options: Normalize biomass dry weight to 1 gram")

    args = parser.parse_args()

    if args.build:
        mode = 'build'
    elif args.annotate:
        mode = 'annotate'
    elif args.curate:
        mode = 'curate'

    if args.cyanobacteria:
        taxa = 'cyanobacteria'
    elif args.archaea:
        taxa = 'archaea'
    else:
        taxa = 'bacteria'

    maincall(mode=mode, inputfile=args.input, outputfile=args.output, biomass=args.biomass,
             biomass_db_path=args.biomass_db, taxa=taxa)


if __name__ == '__main__':
    main()


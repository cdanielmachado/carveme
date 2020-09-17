#!/usr/bin/env python

from carveme import config, project_dir
import argparse
import os
import pandas as pd

from carveme.universe.bigg_download import build_bigg_universe_model, download_model_specific_data, create_gpr_table,\
    download_gene_sequences
from carveme.universe.curation import curate_universe
from carveme.universe.thermodynamics import compute_bigg_gibbs_energy
from carveme.reconstruction.utils import load_biomass_db
from reframed import load_cbmodel


def maincall(mode, noheuristics=False, nothermo=False, allow_unbalanced=False, allow_blocked=False,
         biomass=None, biomass_db_path=None, normalize_biomass=False, taxa=None, outputfile=None):

    if mode == 'draft':

        if outputfile:
            universe_draft = outputfile
            model_specific_data = os.path.splitext(outputfile)[0] + '.csv'
            bigg_gprs = os.path.splitext(outputfile)[0] + '_gprs.csv'
#            fasta_file = os.path.splitext(outputfile)[0] + '.faa'
        else:
            universe_draft = project_dir + config.get('generated', 'universe_draft')
            model_specific_data = project_dir + config.get('generated', 'model_specific_data')
            bigg_gprs = project_dir + config.get('generated', 'bigg_gprs')
#            fasta_file = project_dir + config.get('input', 'fasta_file')

        build_bigg_universe_model(universe_draft)
        data = download_model_specific_data(model_specific_data)
        gprs = create_gpr_table(data, outputfile=bigg_gprs)
#        download_gene_sequences(gprs, outputfile=fasta_file)

    elif mode == 'thermo':
        universe_draft = project_dir + config.get('generated', 'universe_draft')
        equilibrator_compounds = project_dir + config.get('input', 'equilibrator_compounds')

        if outputfile:
            bigg_gibbs = outputfile
        else:
            bigg_gibbs = project_dir + config.get('generated', 'bigg_gibbs')

        compute_bigg_gibbs_energy(universe_draft, equilibrator_compounds, bigg_gibbs)

    elif mode == 'curated':

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

        if nothermo:
            thermodynamics_data = None
            metabolomics_data = None
        else:
            try:
                bigg_gibbs = project_dir + config.get('generated', 'bigg_gibbs')
                thermodynamics_data = pd.read_csv(bigg_gibbs, index_col=0)
            except IOError:
                raise IOError('Thermodynamic data not found. Please run --thermo first to generate thermodynamic data.')

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
                        biomass_eq=biomass_eq,
                        use_heuristics=(not noheuristics),
                        remove_unbalanced=(not allow_unbalanced),
                        remove_blocked=(not allow_blocked))

    else:
        print('Unrecognized option:', mode)


def main():
    parser = argparse.ArgumentParser(description="Generate universal model to use with CarveMe")

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--draft', action='store_true',
                      help='Download all data from BiGG and save uncurated (draft) model to SBML.')
    mode.add_argument('--thermo', action='store_true',
                      help='Compute thermodynamics data for BiGG reactions.')
    mode.add_argument('--curated', action='store_true',
                      help='Generate curated universal model of bacterial metabolism.')

    parser.add_argument('-o', '--output', dest='output', help="Output file")

    parser.add_argument('--nothermo', action='store_true',
                        help="Advanced options: do not use thermodynamics data")
    parser.add_argument('--noheuristics', action='store_true',
                        help="Advanced options: do not apply heuristic reversibility rules")
    parser.add_argument('--unbalanced', action='store_true',
                        help="Advanced options: allow unbalanced reactions")
    parser.add_argument('--blocked', action='store_true',
                        help="Advanced options: allow blocked reactions")

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

    if args.nothermo and (args.draft or args.thermo):
        parser.error('--nothermo cannot be used with --draft or --thermo')

    if args.noheuristics and (args.draft or args.thermo):
        parser.error('--noheuristics cannot be used with --draft or --thermo')

    if args.unbalanced and (args.draft or args.thermo):
        parser.error('--unbalanced cannot be used with --draft or --thermo')

    if args.blocked and (args.draft or args.thermo):
        parser.error('--blocked cannot be used with --draft or --thermo')

    if args.biomass and (args.draft or args.thermo):
        parser.error('--biomass cannot be used with --draft or --thermo')

    if args.cyanobacteria and (args.draft or args.thermo):
        parser.error('--cyanobacteria cannot be used with --draft or --thermo')

    if args.archaea and (args.draft or args.thermo):
        parser.error('--archaea cannot be used with --draft or --thermo')

    if args.draft:
        mode = 'draft'
    elif args.thermo:
        mode = 'thermo'
    elif args.curated:
        mode = 'curated'

    if args.cyanobacteria:
        taxa = 'cyanobacteria'
    elif args.archaea:
        taxa = 'archaea'
    else:
        taxa = 'bacteria'

    maincall(mode=mode,
         nothermo=args.nothermo,
         noheuristics=args.noheuristics,
         allow_unbalanced=args.unbalanced,
         allow_blocked=args.blocked,
         biomass=args.biomass,
         biomass_db_path=args.biomass_db,
         normalize_biomass=args.normalize_biomass,
         taxa=taxa,
         outputfile=args.output)


if __name__ == '__main__':
    main()


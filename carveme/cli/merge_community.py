import argparse
import os
from reframed import load_cbmodel, save_cbmodel, Environment, Community
from carveme import config, project_dir
from carveme.reconstruction.utils import load_media_db


def maincall(inputfiles, flavor=None, split_pool=False, no_biomass=False, init=None, mediadb=None, ext_comp_id=None, outputfile=None):

    if not flavor:
        flavor = config.get('sbml', 'default_flavor')

    if outputfile:
        model_id = os.path.splitext(os.path.basename(outputfile))[0]
    else:
        model_id = 'community'
        outputfile = 'community.xml'

    if ext_comp_id is None:
        ext_comp_id = 'C_e'

    models = [load_cbmodel(inputfile, flavor=flavor) for inputfile in inputfiles]

    community = Community(model_id, models,
                          extracellular_compartment_id=ext_comp_id,
                          merge_extracellular_compartments=(not split_pool),
                          create_biomass=(not no_biomass))

    merged = community.generate_merged_model()

    if init:
        if not mediadb:
            mediadb = project_dir + config.get('input', 'media_library')

        try:
            media_db = load_media_db(mediadb)
        except IOError:
            raise IOError('Failed to load media library:' + mediadb)

        if split_pool:
            fmt_func = lambda x: f"R_EX_M_{x}_e_pool"
        else:
            fmt_func = lambda x: f"R_EX_{x}_e"
        init_env = Environment.from_compounds(media_db[init], fmt_func=fmt_func)
        init_env.apply(merged, inplace=True)

    save_cbmodel(merged, outputfile, flavor=flavor)


def main():
    parser = argparse.ArgumentParser(description="Merge single species models into a microbial community model")

    parser.add_argument('input', metavar='INPUTFILES', nargs='+', help="SBML input files (single species)")

    parser.add_argument('-o', '--output', dest='output', help="SBML output file (community)")

    parser.add_argument('--split-pool', action='store_true', dest='split_pool',
                        help='Keep individual extracellular compartments separated from common metabolite pool.')

    parser.add_argument('--no-community-biomass', action='store_true', dest='no_biomass',
                        help='Do not create a common community biomass equation.')

    parser.add_argument('-i', '--init', dest='init',
                        help="Initialize model with given medium")

    parser.add_argument('--mediadb', help="Media database file")
    parser.add_argument('--ext', help="Extracellular compartment identifier in the models (default 'C_e').")

    sbml = parser.add_mutually_exclusive_group()
    sbml.add_argument('--cobra', action='store_true', help="SBML input/output in old cobra format")
    sbml.add_argument('--fbc2', action='store_true', help="SBML input/output in sbml-fbc2 format")

    args = parser.parse_args()

    if len(args.input) < 2:
        print(args.input)
        parser.error("Please provide two or more single species models as input files.")

    if args.fbc2:
        flavor = 'fbc2'
    elif args.cobra:
        flavor = 'cobra'
    else:
        flavor = None

    maincall(inputfiles=args.input,
         flavor=flavor,
         split_pool=args.split_pool,
         no_biomass=args.no_biomass,
         init=args.init,
         mediadb=args.mediadb,
         ext_comp_id=args.ext,
         outputfile=args.output)


if __name__ == '__main__':
    main()
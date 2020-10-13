import argparse
import os
from reframed import load_cbmodel, save_cbmodel, Environment, Community
from carveme import config, project_dir
from carveme.reconstruction.utils import load_media_db


def maincall(inputfiles, flavor=None, init=None, mediadb=None,  outputfile=None):

    if not flavor:
        flavor = config.get('sbml', 'default_flavor')

    if outputfile:
        model_id = os.path.splitext(os.path.basename(outputfile))[0]
    else:
        model_id = 'community'
        outputfile = 'community.xml'

    models = [load_cbmodel(inputfile, flavor=flavor) for inputfile in inputfiles]
    community = Community(model_id, models)
    model = community.merged_model

    if init:
        if not mediadb:
            mediadb = project_dir + config.get('input', 'media_library')

        try:
            media_db = load_media_db(mediadb)
        except IOError:
            raise IOError('Failed to load media library:' + mediadb)

        init_env = Environment.from_compounds(media_db[init])
        init_env.apply(model, inplace=True)

    save_cbmodel(model, outputfile, flavor=flavor)


def main():
    parser = argparse.ArgumentParser(description="Merge single species models into a microbial community model")

    parser.add_argument('input', metavar='INPUTFILES', nargs='+', help="SBML input files (single species)")

    parser.add_argument('-o', '--output', dest='output', help="SBML output file (community)")

    parser.add_argument('-i', '--init', dest='init',
                        help="Initialize model with given medium")

    parser.add_argument('--mediadb', help="Media database file")

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
        flavor = config.get('sbml', 'default_flavor')

    maincall(inputfiles=args.input,
         flavor=flavor,
         init=args.init,
         mediadb=args.mediadb,
         outputfile=args.output)


if __name__ == '__main__':
    main()
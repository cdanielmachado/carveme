from carveme import config, project_dir
from carveme import __version__ as version
from carveme.reconstruction.carving import carve_model, build_ensemble
from carveme.reconstruction.eggnog import load_eggnog_data
from carveme.reconstruction.gapfilling import multiGapFill
from carveme.reconstruction.utils import load_media_db, load_soft_constraints, load_hard_constraints
from carveme.reconstruction.ncbi_download import load_ncbi_table, download_ncbi_genome
from carveme.reconstruction.scoring import reaction_scoring
from carveme.reconstruction.diamond import run_blast, load_diamond_results
from reframed import load_cbmodel, save_cbmodel, Environment
from reframed.io.sbml import sanitize_id
import argparse
import os
import os.path
import pandas as pd
from multiprocessing import Pool
from glob import glob
import subprocess


def first_run_check():
    diamond_db = project_dir + config.get('generated', 'diamond_db')
    if not os.path.exists(diamond_db):
        print("Running diamond for the first time, please wait while we build the internal database...")
        fasta_file = project_dir + config.get('input', 'fasta_file')
        cmd = ['diamond', 'makedb', '--in', fasta_file, '-d', diamond_db.split('.')[0]]
        try:
            exit_code = subprocess.call(cmd)
        except OSError:
            print('Unable to run diamond (make sure diamond is available in your PATH).')
        else:
            if exit_code != 0:
                print('Failed to run diamond (wrong arguments).')


def build_model_id(name):
    model_id = sanitize_id(name)
    if not model_id[0].isalpha():
        model_id = 'm_' + model_id
    return model_id


def maincall(inputfile, input_type='protein', outputfile=None, diamond_args=None, universe=None, universe_file=None,
         ensemble_size=None, verbose=False, debug=False, flavor=None, gapfill=None, blind_gapfill=False, init=None,
         mediadb=None, default_score=None, uptake_score=None, soft_score=None, soft=None, hard=None, reference=None,
         ref_score=None, recursive_mode=False):

    if recursive_mode:
        model_id = os.path.splitext(os.path.basename(inputfile))[0]

        if outputfile:
            outputfile = '{}/{}.xml'.format(outputfile, model_id)
        else:
            outputfile = os.path.splitext(inputfile)[0] + '.xml'

    else:
        if outputfile:
            model_id = os.path.splitext(os.path.basename(outputfile))[0]
        else:
            model_id = os.path.splitext(os.path.basename(inputfile))[0]
            outputfile = os.path.splitext(inputfile)[0] + '.xml'

    model_id = build_model_id(model_id)

    outputfolder = os.path.abspath(os.path.dirname(outputfile))

    if not os.path.exists(outputfolder):
        try:
            os.makedirs(outputfolder)
        except:
            print('Unable to create output folder:', outputfolder)
            return

    if soft:
        try:
            soft_constraints = load_soft_constraints(soft)
        except IOError:
            raise IOError('Failed to load soft-constraints file:' + soft)
    else:
        soft_constraints = None

    if hard:
        try:
            hard_constraints = load_hard_constraints(hard)
        except IOError:
            raise IOError('Failed to load hard-constraints file:' + hard)
    else:
        hard_constraints = None

    if input_type == 'refseq':

        if verbose:
            print('Downloading genome {} from NCBI...'.format(inputfile))

        ncbi_table = load_ncbi_table(project_dir + config.get('input', 'refseq'))
        inputfile = download_ncbi_genome(inputfile, ncbi_table)

        if not inputfile:
            print('Failed to download genome from NCBI.')
            return

        input_type = 'protein' if inputfile.endswith('.faa.gz') else 'dna'

    if input_type == 'protein' or input_type == 'dna':
        if verbose:
            print('Running diamond...')
        diamond_db = project_dir + config.get('generated', 'diamond_db')
        blast_output = os.path.splitext(inputfile)[0] + '.tsv'
        exit_code = run_blast(inputfile, input_type, blast_output, diamond_db, diamond_args, verbose)

        if exit_code is None:
            print('Unable to run diamond (make sure diamond is available in your PATH).')
            return

        if exit_code != 0:
            print('Failed to run diamond.')
            if diamond_args is not None:
                print('Incorrect diamond args? Please check documentation or use default args.')
            return

        annotations = load_diamond_results(blast_output)
    elif input_type == 'eggnog':
        annotations = load_eggnog_data(inputfile)
    elif input_type == 'diamond':
        annotations = load_diamond_results(inputfile)
    else:
        raise ValueError('Invalid input type: ' + input_type)

    if verbose:
        print('Loading universe model...')

    if not universe_file:
        if universe:
            universe_file = "{}{}universe_{}.xml.gz".format(project_dir, config.get('generated', 'folder'), universe)
        else:
            universe_file = project_dir + config.get('generated', 'default_universe')

    try:
        universe_model = load_cbmodel(universe_file, flavor='cobra')
        universe_model.id = model_id
    except IOError:
        available = '\n'.join(glob("{}{}universe_*.xml.gz".format(project_dir, config.get('generated', 'folder'))))
        raise IOError('Failed to load universe model: {}\nAvailable universe files:\n{}'.format(universe_file, available))

    if reference:
        if verbose:
            print('Loading reference model...')

        try:
            ref_model = load_cbmodel(reference)
        except:
            raise IOError('Failed to load reference model.')
    else:
        ref_model = None

    if gapfill or init:

        if verbose:
            print('Loading media library...')

        if not mediadb:
            mediadb = project_dir + config.get('input', 'media_library')

        try:
            media_db = load_media_db(mediadb)
        except IOError:
            raise IOError('Failed to load media library:' + mediadb)

    if verbose:
        print('Scoring reactions...')

    bigg_gprs = project_dir + config.get('generated', 'bigg_gprs')
    gprs = pd.read_csv(bigg_gprs)
    gprs = gprs[gprs.reaction.isin(universe_model.reactions)]

    debug_output = model_id if debug else None
    scores = reaction_scoring(annotations, gprs, debug_output=debug_output)

    if scores is None:
        print('The input genome did not match sufficient genes/reactions in the database.')
        return

    if not flavor:
        flavor = config.get('sbml', 'default_flavor')

    init_env = None

    if init:
        if init in media_db:
            init_env = Environment.from_compounds(media_db[init])
        else:
            print('Error: medium {} not in media database.'.format(init))

    universe_model.metadata['Description'] = 'This model was built with CarveMe version ' + version

    if ensemble_size is None or ensemble_size <= 1:
        if verbose:
            print('Reconstructing a single model')

        if not gapfill:
            carve_model(universe_model, scores,
                        outputfile=outputfile,
                        flavor=flavor,
                        default_score=default_score,
                        uptake_score=uptake_score,
                        soft_score=soft_score,
                        soft_constraints=soft_constraints,
                        hard_constraints=hard_constraints,
                        ref_model=ref_model,
                        ref_score=ref_score,
                        init_env=init_env,
                        debug_output=debug_output)
        else:
            model = carve_model(universe_model, scores,
                                inplace=False,
                                default_score=default_score,
                                uptake_score=uptake_score,
                                soft_score=soft_score,
                                soft_constraints=soft_constraints,
                                hard_constraints=hard_constraints,
                                ref_model=ref_model,
                                ref_score=ref_score,
                                init_env=init_env,
                                debug_output=debug_output)
    else:
        if verbose:
            print('Building an ensemble of', ensemble_size, 'models')
        build_ensemble(universe_model, scores, ensemble_size, outputfile, flavor, init_env=init_env)

    if gapfill and model is not None:

        media = gapfill.split(',')

        if verbose:
            m1, n1 = len(model.metabolites), len(model.reactions)
            print('Gap filling for {}...'.format(', '.join(media)))

        max_uptake = config.getint('gapfill', 'max_uptake')

        if blind_gapfill:
            scores = None
        else:
            scores = dict(scores[['reaction', 'normalized_score']].values)
        multiGapFill(model, universe_model, media, media_db, scores=scores, max_uptake=max_uptake, inplace=True)

        if verbose:
            m2, n2 = len(model.metabolites), len(model.reactions)
            print('Added {} reactions and {} metabolites'.format((n2 - n1), (m2 - m1)))

        if init_env:  # Initializes environment again as new exchange reactions can be acquired during gap-filling
            init_env.apply(model, inplace=True, warning=False)

        save_cbmodel(model, outputfile, flavor=flavor)

    if verbose:
        print('Done.')


def main():

    parser = argparse.ArgumentParser(description="Reconstruct a metabolic model using CarveMe",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('input', metavar='INPUT', nargs='+',
                        help="Input (protein fasta file by default, see other options for details).\n" +
                             "When used with -r an input pattern with wildcards can also be used.\n" +
                             "When used with --refseq an NCBI RefSeq assembly accession is expected."
                        )

    input_type_args = parser.add_mutually_exclusive_group()
    input_type_args.add_argument('--dna', action='store_true', help="Build from DNA fasta file")
    input_type_args.add_argument('--egg', action='store_true', help="Build from eggNOG-mapper output file")
    input_type_args.add_argument('--diamond', action='store_true', help=argparse.SUPPRESS)
    input_type_args.add_argument('--refseq', action='store_true', help="Download genome from NCBI RefSeq and build")

    parser.add_argument('--diamond-args', help="Additional arguments for running diamond")

    parser.add_argument('-r', '--recursive', action='store_true', dest='recursive',
                        help="Bulk reconstruction from folder with genome files")

    parser.add_argument('-o', '--output', dest='output', help="SBML output file (or output folder if -r is used)")

    univ = parser.add_mutually_exclusive_group()
    univ.add_argument('-u', '--universe', dest='universe', help="Pre-built universe model (default: bacteria)")
    univ.add_argument('--universe-file', dest='universe_file', help="Reaction universe file (SBML format)")

    sbml = parser.add_mutually_exclusive_group()
    sbml.add_argument('--cobra', action='store_true', help="Output SBML in old cobra format")
    sbml.add_argument('--fbc2', action='store_true', help="Output SBML in sbml-fbc2 format")

    parser.add_argument('-n', '--ensemble', type=int, dest='ensemble',
                        help="Build model ensemble with N models")

    parser.add_argument('-g', '--gapfill', dest='gapfill',
                        help="Gap fill model for given media")

    parser.add_argument('-i', '--init', dest='init',
                        help="Initialize model with given medium")

    parser.add_argument('--mediadb', help="Media database file")

    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', help="Switch to verbose mode")
    parser.add_argument('-d', '--debug', action='store_true', dest='debug',
                        help="Debug mode: writes intermediate results into output files")

    parser.add_argument('--soft', help="Soft constraints file")
    parser.add_argument('--hard', help="Hard constraints file")

    parser.add_argument('--reference', help="Manually curated model of a close reference species.")

    parser.add_argument('--default-score', type=float, default=-1.0, help=argparse.SUPPRESS)
    parser.add_argument('--uptake-score', type=float, default=0.0, help=argparse.SUPPRESS)
    parser.add_argument('--soft-score', type=float, default=1.0, help=argparse.SUPPRESS)
    parser.add_argument('--reference-score', type=float, default=0.0, help=argparse.SUPPRESS)

    parser.add_argument('--blind-gapfill', action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.gapfill and args.ensemble:
        parser.error('Gap fill and ensemble generation cannot currently be combined (not implemented yet).')

    if (args.soft or args.hard) and args.ensemble:
        parser.error('Soft/hard constraints and ensemble generation cannot currently be combined (not implemented yet).')

    if args.mediadb and not args.gapfill:
        parser.error('--mediadb can only be used with --gapfill')

    if args.recursive and args.refseq:
        parser.error('-r cannot be combined with --refseq')

    if args.egg:
        input_type = 'eggnog'
    elif args.dna:
        input_type = 'dna'
    elif args.diamond:
        input_type = 'diamond'
    elif args.refseq:
        input_type = 'refseq'
    else:
        input_type = 'protein'

    if args.fbc2:
        flavor = 'fbc2'
    elif args.cobra:
        flavor = 'cobra'
    else:
        flavor = config.get('sbml', 'default_flavor')

    first_run_check()

    if not args.recursive:
        if len(args.input) > 1:
            parser.error('Use -r when specifying more than one input file')

        maincall(
            inputfile=args.input[0],
            input_type=input_type,
            outputfile=args.output,
            diamond_args=args.diamond_args,
            universe=args.universe,
            universe_file=args.universe_file,
            ensemble_size=args.ensemble,
            verbose=args.verbose,
            debug=args.debug,
            flavor=flavor,
            gapfill=args.gapfill,
            blind_gapfill=False,
            init=args.init,
            mediadb=args.mediadb,
            default_score=args.default_score,
            uptake_score=args.uptake_score,
            soft_score=args.soft_score,
            soft=args.soft,
            hard=args.hard,
            reference=args.reference,
            ref_score=args.reference_score
        )

    else:

        def f(x):
            maincall(
                inputfile=x,
                input_type=input_type,
                outputfile=args.output,
                diamond_args=args.diamond_args,
                universe=args.universe,
                universe_file=args.universe_file,
                ensemble_size=args.ensemble,
                verbose=args.verbose,
                flavor=flavor,
                gapfill=args.gapfill,
                blind_gapfill=False,
                init=args.init,
                mediadb=args.mediadb,
                default_score=args.default_score,
                uptake_score=args.uptake_score,
                soft_score=args.soft_score,
                soft=args.soft,
                hard=args.hard,
                reference=args.reference,
                ref_score=args.reference_score,
                recursive_mode=True
            )

        p = Pool()
        p.map(f, args.input)


if __name__ == '__main__':
    main()

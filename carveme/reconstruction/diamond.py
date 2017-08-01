import pandas as pd
from subprocess import call
import os


def load_diamond_results(filename, drop_unused_cols=True):
    """ Load and parse an diamond results file.

    Args:
        filename (str): input file
        drop_unused_cols (bool): remove columns not used for model carving (default: True)

    Returns:
        pandas.DataFrame: eggnog data

    """
    columns = ['query_gene', 'BiGG_gene', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'score']

    data = pd.read_csv(filename, sep='\t', names=columns)

    if drop_unused_cols:
        data = data[['query_gene', 'BiGG_gene', 'score']]

    return data


def run_blast(inputfile, input_type, outputfile, database, args=None, verbose=True):
    """ Run blast aligment with Diamond.

    Args:
        inputfile (str): fasta input file
        input_type (str): sequence type ('protein' or 'dna')
        outputfile (str): output file name
        database (str): path to diamond protein database file
        args (str): additional arguments to be passed to diamond (optional)
        verbose (bool): allow diamond output to stdout (default: True)

    Returns:
        int: diamond exit code

    Notes:
        Default arguments are: --top 10 --more-sensitive

    """

    assert (input_type in ['protein', 'dna']), "Input type must be either 'protein' or 'dna'"

    cmd = ['diamond']

    if input_type == 'protein':
        cmd += ['blastp']
    elif input_type == 'dna':
        cmd += ['blastx']

    cmd += ['-d', database]
    cmd += ['-q', inputfile]
    cmd += ['-o', outputfile]

    if not args:
        args = "--more-sensitive --top 10"

    cmd += args.split()

    if verbose:
        print ' '.join(cmd)

    with open(os.devnull, 'w') as devnull:
        try:
            exit_code = call(cmd, stdout=devnull)
        except OSError:
            exit_code = None

    return exit_code

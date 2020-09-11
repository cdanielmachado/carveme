import urllib.request, urllib.parse, urllib.error
import subprocess
import os
from carveme import project_dir


def main():
    source_url = 'https://github.com/cdanielmachado/carveme/raw/master/carveme/'

    downloads = [
        'data/input/bigg_proteins.faa',
        'data/input/equilibrator_compounds.tsv.gz',
        'data/input/refseq_release_92.tsv.gz',
        'data/input/genbank_release_230.tsv.gz',
        'data/generated/bigg_gprs.csv.gz',
        'data/generated/model_specific_data.csv.gz',
        'data/generated/universe_draft.xml.gz',
        'data/generated/universe_bacteria.xml.gz',
        'data/generated/universe_grampos.xml.gz',
        'data/generated/universe_gramneg.xml.gz',
        'data/generated/universe_archaea.xml.gz',
    ]

    for file_path in downloads:
        print('Downloading', file_path)
        inputfile = source_url + file_path
        outputfile = project_dir + file_path
        urllib.request.urlretrieve(inputfile, outputfile)

    print('Building diamond database...')

    cmd = [
        'diamond',
        'makedb',
        '--in',
        project_dir + 'data/input/bigg_proteins.faa',
        '-d',
        project_dir + 'data/input/bigg_proteins'
    ]

    subprocess.call(cmd)

    with open(os.devnull, 'w') as devnull:
        try:
            exit_code = subprocess.call(cmd, stdout=devnull)
        except OSError:
            exit_code = None

    if exit_code != 0:
        print('Failed to run diamond. Please make sure diamond is installed and add it to your PATH.')
    else:
        print('Done.')


if __name__ == '__main__':
    main()
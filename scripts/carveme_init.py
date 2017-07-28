#!/usr/bin/env python

import urllib
import subprocess

from carveme import project_dir


if __name__ == '__main__':

    downloads = {
        'data/input/bigg_proteins.faa.gz': 'https://oc.embl.de/index.php/s/MeQ2XvaKo3IkEyQ',
        'data/input/equilibrator_compounds.tsv.gz': 'https://oc.embl.de/index.php/s/DJyljWecdgsd0Rm',
        'data/input/refseq_assembly_all.tsv.gz': 'https://oc.embl.de/index.php/s/kEfngj0H5DvjJR4',
        'data/input/genbank_assembly_all.tsv.gz': 'https://oc.embl.de/index.php/s/b5VCK9tfhW9TSgd',
        'data/generated/bigg_gprs.csv.gz': 'https://oc.embl.de/index.php/s/VOJ165bbkfCro2A',
        'data/generated/model_specific_data.csv.gz': 'https://oc.embl.de/index.php/s/83HjrtEIny6tAjy',
        'data/generated/universe_draft.xml.gz': 'https://oc.embl.de/index.php/s/jyAMU4JXB9Rwv9Q',
        'data/generated/universe_bacteria.xml.gz': 'https://oc.embl.de/index.php/s/NEEZzHnWVfhNBWz',
        'data/generated/universe_cyanobacteria.xml.gz': 'https://oc.embl.de/index.php/s/QfehkO5JrGI5epW'
    }

    for file_path, file_url in downloads.items():
        print 'Downloading', file_path
        outputfile = project_dir + file_path
        urllib.urlretrieve(file_url, outputfile)

    print 'Building diamond database...'

    cmd = ['diamond', 'makedb', '--in', 'data/input/bigg_proteins.faa.gz', '-d', 'data/input/bigg_proteins']
    subprocess.call(cmd)

    print 'Done.'

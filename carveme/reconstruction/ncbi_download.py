from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
import pandas as pd
import urllib.request, urllib.parse, urllib.error
import os


def load_ncbi_table(inputfile):
    table = pd.read_csv(inputfile, sep='\t', index_col=0, dtype=str)
    return table


def download_ncbi_genome(accession, refseq_table, prefer_protein=True, overwrite=False):

    if accession not in refseq_table.index:
        print('Invalid accession code')
        return

    entry = refseq_table.loc[accession, :]

    downloaded = False

    if prefer_protein:
        url = 'https://{}/{}_protein.faa.gz'.format(
            entry['ftp_path'][6:], entry['ftp_path'].split('/')[-1])

        outputfile = '{}.faa.gz'.format(accession)

        if os.path.exists(outputfile) and not overwrite:
            print('File exists, skipping.')
            return outputfile

        _, result = urllib.request.urlretrieve(url, outputfile)

        if result.get_content_type() != 'application/x-gzip':
            os.remove(outputfile)
        else:
            downloaded = True

    if not downloaded:
        url = 'https://{}/{}_genomic.fna.gz'.format(
            entry['ftp_path'][6:], entry['ftp_path'].split('/')[-1])

        outputfile = '{}.fna.gz'.format(accession)

        if os.path.exists(outputfile) and not overwrite:
            print('File exists, skipping.')
            return outputfile

        _, result = urllib.request.urlretrieve(url, outputfile)

        if result.get_content_type() != 'application/x-gzip':
            os.remove(outputfile)
        else:
            downloaded = True

    if downloaded:
        return outputfile

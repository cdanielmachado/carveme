import pandas as pd
import urllib
import os


def load_ncbi_table(inputfile):
    table = pd.read_csv(inputfile, sep='\t', index_col=0, dtype=str)
    return table


def download_ncbi_genome(accession, refseq_table, prefer_protein=True):

    if accession not in refseq_table.index:
        print 'Invalid accession code'
        return

    entry = refseq_table.loc[accession, :]

    downloaded = False

    if prefer_protein:
        url = 'https://{}/{}_protein.faa.gz'.format(
            entry['ftp_path'][6:], entry['ftp_path'].split('/')[-1])

        outputfile = '{}.faa.gz'.format(accession)

        _, result = urllib.urlretrieve(url, outputfile)

        if result.gettype() != 'application/x-gzip':
            os.remove(outputfile)
        else:
            downloaded = True

    if not downloaded:
        url = 'https://{}/{}_genomic.fna.gz'.format(
            entry['ftp_path'][6:], entry['ftp_path'].split('/')[-1])

        outputfile = '{}.fna.gz'.format(accession)

        _, result = urllib.urlretrieve(url, outputfile)

        if result.gettype() != 'application/x-gzip':
            os.remove(outputfile)
        else:
            downloaded = True

    if downloaded:
        return outputfile

#!/usr/bin/env python

from carveme import config, project_dir
import argparse
import os

from carveme.universe.download import download_universal_model, download_model_specific_data


def main():
    parser = argparse.ArgumentParser(description="Generate universal model to use with CarveMe")
    parser.add_argument('-o', '--output', dest='output', help="Output file")

    args = parser.parse_args()

    if args.output:
        universe_draft = args.output
        model_specific_data = os.path.splitext(args.output)[0] + '.csv'
        bigg_gprs = os.path.splitext(args.output)[0] + '_gprs.csv'
        fasta_file = os.path.splitext(args.output)[0] + '.faa'
        gene_annotations = os.path.splitext(args.output)[0] + '.tsv'
    else:
        universe_draft = project_dir + config.get('generated', 'bigg_universe')
        model_specific_data = project_dir + config.get('generated', 'model_specific_data')
        bigg_gprs = project_dir + config.get('generated', 'bigg_gprs')
        fasta_file = project_dir + config.get('generated', 'fasta_file')
        gene_annotations = project_dir + config.get('generated', 'gene_annotations')

    cpd_annotation = project_dir + config.get('input', 'mnx_compounds')
    download_universal_model(universe_draft, cpd_annotation)
    download_model_specific_data(model_specific_data, bigg_gprs, fasta_file, gene_annotations)


if __name__ == '__main__':
    main()


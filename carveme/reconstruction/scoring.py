import pandas as pd


def merge_subunits(genes):
    """ Merge list of protein subunit genes into complex

    Args:
        genes (pandas.Series): list of genes

    Returns:
        str: boolean rule

    """
    genes = genes.dropna()

    if len(genes) == 0:
        return None
    else:
        protein = ' and '.join(sorted(genes))
        if len(genes) > 1:
            return '(' + protein + ')'
        else:
            return protein


def merge_subunit_scores(scores):
    """ Merge scores of all genes in a protein complex.
    Calculates the mean score among all subunits.

    Args:
        scores: individual gene scores

    Returns:
        float: merged score

    """

    return scores.fillna(0).mean()


def merge_proteins(proteins):
    """ Merge all isozymes that catalyze a given reaction.
    Automatically removes all isozymes with missing score.

    Args:
        proteins (pandas.Series): list of proteins

    Returns:
        str: boolean rule

    """
    proteins = set(proteins.dropna())
    if not proteins:
        return None

    gpr_str = ' or '.join(sorted(proteins))

    if len(proteins) > 1:
        return '(' + gpr_str + ')'
    else:
        return gpr_str


def merge_protein_scores(scores):
    """ Merge scores of all isozymes that catalyze a given reaction.
    Calculates the maximum score among all isozymes.

    Args:
        scores (pandas.Series): protein scores

    Returns:
        float: merged score

    """
    return scores.max(skipna=True)


def reaction_scoring(annotation, gprs, spontaneous_score=0.0, debug_output=None):
    """ Calculate reaction scores using new eggnog output.

    Args:
        annotation (pandas.DataFrame): gene annotation results
        gprs (pandas.DataFrame): BiGG GPR rules
        spontaneous_score (float): score to give to spontaneous reactions (default: 0.0)

    Returns:
        pandas.DataFrame: reaction scores
    """

    # filter best match for each gene
    gene2gene = annotation.sort_values(by='score', ascending=False) \
                          .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])

    # merge with gpr table
    gprs['BiGG_gene'] = gprs.apply(lambda row: f"{row['model']}.{row['gene'][2:]}", axis=1)
    gene_scores = pd.merge(gene2gene, gprs, how='right')

    # add default scores for spontaneous genes
    spontaneous = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_KPN_SPONT'}
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'score'] = spontaneous_score
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'query_gene'] = 'spontaneous'

    # from gene to protein scores
    protein_scores = gene_scores.groupby(['protein', 'reaction', 'model'], as_index=False) \
        .agg({'query_gene': merge_subunits, 'score': merge_subunit_scores})

    protein_scores.rename(columns={'query_gene': 'GPR'}, inplace=True)

    # from protein to reaction scores
    reaction_scores = protein_scores.groupby(['reaction'], as_index=False) \
        .agg({'GPR': merge_proteins, 'score': merge_protein_scores}).dropna()

    avg_score = reaction_scores['score'].median()

    if avg_score == 0:
        return None, gene2gene

    reaction_scores['normalized_score'] = reaction_scores['score'] / avg_score

    if debug_output:
        gene_scores.to_csv(debug_output + '_gene_scores.tsv', sep='\t', index=False)
        protein_scores.to_csv(debug_output + '_protein_scores.tsv', sep='\t', index=False)
        reaction_scores.to_csv(debug_output + '_reaction_scores.tsv', sep='\t', index=False)

    return reaction_scores, gene2gene

from reframed import CBModel, Compartment, Metabolite, CBReaction, save_cbmodel
from reframed.io.sbml import parse_gpr_rule
from ..reconstruction.utils import to_rdf_annotation
import pandas as pd
import requests
import sys

UNIVERSE_URL = 'http://bigg.ucsd.edu/static/namespace/universal_model.json'
COMPARTMENTS_URL = 'http://bigg.ucsd.edu/api/v2/compartments'
METABOLITES_URL = 'http://bigg.ucsd.edu/api/v2/universal/metabolites/'
REACTIONS_URL = 'http://bigg.ucsd.edu/api/v2/universal/reactions/'
MODELS_URL = 'http://bigg.ucsd.edu/api/v2/models'

MAX_GPR_TOKENS = 50


def progress(i, n):
    p = int((i+1)*100.0 / n)
    sys.stdout.write(f"\r{p}%")
    sys.stdout.flush()


def get_request(url, max_tries=10):
    """ Get JSON data from BiGG RESTful API.

    Args:
        url (str): url request
        max_tries (int): maximum number of communication attempts (default: 10)

    Returns:
        dict: json data

    """
    resp, data = None, None

    for i in range(max_tries):
        try:
            resp = requests.get(url)
        except:
            pass
        if resp is not None:
            data = resp.json()
            break

    if data is None:
        print('max number of attempts exceeded:', max_tries)
        print('try again? [y]/n')
        resp = input()
        if resp.lower() != 'n':
            data = get_request(url)

    return data


def extract_annotation(elem, data):
    annotation = to_rdf_annotation(elem.id, [x[1] for x in data])
    elem.metadata['XMLAnnotation'] = annotation


def load_compartments(model):
    compartments = get_request(COMPARTMENTS_URL)
    compartments = sorted(compartments['compartments'], key=lambda x: x['bigg_id'])

    for entry in compartments:
        c_id = 'C_' + entry['bigg_id']
        comp = Compartment(c_id, entry['name'], external=(c_id == 'C_e'))
        comp.metadata['SBOTerm'] = 'SBO:0000290'
        model.add_compartment(comp)


def load_metabolites(json_model, model, cpds):
    metabolites = sorted(json_model['metabolites'], key=lambda x: x['id'])

    for entry in metabolites:
        c_id = 'C_' + entry['id'].split('_')[-1]
        m_id = 'M_' + entry['id']
        met = Metabolite(m_id, str(entry['name']), c_id)
        met.metadata['SBOTerm'] = 'SBO:0000247'
        extract_annotation(met, entry['annotation'])
        if m_id[2:-2] in cpds.index:
            met.metadata['FORMULA'] = cpds.loc[m_id[2:-2], "formula"]
            met.metadata['CHARGE'] = str(int(cpds.loc[m_id[2:-2], "charge"]))
        model.add_metabolite(met)


def is_pseudo(rxn):
    prefix = rxn.split('_', 1)[0].upper()
    return prefix in ['ATPM', 'BIOMASS', 'DM', 'EX', 'SK']


def load_reactions(json_model, model):
    reactions = sorted(json_model['reactions'], key=lambda x: x['id'])

    for entry in reactions:
        if is_pseudo(entry['id']):
            continue

        r_id = 'R_' + entry['id']
        stoichiometry = {'M_' + met: coeff for met, coeff in entry['metabolites'].items()}
        rxn = CBReaction(r_id, str(entry['name']), stoichiometry=stoichiometry)
        extract_annotation(rxn, entry['annotation'])
        model.add_reaction(rxn)

        if len(model.get_reaction_compartments(r_id)) > 1:
            rxn.metadata['SBOTerm'] = 'SBO:0000185'
        else:
            rxn.metadata['SBOTerm'] = 'SBO:0000176'


def download_universal_model(outputfile, cpd_annotation):
    print("Downloading BiGG universe...")

    cpds = pd.read_csv(cpd_annotation, sep="\t", index_col=0)

    json_model = get_request(UNIVERSE_URL)
    model = CBModel("bigg_universe")
    load_compartments(model)
    load_metabolites(json_model, model, cpds)
    load_reactions(json_model, model)
    save_cbmodel(model, outputfile)


def download_model_specific_data(outputfile, bigg_gprs, fastafile, annotations):

    data = get_request(MODELS_URL)
    models = sorted(data['results'], key=lambda x: x['gene_count'], reverse=True)
    no_seq_data = ['iAM_Pk459', 'iRC1080', 'iYS1720', 'RECON1', 'Recon3D', 'iJB785', 'iAB_RBC_283']
    filtered = {}

    for model in models:
        species = ' '.join(model['organism'].split()[:2])
        if species not in filtered and model['bigg_id'] not in no_seq_data:
            filtered[species] = model

    m = len(filtered)
    model_data = []
    sequences = {}
    gene_links = []

    print("Downloading model-specific data...")

    for i, (species, entry) in enumerate(sorted(filtered.items())):
        print(f'Downloading {species} [{i+1}/{m}]')
        model_id = entry['bigg_id']
        model = get_request(f'{MODELS_URL}/{model_id}/download')

        for rxn in model['reactions']:

            if is_pseudo(rxn['id']):
                continue

            lb, ub = rxn['lower_bound'], rxn['upper_bound']
            subsystem = rxn.get('subsystem', '')
            gpr_tokens = len(rxn['gene_reaction_rule'].split())
            gpr = rxn['gene_reaction_rule'] if gpr_tokens <= MAX_GPR_TOKENS else ''
            rxn_data = (model_id, rxn['id'], lb, ub, subsystem, gpr)
            model_data.append(rxn_data)

        n = len(model['genes'])
        for j, gene in enumerate(model['genes']):
            progress(j, n)
            gene_id = gene["id"]
            gene_data = get_request(f'{MODELS_URL}/{model_id}/genes/{gene_id}')

            seq = gene_data['protein_sequence']
            if seq is not None:
                sequences[(model_id, gene_id)] = seq

            if 'database_links' in gene_data:
                for links in gene_data['database_links'].values():
                    for link in links:
                        gene_links.append((model_id, gene_id, link['link']))

        print('')

    df = pd.DataFrame(model_data, columns=['model', 'reaction', 'lb', 'ub', 'subsystem', 'gpr'])
    df.to_csv(outputfile, index=False)

    df2 = pd.DataFrame(gene_links, columns=['model', 'gene', 'annotation'])
    df2.to_csv(annotations, index=False, sep='\t')

    write_fasta(sequences, fastafile)
    create_gpr_table(df, bigg_gprs)



def write_fasta(sequences, fastafile):

    with open(fastafile, 'w') as f:
        for (model_id, gene), seq in sequences.items():
            f.write(f'>{model_id}.{gene}\n')
            f.write(seq + '\n')


def create_gpr_table(data, outputfile):
    """ Extract GPR associations from data into a relational database format.

    Note:
        The boolean rules are converted to relational format: (Gene, Protein, Reaction). Since GPRs don't
        contain protein identifiers, the protein id is a concatenation of all subunit gene ids.
        Pseudo-genes corresponding to spontaneous reactions are discarded.

    Examples:
        The rule ((G1 and G2) or G3) -> R1, becomes:

        (G1, G1:G2, R1)
        (G2, G1:G2, R1)
        (G3, G3,    R1)

    Args:
        data (pandas.DataFrame): model specific data
        outputfile (str): output TSV file

    Returns:
        pandas.DataFrame: GPR association table
    """
    rows = []
    spontaneous = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_KPN_SPONT'}

    print("Print creating GPR table...")

    n = data.shape[0]

    for i in data.index:
        progress(i, n)

        r_id = 'R_' + data.loc[i, 'reaction']
        model_id = data.loc[i, 'model']
        gpr_str = data.loc[i, 'gpr']

        if pd.notnull(gpr_str) and gpr_str != '':
            gpr = parse_gpr_rule(gpr_str, prefix='G_')
            for protein in gpr.proteins:
                genes = sorted(set(protein.genes) - spontaneous)
                p_id = 'P_' + '+'.join([gene[2:] for gene in genes])
                for gene in genes:
                    rows.append((gene, p_id, r_id, model_id))

    print('\n')

    columns = ['gene', 'protein', 'reaction', 'model']
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(outputfile, index=False)

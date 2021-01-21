from reframed import CBModel, Compartment, Metabolite, CBReaction, save_cbmodel
from reframed.io.sbml import parse_gpr_rule
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
    for key, value in data:
        if key in elem.metadata:
            elem.metadata[key] += (';' + value)
        else:
            elem.metadata[key] = value


def load_compartments(model):
    compartments = get_request(COMPARTMENTS_URL)
    compartments = sorted(compartments['compartments'], key=lambda x: x['bigg_id'])

    for entry in compartments:
        c_id = 'C_' + entry['bigg_id']
        comp = Compartment(c_id, entry['name'], external=(c_id == 'C_e'))
        model.add_compartment(comp)


def load_metabolites(json_model, model):
    metabolites = sorted(json_model['metabolites'], key=lambda x: x['id'])

    for entry in metabolites:
        c_id = 'C_' + entry['id'].split('_')[-1]
        m_id = 'M_' + entry['id']
        met = Metabolite(m_id, str(entry['name']), c_id)
        extract_annotation(met, entry['annotation'])
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


def download_universal_model(outputfile):
    print("Downloading BiGG universe...")

    json_model = get_request(UNIVERSE_URL)
    model = CBModel("bigg_universe")
    load_compartments(model)
    load_metabolites(json_model, model)
    load_reactions(json_model, model)
    save_cbmodel(model, outputfile)


def download_model_specific_data(outputfile, bigg_gprs, fastafile):

    data = get_request(MODELS_URL)
    models = sorted(data['results'], key=lambda x: x['gene_count'], reverse=True)
    no_seq_data = ['iAM_Pk459', 'iRC1080', 'iYS1720', 'RECON1', 'Recon3D', 'iJB785', 'iAB_RBC_283']
    filtered = {}

    for model in models:
        species = ' '.join(model['organism'].split()[:2])
        if species not in filtered and model['bigg_id'] not in no_seq_data:
            filtered[species] = model

    n = len(filtered)
    model_data = []
    sequences = {}
    print("Downloading model-specific data...")

    for i, (species, entry) in enumerate(sorted(filtered.items())):
        progress(i, n)

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

        for gene in model['genes']:
            gene_id = gene["id"]
            gene_data = get_request(f'{MODELS_URL}/{model_id}/genes/{gene_id}')
            seq = gene_data['protein_sequence']

            if seq is not None:
                sequences[(model_id, gene_id)] = seq

    df = pd.DataFrame(model_data, columns=['model', 'reaction', 'lb', 'ub', 'subsystem', 'gpr'])
    df.to_csv(outputfile, index=False)
    write_fasta(sequences, fastafile)

    print('\n')

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

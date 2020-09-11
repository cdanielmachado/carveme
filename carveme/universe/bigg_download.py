import requests
import pandas as pd
import sys

from reframed import CBModel, CBReaction, Metabolite, Compartment
from reframed import save_cbmodel
from reframed.io.sbml import parse_gpr_rule


base_url = 'http://bigg.ucsd.edu/api/v2/'
reactions_url = base_url + 'universal/reactions/'
metabolites_url = base_url + 'universal/metabolites/'
compartments_url = base_url + 'compartments/'
models_url =  base_url + 'models/'


def progress(i, n):
    p = int((i+1)*100.0 / n)
    sys.stdout.write("\r{}%".format(p))
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
    elif i > 0:
        print('failed attempts', i)
    
    return data


def build_reaction(model, bigg_id, ignore_pseudo_reactions=True, only_pseudo_reactions=False):
    """ Retrieve reaction data from BiGG and add reaction to model.

    Args:
        model (CBModel): model to store the reaction
        bigg_id (str): reaction id
        ignore_pseudo_reactions (bool): ignore reactions with pseudoreaction attribute (default: True)

    """

    data = get_request(reactions_url + bigg_id)

    if data['pseudoreaction'] and ignore_pseudo_reactions:
        return

    if not data['pseudoreaction'] and only_pseudo_reactions:
        return

    r_id = 'R_' + str(bigg_id)
    name = str(data['name'])
    reversible = True  # reversibility is not universal
    stoichiometry = []

    for entry in data['metabolites']:
        comp = str(entry['compartment_bigg_id'])
        c_id = 'C_{}'.format(comp)

        if c_id not in model.compartments:
            build_compartment(model, comp)

        met = entry['bigg_id']
        m_id = 'M_{}_{}'.format(met, comp)

        if m_id not in model.metabolites:
            build_metabolite(model, met, comp)

        stoichiometry.append((m_id, entry['stoichiometry']))

    reaction = CBReaction(r_id, name, reversible, stoichiometry)
    model.add_reaction(reaction)
    extract_reaction_metadata(reaction, data)


def build_compartment(model, bigg_id):
    """ Retrieve compartment data from BiGG and add compartment to model.

    Args:
        model (CBModel): model to store the compartment
        bigg_id (str): compartment id

    """

    data = get_request(compartments_url + bigg_id)
    c_id = 'C_' + bigg_id
    name = str(data['name'])
    compartment = Compartment(c_id, name)
    model.add_compartment(compartment)


def build_metabolite(model, bigg_id, compartment):
    """ Retrieve metabolite data from BiGG and add metabolite to model.

    Args:
        model (CBModel): model to store the metabolite
        bigg_id (str): metabolite id
        compartment (str): id of compartment where the metabolite is located

    """
    data = get_request(metabolites_url + bigg_id)
    m_id = 'M_{}_{}'.format(bigg_id, compartment)
    c_id = 'C_' + compartment
    name = str(data['name'])
    metabolite = Metabolite(m_id, name, c_id)
    model.add_metabolite(metabolite, clear_tmp=False)
    extract_metabolite_metadata(metabolite, data, compartment)


def extract_reaction_metadata(reaction, data):
    """ Extract and store metadata information for a given reaction.
    Currently only stores the ids of models containing this reaction.

    Args:
        reaction (Reaction): model reaction
        data (dict): json data from BiGG reaction

    """
    extract_database_links(reaction, data)
    models = data['models_containing_reaction']
    model_ids = [entry['bigg_id'] for entry in models]
    reaction.metadata['BiGG models'] = str(';'.join(model_ids))


def extract_metabolite_metadata(metabolite, data, compartment):
    """ Extract and store metadata information for a given metabolite.
    Currently stores *formula*, *charge*, crossref links and ids of models containing this metabolite.
    Note: only returns the models where metabolite is in the given compartment.

    Args:
        metabolite (Metabolite): model reaction
        data (dict): json data from BiGG reaction
        compartment (str): compartment id

    """

    if data['charges']:
        charge = ';'.join(map(str, data['charges']))
        metabolite.metadata['CHARGE'] = charge

    if data['formulae']:
        formula = str(';'.join(data['formulae']))
        metabolite.metadata['FORMULA'] = formula

    extract_database_links(metabolite, data)

    models = data['compartments_in_models']
    models = [entry['model_bigg_id'] for entry in models
              if entry['bigg_id'] == compartment]
    metabolite.metadata['BiGG models'] = str(';'.join(models))


def extract_database_links(elem, data):
    external = data['database_links']
    for source, entry in external.items():
        key = str(source)
        values = [item['id'] for item in entry]
        elem.metadata[key] = str(';'.join(values))


def build_bigg_universe_model(outputfile=None):
    """ Download the whole BiGG universe database as a CBModel and (optionally) store in SBML.

    Args:
        outputfile (str): SBML output file (optional)

    Returns:
        CBModel: universe model
    """

    print('Downloading universal data from BiGG...')
    model = CBModel('bigg_universe')
    bigg_rxns = get_request(reactions_url)

    n = len(bigg_rxns['results'])
    for i, entry in enumerate(bigg_rxns['results']):
        build_reaction(model, entry['bigg_id'])
        progress(i, n)

    print('\n')

    if outputfile:
        save_cbmodel(model, outputfile)

    return model


def download_model_specific_data(outputfile=None):
    """ Download model specific data from BiGG (i.e. all data that is not in universal section of BiGG).
    Currently includes reaction flux bounds, subsystems and GPR associations.

    Args:
        outputfile (str): output CSV file (optional)

    Returns:
        pandas.DataFrame: model specific data

    Notes:
        It currently discards GPRs longer than 1000 characters. This avoids some troublemakers from RECON1, which
        are too complex to be parsed by sympy within reasonable time.

    """

    print('Downloading model-specific data from BiGG...')

    data = get_request(models_url)
    rows = []

    n = len(data['results'])

    for i, entry in enumerate(data['results']):
        model_id = entry['bigg_id']
        model = get_request('{}{}/download'.format(models_url, model_id))

        for reaction in model['reactions']:
            r_id = reaction['id']
            if r_id[:-1].endswith('_copy'):
                r_id = r_id[:-6]
            lb = reaction.pop('lower_bound', None)
            ub = reaction.pop('upper_bound', None)
            subsystem = reaction.pop('subsystem', None)
            gpr = reaction.pop('gene_reaction_rule', None)

            if len(gpr) > 1000:
                gpr = None  # avoid troublemakers (mostly from RECON1)

            rows.append((r_id, model_id, lb, ub, subsystem, gpr))

        progress(i, n)
    print('\n')

    columns = ['reaction', 'model',  'lower_bound', 'upper_bound', 'subsystem', 'GPR']
    df = pd.DataFrame(rows, columns=columns)
    df.sort_values(by='reaction', inplace=True)

    if outputfile:
        df.to_csv(outputfile, index=False)

    return df


def create_gpr_table(model_specific_data, reactions=None, outputfile=None):
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
        model_specific_data (pandas.DataFrame): model specific data
        reactions (list): only extract data for given reactions (optional)
        outputfile (str): output CSV file (optional)

    Returns:
        pandas.DataFrame: GPR association table
    """
    rows = []
#    spontaneous = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_KPN_SPONT'}

    for i, row in model_specific_data.iterrows():
        rxn, model_id, _, _, _, gpr_str = row
        r_id = 'R_' + rxn
        if (reactions is None or r_id in reactions) and pd.notnull(gpr_str):
            gpr = parse_gpr_rule(gpr_str)
            for protein in gpr.proteins:
                genes = sorted(set(protein.genes))
#                genes = sorted(set(protein.genes) - spontaneous)
                p_id = 'P_' + ':'.join([gene[2:] for gene in genes])
                for gene in genes:
                    rows.append((gene, p_id, r_id, model_id))

    columns = ['gene', 'protein', 'reaction', 'model']
    df = pd.DataFrame(rows, columns=columns)

    if outputfile:
        df.to_csv(outputfile, index=False)

    return df

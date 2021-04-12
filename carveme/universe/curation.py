from carveme.reconstruction.utils import create_exchange_reactions, add_biomass_equation, add_maintenance_atp, create_sink_reactions
from reframed import save_cbmodel, simplify
from reframed.core.transformation import disconnected_metabolites
from reframed.core.elements import parse_formula, molecular_weight
from .annotate import compute_missing_formulae, fix_protons_and_charge, fix_hydrogen_stoichiometry

import pandas as pd


def remove_compartments(model, taxa):

    print(f'Removing compartments that do not belong to {taxa}..')

    if taxa == 'cyanobacteria':
        compartments = {'C_c', 'C_p', 'C_u', 'C_e'}
    else:
        compartments = {'C_c', 'C_p', 'C_e'}

    rxns_to_remove = []
    for r_id in model.reactions:
        rxn_comps = model.get_reaction_compartments(r_id)
        if len(set(rxn_comps) - compartments) > 0:
            rxns_to_remove.append(r_id)

    model.remove_reactions(rxns_to_remove)

    to_remove = set(model.compartments.keys()) - compartments
    model.remove_compartments(to_remove)

    print(f'Current model size: {len(model.metabolites)} x {len(model.reactions)}')


def filter_reactions_by_kingdom(model, taxa, model_specific_data, bigg_models):

    print(f'Removing reactions that do not belong to {taxa}..')

    df = pd.merge(model_specific_data, bigg_models)
    df = df.groupby('kingdom').agg({'reaction': set})

    if taxa == 'archaea':
        valid_rxns = df.loc['Archaea', 'reaction'] | df.loc['Bacteria', 'reaction']
    else:
        valid_rxns = df.loc['Bacteria', 'reaction']

    rxns_to_remove = set(model.reactions) - valid_rxns
    model.remove_reactions(rxns_to_remove)

    mets_to_remove = disconnected_metabolites(model)
    model.remove_metabolites(mets_to_remove)

    print(f'Current model size: {len(model.metabolites)} x {len(model.reactions)}')


def clean_up_atp_synthases(model):
    rxns = set(model.search_reactions("R_ATPS")) - {"R_ATPS4rpp"}
    model.remove_reactions(rxns)


def elemental_balance(formulas, tol=1e-6):
    """ Calculate elemental balance for a given reaction.

    Args:
        formulas (list): formulas of compounds and respective stoichiometric coefficients (as tuples)
        tol (float): absolute tolerance to consider an element to be balanced (default: 1e-6)

    Returns:
        dict: elements and respective balance (only if element is unbalanced)
    """

    balance = {}
    for formula, coeff in formulas:
        elements = parse_formula(formula)
        for element, number in elements.items():
            if element in balance:
                balance[element] += coeff * number
            else:
                balance[element] = coeff * number

    balance_clean = {}
    for element, coeff in balance.items():
        if abs(coeff) > tol:
            balance_clean[element] = coeff

    return balance_clean


def check_elemental_balance(model, r_id, allow_groups=False, proton_relax=False):
    balanced = False
    formulas = []

    for m_id, coeff in model.reactions[r_id].stoichiometry.items():
        metabolite = model.metabolites[m_id]

        formula = metabolite.metadata.get('FORMULA', '')

        if formula == '':
            break

        if '*' in formula:
            if allow_groups:
                formula = formula.replace('*', '')
            else:
                break

        formulas.append((formula, coeff))

    else:
        balance = elemental_balance(formulas)

        if balance == {}:
            balanced = True
        elif proton_relax:
            balanced = set(balance.keys()) == {'H'}

    return balanced


def remove_unbalanced_reactions(model, proton_relax=False):

    print('Removing unbalanced reactions..')

    unbalanced = [r_id for r_id in model.reactions
                  if not check_elemental_balance(model, r_id, proton_relax=proton_relax)]
    print(f'found {len(unbalanced)} reactions')
    model.remove_reactions(unbalanced)

    mets_to_remove = disconnected_metabolites(model)
    model.remove_metabolites(mets_to_remove)

    print(f'Current model size: {len(model.metabolites)} x {len(model.reactions)}')


def curate_transport_by_weight(model, threshold=1500):
    """
    Removes the unrealistic transport of very large molecules.
    By default >1500 Da, this allows for transport of cobalamin (1330 Da) but not larger macromolecules.

    """
    rxns_to_delete = set()
    mets_to_delete = []

    for m_id, met in model.metabolites.items():

        if not m_id.endswith('_e'):
            continue

        mass = molecular_weight(met.metadata.get('FORMULA', ''))

        if mass < threshold:
            continue

        mets_to_delete.append(m_id)
        rxns_to_delete |= set(model.get_metabolite_reactions(m_id))

    model.remove_reactions(rxns_to_delete)
    model.remove_metabolites(mets_to_delete)


def constrain_reversibility(model, model_specific_data, manually_curated=None):

    if manually_curated is None:
        manually_curated = {}

    def lb_consensus(lbs, default_lb=-1000, abstol=1e-6):
        n_open = len([lb for lb in lbs if lb < -abstol])
        n_closed = len([lb for lb in lbs if lb >= -abstol])
        if n_closed > n_open:
            return 0
        else:
            return default_lb

    def ub_consensus(ubs, default_ub=1000, abstol=1e-6):
        n_open = len([ub for ub in ubs if ub > abstol])
        n_closed = len([ub for ub in ubs if ub <= abstol])
        if n_closed > n_open:
            return 0
        else:
            return default_ub

    df = model_specific_data.groupby("reaction").agg({"lb": lb_consensus, "ub": ub_consensus})

    for r_id, rxn in model.reactions.items():
        if r_id in manually_curated.index:
            rxn.lb = manually_curated.loc[r_id, 'lb']
            rxn.ub = manually_curated.loc[r_id, 'ub']
            rxn.reversible == rxn.lb < 0
        elif r_id in df.index:
            rxn.lb = df.loc[r_id, 'lb']
            rxn.ub = df.loc[r_id, 'ub']
            rxn.reversible == rxn.lb < 0


def reversibility_heuristics(model):
    """ Apply heuristic rules to constrain reaction reversibilities. """

    # ABC transporters

    for r_id, rxn in model.reactions.items():
        if 'M_atp_c' in rxn.get_substrates() and len(model.get_reaction_compartments(r_id)) == 2:
            model.set_flux_bounds(r_id, 0, 1000)

    # proton pumps

    h_dir = {}
    pump_rxns = {}
    other = set()

    for r_id, rxn in model.reactions.items():
        substrates = {m_id[2:-2]: m_id[-1] for m_id in rxn.get_substrates()}
        products = {m_id[2:-2]: m_id[-1] for m_id in rxn.get_products()}

        if set(substrates) != set(products):
            continue

        if 'h' in substrates:
            pump_rxns[r_id] = set(substrates)
            h_dir[r_id] = (substrates['h'], products['h'])
        else:
            other.update(substrates)

    for r_id, cpds in pump_rxns.items():

        if len(cpds & other) == 0: # skip compounds that do not have other transport mechanisms
            continue

        if h_dir[r_id] in {('c', 'p'), ('c', 'e'), ('p', 'e')}:
            model.set_flux_bounds(r_id, -1000, 0)
        if h_dir[r_id] in {('p', 'c'), ('e', 'c'), ('e', 'p')}:
            model.set_flux_bounds(r_id, 0, 1000)


def curate_universe(model, outputfile, model_specific_data, bigg_models, taxa, biomass_eq,
                    manually_curated=None, unbalanced_metabolites=None):

    """ Curate universal reaction database from initial database dump.

    Args:
        model (CBModel): universal model
        outputfile (str): output SBML file (optional)
        model_specific_data (pandas.DataFrame): model specific data downloaded from BiGG
        bigg_models (pandas.DataFrame): Additional information on BiGG models
        taxa (str): filter by taxa (optional)
        biomass_eq (str): default biomass equation

    """

    print(f'Curating {taxa} universe...')
    print(f'Initial model size: {len(model.metabolites)} x {len(model.reactions)}')

    remove_compartments(model, taxa)

    model_specific_data['reaction'] = model_specific_data['reaction'].apply(lambda x: 'R_' + x)
    filter_reactions_by_kingdom(model, taxa, model_specific_data, bigg_models)

    clean_up_atp_synthases(model)

    compute_missing_formulae(model)

    curate_transport_by_weight(model)

    remove_unbalanced_reactions(model, proton_relax=True)

    fix_protons_and_charge(model)

    fix_hydrogen_stoichiometry(model)

    constrain_reversibility(model, model_specific_data, manually_curated)

    reversibility_heuristics(model)

    if unbalanced_metabolites is not None:
        create_sink_reactions(model, unbalanced_metabolites)

    create_exchange_reactions(model, default_lb=-1000, default_ub=1000)

    add_biomass_equation(model, biomass_eq)

    add_maintenance_atp(model)

    print('Removing blocked reactions and dead-end metabolites...')
    simplify(model)
    print(f'Final model size: {len(model.metabolites)} x {len(model.reactions)}')

    save_cbmodel(model, outputfile)
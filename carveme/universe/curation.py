import numpy as np
from framed.model.transformation import disconnected_metabolites
from carveme.reconstruction.utils import create_exchange_reactions, add_biomass_equation, create_sink_reactions, \
    add_maintenance_atp
from framed import simplify
from framed import save_cbmodel
from carveme.universe.thermodynamics import compute_flux_bounds
from framed.experimental.elements import parse_formula


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


def check_elemental_balance(model, r_id, tol=1e-6):
    """ Check elemental balance of a reaction in a model.

    Notes:
        Since some metabolites can have multiple formulas associated, it will check if at least one possible
        combination of all formulas results in a balanced reaction.

    Args:
        model (CBModel): metabolic model
        r_id (str): reaction id
        tol (float): absolute tolerance to consider an element to be balanced (default: 1e-6)

    Returns:
        bool: balance status (False if unbalanced or unable to calculate)
        str: error message (if reaction is not balanced)

    """
    reaction = model.reactions[r_id]
    combinations = [[]]
    error = None

    for m_id, coeff in reaction.stoichiometry.items():
        metabolite = model.metabolites[m_id]
        if 'FORMULA' in metabolite.metadata:
            formulae = metabolite.metadata['FORMULA'].split(';')
            combinations = [combination + [(formula, coeff)]
                            for formula in formulae
                            for combination in combinations]
        else:
            error = 'missing formula'
            break

    if error:
        return False, error

    balance_all = [elemental_balance(combination, tol=tol) for combination in combinations]
    balanced = [balance for balance in balance_all if len(balance) == 0]

    if len(balanced) > 0:  # at least one combination must be valid
        return True, None
    else:
        return False, 'unbalanced'


def check_charge_balance(model, r_id, tol=1e-6):
    """ Check charge balance of a reaction in a model.

    Notes:
        Since some metabolites can have multiple charges associated, it will check if at least one possible
        combination of all charges results in a balanced reaction.

    Args:
        model (CBModel): metabolic model
        r_id (str): reaction id
        tol (float): absolute tolerance to consider the charge to be balanced (default: 1e-6)

    Returns:
        bool: balance status (False if unbalanced or unable to calculate)
        str: error message (if reaction is not balanced)

    """

    reaction = model.reactions[r_id]
    combinations = [0]
    error = None

    for m_id, coeff in reaction.stoichiometry.items():
        metabolite = model.metabolites[m_id]
        if 'CHARGE' in metabolite.metadata:
            charges = metabolite.metadata['CHARGE'].split(';')
            combinations = [combination + int(charge)*coeff
                            for charge in charges
                            for combination in combinations]
        else:
            error = 'missing charge'
            break

    if error:
        return False, error

    balanced = [balance for balance in combinations if abs(balance) < tol]

    if len(balanced) >= 1:  # at least one combination must be valid
        return True, None
    else:
        return False, 'unbalanced'


def remove_unbalanced_reactions(model, ignore_list=None, test_formula=True, test_charge=False):
    """ Remove elemental and/or charge unbalanced reactions from model.

    Args:
        model (CBModel): model
        ignore_list (list): list of reactions not to be tested (optional)
        test_formula (bool): test for elemental balance
        test_charge (bool): test for charge balance

    """
    unbalanced = []

    if not ignore_list:
        ignore_list = []

    mass_balance = True
    charge_balance = True

    for r_id in model.reactions:
        if r_id not in ignore_list:
            if test_formula:
                mass_balance, error_msg = check_elemental_balance(model, r_id)
            if test_charge:
                charge_balance, error_msg = check_charge_balance(model, r_id)
            if not mass_balance or not charge_balance:
                unbalanced.append(r_id)

    model.remove_reactions(unbalanced)


def add_bounds_from_extracted_data(model, data, trusted_models=None):
    """ Determine flux bounds from model-specific data extracted from BiGG.

    Notes:
        If there are multiple values, it will use the most common case.
        A list of **trusted** models can be given. These will be considered first.

    Args:
        model (CBModel): model
        data (pandas.DataFrame): model specific data
        trusted_models (list): list of trusted model ids (optional)

    Returns:

    """

    mapper = {'lower_bound': lambda x: -1000 if np.median(x) <= -500 else 0,
              'upper_bound': lambda x: 1000 if np.median(x) >= 500 else 0}

    data.loc[data['upper_bound'] > 0, 'upper_bound'] = 1000
    data.loc[data['lower_bound'] < 0, 'lower_bound'] = -1000

    merged_all = data.groupby('reaction').agg(mapper)

    if trusted_models:
        trusted_data = data.query('model in {}'.format(trusted_models))
        merged_trusted = trusted_data.groupby('reaction').agg(mapper)
    else:
        merged_trusted = merged_all

    for r_id in model.reactions.keys():
        if r_id[2:] in merged_trusted.index:
            lb = merged_trusted.loc[r_id[2:], 'lower_bound']
            ub = merged_trusted.loc[r_id[2:], 'upper_bound']
            model.reactions[r_id].trusted = True
        else:
            lb = merged_all.loc[r_id[2:], 'lower_bound']
            ub = merged_all.loc[r_id[2:], 'upper_bound']
            model.reactions[r_id].trusted = False

        model.set_flux_bounds(r_id, lb, ub)
        model.reactions[r_id].reversible = bool(lb < 0)


def filter_reactions_by_kingdoms(model, kingdoms, kingdom_map, exclusive=False, inplace=False):
    """ Filter reactions in model by Kingdoms.
    A reaction is considered valid for a given Kingdom if it is present in at least one model from that Kingdom.

    Args:
        model (CBModel): model
        kingdoms (set): Kingdom names
        kingdom_map (dict): mapping between model ids and kingdoms
        exclusive (bool): only accept reactions *exclusive* to given Kingdom (default: False)
        inplace (bool): automatically remove invalid reactions from model (default: False)

    Returns:
        list: list of reactions that belong to the given Kingdom

    """

    valid = []
    invalid = []

    for r_id, rxn in model.reactions.items():
        model_ids = rxn.metadata['BiGG models'].split(';')
        model_kingdoms = set([kingdom_map[model_id] for model_id in model_ids])

        if exclusive and model_kingdoms == kingdoms:
            valid.append(r_id)
        elif not exclusive and len(model_kingdoms & kingdoms) > 0:
            valid.append(r_id)
        else:
            invalid.append(r_id)

    if inplace:
        model.remove_reactions(invalid)
        model.remove_metabolites(disconnected_metabolites(model))

    else:
        return valid


def reversibility_heuristics(model, no_reverse_atp=True, no_proton_pumps=True, override_trusted=False):
    """ Apply heuristic rules to constrain reaction reversibilities.

    Args:
        model (CBModel): model
        no_reverse_atp (bool): no reversible ATP production from ATP consumers, except for trusted reactions (default: True)
        no_proton_pumps (bool): no pumping of electrons from cytosol to periplasm or extracellular space
    """

    if no_reverse_atp:
        for r_id, rxn in model.reactions.items():
            if rxn.trusted and not override_trusted:
                continue

            if 'M_atp_c' in rxn.get_substrates():
                model.set_flux_bounds(r_id, 0, 1000)
                rxn.reversible = False

    if no_proton_pumps:

        protons = {'M_h_c', 'M_h_p', 'M_h_e'}

        for r_id, rxn in model.reactions.items():
            if rxn.trusted and not override_trusted:
                continue

            substrates = set(rxn.get_substrates())
            products = set(rxn.get_products())
            if len(substrates) == len(products) == 2:
                if ('M_h_p' in substrates or 'M_h_e' in substrates) and 'M_h_c' in products:
                    substrate = (substrates - protons).pop()
                    product = (products - protons).pop()
                    if substrate[:-2] == product[:-2]:
                        model.set_flux_bounds(r_id, 0, 1000)
                        rxn.reversible = False
                if ('M_h_p' in products or 'M_h_e' in products) and 'M_h_c' in substrates:
                    substrate = (substrates - protons).pop()
                    product = (products - protons).pop()
                    if substrate[:-2] == product[:-2]:
                        model.set_flux_bounds(r_id, -1000, 0)


def curate_universe(model, model_specific_data, bigg_models, biomass_eq, taxa=None,
                    thermodynamics_data=None, metabolomics_data=None, thermodynamics_method=None,
                    manually_curated=None, unbalanced_metabolites=None, use_heuristics=True,
                    remove_unbalanced=True, remove_blocked=True, outputfile=None):

    """ Curate universal reaction database from initial database dump.

    Args:
        model (CBModel): universal model
        model_specific_data (pandas.DataFrame): model specific data downloaded from BiGG
        bigg_models (pandas.DataFrame): Additional information on BiGG models
        biomass_eq (str): default biomass equation
        taxa (str): filter by taxa (optional)
        thermodynamics_data (pandas.DataFrame): used for reversibility estimation (optional)
        metabolomics_data (pandas.DataFrame): used for reversibility estimation (optional)
        thermodynamics_method (str): thermodynamics method to use (optional)
        manually_curated (pandas.DataFrame): manually curated reaction bounds (optional)
        unbalanced_metabolites (list): unbalanced metabolites that require sink reactions (optional)
        use_heuristics (bool): apply heuristic rules (no proton pumps, no reversible ATP consumers) (default: True)
        remove_unbalanced (bool): remove unbalanced reactions from model (default: True)
        remove_blocked (bool): remove blocked reactions and dead-end metabolites (default: True)
        outputfile (str): output SBML file (optional)

    Returns:
        CBModel: curated universal model

    Notes:
        Combines thermodynamics and heuristic rules to determine reaction reversibility.
        Adds exchange reactions for all extracellular metabolites.
        Adds sinks reactions for a list of known unbalanced compounds.
        Adds biomass equations from local biomass database (avoids discarding biomass precursors and other
        essencial reactions/metabolites that would otherwise be structurally blocked).
    """

    print('Starting universe curation...')
    print('(initial size: {} x {})\n'.format(len(model.metabolites), len(model.reactions)))

    trusted_models = bigg_models.query('trusted == True').index.tolist()

    add_bounds_from_extracted_data(model, model_specific_data, trusted_models)

    if taxa:
        print('Filtering by taxa:', taxa)
        kingdom_map = bigg_models['domain'].to_dict()

        if taxa in {'cyanobacteria', 'bacteria'}:
            kingdoms = {'Bacteria'}
        elif taxa == 'archaea':
            kingdoms = {'Archaea', 'Bacteria'}
        else:
            raise ValueError('Unsupported taxa:' + taxa)

        filter_reactions_by_kingdoms(model, kingdoms, kingdom_map, inplace=True)

        if taxa in {'bacteria', 'archaea'}:
            valid_compartments = {'C_c', 'C_p', 'C_e'}
        elif taxa == 'cyanobacteria':
            valid_compartments = {'C_c', 'C_p', 'C_e', 'C_u'}

        other_compartments = set(model.compartments.keys()) - valid_compartments
        model.remove_compartments(other_compartments, delete_metabolites=True, delete_reactions=True)

        print('(size: {} x {})\n'.format(len(model.metabolites), len(model.reactions)))

    if thermodynamics_data is not None:
        print('Computing thermodynamics...', end=' ')

        dG0 = thermodynamics_data['dG0'].to_dict()
        sdG0 = thermodynamics_data['sdG0'].to_dict()

        if metabolomics_data is not None:
            x0 = metabolomics_data.median(axis=1).to_dict()
        else:
            x0 = None

        compute_flux_bounds(model, dG0, sdG0, x0, method=thermodynamics_method, inplace=True, override_trusted=False)
        print('done\n')


    print('Applying manual curation rules...', end=' ')

    if use_heuristics:
        reversibility_heuristics(model, no_reverse_atp=True, no_proton_pumps=False, override_trusted=False)

    # manually curated reactions
    if manually_curated is not None:
        for r_id, (lb, ub) in manually_curated.iterrows():
            if r_id in model.reactions:
                model.set_flux_bounds(r_id, lb, ub)

    print('done\n')

    if remove_unbalanced:

        # remove arbitrary 'Z' formula from photons
        if taxa == 'cyanobacteria':
            for m_id in ['M_photon_e', 'M_photon_p', 'M_photon_c']:
                model.metabolites[m_id].metadata['FORMULA'] = ''

        print('Removing unbalanced reactions...')
        remove_unbalanced_reactions(model)
        print('(size: {} x {})\n'.format(len(model.metabolites), len(model.reactions)))

    print('Creating pseudo-reactions...')

    create_exchange_reactions(model, default_lb=-1000, default_ub=1000)

    if unbalanced_metabolites:
        create_sink_reactions(model, unbalanced_metabolites)

    add_biomass_equation(model, biomass_eq)

    add_maintenance_atp(model)

    print('(size: {} x {})\n'.format(len(model.metabolites), len(model.reactions)))

    if remove_blocked:
        print('Removing blocked reactions and dead-end metabolites...')
        simplify(model)
        print('(size: {} x {})\n'.format(len(model.metabolites), len(model.reactions)))

    if outputfile:
        save_cbmodel(model, outputfile)

    print('Done.')

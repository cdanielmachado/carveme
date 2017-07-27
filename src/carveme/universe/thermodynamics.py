import pandas as pd
import numpy as np
from framed import load_cbmodel
from framed import NET, TVA
from math import sqrt

try:
    from component_contribution.kegg_model import KeggModel
    from component_contribution.component_contribution_trainer import ComponentContribution as CC
    from component_contribution.thermodynamic_constants import default_I, default_T, default_pH, R
except ImportError:
    R = 8.31e-3
    default_T = 298.15
    default_I = 0.25
    default_pH = 7.0


from carveme import config

measured_fold_change = config.getfloat('thermodynamics', 'measured_fold_change')
concentration_min = config.getfloat('thermodynamics', 'concentration_min')
concentration_max = config.getfloat('thermodynamics', 'concentration_max')


def BiGG_to_KEGG_reaction(model, r_id, kegg_compounds=None):
    """ Convert reaction from BiGG to KEGG notation to be used by eQuilibrator.

    Args:
        model (CBModel): model
        r_id (str): reaction id
        kegg_compounds (list): KEGG compounds accepted by eQuilibrator

    Returns:
        str: KEGG reaction

    Notes:
        If there are multiple KEGG compounds associated with a BiGG id it will use the first match. (As a rule of thumb,
         the first always seems to be the correct one).
        No result is returned if unable to translate any compound.
    """
    reaction = model.reactions[r_id]
    coeffs = []

    for m_id, coeff in reaction.stoichiometry.items():
        metabolite = model.metabolites[m_id]
        if 'KEGG Compound' in metabolite.metadata:
            match = False
            for kegg_id in metabolite.metadata['KEGG Compound'].split(';'):
                if kegg_compounds and kegg_id in kegg_compounds or kegg_compounds is None:
                    coeffs.append((kegg_id, coeff))
                    match = True
                    break
            if not match:
                return None
        else:
            return None

    substrates = [kegg_id if coeff == -1 else '{} {}'.format(-coeff, kegg_id)
                  for kegg_id, coeff in coeffs if coeff < 0]
    products = [kegg_id if coeff == 1 else '{} {}'.format(coeff, kegg_id)
                for kegg_id, coeff in coeffs if coeff > 0]
    rxn_str = ' + '.join(substrates) + ' <=> ' + ' + '.join(products)

    return rxn_str


def build_kegg_reactions(model, kegg_compounds):
    """ Translate reactions in model to KEGG reaction format accepted by eQuilibrator.

    Args:
        model (CBModel): model
        kegg_compounds (list): KEGG compounds accepted by eQuilibrator

    Returns:
        dict: mapping from reaction id to KEGG reaction strings
    """
    kegg_rxns = {}

    for r_id in model.reactions:
        rxn_str = BiGG_to_KEGG_reaction(model, r_id, kegg_compounds)
        if rxn_str:
            kegg_rxns[r_id] = rxn_str

    return kegg_rxns


def calculate_deltaG0s(model, kegg_compounds, pH=default_pH, I=default_I, T=default_T):
    """ Calculate standard Gibbs Energy for reactions in model (as many as possible) using eQuilibrator.

    Args:
        model (CBModel): model
        kegg_compounds (list): KEGG compounds accepted by eQuilibrator
        pH (float): pH (default: 7.0)
        I (float): ionic strenght (default: 0.25)
        T (float): temperature (default: 298.15 K (25.0 C))

    Returns:
        dict: standard Gibbs Energies indexed by reaction ids
        dict: estimation error indexed by reaction ids

    """
    kegg_rxns = build_kegg_reactions(model, kegg_compounds)
    kmodel = KeggModel.from_formulas(kegg_rxns.values(), raise_exception=True)
    kmodel.add_thermo(CC.init())
    dG0, sdG0, _ = kmodel.get_transformed_dG0(pH, I, T)

    dG0 = dict(zip(kegg_rxns.keys(), dG0.A1))
    sdG0 = dict(zip(kegg_rxns.keys(), sdG0.A1))

    return dG0, sdG0


def compute_bigg_gibbs_energy(modelfile, equi_cmpds_file, outputfile=None):
    """ Calculate standard Gibbs Energy for reactions in a model (as many as possible) using eQuilibrator.

    Args:
        modelfile (str): SBML file
        equi_cmpds_file (str): file containing KEGG compounds accepted by eQuilibrator
        outputfile (str): output CSV file (optional)

    Returns:
        dict: standard Gibbs Energies indexed by reaction ids
        dict: estimation error indexed by reaction ids
    """

    model = load_cbmodel(modelfile)

    kegg_compounds = pd.read_csv(equi_cmpds_file, sep='\t')
    kegg_compounds = set(kegg_compounds['compound_id'])

    dG0, sdG0 = calculate_deltaG0s(model, kegg_compounds, pH=default_pH, I=default_I, T=default_T)

    if outputfile:
        data = pd.DataFrame({'dG0': dG0, 'sdG0': sdG0})
        data.to_csv(outputfile)
    else:
        return dG0, sdG0


def dG_bounds(model, r_id, dG0, sdG0=None, x0=None, excluded=None, T=default_T):
    """ Estimate lower and upper bounds for standard Gibbs energy of a reaction.

    Args:
        model (CBModel): model
        r_id (str): reaction id
        dG0 (dict): standard Gibbs energies of reactions
        sdG0 (dict): estimation errors (optional)
        x0 (dict): reference metabolite concentrations (optional)
        excluded (list): compounds to exclude from calculation (waters, protons, etc.)
        T (float): temperature (default: 298.15 K (25.0 C))

    Returns:

    """

    if not x0:
        x0 = {}

    if not sdG0:
        sdG0 = {r_id: 0 for r_id in dG0}

    if not excluded:
        excluded = []

    reac_min, reac_max = [], []
    prod_min, prod_max = [], []

    for m_id, coeff in model.reactions[r_id].stoichiometry.items():

        if m_id in excluded:
            continue
        elif m_id in x0:
            x_min = x0[m_id] / sqrt(measured_fold_change)
            x_max = x0[m_id] * sqrt(measured_fold_change)
        else:
            x_min = concentration_min
            x_max = concentration_max

        if coeff < 0:
            reac_min.append(x_min ** abs(coeff))
            reac_max.append(x_max ** abs(coeff))
        else:
            prod_min.append(x_min ** coeff)
            prod_max.append(x_max ** coeff)

    dG_min = dG0[r_id] - sdG0[r_id] + R * T * np.log(np.prod(prod_min) / np.prod(reac_max))
    dG_max = dG0[r_id] + sdG0[r_id] + R * T * np.log(np.prod(prod_max) / np.prod(reac_min))

    return dG_min, dG_max


def dG_to_flux_bounds(dG_min, dG_max, infinity=1000, abstol=1e-6):
    """ Convert standard Gibbs energy range to reaction flux bounds.

    Args:
        dG_min (float): minimum standard Gibbs energy
        dG_max (float): maximum standard Gibbs energy
        infinity (float): value to represent infinity (default: 1000)
        abstol (float): absolute tolerance to consider a value larger than zero (default: 1e-6)

    Returns:

    """

    if dG_min > abstol:
        lb, ub = -infinity, 0
    elif dG_max < -abstol:
        lb, ub = 0, infinity
    else:
        lb, ub = -infinity, infinity

    return lb, ub


def compute_flux_bounds(model, dG0, sdG0=None, x0=None, method=None, inplace=False, override_trusted=False, conservative=True):
    """ Compute flux bounds for model using thermodynamics

    Args:
        model (CBModel): model
        dG0 (dict): standard Gibbs energies of reactions
        sdG0 (dict): estimation errors (optional)
        x0 (dict): reference metabolite concentrations (optional)
        method (str): thermodynamic analysis method (currently not available) # TODO: maybe remove this?
        inplace (bool): update model (default: False)
        override_trusted (bool): override bounds of reactions tagged as **trusted** (default: False)
        conservative (bool): allow constraining reference model bounds but not relaxing (default: True)

    Returns:
        dict: flux bounds

    """

    bounds = {}

    if not x0:
        x0 = {}

    if not sdG0:
        sdG0 = {r_id: 0 for r_id in dG0}

    dG0_new = {}
    sdG0_new = {}

    excluded = ['h2o', 'h']
    excluded = [m_id for m_id in model.metabolites if m_id[2:-2] in excluded]

    for r_id in model.reactions:
        if r_id not in dG0:
            continue

        if not 0 < sdG0[r_id] < 1000:
            continue

        if len(model.get_reaction_compartments(r_id)) > 1:
            continue

        dG0_new[r_id], sdG0_new[r_id] = dG0[r_id], sdG0[r_id]

    dG0, sdG0 = dG0_new, sdG0_new

    if not method:
        for r_id in dG0:
            if 'M_atp_c' in model.reactions[r_id].get_substrates(): # TODO: CHECK IF THIS IS REALLY NECESSARY!
                model.set_flux_bounds(r_id, 0, 1000)
            else:
                dG_min, dG_max = dG_bounds(model, r_id, dG0, sdG0, x0, excluded)
                bounds[r_id] = dG_to_flux_bounds(dG_min, dG_max)

    # elif method == 'NET':
    #
    #     directions = {}
    #
    #     for r_id in dG0:
    #         if 'M_atp_c' in model.reactions[r_id].get_substrates():
    #             directions[r_id] = 1
    #
    #     dGbounds = NET(model, dG0, sdG0, x0,
    #                    excluded=excluded,
    #                    reaction_directions=directions,
    #                    measured_fold_change=measured_fold_change,
    #                    concentration_min=concentration_min,
    #                    concentration_max=concentration_max)
    #
    #     for r_id, (dG_min, dG_max) in dGbounds.items():
    #         bounds[r_id] = dG_to_flux_bounds(dG_min, dG_max)
    #
    # elif method == 'TVA':
    #
    #     for r_id in dG0:
    #         if 'M_atp_c' in model.reactions[r_id].get_substrates():
    #             model.set_flux_bounds(r_id, 0, 1000)
    #         else:
    #             model.set_flux_bounds(r_id, -1000, 1000)
    #
    #     flux_bounds = TVA(model, dG0, sdG0, x0,
    #                       excluded=excluded,
    #                       concentration_max=concentration_max,
    #                       concentration_min=concentration_min,
    #                       measured_fold_change=measured_fold_change,
    #                       ignore_model_bounds=False,
    #                       reactions=dG0.keys())
    #
    #     for r_id in dG0:
    #         lb, ub = flux_bounds[r_id]
    #         lb = -1000 if lb < -1e-9 else 0
    #         ub = 1000 if ub > 1e-9 else 0
    #         bounds[r_id] = (lb, ub)

    if inplace:
        for r_id, (lb, ub) in bounds.items():
            if not model.reactions[r_id].trusted or override_trusted:
                if conservative:
                    lb = max(lb, model.reactions[r_id].lb)
                    ub = min(ub, model.reactions[r_id].ub)
                model.set_flux_bounds(r_id, lb, ub)
                model.reactions[r_id].reversible = bool(lb < 0)
    else:
        return bounds

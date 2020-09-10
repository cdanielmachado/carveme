import numpy as np
import warnings
import pandas as pd

from framed.cobra.ensemble import EnsembleModel, save_ensemble
from framed.io.sbml import parse_gpr_rule, save_cbmodel
from framed.model.transformation import disconnected_metabolites
from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status


def inactive_reactions(model, solution):
    inactive = []

    internal = [r_id for r_id in model.reactions if not r_id.startswith('R_EX')]
    external = [r_id for r_id in model.reactions if r_id.startswith('R_EX')]

    for r_id in internal:
        if ((abs(solution.values[r_id]) < 1e-6)
                and (solution.values.get('yf_' + r_id, 0) < 0.5)
                and (solution.values.get('yr_' + r_id, 0) < 0.5)):

            inactive.append(r_id)

    m_r_lookup = model.metabolite_reaction_lookup()
    inactive_ext = []

    for r_id in external:
        m_id = model.reactions[r_id].get_substrates()[0]
        neighbors = m_r_lookup[m_id]
        if len(set(neighbors) - set(inactive)) == 1:
            inactive_ext.append(r_id)

    return inactive + inactive_ext


def minmax_reduction(model, scores, min_growth=0.1, min_atpm=0.1, eps=1e-3, bigM=1e3, default_score=-1.0,
                     uptake_score=0.0, soft_score=1.0, soft_constraints=None, hard_constraints=None,
                     ref_reactions=None, ref_score=0.0, solver=None, debug_output=None):
    """ Apply minmax reduction algorithm (MILP).

    Computes a binary reaction vector that optimizes the agreement with reaction scores (maximizes positive scores,
    and minimizes negative scores). It generates a fully connected reaction network (i.e. all reactions must be able
    to carry some flux).

    Args:
        model (CBModel): universal model
        scores (dict): reaction scores
        min_growth (float): minimal growth constraint
        min_atpm (float): minimal maintenance ATP constraint
        eps (float): minimal flux required to consider leaving the reaction in the model
        bigM (float): maximal reaction flux
        default_score (float): penalty score for reactions without an annotation score (default: -1.0).
        uptake_score (float): penalty score for using uptake reactions (default: 0.0).
        soft_score (float): score for soft constraints (default: 1.0)
        soft_constraints (dict): dictionary from reaction id to expected flux direction (-1, 1, 0)
        hard_constraints (dict): dictionary of flux bounds
        solver (Solver): solver instance (optional)

    Returns:
        Solution: optimization result
    """

    if not solver:
        solver = solver_instance(model)

    objective = {}

    scores = scores.copy()
    reactions = list(scores.keys())

    if soft_constraints:
        reactions += [r_id for r_id in soft_constraints if r_id not in reactions]
    else:
        soft_constraints = {}

    if not ref_reactions:
        ref_reactions = {}

    if hard_constraints:
        solver.set_bounds(hard_constraints)

    if default_score != 0:
        for r_id in model.reactions:
            if r_id not in reactions and r_id not in ref_reactions and not r_id.startswith('R_EX') and r_id != 'R_ATPM':
                scores[r_id] = default_score
                reactions.append(r_id)

    if ref_score != 0:
        for r_id in ref_reactions:
            if r_id not in reactions and r_id != 'R_ATPM':
                scores[r_id] = ref_score
                reactions.append(r_id)

    if not hasattr(solver, '_carveme_flag'):
        solver._carveme_flag = True

        biomass = model.biomass_reaction
        solver.add_constraint('min_growth', {biomass: 1}, '>', min_growth, update_problem=False)
        solver.add_constraint('min_atpm', {'R_ATPM': 1}, '>', min_atpm, update_problem=False)

        solver.neg_vars = []
        solver.pos_vars = []

        for r_id in reactions:
            if model.reactions[r_id].lb is None or model.reactions[r_id].lb < 0:
                y_r = 'yr_' + r_id
                solver.add_variable(y_r, 0, 1, vartype=VarType.BINARY, update_problem=False)
                solver.neg_vars.append(y_r)
            if model.reactions[r_id].ub is None or model.reactions[r_id].ub > 0:
                y_f = 'yf_' + r_id
                solver.add_variable(y_f, 0, 1, vartype=VarType.BINARY, update_problem=False)
                solver.pos_vars.append(y_f)

        if uptake_score != 0:
            for r_id in model.reactions:
                if r_id.startswith('R_EX'):
                    solver.add_variable('y_' + r_id, 0, 1, vartype=VarType.BINARY, update_problem=False)

        solver.update()

        for r_id in reactions:
            y_r, y_f = 'yr_' + r_id, 'yf_' + r_id
            if y_r in solver.neg_vars and y_f in solver.pos_vars:
                solver.add_constraint('lb_' + r_id, {r_id: 1, y_f: -eps, y_r: bigM}, '>', 0, update_problem=False)
                solver.add_constraint('ub_' + r_id, {r_id: 1, y_f: -bigM, y_r: eps}, '<', 0, update_problem=False)
                solver.add_constraint('rev_' + r_id, {y_f: 1, y_r: 1}, '<', 1, update_problem=False)
            elif y_f in solver.pos_vars:
                solver.add_constraint('lb_' + r_id, {r_id: 1, y_f: -eps}, '>', 0, update_problem=False)
                solver.add_constraint('ub_' + r_id, {r_id: 1, y_f: -bigM}, '<', 0, update_problem=False)
            elif y_r in solver.neg_vars:
                solver.add_constraint('lb_' + r_id, {r_id: 1, y_r: bigM}, '>', 0, update_problem=False)
                solver.add_constraint('ub_' + r_id, {r_id: 1, y_r: eps}, '<', 0, update_problem=False)

        if uptake_score != 0:
            for r_id in model.reactions:
                if r_id.startswith('R_EX'):
                    solver.add_constraint('lb_' + r_id, {r_id: 1, 'y_' + r_id: bigM}, '>', 0, update_problem=False)

        solver.update()

    for r_id in reactions:
        y_r, y_f = 'yr_' + r_id, 'yf_' + r_id

        if r_id in soft_constraints:
            sign = soft_constraints[r_id]

            if sign > 0:
                w_f, w_r = soft_score, 0
            elif sign < 0:
                w_f, w_r = 0, soft_score
            else:
                w_f, w_r = -soft_score, -soft_score

        if y_f in solver.pos_vars:
            if r_id in soft_constraints:
                 objective[y_f] = w_f
            elif r_id in ref_reactions:
                objective[y_f] = 2 * scores[r_id] + ref_score
            else:
                objective[y_f] = scores[r_id]

        if y_r in solver.neg_vars:
            if r_id in soft_constraints:
                objective[y_r] = w_r
            elif r_id in ref_reactions:
                objective[y_r] = 2 * scores[r_id] + ref_score
            else:
                objective[y_r] = scores[r_id]

    if uptake_score != 0:
        for r_id in model.reactions:
            if r_id.startswith('R_EX') and r_id not in soft_constraints:
                objective['y_' + r_id] = uptake_score

    solver.set_objective(linear=objective, minimize=False)

    if debug_output:
        solver.write_to_file(debug_output + "_milp_problem.lp")

    solution = solver.solve()

    return solution


def carve_model(model, reaction_scores, outputfile=None, flavor=None, inplace=True,
                default_score=-1.0, uptake_score=0.0, soft_score=1.0,
                soft_constraints=None, hard_constraints=None, ref_model=None, ref_score=0.0,
                init_env=None, debug_output=None):
    """ Reconstruct a metabolic model using the CarveMe approach.

    Args:
        model (CBModel): universal model
        reaction_scores (pandas.DataFrame): reaction scores
        outputfile (str): write model to SBML file (optional)
        flavor (str): SBML flavor ('cobra' or 'fbc2', optional)
        inplace (bool): Change model in place (default: True)
        default_score (float): penalty for non-annotated intracellular reactions (default: -1.0)
        uptake_score (float): penalty for utilization of extracellular compounds (default: 0.0)
        soft_score (float): score for soft constraints (default: 1.0)
        soft_constraints (dict): dictionary from reaction id to expected flux direction (-1, 1, 0)
        hard_constraints (dict): dictionary of flux bounds
        init_env (Environment): initialize final model with given Environment (optional)

    Returns:
        CBModel: reconstructed model
    """

    if not inplace:
        model = model.copy()

    scores = dict(reaction_scores[['reaction', 'normalized_score']].values)

    if soft_constraints:
        not_in_model = set(soft_constraints) - set(model.reactions)
        if not_in_model:
            soft_constraints = {r_id: val for r_id, val in soft_constraints.items() if r_id in model.reactions}
            warnings.warn("Soft constraints contain reactions not in the model:\n" + "\n".join(not_in_model))

    if hard_constraints:
        not_in_model = set(hard_constraints) - set(model.reactions)
        if not_in_model:
            hard_constraints = {r_id: (lb, ub) for r_id, (lb, ub) in hard_constraints.items() if r_id in model.reactions}
            warnings.warn("Hard constraints contain reactions not in the model:\n" + "\n".join(not_in_model))

    if ref_model:
        ref_reactions = set(model.reactions) & set(ref_model.reactions)
    else:
        ref_reactions = None

    sol = minmax_reduction(model, scores, default_score=default_score, uptake_score=uptake_score, soft_score=soft_score,
                           soft_constraints=soft_constraints, hard_constraints=hard_constraints,
                           ref_reactions=ref_reactions, ref_score=ref_score, debug_output=debug_output)

    if sol.status == Status.OPTIMAL:
        inactive = inactive_reactions(model, sol)
    else:
        print("MILP solver failed: {}".format(sol.message))
        return

    if debug_output:
        pd.DataFrame.from_dict(sol.values, orient='index').to_csv(debug_output + '_milp_solution.tsv',
                                                                  sep='\t', header=False)

    model.remove_reactions(inactive)

    del_metabolites = disconnected_metabolites(model)
    model.remove_metabolites(del_metabolites)

    for i, row in reaction_scores.iterrows():
        r_id = row['reaction']
        if r_id in model.reactions:
            try:
                gpr = parse_gpr_rule(row['GPR'], prefix='G_')
                model.set_gpr_association(r_id, gpr, add_genes=True)
            except:
                print('Failed to parse:', row['GPR'])

    cleanup_metadata(model)

    if init_env:
        init_env.apply(model, inplace=True, warning=False)

    if outputfile:
        save_cbmodel(model, outputfile, flavor=flavor)

    return model


def build_ensemble(model, reaction_scores, size, outputfile=None, flavor=None, init_env=None):
    """ Reconstruct a model ensemble using the CarveMe approach.

    Args:
        model (CBModel): universal model
        reaction_scores (dict): reaction scores
        size (int): ensemble size
        outputfile (str): write model to SBML file (optional)
        flavor (str): SBML flavor ('cobra' or 'fbc2', optional)
        init_env (Environment): initialize final model with given Environment (optional)

    Returns:
        EnsembleModel: reconstructed ensemble
    """

    scores = dict(reaction_scores[['reaction', 'normalized_score']].values)
    unscored = [r_id for r_id in model.reactions if r_id not in scores and not r_id.startswith('R_EX')]
    logstd = np.std(np.log([x for x in scores.values() if x > 0]))

    reaction_status = {r_id: [] for r_id in model.reactions}
    solver = solver_instance(model)
    failed = 0

    for i in range(size):
        random_scores = -np.exp(logstd * np.random.randn(len(unscored)))
        all_scores = dict(zip(unscored, random_scores))
        all_scores.update(scores)

        sol = minmax_reduction(model, all_scores, solver=solver)

        if sol.status == Status.OPTIMAL:
            for r_id in model.reactions:
                active = (abs(sol.values[r_id]) >= 1e-6
                          or (sol.values.get('yf_' + r_id, 0) > 0.5)
                          or (sol.values.get('yr_' + r_id, 0) > 0.5))
                reaction_status[r_id].append(active)
        else:
            failed += 1

    ensemble_size = size - failed
    ensemble = EnsembleModel(model, ensemble_size, reaction_status)
    ensemble.simplify()

    for i, row in reaction_scores.iterrows():
        r_id = row['reaction']
        if r_id in ensemble.model.reactions:
            gpr = parse_gpr_rule(row['GPR'])
            ensemble.model.reactions[r_id].set_gpr_association(gpr)

    if init_env:
        init_env.apply(ensemble.model, inplace=True, warning=False)

    if outputfile:
        cleanup_metadata(ensemble.model)
        save_ensemble(ensemble, outputfile, flavor=flavor)


def cleanup_metadata(model):

    for met in model.metabolites.values():
        if 'BiGG models' in met.metadata:
            del met.metadata['BiGG models']
        if 'CHEBI' in met.metadata:
            del met.metadata['CHEBI']
        if 'CHARGE' in met.metadata:
            del met.metadata['CHARGE']

    for rxn in model.reactions.values():
        if 'BiGG models' in rxn.metadata:
            del rxn.metadata['BiGG models']
        if 'Reactome' in rxn.metadata:
            del rxn.metadata['Reactome']
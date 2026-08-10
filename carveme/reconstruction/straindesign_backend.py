"""StrainDesign backend for gap-filling.

Gap-filling asks for the cheapest set of reactions to add from a universe so that the model
grows. That is the same problem StrainDesign solves as a *knock-in* design: candidate
reactions are absent until bought, each carries a cost, and a PROTECT module states the
behaviour that must remain feasible.

    reframed / default backend            StrainDesign backend
    ----------------------------------    --------------------------------------
    y_r in {0,1} per candidate            knock-in candidate (inverted binary)
    v_r linked to y_r by bigM (1e3)       indicator constraint
    biomass >= min_growth                 PROTECT module
    min sum 1/(1+score_r) * y_r           ki_cost[r] = 1/(1+score_r)

Two differences are the reason to try it. StrainDesign links binaries with **indicator
constraints** rather than a big-M, which removes the 1e3 bound as a source of numerical
error; and it **compresses** the network first, which the default path does not do at all --
on a model merged with a universe that is where the work is.
"""

def _to_cobra(model, constraints=None):
    """reframed's own converter, plus CarveMe's medium constraints.

    `to_cobrapy` already carries stoichiometry, bounds, GPRs and the objective across, so
    there is no reason to reimplement it -- using theirs also means we track their model
    semantics instead of drifting from them.
    """
    from reframed import to_cobrapy

    # Note: this builds an optlang problem that StrainDesign never touches -- its FVA and
    # compression construct their own MILP_LP straight from the stoichiometry, and it keeps
    # cobra's solver suppressed throughout. On a universe-sized model the conversion costs
    # ~3.5 s for 5532 reactions, so it is worth revisiting if it ever dominates.
    out = to_cobrapy(model)
    for r_id, bounds in (constraints or {}).items():
        if r_id not in out.reactions:
            continue
        rxn = out.reactions.get_by_id(r_id)
        if isinstance(bounds, (tuple, list)):
            rxn.bounds = (float(bounds[0]), float(bounds[1]))
        else:                                        # a single value fixes the flux
            rxn.bounds = (float(bounds), float(bounds))
    return out


def gapfill_straindesign(model, new_reactions, scores=None, min_growth=0.1, constraints=None,
                         solver=None, time_limit=None):
    """Cheapest set of candidates to add so the model grows.

    Args:
        model (CBModel): merged model (original plus universe)
        new_reactions (iterable): candidate reaction ids, absent until bought
        scores (dict): reaction scores; penalty is 1/(1+score), as in the default backend
        min_growth (float): growth the gap-filled model must reach
        constraints (dict): medium and other bound overrides
        solver (str): MILP solver for StrainDesign ('gurobi', 'cplex', 'scip', 'glpk')
        time_limit (float): seconds, optional

    Returns:
        set: the candidates that must be ADDED. Everything else in `new_reactions` is
        inactive and gets removed by the caller.
    """
    from straindesign import SDModule, compute_strain_designs

    scores = scores or {}
    new_reactions = list(new_reactions)
    cobra_model = _to_cobra(model, constraints)

    biomass = model.biomass_reaction
    ki_cost = {r_id: 1.0 / (1.0 + scores.get(r_id, 0.0)) for r_id in new_reactions}

    kwargs = dict(compress=True, solution_approach='best', max_solutions=1,
                  ki_cost=ki_cost)
    if solver:
        kwargs['solver'] = solver
    if time_limit:
        kwargs['time_limit'] = time_limit

    module = SDModule(cobra_model, 'protect', constraints=[f'{biomass} >= {min_growth}'])
    solution = compute_strain_designs(cobra_model, sd_modules=[module], **kwargs)

    designs = solution.get_reaction_sd() if hasattr(solution, 'get_reaction_sd') else solution.reaction_sd
    if not designs:
        raise RuntimeError('StrainDesign found no gap-filling solution '
                           f'(status {getattr(solution, "status", "unknown")})')

    # a knock-in carries a positive value; anything else was not bought
    return {r_id for r_id, v in designs[0].items() if v > 0 and r_id in ki_cost}

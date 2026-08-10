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

from math import inf


def _to_cobra(model, constraints=None):
    """reframed's own converter, plus CarveMe's medium constraints.

    `to_cobrapy` already carries stoichiometry, bounds, GPRs and the objective across, so
    there is no reason to reimplement it -- using theirs also means we track their model
    semantics instead of drifting from them.
    """
    import cobra
    from reframed import to_cobrapy
    from straindesign.networktools import suppress_lp_context, set_suppressed_objective

    # Built inside StrainDesign's suppression context: cobra otherwise updates its optlang
    # problem on every metabolite, reaction and bound change, and nothing downstream reads it
    # -- StrainDesign's FVA and compression build their own MILP_LP from the stoichiometry.
    # Gene rules are dropped for the conversion: gap-filling buys reactions, nothing here ever
    # looks at a gene, and cobra parses every rule into a GPR tree on assignment -- 0.34 s to
    # build them and another 0.25 s every time the model is copied. They are put back on the
    # caller's model immediately, so this is invisible from the outside; blanking them first is
    # what lets us keep using reframed's own converter instead of maintaining a second one.
    gprs = {r_id: rxn.gpr for r_id, rxn in model.reactions.items()}
    with suppress_lp_context(cobra.Model('shell')):
        try:
            for rxn in model.reactions.values():
                rxn.gpr = None
            out = to_cobrapy(model)
        finally:
            for r_id, gpr in gprs.items():
                model.reactions[r_id].gpr = gpr
        # to_cobrapy's `cb_model.objective = ...` writes into the optlang problem, which is
        # suppressed here and has no variables to write to, so the objective would be dropped.
        # Record it where StrainDesign reads it instead.
        set_suppressed_objective(out, model.get_objective())
        for r_id, bounds in (constraints or {}).items():
            if r_id not in out.reactions:
                continue
            rxn = out.reactions.get_by_id(r_id)
            if isinstance(bounds, (tuple, list)):
                # multiGapFill writes bounds like (-max_uptake, None); None means unbounded
                lo, hi = bounds
                rxn.bounds = (-inf if lo is None else float(lo),
                              inf if hi is None else float(hi))
            else:                                    # a single value fixes the flux
                rxn.bounds = (float(bounds), float(bounds))
    return out


def gapfill_straindesign(model, new_reactions, scores=None, min_growth=0.1, constraints=None,
                         solver=None, time_limit=None, compress=False, skip_fvas=True):
    """Cheapest set of candidates to add so the model grows.

    Args:
        model (CBModel): merged model (original plus universe)
        new_reactions (iterable): candidate reaction ids, absent until bought
        scores (dict): reaction scores; penalty is 1/(1+score), as in the default backend
        min_growth (float): growth the gap-filled model must reach
        constraints (dict): medium and other bound overrides
        solver (str): MILP solver for StrainDesign ('gurobi', 'cplex', 'scip', 'glpk')
        time_limit (float): seconds, optional
        compress (bool): run StrainDesign's network compression (default False)
        skip_fvas (bool): skip StrainDesign's preprocessing FVAs (default True) -- they
            narrow the problem for the MILP, which gap-filling does not need

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

    # compression is off by default: on a model merged with a universe it costs far more
    # than it saves -- measured 308 s against 0.9 s for the default backend on a gap-fill
    # that adds nothing at all
    kwargs = dict(compress=compress, solution_approach='best', max_solutions=1,
                  ki_cost=ki_cost, skip_preprocessing_fvas=skip_fvas)
    if solver:
        kwargs['solver'] = solver
    if time_limit:
        kwargs['time_limit'] = time_limit

    # skip_checks: the module's feasibility FBA reads objective coefficients one reaction at
    # a time through cobra's optlang problem, which is the only thing here that needs a live
    # cobra solver. Gap-filling has already established the model is the one we built.
    module = SDModule(cobra_model, 'protect', constraints=[f'{biomass} >= {min_growth}'],
                      skip_checks=True)
    solution = compute_strain_designs(cobra_model, sd_modules=[module], **kwargs)

    designs = solution.get_reaction_sd() if hasattr(solution, 'get_reaction_sd') else solution.reaction_sd
    if not designs:
        raise RuntimeError('StrainDesign found no gap-filling solution '
                           f'(status {getattr(solution, "status", "unknown")})')

    # a knock-in carries a positive value; anything else was not bought
    return {r_id for r_id, v in designs[0].items() if v > 0 and r_id in ki_cost}

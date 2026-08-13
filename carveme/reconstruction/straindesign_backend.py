"""StrainDesign backends for gap-filling and for carving.

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

The second entry point, :func:`complete_model`, replaces carving itself rather than the
gap-filling that follows it. See its docstring for what changes.
"""

import warnings
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


# ---------------------------------------------------------------------------------------------
# Carving, posed as a completion
# ---------------------------------------------------------------------------------------------

SINK_COST = 50.0        # an artificial sink is a hole in the mass balance, not a reaction
SPONT_COST = 0.1        # spontaneous chemistry needs no gene, so the genome cannot argue against it
EXCHANGE_COST = 0.1     # an exchange needs no gene either, but should not be free (see below)
HET_COST = 1.0          # a reaction with no annotation support at all
MIN_SCORE = 0.1         # floor on the reward, so a barely-supported hit is not free
MIN_ATPM = 0.1          # maintenance demand the network must be able to meet
EXTERNAL = 'C_e'


def _spontaneous(gprs, reactions):
    """Reactions whose only 'gene' is a spontaneity placeholder."""
    if gprs is None:
        return set()
    is_spont = gprs.gene.astype(str).str.contains('s0001|SPONT', case=False, regex=True)
    spont, catalysed = set(gprs[is_spont].reaction), set(gprs[~is_spont].reaction)
    return (spont - catalysed) & set(reactions)


def _merge_reverse_duplicates(model, annotated, score_of, gpr_of):
    """Reactions that are the same chemistry written backwards.

    BiGG carries many of these, and left alone they inflate the count of things the completion
    can buy while adding no chemistry. Merging them is only safe when the two do not carry
    *disjoint* enzyme evidence: two genuinely different enzymes running the same conversion in
    opposite directions are two reactions, and collapsing them would erase one gene's support.
    """
    def signature(stoich, flip=False):
        sign = -1 if flip else 1
        return tuple(sorted((m, sign * v) for m, v in stoich.items()))

    def genes_of(r_id):
        import re
        rule = gpr_of.get(r_id) or ''
        return set(re.findall(r'[A-Za-z_][A-Za-z0-9_.\-]*', rule)) - {'and', 'or'}

    forward = {}
    for r_id, rxn in model.reactions.items():
        if len(rxn.stoichiometry) >= 2:
            forward.setdefault(signature(rxn.stoichiometry), []).append(r_id)

    dropped, kept_apart = set(), 0
    for r_id, rxn in model.reactions.items():
        if len(rxn.stoichiometry) < 2 or r_id in dropped:
            continue
        for other in forward.get(signature(rxn.stoichiometry, flip=True), []):
            if other == r_id or other in dropped:
                continue
            genes_a, genes_b = genes_of(r_id), genes_of(other)
            if genes_a and genes_b and not (genes_a <= genes_b or genes_b <= genes_a):
                kept_apart += 1
                continue
            # keep the annotated one, then the reversible one, so no evidence is lost
            keep, gone = sorted((r_id, other),
                                key=lambda x: (x in annotated, model.reactions[x].lb < 0, x),
                                reverse=True)
            dropped.add(gone)
            if gone in annotated:
                annotated.add(keep)
                score_of[keep] = max(score_of.get(keep, 0.0), score_of.get(gone, 0.0))
    return dropped, kept_apart


def complete_model(model, reaction_scores, gprs=None, min_growth=0.1, min_atpm=MIN_ATPM,
                   constraints=None, extra_conditions=None, solver=None, threads=None,
                   time_limit=None, thermodynamic='loopless', verbose=False, compress=None,
                   skip_fvas=None):
    """Which reactions to leave out of the universe, chosen by completion instead of carving.

    Carving picks a threshold on annotation score and deletes below it, then repairs whatever
    broke. Nothing in that formulation ties a reaction's presence to its ability to carry flux,
    so a draft model can contain reactions that can never run and can lose reactions the genome
    plainly supports. Completion states both as constraints instead: an annotated reaction that
    survives must demonstrably carry flux, and a reaction that was not bought carries none. The
    result therefore has no blocked reactions, and the annotated reactions that could not be
    connected at any price are named rather than silently dropped.

    Costs follow CarveMe's own economics, so the two are comparable: an annotated reaction is
    rewarded by its normalized score, an unannotated one costs 1. Two exceptions, both measured:
    spontaneous reactions cost almost nothing, because no genome can be evidence against
    chemistry that needs no enzyme; and artificial sinks cost 50, because a sink is a hole in the
    mass balance that lets the network meet any demand without producing it. Pricing sinks like
    ordinary reactions is what lets a completion satisfy biomass through a shortcut.

    Args:
        model (CBModel): universal model
        reaction_scores (pandas.DataFrame): reaction scores, as `reaction_scoring` returns them
        gprs (pandas.DataFrame): the BiGG GPR table, used to identify spontaneous reactions
        min_growth (float): growth the completed model must reach
        min_atpm (float): maintenance flux the completed model must be able to carry
        constraints (dict): medium and other bound overrides
        extra_conditions (list of dict): further media the model must also grow on. Each becomes
            a PROTECT module over the same binaries, so every condition gets its own flux system
            inside one MILP and a reaction shared by two of them is paid for once -- which
            iterated single-condition gap-filling does not achieve.
        solver (str): MILP solver ('gurobi', 'cplex', 'scip', 'glpk')
        threads (int), time_limit (float), verbose (bool)
        thermodynamic (str): 'loopless' (default) forbids core reactions from satisfying their
            must-run condition inside a thermodynamically infeasible cycle; None omits it.
        compress, skip_fvas: StrainDesign pipeline options. Both default to off for a CarveMe
            module -- on a universe there is little to compress and most of what an FVA scans
            will not be bought. Note that coupled lumping sums knock-in costs, so
            compress='coupled' is the harsh option here, not the gentle one.

    Returns:
        set: reaction ids to REMOVE from the universe.
    """
    import straindesign as sd
    from straindesign.names import CARVEME, PROTECT, BEST

    score_of = dict(zip(reaction_scores['reaction'], reaction_scores['normalized_score']))
    gpr_of = dict(zip(reaction_scores['reaction'], reaction_scores.get('GPR', '')))
    annotated = set(reaction_scores['reaction']) & set(model.reactions)

    duplicates, kept_apart = _merge_reverse_duplicates(model, annotated, score_of, gpr_of)
    spontaneous = _spontaneous(gprs, model.reactions)

    cobra_model = _to_cobra(model, constraints)
    biomass = model.biomass_reaction
    atpm = next((r for r in ('R_ATPM', 'ATPM') if r in cobra_model.reactions), None)

    # A merged duplicate is not something the completion may buy -- it is not in the universe.
    # Its bounds are closed rather than the reaction deleted: this cobra model is built with a
    # suppressed optlang problem, which cobra's remove_reactions reaches into and ours has no
    # variables for. Excluded from the candidates below, it is inert either way.
    from straindesign.networktools import suppress_lp_context
    with suppress_lp_context(cobra_model):
        for r_id in duplicates:
            if r_id in cobra_model.reactions:
                cobra_model.reactions.get_by_id(r_id).bounds = (0.0, 0.0)

    exchange, sinks = set(), set()
    for rxn in cobra_model.reactions:
        if len(rxn.metabolites) != 1:
            continue
        (exchange if list(rxn.metabolites)[0].compartment == EXTERNAL else sinks).add(rxn.id)

    # Exchanges are candidates too, at a small price. Left out of the candidate set they can never
    # be dropped, and the result then carries a transporter for every metabolite in the universe --
    # 645 of them, nearly all unusable by the chemistry that was actually bought, and every one
    # counted as blocked. Carving drops them for the same reason, so leaving them in would also
    # make the two backends incomparable. The price has to be positive rather than zero: at zero
    # an unusable exchange costs nothing to buy and the choice is arbitrary.
    fixed = duplicates | {biomass} | ({atpm} if atpm else set())
    candidates = [r.id for r in cobra_model.reactions if r.id not in fixed]

    def cost_of(r):
        if r in sinks:
            return SINK_COST
        if r in annotated:
            return -max(score_of.get(r, 1.0), MIN_SCORE)
        if r in exchange:
            return EXCHANGE_COST
        if r in spontaneous:
            return SPONT_COST
        return HET_COST

    ki_cost = {r: cost_of(r) for r in candidates}
    core = [r for r in candidates if r in annotated]

    demands = ['%s >= %g' % (biomass, min_growth)]
    if atpm:
        demands.append('%s >= %g' % (atpm, min_atpm))

    if verbose:
        print(f'Completion: {len(candidates)} candidates, {len(core)} annotated '
              f'({len(duplicates)} reverse duplicates merged, {kept_apart} kept apart on '
              f'disjoint enzymes)')

    module = sd.SDModule(cobra_model, CARVEME, constraints=demands, core_reactions=core,
                         thermodynamic=thermodynamic, skip_checks=True)
    # compression is off by default: on a universe it costs far more than it saves. It is a plain
    # pipeline option here, not something this module type bypasses -- pass compress=True to use it.
    kwargs = dict(ki_cost=ki_cost, solution_approach=BEST, max_solutions=1)
    if compress is not None:          # otherwise StrainDesign's CarveMe defaults apply
        kwargs['compress'] = compress
    if skip_fvas is not None:
        kwargs['skip_preprocessing_fvas'] = skip_fvas
    if solver:
        kwargs['solver'] = solver
    if threads:
        kwargs['milp_threads'] = threads
    if time_limit:
        kwargs['time_limit'] = time_limit
    modules = [module]
    for condition in (extra_conditions or []):
        modules.append(sd.SDModule(cobra_model, PROTECT, skip_checks=True,
                                   constraints=demands + _bounds_to_constraints(cobra_model,
                                                                                condition)))

    solution = sd.compute_strain_designs(cobra_model, sd_modules=modules, **kwargs)
    designs = solution.reaction_sd
    if not designs:
        raise RuntimeError('StrainDesign found no completion '
                           f'(status {getattr(solution, "status", "unknown")})')

    status = getattr(solution, 'status', None)
    if status != 'optimal':
        warnings.warn(f'The completion did not prove optimality (status {status}). The model is '
                      f'valid -- every constraint holds -- but a cheaper one may exist. Raise '
                      f'time_limit to close the gap.')

    # read against the candidate list, not the design's keys: a candidate whose binary the MILP
    # fixed to zero is reported in no design at all
    bought = {r for r in candidates if designs[0].get(r)}
    inactive = (set(candidates) - bought) | duplicates
    if verbose:
        print(f'Completion kept {len(bought & annotated)}/{len(core)} annotated and bought '
              f'{len(bought - annotated)} unannotated reactions')
    return inactive


def _bounds_to_constraints(cobra_model, bounds):
    """CarveMe writes media as {reaction: (lb, ub)}; a constraint block wants expressions."""
    out = []
    for r_id, bnd in (bounds or {}).items():
        if r_id not in cobra_model.reactions:
            continue
        lo, hi = bnd if isinstance(bnd, (tuple, list)) else (bnd, bnd)
        if lo is not None:
            out.append('%s >= %g' % (r_id, lo))
        if hi is not None:
            out.append('%s <= %g' % (r_id, hi))
    return out

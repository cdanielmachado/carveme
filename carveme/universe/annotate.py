from reframed.core.elements import parse_formula, unparse_formula
from reframed import solver_instance
from reframed.solvers.solver import VarType


def determine_elem_composition(model, elem, known, unknown, incomplete):
    solver = solver_instance()

    reactions = []
    objective = {}
    all_unknown = unknown | set(incomplete)

    for m_id in unknown:
        solver.add_variable(m_id, lb=0, vartype=VarType.INTEGER, update=False)

    for m_id, formula in incomplete.items():
        solver.add_variable(m_id, lb=formula.get(elem, 0), update=False)

    for r_id, rxn in model.reactions.items():

        if len(rxn.stoichiometry) == 1:
            continue  # exclude sink and exchange reactions

        if len(set(rxn.stoichiometry) & all_unknown) == 0:
            continue  # no unknowns for this reaction

        reactions.append(r_id)
        solver.add_variable(f'{r_id}_e+', lb=0, update=False)
        solver.add_variable(f'{r_id}_e-', lb=0, update=False)
        objective[f'{r_id}_e+'] = 1
        objective[f'{r_id}_e-'] = 1

    solver.update()

    for r_id in reactions:
        stoichiometry = model.reactions[r_id].stoichiometry

        balance = sum(coeff * known[m_id].get(elem, 0)
                      for m_id, coeff in stoichiometry.items() if m_id in known)

        lhs = {m_id: coeff for m_id, coeff in stoichiometry.items()
               if m_id in all_unknown}

        lhs[f'{r_id}_e+'] = 1
        lhs[f'{r_id}_e-'] = -1

        solver.add_constraint(r_id, lhs, '=', -balance, update=False)

    solver.update()
    sol = solver.solve(objective, minimize=True)

    solver.add_constraint('step1', objective, '<', sol.fobj * (1 + 1e-6))
    objective2 = {m_id: 1 for m_id in all_unknown}
    sol2 = solver.solve(objective2, minimize=True)

    return sol2


def compute_missing_formulae(model):

    print('Computing missing formulae...')

    known = {}
    unknown = set()
    incomplete = {}

    for m_id, met in model.metabolites.items():

        if 'FORMULA' in met.metadata:
            formula = met.metadata['FORMULA']
            if '*' in formula:
                formula = formula.replace('*', '')
                incomplete[m_id] = parse_formula(formula)
            else:
                known[m_id] = parse_formula(formula)
        else:
            unknown.add(m_id)

    elements = {elem for formula in known.values() for elem in formula}

    solutions = {}
    for elem in elements:
        solutions[elem] = determine_elem_composition(model, elem, known, unknown, incomplete)

    for m_id in unknown | set(incomplete):
        elems_dict = {elem: solutions[elem].values[m_id] for elem in elements}
        formula = unparse_formula(elems_dict)
        model.metabolites[m_id].metadata['FORMULA'] = formula


def determine_proton_composition(model, alpha=0.5):
    solver = solver_instance()

    reactions = [r_id for r_id, rxn in model.reactions.items() if len(rxn.stoichiometry) > 1]

    cpds_nH = {m_id[2:-2]: parse_formula(met.metadata['FORMULA']).get('H', 0)
               for m_id, met in model.metabolites.items()}

    objective = {}

    for cpd, nH in cpds_nH.items():
        solver.add_variable(f'{cpd}_d+', lb=0, vartype=VarType.INTEGER, update=False)
        solver.add_variable(f'{cpd}_d-', lb=0, ub=nH, vartype=VarType.INTEGER, update=False)
        objective[f'{cpd}_d+'] = alpha / len(cpds_nH)
        objective[f'{cpd}_d-'] = alpha / len(cpds_nH)

    for r_id in reactions:
        solver.add_variable(f'{r_id}_e+', lb=0, vartype=VarType.INTEGER, update=False)
        solver.add_variable(f'{r_id}_e-', lb=0, vartype=VarType.INTEGER, update=False)
        objective[f'{r_id}_e+'] = (1 - alpha) / len(reactions)
        objective[f'{r_id}_e-'] = (1 - alpha) / len(reactions)

    solver.update()

    for r_id in reactions:
        stoich = model.reactions[r_id].stoichiometry
        lhs = {f'{r_id}_e+': 1, f'{r_id}_e-': -1}
        lhs.update({f'{m_id[2:-2]}_d+': coeff for m_id, coeff in stoich.items()})
        lhs.update({f'{m_id[2:-2]}_d-': -coeff for m_id, coeff in stoich.items()})
        rhs = -sum(cpds_nH[m_id[2:-2]] * coeff for m_id, coeff in stoich.items())
        solver.add_constraint(r_id, lhs, '=', rhs, update=False)

    solver.update()

    sol = solver.solve(objective, minimize=True)

    return sol


def fix_protons_and_charge(model):

    print('Trying to correct proton and charge balance...')

    sol = determine_proton_composition(model, alpha=0.1)

    for m_id, met in model.metabolites.items():
        diff = sol.values[f'{m_id[2:-2]}_d+'] - sol.values[f'{m_id[2:-2]}_d-']

        if abs(diff) > 0.5:
            coeffs = parse_formula(met.metadata['FORMULA'])
            coeffs['H'] = coeffs.get('H', 0) + int(diff)
            met.metadata['FORMULA'] = unparse_formula(coeffs)

            charge = int(met.metadata.get('CHARGE', 0)) + diff
            met.metadata['CHARGE'] = str(charge)


def fix_hydrogen_stoichiometry(model):

    print('Trying to fix hydrogen stoichiometry...')

    nH = {m_id: parse_formula(met.metadata['FORMULA']).get('H', 0)
          for m_id, met in model.metabolites.items()}

    for r_id, rxn in model.reactions.items():

        if len(rxn.stoichiometry) == 1:
            continue

        balance = int(sum(nH[m_id] * coeff for m_id, coeff in rxn.stoichiometry.items()))

        if abs(balance) == 0:
            continue

        comps = model.get_reaction_compartments(r_id)
        if len(comps) > 1:
            continue

        proton = 'M_h_' + list(comps)[0][-1]

        rxn.stoichiometry[proton] = rxn.stoichiometry.get(proton, 0) - balance
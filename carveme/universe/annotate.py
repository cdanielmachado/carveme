from reframed.core.elements import parse_formula
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

    print('Computing missing formulae in model...')

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
        formula = ''
        for elem in sorted(elements):
            coeff = solutions[elem].values[m_id]
            if coeff == 0:
                continue
            elif coeff == 1:
                formula += elem
            else:
                formula += elem + str(int(coeff))

        if formula != '':
            model.metabolites[m_id].metadata['FORMULA'] = formula


def annotate_charges(model):
    # procedure similar to the one above for charges is too computationally expensive
    # will simply create a default value to annotate the SBML file

    for m_id, met in model.metabolites.items():
        if 'CHARGE' not in met.metadata:
            met.metadata['CHARGE'] = '0'

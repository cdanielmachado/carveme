from reframed import FBA, solver_instance, Environment, essential_genes
from reframed.solvers.solution import Status


def benchmark_biolog(model, medium, data, min_growth=0.1, max_uptake=10):

    env = Environment.from_compounds(medium)
    constraints = env.apply(model, inplace=False, warning=False)
    solver = solver_instance(model)
    data = data[["bigg_id", "growth"]].dropna()
    result = {}

    for _, row in data.iterrows():
        met = row["bigg_id"]
        in_vivo_growth = row["growth"] in {'++', '+'}
        r_id = f"R_EX_{met}_e"

        if r_id in model.reactions:
            tmp = constraints[r_id] if r_id in constraints else (0, 0)
            constraints[r_id] = (-max_uptake, 0)
            sol = FBA(model, constraints=constraints, solver=solver)
            in_silico_growth = Status.OPTIMAL and sol.fobj > min_growth
            constraints[r_id] = tmp
        else:
            in_silico_growth = False

        if in_silico_growth:
            result[met] = 'TP' if in_vivo_growth else 'FP'
        else:
            result[met] = 'FN' if in_vivo_growth else 'TN'

    return result


def benchmark_essentiality(model, medium, in_vivo):

    if medium is not None:
        env = Environment.from_compounds(medium)
    else:
        env = Environment.complete(model)
        
    constraints = env.apply(model, inplace=False, warning=False)
    in_silico = essential_genes(model, constraints=constraints, min_growth=0.1)

    result = {}
    for gene, is_essential in in_vivo.items():
        if is_essential:
            if gene in in_silico:
                result[gene] = 'TP'
            else:
                result[gene] = 'FN'
        else:
            if gene in in_silico:
                result[gene] = 'FP'
            else:
                result[gene] = 'TN'
    
    return result


def mcc(data):
    data = list(data['value'])
    tp = data.count('TP')
    fp = data.count('FP')
    fn = data.count('FN')
    tn = data.count('TN')

    num = (tp * tn - fp * fn)
    dom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    return num / dom


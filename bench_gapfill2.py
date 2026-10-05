"""Gap-filling with real gaps: knock reactions out of a working model, then restore growth.

The previous run was vacuous -- the model already grew, so nothing was added and the MILP was
trivial. Removing k reactions that the model needs creates a controlled difficulty knob and
forces the solver to actually search.
"""
import warnings, logging, time, sys, random
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
from reframed import load_cbmodel, FBA
from reframed.core.environment import Environment
from carveme.reconstruction.gapfilling import gapFill

ORG = sys.argv[1] if len(sys.argv) > 1 else 'ecol'
KS = [int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else '5,20').split(',')]
UNI = 'carveme/data/generated/universe_bacteria.xml.gz'
COMPOUNDS = ['glc__D', 'o2', 'nh4', 'pi', 'so4', 'h2o', 'h', 'k', 'mg2', 'fe2', 'ca2', 'cl', 'na1']
CONSTR = dict(Environment.from_compounds(COMPOUNDS, max_uptake=10))

def load(): return load_cbmodel(f'carveme/data/benchmark/models/{ORG}.xml', flavor='bigg')

base = load()
sol = FBA(base, constraints=CONSTR)
print(f'MODEL {ORG}: {len(base.reactions)} rxns, wild-type growth={sol.fobj:.4f}', flush=True)

# pick reactions whose removal actually kills growth, so the gap is real
rng = random.Random(1)
carrying = [r for r, v in sol.values.items() if abs(v) > 1e-9 and not r.startswith('R_EX')]
rng.shuffle(carrying)

for k in KS:
    removed = carrying[:k]
    for backend in ('reframed', 'straindesign'):
        model = load()
        model.remove_reactions(removed)
        g = FBA(model, constraints=CONSTR)
        pre = g.fobj if g.fobj else 0.0
        uni = load_cbmodel(UNI, flavor='bigg')
        t0 = time.time()
        try:
            out = gapFill(model, uni, constraints=CONSTR, min_growth=0.05,
                          inplace=False, backend=backend)
            wall = time.time() - t0
            added = sorted(set(out.reactions) - (set(load().reactions) - set(removed)))
            print(f'GAP k={k:3d} {backend:13s} pre_growth={pre:.4f} wall={wall:8.1f}s '
                  f'added={len(added)}', flush=True)
        except Exception as e:
            print(f'GAP k={k:3d} {backend:13s} FAILED {time.time()-t0:.1f}s '
                  f'{type(e).__name__}: {str(e)[:120]}', flush=True)

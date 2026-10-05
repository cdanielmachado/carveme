"""Real-scale gap-fill: a CarveMe model against the bacterial universe, both backends.

The toy case measured process warm-up, not the algorithm. Here the merged model carries the
full universe, which is the regime where compression could matter -- or not, if the MILP is
so easy that no backend helps.
"""
import warnings, logging, time, sys
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
from reframed import load_cbmodel
from reframed.core.environment import Environment
from carveme.reconstruction.gapfilling import gapFill

ORG = sys.argv[1] if len(sys.argv) > 1 else 'ecol'
UNI = 'carveme/data/generated/universe_bacteria.xml.gz'

def fresh():
    m = load_cbmodel(f'carveme/data/benchmark/models/{ORG}.xml', flavor='bigg')
    u = load_cbmodel(UNI, flavor='bigg')
    return m, u

m0, u0 = fresh()
print(f'MODEL {ORG}: {len(m0.reactions)} rxns | universe {len(u0.reactions)} rxns | '
      f'candidates {len(set(u0.reactions) - set(m0.reactions))}', flush=True)

# a minimal-ish medium: glucose aerobic, plus the model's own exchange set closed
compounds = ['glc__D', 'o2', 'nh4', 'pi', 'so4', 'h2o', 'h', 'k', 'mg2', 'fe2', 'ca2', 'cl', 'na1']
results = {}
for backend in ('reframed', 'straindesign'):
    model, uni = fresh()
    env = Environment.from_compounds(compounds, max_uptake=10)
    merged_probe = model                      # constraints are applied to the merged model
    constraints = {r: b for r, b in dict(env).items()}
    t0 = time.time()
    try:
        out = gapFill(model, uni, constraints=constraints, min_growth=0.05,
                      inplace=False, backend=backend)
        wall = time.time() - t0
        added = sorted(set(out.reactions) - set(m0.reactions))
        results[backend] = (wall, added)
        print(f'GF {backend:13s} wall={wall:8.1f}s added={len(added)} reactions', flush=True)
    except Exception as e:
        print(f'GF {backend:13s} FAILED after {time.time()-t0:.1f}s: {type(e).__name__}: {str(e)[:160]}', flush=True)

if len(results) == 2:
    a, b = results['reframed'][1], results['straindesign'][1]
    print(f'CMP identical={a == b} | reframed {len(a)} vs SD {len(b)} | '
          f'only-reframed {len(set(a)-set(b))} only-SD {len(set(b)-set(a))}', flush=True)
    print(f'CMP speed: SD x{results["straindesign"][0]/max(0.01,results["reframed"][0]):.2f} vs reframed', flush=True)

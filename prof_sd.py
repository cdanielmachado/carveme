import warnings, logging, time
warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(message)s', force=True)
from reframed import load_cbmodel
from reframed.core.environment import Environment
from carveme.reconstruction.gapfilling import merge_models
from carveme.reconstruction.straindesign_backend import gapfill_straindesign

m = load_cbmodel('carveme/data/benchmark/models/ecol.xml', flavor='bigg')
u = load_cbmodel('carveme/data/generated/universe_bacteria.xml.gz', flavor='bigg')
new = set(u.reactions) - set(m.reactions)
merged = merge_models(m, u, inplace=False)
constr = dict(Environment.from_compounds(
    ['glc__D','o2','nh4','pi','so4','h2o','h','k','mg2','fe2','ca2','cl','na1'], max_uptake=10))
t0 = time.time()
added = gapfill_straindesign(merged, new, scores={}, min_growth=0.05, constraints=constr, compress=False)
print(f'TOTAL {time.time()-t0:.1f}s added={len(added)}', flush=True)

from equilibrator_api import ComponentContribution, Reaction
from equilibrator_cache.exceptions import MissingDissociationConstantsException
from reframed import ReactionType


def annotate_with_eQuilibrator(model):
    cc = ComponentContribution()

    cc_cpds = {}

    for m_id, met in model.metabolites.items():
        cpd = m_id[2:-2]

        if cpd not in cc_cpds:
            cc_cpds[cpd] = cc.get_compound("bigg.metabolite:" + cpd)

        cc_cpd = cc_cpds[cpd]

        if cc_cpd is None:
            continue

        dGf = cc.standard_dg_formation(cc_cpd)

        if dGf[0] is not None:
            met.metadata['dGf'] = f'{dGf[0]:.6f}'

    for r_id, rxn in model.reactions.items():
        if rxn.reaction_type in {ReactionType.TRANSPORT, ReactionType.EXCHANGE}:
            continue

        stoichiometry = {cc_cpds[m_id[2:-2]]: coeff for m_id, coeff in rxn.stoichiometry.items()}

        if None in stoichiometry:
            continue

        cc_rxn = Reaction(sparse=stoichiometry, rid=r_id)

        try:
            dG0r = cc.standard_dg_prime(cc_rxn)
        except MissingDissociationConstantsException:
            continue

        if dG0r.s > 100:  # exclude erroneous estimations
            continue

        rxn.metadata['dG0r'] = f'{dG0r.n:.6f}'
        rxn.metadata['dG0r_std'] = f'{dG0r.s:.6f}'


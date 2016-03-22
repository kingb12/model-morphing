import Plotter
import Helpers
import Model, Reaction
import math
from Morph import *
def reaction_analysis(model_identities):
    models = list()
    for m in model_identities:
        model = Helpers.get_object(m[0], m[1])
        models.append((m[2], frozenset([Helpers.get_rxn_id(r) for r in Model.get_reactions(model['data'])])))
    all_models = set(models)
    rxn_sets = dict()
    for model_set in powerset(models):
        if len(model_set) > 0:
            rxn_set = model_set[0][1]
            keys = set(model_set)
            skips = all_models - keys
            key = set()
            for m in model_set:
                key.add(m[0])
            for m in model_set:
                rxn_set = rxn_set.intersection(m[1])
            for m in skips:
                rxn_set = rxn_set - m[1]
            rxn_sets[frozenset(key)] = rxn_set
    return rxn_sets
def gene_percents(model, rxn_analysis):
    #TODO: Find a creative way to avoid explicitly calling the model
    result = dict()
    gprs = dict()
    for r in model['modelreactions']:
        gprs[Reaction.get_rxn_id(r)] = Gpr(r)
    for group in rxn_analysis:
        i = 0
        for r in rxn_analysis[group]:
            try:
                if gprs[r].gpr_type == 'no-gene':
                    i += 1
            except KeyError:
                i += 1
        if len(rxn_analysis[group]) > 0:
            percentage = float(len(rxn_analysis[group]) - i) / float(len(rxn_analysis[group])) * 100.0
        else:
            percentage = 0.0
        result[group] = percentage
    return result
def powerset(lst):
    result = [[]]
    for x in lst:
        result.extend([subset + [x] for subset in result])
    return result
def three_model_venn(model_identities, filename, genes=False, morphed_model=None, title=None):
    rxn_sets=[]
    for m in model_identities:
        rxn_sets.append(set([Reaction.get_rxn_id(r) for r in Model.get_reactions(Helpers.get_object(m[0], m[1])['data'])]))
    gene_dict = None
    if genes:
        rxn_analysis = reaction_analysis(model_identities)
        genes = gene_percents(morphed_model, rxn_analysis)
        gene_dict = dict()
        for model_set in genes:
            key = 0
            for i in range(0, 3):
                if model_identities[i][2] in model_set:
                   key += int(math.pow(10, i))
            key = str(key)
            while len(key) < 3:
                key = '0' + key
            value = '{0:.3g}'.format(genes[model_set])
            gene_dict[key] = str(value) + '%'



    Plotter.venn3(rxn_sets, title, 'Reaction', filename, set_labels=[m[2] for m in model_identities], annotation=gene_dict)

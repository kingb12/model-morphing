from operator import attrgetter

import Plotter
import Helpers
import math
from objects import *
from Morph import *
def common_reaction_analysis(model_identities):
    models = list()
    for m in model_identities:
        model = FBAModel(m[0], m[1])
        models.append((m[2], frozenset([r.rxn_id() for r in model.get_reactions()])))
    return venn_analysis(models)

def venn_analysis(list_of_set_pairs):
    """
    Performs a venn analysis on a list of sets (as tuples with labels at tup[0])

    Imagine this as a venn diagram, but any dimensions are possible. Returns a dictionary with the elements belonging
    to each section of the supposed diagram and a label identifying which sets contribute to this area.

    :param list_of_set_pairs: A list of tuples, where each tuple has a set object at tup[1] and a name for the frozenset
     at tup[0]
    :return: a dictionary with a frozen set of names as a key to the elements from the sets in these keys and not others
    """
    all_sets = set(list_of_set_pairs)
    rxn_sets = dict()
    # for each subset of the power set of the sets passed as arguments...
    for powsub_set in powerset(list_of_set_pairs):
        if len(powsub_set) > 0:
            rxn_set = powsub_set[0][1]
            keys = set(powsub_set)
            skips = all_sets - keys
            key = set()
            for m in powsub_set:
                key.add(m[0])
            for m in powsub_set:
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
        gprs[ModelReaction(r).rxn_id()] = Gpr(r)
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
    rxn_sets = []
    for m in model_identities:
        rxn_sets.append(set([r.rxn_id() for r in FBAModel(m[0], m[1]).get_reactions()]))
    gene_dict = None
    if genes:
        rxn_analysis = common_reaction_analysis(model_identities)
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

def reaction_sets(list_of_tuples):
    result = dict()
    for m in list_of_tuples:
        model = FBAModel(m[0], m[1])
        result[m[2]] = set([r.rxn_id() for r in model.get_reactions()])


def pairwise_venn_analysis(model_indentities, comparison_tup):
    result = dict()
    for m in model_indentities:
        result[(comparison_tup[2], m[2])] = reaction_analysis([comparison_tup, m])
    return result

def reaction_analysis(reactions, model, media, rxn_labels):
    # reaction, genes?, label?, subsystem?
    # header = ('Reaction', 'genes in morph?', '# of genes', 'morph labels', 'class', 'sublclass', 'subsystem')
    header = ('Reaction', 'genes in morph?', '# of genes', 'morph labels', 'formula')
    data = Plotter.SimpleTable(header)
    model_rxns = dict([(r.rxn_id(), r) for r in model.get_reactions()])
    model_info, reaction_info = service.model_info(model)
    for r in reactions:
        reaction = model_rxns[r]
        genes = len(reaction.gpr.features()) > 0
        num_genes = len(reaction.gpr.features())
        labels = set()
        for key in rxn_labels:
            if r in rxn_labels[key]:
                labels.add(key)
        label_str = ''
        if len(labels) > 1:
            for l in labels:
                label_str += l + ', '
        else:
            label_str = list(labels)[0]
        rxn_class_str = None
        subclass_str = None
        subsystem_str = None
        try:
            info = reaction_info[r.split('_')[0]]
            rxn_class = set(info['primclass'].split(';'))
            rxn_class_str = ''
            for r1 in rxn_class:
                rxn_class_str += str(r1) + ', '
            subclass = set(info['subclass'].split(';'))
            subclass_str = ''
            for r1 in subclass:
                subclass_str += str(r1) + ', '
            subsystem = set(info['subsytem'].split(';'))  # MISSPELLED ON PURPOSE, KBase's key is spelled this way
            subsystem_str = ''
            for r1 in subsystem:
                subsystem_str += str(r1) + ', '
        except KeyError:
            pass
        formula = reaction.get_equation()
        formula.sort(key=attrgetter('coeff'))
        print formula
        formula_str = ''
        for c in formula:
            formula_str += str(c.coeff) + '\*' + str(c.name()) + ', '

       # data.add((r, genes, num_genes, label_str, rxn_class_str, subclass_str, subsystem_str))
        data.add((r, genes, num_genes, label_str, formula_str))

    return data

def fva_change_analysis(models, media, workspace=None):
    """
    runs fva on all models and returns a dictionary od reactions which are blocked in one, the other, or both
    :param models:
    :return:
    """
    if workspace is None:
        workspace = models[0].workspace_id
    blocked_sets = []
    for m in models:
        info = service.runfva(m, media, workspace=workspace)
        fva = FBA(info[0], info[1])
        blocked_sets.append((m.name, frozenset(fva.blocked_reactions())))
    return venn_analysis(blocked_sets)





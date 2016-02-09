import Helpers
import Client

def gpr_set(reaction):
    '''
    deprecate - gpr is now a type with its own methods
    '''
    rxn_proteins = reaction['modelreactionproteins']
    prots = set()
    for i in range(0, len(rxn_proteins)):
        prot = set()
        subunits = rxn_proteins[i]['modelreactionproteinsubunits']
        for j in range(0, len(subunits)):
            unit = subunits[j]
            ftrs = [f.split('/')[-1] for f in unit['feature_refs']]
            if len(ftrs) > 0:
                prot.add(frozenset(ftrs))
        if len(prot) > 0:
            prots.add(frozenset(prot))
    if len(prots) > 0:
        return frozenset(prots)
    return none

def ftr_set(reaction):
    features = set()
    for protein in reaction['modelreactionproteins']:
        for sub in protein['modelreactionproteinsubunits']:
            for f in [i.split('/')[-1] for i in sub['feature_refs']]:
                features.add(f)
    return features
def get_rxn_id(modelreaction_obj):
    '''
    returns the reaction id of a model reaction object (i.e. rxn34565_c0)
    '''
    #rxn_id = modelreaction_obj['id']
    rxn_id = modelreaction_obj['reaction_ref'].split('/')[-1] + '_' + modelreaction_obj['modelcompartment_ref'].split('/')[-1]
    if rxn_id.startswith('rxn00000'):
        return modelreaction_obj['id']
    return rxn_id
def get_equation(rxn_object):
    '''
    takes a modelreaction object (type dict) and returns it's equation
    '''
    eq = ''
    for cpd in rxn_object['modelReactionReagents']:
        if cpd['coefficient'] < 0:
            if cpd['coefficient'] != -1:
                eq += '(' + str(-1*cpd['coefficient']) + ')'
            eq += cpd['modelcompound_ref'].split('/')[-1] + ' + '
    if eq.endswith('+ '):
        eq = eq[:-2]
    eq += rxn_object['direction'] + ' '
    for cpd in rxn_object['modelReactionReagents']:
        if cpd['coefficient'] > 0:
            if cpd['coefficient'] != 1:
                eq += '(' + str(cpd['coefficient']) + ')'
            eq += cpd['modelcompound_ref'].split('/')[-1] + ' + '
    if eq.endswith('+ '):
        eq = eq[:-2]
    return eq

def get_direction(reaction_obj):
    return reaction_obj['direction']
def get_comp_ref(reaction_obj):
    return str(reaction_obj['modelcompartment_ref'].split('/')[-1][0])
def get_rxn_ref(reaction_obj):
    '''
    gets the reaction object reference if it is associated with a KBase Biochemistry
    '''
    return reaction_obj['reaction_ref'].split('/')[-1]

def get_removal_id(reaction_obj):
    return reaction_obj['id']

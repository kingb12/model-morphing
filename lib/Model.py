import Helpers
import Morph
import Reaction
import re

'''
A Module for performing functions on Model Objects, as they are represented in KBase or otherwise
in the morphing algorithm
'''

def get_reactions(model_obj):
    return [r for r in model_obj['modelreactions']]
def find_reaction_id_disagreements(model_obj):
    '''
    finds reactions in the model that seem to represent KBase reactions but don't truly
    '''
    result = list()
    reactions = get_reactions(model_obj)
    for r in reactions:
        ref = r['reaction_ref'].split('/')[-1]
        if ref == 'rxn00000':
            if re.search('rxn[0-9]{5}', r['id']) is not None:
                # you used a KBase style ID but KBase didnt recognize it
                result.append(r)
    return Reaction.find_biochem_disagreements(result)

def find_compound_uses(model_obj, compound_id):
    '''
    finds all reactions which have compound_id as a reagent

    returns the reaction_obj's as a list
    '''
    reactions = get_reactions(model_obj)
    result = list()
    for r in reactions:
        if Reaction.has_compound(r, compound_id):
            result.append(r)
    return result

def find_disagreement_impact(model_obj, disagreements):
    '''
    takes the return list from find_reaction_id_disagreements as an argument...

    returns a dictionary of the prevalence of each compound disagreement in the model,
    '''
    result = dict()
    for r in disagreements:
        d = disagreements[r]
        if d['type'] == 'compounds' and d['guess'] == 'replacement':
            me = d['difference']['me']
#            kb = d['difference']['kb']
            for cpd in me:
                info = {'contained' : True, 'reactions': dict(), 'uncontainedReactions': dict()}
                uses = find_compound_uses(model_obj, cpd)
                for rxn in uses:
                    rxn_id = Reaction.get_rxn_id(rxn)
                    if rxn_id not in disagreements:
                        info['contained'] = False
                        info['uncontainedReactions'][rxn_id] = rxn
                    info['reactions'][rxn_id] = rxn
            result[cpd] = info
    return result
def has_reaction(model_obj, reaction_id): # NOT PRAGMATIC
    reactions = [Reaction.get_rxn_id(r) for r in get_reactions(model_obj)]
    return reaction_id in reactions





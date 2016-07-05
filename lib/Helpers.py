import service
from Morph import Morph
import json
import pickle
import service
import random
import os
from lib.objects import *


def get_object(objid, wsid, name=None):
    '''
    Returns an object and it's associated KBase information

    Returns an ObjectData (dictionary) like what is returned in the workspace service 'get_objects' function:

	/* The data and supplemental info for an object.

		UnspecifiedObject data - the object's data or subset data.
		object_info info - information about the object.
		list<ProvenanceAction> provenance - the object's provenance.
		username creator - the user that first saved the object to the
			workspace.
		timestamp created - the date the object was first saved to the
			workspace.
		list<obj_ref> - the references contained within the object.
		obj_ref copied - the reference of the source object if this object is
			a copy and the copy source exists and is accessible.
			null otherwise.
		boolean copy_source_inaccessible - true if the object was copied from
			another object, but that object is no longer accessible to the
			user. False otherwise.
		mapping<id_type, list<extracted_id>> extracted_ids - any ids extracted
			from the object.
		string handle_error - if an error occurs while setting ACLs on
			embedded handle IDs, it will be reported here.
		string handle_stacktrace - the stacktrace for handle_error.

	*/
	typedef structure {
		UnspecifiedObject data;
		object_info info;
		list<ProvenanceAction> provenance;
		username creator;
		timestamp created;
		list<obj_ref> refs;
		obj_ref copied;
		boolean copy_source_inaccessible;
		mapping<id_type, list<extracted_id>> extracted_ids;
		string handle_error;
		string handle_stacktrace;
	} ObjectData;

    '''
    obj = service.get_object(objid, wsid, name=name)
    return {'data': obj[0], 'info': obj[1]}

def dump(object, filepath):
    with open(filepath, 'w') as f:
        pickle.dump(object, f)

def load(filepath):
    with open(filepath, 'rb') as f:
        return pickle.load(f)


def save_object(data, type, wsid, objid=None, name=None):
    return service.save_object(data,type, wsid, objid=objid, name=name)


def runfba(morph):
    fba_formulation = {'media': morph.media, 'media_workspace': morph.mediaws}
    fba_params = {'workspace': morph.ws_id, 'model' : morph.model, 'model_workspace':morph.ws_id,
                  'formulation': fba_formulation, 'fva': True}
    fbaMeta = service.fba_client.runfba(fba_params)
    return fbaMeta


def runmodelfba(model, wsid, media, mediaws):
    fba_formulation = {'media': media, 'media_workspace': mediaws}
    fba_params = {'workspace': wsid, 'model' : model, 'model_workspace':wsid,
                  'formulation': fba_formulation}
    fbaMeta = service.fba_client.runfba(fba_params)
    return fbaMeta


def copy_object(objid, wsid, new_ws, name=None):
    obj = get_object(objid, wsid, name=name)
    return save_object(obj['data'], obj['info'][2], new_ws, name=obj['info'][1])[0]


def get_rxn_id(modelreaction_obj):
    '''
    returns the reaction id of a model reaction object (i.e. rxn34565_c0)
    '''
    rxn_id = modelreaction_obj['reaction_ref'].split('/')[-1] + '_' + modelreaction_obj['modelcompartment_ref'].split('/')[-1]
    if rxn_id.startswith('rxn00000'):
        return modelreaction_obj['id']
    return rxn_id


def classify_gprs(gpr_comparison_set):
    gprs = gpr_comparison_set
    for rxn in gprs:
        print gprs[rxn][0] - gprs[rxn][1]


        # differences = ['','',''] #prot, sub, ftr
        # for prot in gprs[0]: # iterate over the set of proteins of first model

            # if prot not in gprs[1]:
                # they differ by this protein somehow
                # get the list of proteins in gprs[1] that contain some ftr in
                # prot
                # if the length of this list is 0, they differ by protein OR by
                # feature


def same_gpr(rxn1, rxn2):
    rxn1_set = gpr_set(rxn1)
    rxn2_set = gpr_set(rxn2)
    return rxn1_set == rxn2_set


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
    return None


def ftr_set(reaction):
    features = set()
    for protein in reaction['modelreactionproteins']:
        for sub in protein['modelreactionproteinsubunits']:
            for f in [i.split('/')[-1] for i in sub['feature_refs']]:
                features.add(f)
    return features


def gpr_exclude_features(reaction, features):
    '''
    returns a gpr set like gpr_set but excludes features in features
    '''
    rxn_proteins = reaction['modelReactionProteins']
    prots = set()
    for i in range(0, len(rxn_proteins)):
        prot = set()
        subunits = rxn_proteins[i]['modelReactionProteinSubunits']
        for j in range(0, len(subunits)):
            unit = subunits[j]
            ftrs = [f.split('/')[-1] for f in unit['feature_refs'] if f.split('/')[-1] not in features]
            if len(ftrs) > 0:
                prot.add(frozenset(ftrs))
        if len(prot) > 0:
            prots.add(frozenset(prot))
    if len(prots) > 0:
        return frozenset(prots)
    return None


def get_reaction_dict(modelobj):
    '''takes a model object and returns a set of the reactions in the model
    in rxn12345_c0 form'''
    list = modelobj['data']['modelreactions']
    rxn_dict = dict()
    for rxn in list:
        reaction = rxn['reaction_ref'].split('/')[-1] + '_' + rxn['modelcompartment_ref'].split('/')[-1]
        rxn_dict[reaction] = rxn
    return rxn_dict

# deletes all objects except for the narrative object in the test_space or
# specified workspade

def _unnest_sets(nested_set):
    '''
    Takes anything in a nested set form and returns a single set of it's non-set elements

    i.e. {{{a}}{{b, c}}{{d}{e}}} -->  {a, b, c, d, e,}
    '''
    single_set = set()
    for item in nested_set:
        if type(item) is set or type(item) is frozenset:
            single_set |= _unnest_sets(item)
        else:
            single_set.add(item)
    return single_set


def gpr_tostring(gpr_set):
    proteins = list(gpr_set)
    proteins_str = ""
    for i in range(0, len(proteins)):
        units_str = ""
        subunits = list(proteins[i])
        for j in range(0, len(subunits)):
            unit = list(subunits[j])
            features = ""
            for k in range(0, len(unit)):
                feature = unit[k]
                if (k > 0):
                    feature = " or " + feature
                features += feature
            unit_str = "(" + features + ")"
            if j > 0:
                unit_str = " and " + unit_str
            units_str += unit_str
        protein = "(" + units_str + ")"
        if (i > 0):
            protein = " or " + protein
        proteins_str += protein
    gpr = "(" + proteins_str + ")"
    return gpr


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


def print_dict(dic, f=None, g=None):
    for key in dic:
        if f is None:
            if g is None:
                print str(key) + ': ' + str(dic[key])
            else:
                if g(dic[key]):
                    print str(key) + ': ' + str(dic[key])
        else:
            if g is not None:
                if g(dic[key]):
                    print str(key) + ': ' + str(f(dic[key]))
            else:
                print str(key) + ': ' + str(f(dic[key]))


def biomass_additions(model, wsid, source_biomass, media, mediaws):
    #parse a list of compounds and coeffs from source_model that need to be added
    # for loop adjust biomass. run fba after every oteration and make sure flux > 0
    a = service.ws_client.copy_object({'from': {'objid': model, 'wsid': wsid}, 'to':{'name': 'model_copy', 'wsid': wsid}})
    model = a[0]
    model_obj = get_object(model, wsid)
    bio1 = model_obj['data']['biomasses'][0]
    compounds = list()
    coeffs = list()
    for c in source_biomass['biomasscompounds']:
        cpd = c['modelcompound_ref'].split('/')[-1]
        coeff = c['coefficient']
        compounds.append(cpd)
        coeffs.append(coeff)
    bio1_cpds = set([c['modelcompound_ref'].split('/')[-1] for c in bio1['biomasscompounds']])
    curr_model = model
    assert len(compounds) == len(coeffs)
    added_cpds = set()
    non_added = set()
    for i in range(0, len(compounds)):
        if compounds[i] not in bio1_cpds:
            a = service.ws_client.copy_object({'from': {'objid': curr_model, 'wsid': wsid}, 'to':{'name': 'biomass' + str(i), 'wsid': wsid}})
            info = service.fba_client.adjust_biomass_reaction({'model': a[0], 'workspace': wsid, 'compounds': [compounds[i]], 'coefficients':[coeffs[i]]})
            new_model = info[0]
            fba = runmodelfba(new_model, wsid, media, mediaws)
            if(fba[-1]['Objective'] > 0):
                curr_model = new_model
                added_cpds.add(compounds[i])
            else:
                non_added.add(compounds[i])
    return curr_model


def ace_to_bark(ws_id=None):
    args = dict()
    args['genome'] = Genome(3, 9145)
    args['src_model'] = FBAModel(5, 9145)
    args['probanno'] = ReactionProbabilities(15, 9145)
    args['protcomp'] = ProteomeComparison(6, 9145)
    args['media'] = Media(24, 9145)
    args['ws_id'] = ws_id

    return Morph(args)


def mari_to_janna(ws_id=None):
    args = dict()
    args['genome'] = Genome(36, 9145)
    args['src_model'] = FBAModel(58, 9145)
    args['probanno'] = ReactionProbabilities(31, 9145)
    args['protcomp'] = ProteomeComparison(27, 9145)
    args['media'] = Media(24, 9145)
    args['ws_id'] = ws_id
    return Morph(args)


def mari_to_mari(ws_id=None):
    args = dict()
    args['genome'] = Genome(26, 9145)
    args['src_model'] = FBAModel(58, 9145)
    args['probanno'] = ReactionProbabilities(30, 9145)
    args['protcomp'] = ProteomeComparison(50, 9145)
    args['media'] = Media(24, 9145)
    args['ws_id'] = ws_id
    return Morph(args)


def mari_to_bark(ws_id=None):
    args = dict()
    args['genome'] = Genome(3, 9145)
    args['src_model'] = FBAModel(58, 9145)
    args['probanno'] = ReactionProbabilities(15, 9145)
    args['protcomp'] = ProteomeComparison(51, 9145)
    args['media'] = Media(24, 9145)
    args['ws_id'] = ws_id
    return Morph(args)


def mari_to_stadt(ws_id=None):
    args = dict()
    args['genome'] = Genome(35, 9145)
    args['src_model'] = FBAModel(58, 9145)
    args['probanno'] = ReactionProbabilities(33, 9145)
    args['protcomp'] = ProteomeComparison(28, 9145)
    args['media'] = Media(24, 9145)
    args['ws_id'] = ws_id
    return Morph(args)


def _get_gpr(reaction):
    rxn_proteins = reaction['modelReactionProteins']
    proteins_str = ""
    for i in range(0, len(rxn_proteins)):
        units_str = ""
        subunits = rxn_proteins[i]['modelReactionProteinSubunits']
        for j in range(0, len(subunits)):
            unit = subunits[j]
            ftrs = unit['feature_refs']
            features = ""
            if (len(ftrs) == 0):
                if unit['note'] != "":
                    features = "(Unknown)"
                else:
                   continue
            else:
                for k in range(0, len(ftrs)):
                    feature = ftrs[k].split('/')[-1]
                    if (k > 0):
                        feature = " or " + feature
                    features += feature
            unit_str = "(" + features + ")"
            if j > 0:
                unit_str = " and " + unit_str
            units_str += unit_str
        protein = "(" + units_str + ")"
        if (i > 0):
            protein = " or " + protein
        proteins_str += protein
    gpr = "(" + proteins_str + ")"
    return gpr


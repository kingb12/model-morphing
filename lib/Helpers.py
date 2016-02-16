from Morph import Morph
from Morph import Gpr
import copy
import json
import pickle
import random
import Client
from firebase import firebase
import os

def make_morph(ws_id=None):
    args = dict()
    args['genome'] = '3'
    args['src_model'] = '5'
    args['probanno'] = '15'
    args['protcomp'] = '6'
    args['genomews'] = '9145'
    args['src_modelws'] = '9145'
    args['probannows'] = '9145'
    args['protcompws'] = '9145'
    args['mediaws'] = '9145'
    args['media'] = '24'
    args['ws_id'] = ws_id
    return Morph(args)
def modelargs(morph):
    args = dict()
    args['model'] = morph.model
    args['model_workspace'] = morph.ws_id
    return args

def essential_test(morph, rxn_list):
    results = dict()
    for rxn in rxn_list:
        print rxn
        a,b,c = Client.find_alternative(rxn, morph=morph)
        analysis = Client.analyze_alternative(morph, rxn, c)
        results[rxn] = (a, b, c, analysis)
        filename = "essential_results/" + str(rxn) + ".pkl"
        pickle.dump(results, open(filename, "wb"))
        print "saved"
# CURRENT USE SPECIFIC INITIALIZATION CODE

def unpack_essen_results():
    path = str(os.getcwd()) + '/../data/essential_results'
    reactions = dict()
    for filename in os.listdir(path):
        filename = path + str('/') + filename
        with open(filename, "rb") as f:
            print str(filename)
            a = pickle.load(f)
            reactions[filename.split('/')[-1].split('.')[0]] = a
    filename = str(os.getcwd()) + '/../data/reactions.pkl'
    return reactions

def _count_essential(reactions):
    rxns = reactions
    reactions = dict()
    for rxn in rxns:
        data = rxns[rxn]
        for r in data[rxn][2]:
            if r not in reactions:
                reactions[r] = [rxn]
            else:
                reactions[r].append(rxn)
    return reactions

def pp_eval_rxns(reactions):
    '''
    prints the essential reactions added  fro an essential reaction test (essential_test)

    Requires: The input to the function should be the results of the unpack_essen_results() function. These are separated to save computational time in some circumstances
    '''
    count_rxns = _count_essential(reactions)
    added_rxns = Client.fba_client.get_reactions({'reactions': [r.split('_')[0] for r in count_rxns.keys()]})
    with open('../analysis/fill_frequency.txt', 'w') as f:
        for r in added_rxns:
            f.write(str(r['name']) + '\n')
            f.write(str(r['definition']) + '\n')
            try:
                f.write('Frequency in which it was added in gapfilling: ' + str(float(len(count_rxns[str(r['id'] + '_c0')])) / float(len(count_rxns.keys()))) + '\n' + '\n' + '\n')
            except KeyError:
                f.write('Frequency in which it was added in gapfilling: ' + str(float(len(count_rxns[str(r['id'] + '_e0')])) / float(len(count_rxns.keys()))) + '\n' + '\n' + '\n')
    with open('../analysis/essen_rxn_info.txt', 'w') as f:
         for rxn in reactions:
             print type(rxn)
    return reactions

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
    if (name is None):
        return Client.ws_client.get_objects([{'objid':objid, 'wsid':wsid}])[0]
    else:
        return Client.ws_client.get_objects([{'name':name, 'wsid':wsid}])[0]
def get_info(objid, wsid, name=None):
    if (name is None):
        return Client.ws_client.get_object_info_new({'objects': [{'objid':objid, 'wsid':wsid}]})[0]
    else:
        return Client.ws_client.get_objects_info_new({'objects': [{'name':name, 'wsid':wsid}]})[0]

def dump(object, filepath):
    with open(filepath, 'w') as f:
        pickle.dump(object, f)

def load(filepath):
    with open(filepath, 'rb') as f:
        return pickle.load(f)
def save_object(data, type, wsid, objid=None, name=None):
    sv = {u'data':data, u'type':type, u'name':name}
    if objid is not None:
        sv[u'objid'] = objid
    if name is not None:
        sv[u'name'] = name
    return Client.ws_client.save_objects({u'id':wsid, u'objects':[sv]})
def load_morph():
    """
    Loads a morph that has completed
    """
    data_unicode = fb.get('standard_morph', 'acebar')
    data = json.loads(data_unicode, object_hook=Morph)
    return data

def put(url, key, value):
    """
    Puts value at the givenkey at url in the Firebase

    (the url is just the subdirectory within firebase, ie. to store "Morty" at 'https://my_base/firebaseio.com/BestCats'
    call put('/BestCats', "Morty")
    """
    try:
        fb.put(url, key, value.to_JSON())
    except AttributeError:
        fb.put(url, key, value)

def check_for_duplicates(model_id, ws_id):
    object = get_object(model_id, ws_id)
    model = object['data']
    cpds = [(c['compound_ref'], c['modelcompartment_ref']) for c in model['modelcompounds']]
    compound_dict = dict()
    for cpd in cpds:
        cpd_id = cpd[0].split('/')[-1]
        if (cpd_id not in compound_dict):
            compound_dict[cpd_id] = [cpd]
        else:
            compound_dict[cpd_id].append(cpd)
    dups = list()
    for i in compound_dict:
        compartments = [j[1] for j in compound_dict[i]]
        # Sets can only have unique values, so if there is a duplicate
        # compartment, it woll cause a change in length from list -> set
        compset = set(compartments)
        if(len(compset) != len(compartments)):
            dups.append(i)
    duplicates = dict()
    for i in dups:
        duplicates[i] = compound_dict[i]
    mdl = model
    return duplicates, mdl
def process_iterated(morph):
    gnm = Client.removal_list(morph.rxn_labels['gene-no-match'])
    for i in gnm:
        (morph, a, b) = Client.process_reactions(morph, rxn_list=[i])
        dump(morph, '../data/morph-Hsmedia.pkl')
    ng = Client.removal_list(morph.rxn_labels['no-gene'])
    for i in ng:
        (morph, a, b) = Client.process_reactions(morph, rxn_list=[i])
        dump(morph, '../data/morph-Hsmedia.pkl')
    return morph

def runfba(morph):
    fba_formulation = {'media': morph.media, 'media_workspace': morph.mediaws}
    fba_params = {'workspace': morph.ws_id, 'model' : morph.model, 'model_workspace':morph.ws_id,
                  'formulation': fba_formulation}
    fbaMeta = Client.fba_client.runfba(fba_params)
    return fbaMeta

def runmodelfba(model, wsid, media, mediaws):
    fba_formulation = {'media': media, 'media_workspace': mediaws}
    fba_params = {'workspace': wsid, 'model' : model, 'model_workspace':wsid,
                  'formulation': fba_formulation}
    fbaMeta = Client.fba_client.runfba(fba_params)
    return fbaMeta

def clone_morph(morph):
    morph = copy.deepcopy(morph)
    info = Client.ws_client.clone_workspace({'wsi': {'id': morph.ws_id}, 'workspace':str('morphclone' + str(morph.ws_id) + str(random.randint(0, 1000000)))})
    morph.ws_id = info[0]
    return morph

def copy_object(objid, wsid, new_ws, name=None):
    obj = get_object(objid, wsid, name=name)
    return save_object(obj['data'], obj['info'][2], new_ws, name=obj['info'][1])[0]

def empty_ws(ws_id):
    # get all object ids in {'id': x} form excluding the narrative object
    object_ids = [{'objid': info[0], 'wsid': ws_id} for info in
                 Client.ws_client.list_objects({'ids': [ws_id]}) if not info[2].startswith('KBaseNarrative')]
    if len(object_ids) > 0:
        Client.ws_client.delete_objects(object_ids)

def rename_object(obj_id, ws_id, new_name):
    Client.ws_client.rename_object({'obj': {'wsid': ws_id, 'objid': obj_id}, 'new_name': new_name})

fb = firebase.FirebaseApplication('https://fiery-fire-3234.firebaseio.com', None)

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
    return none

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




def mari_to_janna(ws_id=None):
    args = dict()
    args['genome'] = '36'
    args['src_model'] = '37'
    args['probanno'] = '31'
    args['protcomp'] = '27'
    args['genomews'] = '9145'
    args['src_modelws'] = '9145'
    args['probannows'] = '9145'
    args['protcompws'] = '9145'
    args['mediaws'] = '9145'
    args['media'] = '24'
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


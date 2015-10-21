from Morph import Morph
import pickle
import Client
import os
def make_morph():
    args = dict()
    args['genome'] = '3'
    args['src_model'] = '19'
    args['probanno'] = '15'
    args['protcomp'] = '6'
    args['genomews'] = '9145'
    args['src_modelws'] = '9145'
    args['probannows'] = '9145'
    args['protcompws'] = '9145'
    args['mediaws'] = '9145'
    args['media'] = '18'
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
def dump(object, filepath):
    with open(filepath, 'w') as f:
        pickle.dump(object, f)

def load(filepath):
    with open(filepath, 'rb') as f:
        return pickle.load(f)
def save_object(data, type, wsid, objid=None):
    sv = {u'data':data, u'type':type}
    if objid is not None:
        sv[u'objid'] = objid
    return Client.ws_client.save_objects({u'id':wsid, u'objects':[sv]})
def load_morph():
    """
    Loads a morph that has completed
    """
    return load('../data/morph_data.pkl')
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

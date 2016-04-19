import Helpers
import Reaction
import Model
import GrowthConditions
import itertools
from Helpers import *
import Client
import Plotter
import Analysis

from biokbase.fbaModelServices.Client import ServerError
# This is a module of functions that perform an unorganized collection of runs
# and routines, not for production usage


def alternative_shotgun():
    morph = Helpers.load('../data/morph-Hsmedia.pkl')
    models = list()
    stats = list()
    morphs = list()

    for n in range(1, 40, 3):
        m, a, b, c = Client.find_alternative(morph, morph.essential_ids.items()[n])
        models.append(m.model)
        morphs.append(m)
        probs = list()
        for r in b:
            try:
                probs.append(morph.probhash[r.split('_')[0]])
            except KeyError:
                probs.append(-1.0)
        stats.append({'new_reactions': b, 'probabilities': probs })
        filename = '../data/altshot/altshot-' + str(n) + '.pkl'
        dump({'morph': morph, 'morphs': morphs, 'models': models, 'stats': stats}, filename)
        dump(models, str(filename.split('.')[0] + '-model.pkl'))

def get_altshot_results():
    path = str(os.getcwd()) + '/../data/altshot'
    results = list()
    for filename in os.listdir(path):
        filename = path + str('/') + filename
        with open(filename, "rb") as f:
            print str(filename)
            a = pickle.load(f)
            results.append(a)
    filename = str(os.getcwd()) + '/../data/reactions.pkl'
    return results

def print_rxnsprobs(result):
    s = result
    rxns = Client.fba_client.get_reactions({'reactions': [a.split('_')[0] for a in s['new_reactions']]})
    for r in rxns:
        print r['definition']
    for r in rxns:
        print r['id']
    for r in s['new_reactions']:
        print r
    for r in s['probabilities']:
        print r

def check_result_fbas(result):
    for m in result['morphs']:
        a = runfba(m)
        print a[10]['Objective']

def move_models(result, ws):
    for m in result['models']:
        obj = get_object(m, result['morph'].ws_id)
        save_object(obj['data'], obj['info'][2], ws, name=obj['info'][1])

def check_essentiality(result, index):
    for s in result['stats'][index]['new_reactions']:
        m = copy.deepcopy(result['morphs'][index])
        m.model = Client.fba_client.remove_reactions({'model': result['models'][index], 'model_workspace': result['morphs'][index].ws_id, 'output_id': 'essncheck', 'workspace':result['morphs'][index].ws_id, 'reactions':[s]})[0]
        a = runfba(m)
        print a[10]['Objective']
def supermodel_check(morph):
    super_rxns = get_object(morph.model, morph.ws_id)['data']['modelreactions']
    recon_rxns = morph.objects['recon_model']['modelreactions']
    super_hash = dict()
    for i in range(0, len(super_rxns)):
        mdlrxn = super_rxns[i]
        rxn_id = mdlrxn['reaction_ref'].split('/')[-1] + '_' + str(mdlrxn['modelcompartment_ref'].split('/')[-1]) # -1 index gets the last in the list
        super_hash[rxn_id] = i
    for rxn in morph.rxn_labels['recon']:
        index1 =  morph.rxn_labels['recon'][rxn][0]
        index2 = super_hash[rxn]
        recon_proteins = recon_rxns[index1]['modelReactionProteins']
        super_proteins = super_rxns[index2]['modelReactionProteins']
        assert len(recon_proteins) == len(super_proteins)
        for i in range(0, len(recon_proteins)):
            recon_subs = recon_proteins[i]['modelReactionProteinSubunits']
            super_subs = super_proteins[i]['modelReactionProteinSubunits']
            assert len(recon_subs) == len(super_subs)
            for j in range(0, len(recon_subs)):
                recon_unit = recon_subs[j]
                super_unit = super_subs[j]
                assert recon_unit['feature_refs'] == super_unit['feature_refs']

def five_model_morpher(media, mediaws, source_biomass=None,src_morph=None):
    workspaces = (11782, 11783)
    #one by one
    morph = make_morph(ws_id=workspaces[0])
    if src_morph is None:
        morph = make_morph(ws_id=workspaces[0])
    else:
        morph = src_morph
        morph.ws_id = workspaces[0]
    #update media
    morph.media = media
    morph.mediaws = mediaws
    media_name = Helpers.get_object(media, mediaws)['info'][1]
    redo = 0
    while(redo < 5):
        try:
            morph = Client.prepare_supermodel(morph)
            if source_biomass is not None:
                morph.model = biomass_additions(morph.model, morph.ws_id, source_biomass, morph.media, morph.mediaws)
            redo = 9
        except ServerError:
            redo += 1
    reaction_list = Client.removal_list(morph.rxn_labels['gene-no-match'])
    reaction_list = [r for r in reaction_list if r[0] not in morph.rxn_labels['common']]
    morph = Client.process_reactions(morph, rxn_list=reaction_list)
    reaction_list = Client.removal_list(morph.rxn_labels['no-gene'])
    reaction_list = [r for r in reaction_list if r[0] not in morph.rxn_labels['common']]
    morph = Client.process_reactions(morph, rxn_list=reaction_list)
    dump(morph, '../data/morph-' + media_name + '.pkl')

    #all at once removal
    if src_morph is None:
        morph = make_morph(ws_id=workspaces[1])
    else:
        morph = src_morph
        morph.ws_id = workspaces[1]
    morph.media = media
    morph.mediaws = mediaws
    while(redo < 5):
        try:
            morph = Client.prepare_supermodel(morph)
            if source_biomass is not None:
                morph.model = biomass_additions(morph.model, morph.ws_id, source_biomass, morph.media, morph.mediaws)
            redo = 9
        except ServerError:
            redo += 1
    gnm_list = Client.removal_list(morph.rxn_labels['gene-no-match'])
    gnm_list = [r[0] for r in gnm_list if r[0] not in morph.rxn_labels['common']]
    ng_list = Client.removal_list(morph.rxn_labels['no-gene'])
    ng_list = [r[0] for r in ng_list if r[0] not in morph.rxn_labels['common']]
    morph = Client.remove_reactions(morph, gnm_list, output_id =media_name +  '-gnm-removed')
    morph_ng = Client.remove_reactions(morph, ng_list, output_id =media_name +  '-gnm-and-ng-removed')

    #various gapfillings
    redo = 0
    while(redo < 5):
        try:
            m2 = Client.probanno_fill(morph, name=media_name + '-gnm-rmv-prob-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m3 = Client.probanno_fill(morph_ng, name=media_name + '-gnm-and-ng-rmv-prob-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m4 = Client.parse_fill(morph, name=media_name + '-gnm-rmv-parse-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m5 = Client.parse_fill(morph_ng, name=media_name + '-gnm-and-ng-rmv-parse-fill')
            redo = 9
        except ServerError:
            redo += 1

    dump(m2, '../data/morph-' + media_name + '-all1-prob.pkl')
    dump(m3, '../data/morph-' + media_name + '-all2-prob.pkl')
    dump(m4, '../data/morph-' + media_name + '-all1-parse.pkl')
    dump(m5, '../data/morph-' + media_name + '-all2-parse.pkl')

def all_removal(media, mediaws, source_biomass=None, src_morph=None):

    #all at once removal
    redo = 0
    if src_morph is None:
        morph = make_morph(ws_id=11783)
    else:
        morph = src_morph
        morph.ws_id = 11783
    #all at once removal
    morph.media = media
    morph.mediaws = mediaws
    media_name = Helpers.get_object(media, mediaws)['info'][1]
    while(redo < 5):
        try:
            morph = Client.prepare_supermodel(morph)
            if source_biomass is not None:
                morph.model = biomass_additions(morph.model, morph.ws_id, source_biomass, morph.media, morph.mediaws)
            redo = 9
        except ServerError:
            redo += 1
    gnm_list = Client.removal_list(morph.rxn_labels['gene-no-match'])
    gnm_list = [r[0] for r in gnm_list if r[0] not in morph.rxn_labels['common']]
    ng_list = Client.removal_list(morph.rxn_labels['no-gene'])
    ng_list = [r[0] for r in ng_list if r[0] not in morph.rxn_labels['common']]
    morph = Client.remove_reactions(morph, gnm_list, output_id =media_name +  '-gnm-removed')
    morph_ng = Client.remove_reactions(morph, ng_list, output_id =media_name +  '-gnm-and-ng-removed')

    #various gapfillings
    redo = 0
    while(redo < 5):
        try:
            m2 = Client.probanno_fill(morph, name=media_name + '-gnm-rmv-prob-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m3 = Client.probanno_fill(morph_ng, name=media_name + '-gnm-and-ng-rmv-prob-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m4 = Client.parse_fill(morph, name=media_name + '-gnm-rmv-parse-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m5 = Client.parse_fill(morph_ng, name=media_name + '-gnm-and-ng-rmv-parse-fill')
            redo = 9
        except ServerError:
            redo += 1

    dump(m2, '../data/morph-' + media_name + '-all1-prob.pkl')
    dump(m3, '../data/morph-' + media_name + '-all2-prob.pkl')
    dump(m4, '../data/morph-' + media_name + '-all1-parse.pkl')
    dump(m5, '../data/morph-' + media_name + '-all2-parse.pkl')

def find_kbaliases(filename, genome_id, ws_id):
    """
    given a list of genes, finds their KBase ID aliases and return them as a list
    """
    genes = dict()
    with open(filename, "r") as file:
        gene_list = [g.split('\n')[0] for g in file.readlines()]
    genome = get_object(genome_id, ws_id)['data']
    # create a hash of a set of aliases to their unique KBase ID
    # Consider making the hash a morph attribute TODO
    gene_hash = dict()
    for feat in genome['features']:
        kb_id = feat['id']
        if 'aliases' in feat.keys():
            for alias in feat['aliases']:
                gene_hash[alias] = kb_id
    # find the appropriate genes from the hash
    for i in range(0, len(gene_list)):
        genes[gene_list[i]] = gene_hash[gene_list[i].split('\n')[0]]
    return genes

def simulate_all_models(ws_id, phenotype_id):
    #Get all the objects in ws_id
    models = list() # A list of tuples of the form (obj_id, ws_id, obj_name)
    objects = Client.ws_client.list_objects({'ids' : [ws_id]})
    for obj in objects:
        type = obj[2]
        if 'FBAModel' in type: #is it a model
            if 'iAF' not in obj[1]:
                models.append((obj[0], obj[6], obj[1]))
        if phenotype_id is None and 'PhenotypeSet' in type:
            phenotype_id = obj[0]
    for model in models:
        args = dict()
        args['model'] = model[0]
        args['model_workspace'] = model[1]
        args['workspace'] = ws_id
        args['phenotypeSet'] = phenotype_id
        args['phenotypeSet_workspace'] = ws_id
        args['phenotypeSimulationSet'] = str(model[2]) + '-pheno_data'
        args['overwrite'] = False
        print args
        info = Client.fba_client.simulate_phenotypes(args)
        Helpers.rename_object(info[0], ws_id, args['phenotypeSimulationSet'])

def try_FBA_geneKO(model, media, ws_id):
    formulation = {'media': media, 'media_workspace': ws_id, 'rxnko': ['rxn03079_c0']}
    Client.fba_client.runfba({'model' : model, 'model_workspace': ws_id, 'formulation':formulation, 'fba':'testFBAgeneKO', 'workspace':ws_id})
def compare_all_models(ws_id):
    models = _allmodelsinws(ws_id)
    args = {'models': [m[0] for m in models], 'workspaces': [m[1] for m in models]}
    return Client.fba_client.compare_models(args)

def all_model_stats(ws_id):
    models = _allmodelsinws(ws_id)
    results = dict()
    names = _names_dict(ws_id)
    for m in models:
        stats = Client.fba_client.generate_model_stats({'model': m[0], 'model_workspace': m[1]})
        results[names[m[0]]] = stats
    return results

def _allmodelsinws(ws_id):
    models = list() # A list of tuples of the form (obj_id, ws_id, obj_name)
    objects = Client.ws_client.list_objects({'ids' : [ws_id]})
    for obj in objects:
        type = obj[2]
        if 'FBAModel' in type: #is it a model
            models.append((obj[0], obj[6], obj[1]))
    return models
def reaction_in_model(model_obj, reaction):
    compartment = None
    if len(reaction.split('_')) == 2:
        compartment = reaction.split('_')[1]
        reaction = reaction.split('_')[0]
    for rxn in model_obj['modelreactions']:
        if rxn['reaction_ref'].split('/')[-1] == reaction:
            if compartment is None or rxn['modelcompartment_ref'].split('/')[-1] == compartment:
                return True
    return False


def all_models_counts(ws_id, model_names):
    names = _names_dict(ws_id)
    comp = compare_all_models(ws_id)
    results = dict()
    for rxn in comp['reaction_comparisons']:
        val = dict()
        val['number_models'] = rxn['number_models']
        # add a key for all the model_names
        for name in model_names:
            val[name] = (0,0)
        for model_ref in rxn['model_features']:
            numgenes = len(rxn['model_features'][model_ref])
            name = names[int(model_ref.split('/')[-1])]
            if(numgenes == 1):
                if rxn['model_features'][model_ref][0][0:2] != u'kb':
                    numgenes = numgenes - 1
            val[name] = (1, numgenes)
        key = rxn['reaction'] + '_' + rxn['compartment'] + str(0)
        results[key] = val
    return results

def make_tab_file(reaction_info, models, filename):
    with open(filename, 'w') as f:
        for rxn in reaction_info:
            f.write(str(rxn) + '\t')
            for m in models:
                f.write(str(reaction_info[rxn][m][0]) + '\t' + str(reaction_info[rxn][m][1]) + '\t')
            f.write('\n')


def _names_dict(ws_id):
    models = dict() # A list of tuples of the form (obj_id, ws_id, obj_name)
    objects = Client.ws_client.list_objects({'ids' : [ws_id]})
    for i in objects:
        models[i[0]] = i[1]
    return models
def venn_analysis(reaction_info, models):
    venn = dict()
    letters = ['A', 'B', 'C', 'D', 'E', 'F']
    for rxn in reaction_info:
        key = ""
        for i in range(0, len(models)):
            if reaction_info[rxn][models[i]][0] > 0:
                key = key + letters[i]
                print key
        try:
            venn[key].append(rxn)
        except KeyError:
            venn[key] = [rxn]
    return venn

def get_rxninfo_sublist(reaction_info, names, no_gene=True):
    result = dict()
    for name in names:
        result[name] = []
    for key in reaction_info.keys():
        rxn = reaction_info[key]
        for k in rxn.keys():
            if k in result: #if its a name we are concerned with
                if (no_gene and rxn[k][0] == 1) or rxn[k][1] > 0:
                    result[k].append(key)
    return result

def get_comp_info(reactions, comparison):
    comp_dict = dict()
    result = list()
    for rxn in comparison:
        key = rxn['reaction'] +  '_' + rxn['compartment'] + '0'
        comp_dict[key] = rxn
    for rxn in reactions:
        result.append(comp_dict[rxn])
    return result

def find_gene_relationships(ws_id=None):
    # Make a morph, a translation, and a reconstruction
    m = ace_to_bark(ws_id)
    m = Client.translate_features(m)
    m = Client.reconstruct_genome(m)
    # get the model objects
    recon_model = get_object(m.recon_model, m.ws_id)
    trans_model = get_object(m.trans_model, m.ws_id)
    #return a set of all disagreeing GPR's
    result = dict() # rxn_id -> GPR (maybe get_gpr?)
    recon_rxns = get_reaction_dict(recon_model)
    trans_rxns = get_reaction_dict(trans_model)
    recon_set = set(recon_rxns.keys())
    trans_set = set(trans_rxns.keys())
    common_set = recon_set & trans_set
    for rxn in common_set:
        if  not same_gpr(recon_rxns[rxn], trans_rxns[rxn]):
            result[rxn] = (gpr_set(trans_rxns[rxn]), gpr_set(recon_rxns[rxn]))
    classes = simple_classify(result.keys(), recon_rxns, trans_rxns)
    return result, common_set, classes

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

def simple_classify(common_set, recon_rxns, trans_rxns):
    result = {'recon_ftrs': dict(), 'trans_ftrs': dict(), 'both': dict()}
    for rxn in common_set:
        if len(ftr_set(recon_rxns[rxn]) - ftr_set(trans_rxns[rxn])) > 0:
            if len(ftr_set(trans_rxns[rxn]) - ftr_set(recon_rxns[rxn])) > 0:
                result['both'][rxn] = (recon_rxns[rxn], trans_rxns[rxn])
            else:
                result['recon_ftrs'][rxn] = (recon_rxns[rxn], trans_rxns[rxn])
        else:
            result['trans_ftrs'][rxn] = (recon_rxns[rxn], trans_rxns[rxn])

    return result

def genes_to_reactions(model_object, list_of_genejs):
    '''
    returns a dictionary of reactions associateed with the listed genes. The genes are keys
    '''
    ftr_hash = feature_hash(model_object)
    result = dict()
    for gene in list_of_genes:
        try:
            result[gene] = ftr_hash[gene]
        except KeyError:
            print str(gene) + ' not found in model'
    return result

def feature_hash(model_object):
    '''
    RETURNS A DICTIONARY OF ALL GENE FEATURES TO THEIR ASSOCIATED REACTIONS
    '''
    result = dict()
    for rxn in model_object['modelreactions']:
        features = set()
        for protein in rxn['modelReactionProteins']:
            for sub in protein['modelReactionProteinSubunits']:
                features |= set([i.split('/')[-1] for i in sub['feature_refs']])
        for ftr in features:
           result[ftr] = Helpers.get_rxn_id(rxn)
    return result

def isInModel(model_obj, compound_id):
    for rxn in model_obj['data']['modelreactions']:
        for cpd in rxn['modelReactionReagents']:
            if(cpd['modelcompound_ref'].find(compound_id) >= 0):
                return True
    return False

def biomass_additions(model, wsid, source_biomass, media, mediaws):
    #parse a list of compounds and coeffs from source_model that need to be added
    # for loop adjust biomass. run fba after every oteration and make sure flux > 0
    a = Client.ws_client.copy_object({'from': {'objid': model, 'wsid': wsid}, 'to':{'name': 'model_copy', 'wsid': wsid}})
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
            a = Client.ws_client.copy_object({'from': {'objid': curr_model, 'wsid': wsid}, 'to':{'name': 'biomass' + str(i), 'wsid': wsid}})
            info = Client.fba_client.adjust_biomass_reaction({'model': a[0], 'workspace': wsid, 'compounds': [compounds[i]], 'coefficients':[coeffs[i]]})
            new_model = info[0]
            fba = Helpers.runmodelfba(new_model, wsid, media, mediaws)
            if(fba[-1]['Objective'] > 0):
                curr_model = new_model
                added_cpds.add(compounds[i])
            else:
                non_added.add(compounds[i])
    return curr_model

def compare_biomasses(bio1, bio2):
    bio1_cpds = set([c['modelcompound_ref'].split('/')[-1] for c in bio1['biomasscompounds']])
    bio2_cpds = set([c['modelcompound_ref'].split('/')[-1] for c in bio2['biomasscompounds']])
    both = set()
    firstonly = set()
    secondonly = set()
    both = bio1_cpds & bio2_cpds
    firstonly = bio1_cpds - bio2_cpds
    secondonly = bio2_cpds - bio1_cpds
    return both, firstonly, secondonly

def reconcile_models(good_model, working_model, ws_id):
    '''
    Takes a working model and tries to make it as close to good model as possible. requires in same ws
    '''
    # Classify reactions in models
    good = get_object(good_model, ws_id)
    working = get_object(working_model, ws_id)
    good_rxns = dict()
    working_rxns = dict()
    common = set()
    good_only = set()
    working_only = set()
    for r in good['data']['modelreactions']:
        good_rxns[get_rxn_id(r)] = r
    for r in working['data']['modelreactions']:
        working_rxns[get_rxn_id(r)] = r
        if get_rxn_id(r) in good_rxns:
            common.add(get_rxn_id(r))
        else:
            working_only.add(get_rxn_id(r))
    for rxn_id in good_rxns:
        if rxn_id not in common:
            good_only.add(rxn_id)
    # Add good reactions
    new_model = None

    # Adjust common reactions
    # Remove unnecessary reactions
def model_checker(model_rxns, super_rxns):
    results = list()
    for r in super_rxns:
        if r in model_rxns:
            if super_rxns[r]['modelReactionReagents'] != model_rxns[r]['modelReactionReagents']:
                for c in super_rxns[r]['modelReactionReagents']:
                    if c not in model_rxns[r]['modelReactionReagents']:
                        results.append(r)
                for c in model_rxns[r]['modelReactionReagents']:
                    if c not in super_rxns[r]['modelReactionReagents']:
                        results.append(r)
    return results

# save functions as variables in Runs
same_gpr = Helpers.same_gpr
gpr_set = Helpers.gpr_set
ftr_set = Helpers.ftr_set

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

def morph_comparison(morph, super_model):
    comp = Client.fba_client.compare_models({'models': [morph.model, morph.src_model, morph.recon_model, morph.trans_model], 'workspaces': [morph.ws_id, morph.src_modelws, morph.ws_id, morph.ws_id]})
    # General Comparisons
    a = Client.fba_client.generate_model_stats({'model': morph.model, 'model_workspace': morph.ws_id})
    b = Client.fba_client.generate_model_stats({'model': morph.src_model, 'model_workspace': morph.src_modelws})
    c = Client.fba_client.generate_model_stats({'model': morph.trans_model, 'model_workspace': morph.ws_id})
    d = Client.fba_client.generate_model_stats({'model': morph.recon_model, 'model_workspace': morph.ws_id})
    # Gene comparisons
    return [a,b,c,d]
    # Probanno comparisons

def test(m):
    a = list()
    supm = get_object(m.model, m.ws_id)['data']
    src = get_object(m.src_model, m.ws_id)['data']
    b = dict()
    c = dict()
    for r in supm['modelreactions']:
        b[get_rxn_id(r)] = r
    for r in src['modelreactions']:
        c[get_rxn_id(r)] = r
        k = None
        for key in m.rxn_labels:
            if get_rxn_id(r) in m.rxn_labels[key]:
                k = key
    for r in c:
        if r not in b:
            for key in m.rxn_labels:
                if r in m.rxn_labels[key]:
                    k = key
            print key
            a.append((c[r], key))
    return a


def mari_to_janna_morph(ws_id, newmedia=None, newmediaws=None):
    m = mari_to_janna(ws_id=ws_id)
    m = Client.prepare_supermodel(m, fill_src=False)
    if newmedia is not None:
        m2 = Client.translate_media(m, newmedia, newmediaws)
    else:
        m2 = copy.deepcopy(m)
    m3 = Client.process_reactions(m2)
    dump(m3, '../data/mari_to_janna_morph.pkl')

def mari_to_bark_morph(ws_id, newmedia=None, newmediaws=None, filename='../data/mari_to_bark_morph.pkl'):
    m = mari_to_bark(ws_id=ws_id)
    m = Client.prepare_supermodel(m, fill_src=False)
    if newmedia is not None:
        m2 = Client.translate_media(m, newmedia, newmediaws)
    else:
        m2 = copy.deepcopy(m)
    m3 = Client.process_reactions(m2)
    dump(m3, filename)

def mari_to_bark3_morph(ws_id, newmedia=None, newmediaws=None, filename='../data/mari_to_bark_morph.pkl'):
    m = mari_to_bark(ws_id=ws_id)
    m = Client.prepare_supermodel(m, fill_src=False)
    media = [(24, 9145), (7, 12091), (8, 12091)]
    for med in media:
        m = Client.translate_media(m, med[0], med[1])
    else:
        m = copy.deepcopy(m)
    m3 = Client.process_reactions(m, growth_condition= GrowthConditions.AllMedia(media))
    dump(m3, filename)

def mari_to_stadt_morph(ws_id, newmedia=None, newmediaws=None):
    m = mari_to_stadt(ws_id=ws_id)
    m = Client.prepare_supermodel(m, fill_src=False)
    if newmedia is not None:
        m2 = Client.translate_media(m, newmedia, newmediaws)
    else:
        m2 = copy.deepcopy(m)
    m3 = Client.process_reactions(m2)
    dump(m3, '../data/mari_to_stadt_morph.pkl')

def mari_to_mari_morph(ws_id):
    m = mari_to_mari(ws_id=ws_id)
    m = Client.prepare_supermodel(m, fill_src=False)
    m2 = copy.deepcopy(m)
    m3 = Client.process_reactions(m2)
    dump(m3, '../data/mari_to_mari_morph.pkl')

def gene_analysis(model, ws_id):
    model = Helpers.get_object(model, ws_id)['data']
    gprs = dict()
    genes = set()
    spontaneous = 0
    reactions_with_genes = 0
    for r in model['modelreactions']:
        gpr = Gpr(reaction=r)
        if gpr.gpr_type is 'genes':
            reactions_with_genes += 1
            genes |= gpr.ftrs
        if gpr.gpr_type == 'spontaneous' or gpr.gpr_type == 'universal':
            spontaneous += 1
        gprs[Reaction.get_rxn_id(r)] = gpr
    result = {'genes': len(genes), 'reactions_with_genes': reactions_with_genes, 'gprs': gprs, 'spontaneous': spontaneous}
    return result


def reaction_analysis(model_identities):
    models = list()
    for m in model_identities:
        model = Helpers.get_object(m[0], m[1])
        models.append((model['info'][1], frozenset([Helpers.get_rxn_id(r) for r in Model.get_reactions(model['data'])])))
    all_models = set(models)
    rxn_sets = dict()
    for model_set in powerset(models):
        if len(model_set) > 0:
            rxn_set = model_set[0][1]
            keys = set(model_set)
            skips = all_models - keys
            key = ""
            for m in model_set:
                key += m[0] + " + "
            key = key[0:-3]
            for m in model_set:
                rxn_set = rxn_set.intersection(m[1])
            print [m[0] for m in skips]
            for m in skips:
                rxn_set = rxn_set - m[1]


            rxn_sets[key] = rxn_set
    return rxn_sets

def gene_percents(model, rxn_analysis):
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

def reaction_sources(morph):
    result = dict()
    for key in morph.rxn_labels:
        result[key] = []
    model = Helpers.get_object(morph.model, morph.ws_id)['data']
    reactions = set([get_rxn_id(r) for r in model['modelreactions']])
    for r in reactions:
        for key in morph.rxn_labels:
            if r in morph.rxn_labels[key]:
                result[key].append(r)
    return result

def reaction_analysis_info(rxn_analysis):
    result = dict()
    for key in rxn_analysis:
        rxn_list = list()
        for r in rxn_analysis[key]:
            rxn_list.append(r.split('_')[0])
        result[key] = Client.fba_client.get_reactions({'reactions' : rxn_list})
    return result

def feature_sources(morph):
    result = dict()
    translation = set()
    reconstruction = set()
    for r in Model.get_reactions(morph.objects['trans_model']):
        gpr = Gpr(reaction=r)
        translation |= gpr.ftrs
    for r in Model.get_reactions(morph.objects['recon_model']):
        gpr = Gpr(reaction=r)
        reconstruction |= gpr.ftrs
    result['both'] = translation.intersection(reconstruction)
    result['recon_only'] = reconstruction - translation
    result['trans_only'] = translation - reconstruction
    return result

def gene_reactions(models, rxn_list):
    for m in models:
        model = get_object(m[0], m[1])
        reactions = dict()
        for r in Model.get_reactions(model['data']):
            reactions[Reaction.get_rxn_id(r)] = r
        for r in rxn_list:
            gpr = Gpr(reactions[r])
            if gpr.gpr_type == 'no-gene':
                rxn_list.remove(r)
            print gpr.gpr_type


def find_reactions_by_compounds(models, compound_list):
    result = set()
    for m in models:
        model = Helpers.get_object(m[0], m[1])['data']
        result |= set([Reaction.get_rxn_id(r) for r in Model.find_compound_uses(model, compound_list[0])])
    for c in compound_list[1:len(compound_list)]:
        rxn_set = set()
        for m in models:
            model = Helpers.get_object(m[0], m[1])['data']
            rxn_set |= set([Reaction.get_rxn_id(r) for r in Model.find_compound_uses(model, c)])
        result = result.intersection(rxn_set)
    return result

def find_reaction_uses(model_list, reaction_list):
    '''

    :param model_list: model tuples
    :param reaction_list: reaction_ids
    :return: a dictionary where keys are elements of reaction_list and values are the models which have this reaction
    '''
    result = dict()
    for r in reaction_list:
        models = list()
        for m in model_list:
            model = get_object(m[0], m[1])['data']
            if Model.has_reaction(model, r):
                models.append(m[2])
        result[r] = models
    return result

def find_missing_reactions(model, filename):
    '''

    :param model: model_object
    :param filename: tab separated file for KBase Import
    :return: list of reactions not correctly added
    '''
    reactions = list()
    result = list()
    with open(filename) as f:
        # SKIPPING HEADER AND BIOMASS
        for line in f.readlines()[2:]:
            elements = line.split('\t')
            reaction = elements[0].split('[')
            try:
                reaction = reaction[0] + '_' + reaction[1][0:2]
            except:
                print reaction
                raise RuntimeError
            print reaction
            reactions.append(reaction)
    for r in reactions:
        if not Model.has_reaction(model, r) and not r.startswith('EX'):
            result.append(r)
    return result

def find_missing_genes(model, filename):
    '''

    :param model: model_object
    :param filename: tab separated file for KBase Import
    :param reactions: list of reactions with no gene in model
    :return: dict of rxns to GPRS in file
    '''
    genes = set()
    file = set()
    for r in Model.get_reactions(model):
        gpr = Gpr(r)
        genes |= gpr.features()
    with open(filename) as f:
        lines = f.readlines()
        file = set([g.split('\'')[1] for g in lines])
    return file - genes

def value_added_analysis():
    """
    Runs an analysis in the value added by the morphing process
    - Runs a venn analysis on source, morph(s), target_recon
    - compiles information on this in the source and morph and not in recon (value added)
    :return: A SimpleTable of information as well as result in dict form
    """
    models = _allmodelsinws(13638)
    my_models = []
    for m in models:
        if m[2] in {'MaripaludisModel', 'Mari_to_Janna_Morph', 'JannashciiReconstruction'}:
            my_models.append(m)
        if m[2] == 'Mari_to_Janna_Morph':
            model = FBAModel(m[0], m[1])
    assert len(my_models) == 3, str(models) + '\n\n' + str(my_models)
    reactions = Analysis.common_reaction_analysis(my_models)[frozenset({'MaripaludisModel', 'Mari_to_Janna_Morph'})]
    rxn_labels = Helpers.load('../data/mari_to_janna_morph.pkl').rxn_labels
    data = Analysis.reaction_analysis(reactions, model, None, rxn_labels)
    total_genes = sum([r[2] for r in data.rows])
    percent_gene_based = format((float(sum([r[1] for r in data.rows])) / len(data.rows)) * 100, '0.2f')
    percent_genes = Plotter.SimpleTable((('Reactions with a genetic basis: ', percent_gene_based)))

    labels = Plotter.SimpleTable(('label', 'percentage'))
    mislabels = Plotter.SimpleTable(data.header)
    label_keys = dict()
    for r in data.rows:
        if r[2] > 0 and r[3] != 'gene-match':
            mislabels.add(r)
        if r[3] not in label_keys:
            label_keys[r[3]] = 1
        else:
            label_keys[r[3]] += 1
    for l in label_keys:
        percentage = format((float(label_keys[l]) / len(data.rows)) * 100, '0.2f')
        labels.add((l, percentage))
    sub_data = dict()
    for l in label_keys:
        sub_data[l] = Plotter.SimpleTable(data.header)
    for r in data.rows:
        sub_data[r[3]].add(r)



    with open('../data/value_added.md', 'w') as f:
        f.write('# Value Added Analysis')
        f.write('This is an autogenerated analysis looking at potential value added in the morphing process. Current' +
                ' implementation looks at genes added, relative labels, reaction classes and subsystems. Future may' +
                ' try to link to online databases and potential literature sources, allowing automated or convenient' +
                ' evaluation of the process')
        f. write(percent_genes.markdown())
        f.write('\n\n')
        f. write(labels.markdown())
        f.write('\n\n')
        f.write('## Mis-labels: \n' + mislabels.markdown())
        f.write('\n\n')
        for d in sub_data:
            my_data = sub_data[d]
            percentage = format((float(label_keys[d]) / len(data.rows)) * 100, '0.2f')
            f.write('\n### ' + d + ' (' + str(percentage) + ')')
            percent_genes = format((float(sum([r[1] for r in my_data.rows])) / len(data.rows)) * 100, '0.2f')
            f.write('##### Percent with genetic basis: ' + str(percent_genes))
            f.write('\n\n' + my_data.markdown())

        f.write(data.markdown())
        




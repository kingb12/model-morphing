import Helpers
from Helpers import *
import Client
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

def five_model_morpher():
    workspaces = (11782, 11783)
    #one by one
    morph = make_morph(ws_id=workspaces[0])
    redo = 0
    while(redo < 5):
        try:
            morph = Client.prepare_supermodel(morph)
            redo = 9
        except ServerError:
            redo += 1
    reaction_list = Client.removal_list(morph.rxn_labels['gene-no-match'])
    morph = Client.process_reactions(morph, rxn_list=reaction_list)
    reaction_list = Client.removal_list(morph.rxn_labels['no-gene'])
    morph = Client.process_reactions(morph, rxn_list=reaction_list)
    dump(morph, '../data/morph-Hsmedianew.pkl')

    #all at once removal
    morph = make_morph(ws_id=workspaces[1])
    while(redo < 5):
        try:
            morph = Client.prepare_supermodel(morph)
            redo = 9
        except ServerError:
            redo += 1
    morph = remove_rxns_by_dict(morph, morph.rxn_labels['gene-no-match'], output_id = 'gnm-removed')
    morph_ng = remove_rxns_by_dict(morph, morph.rxn_labels['no-gene'], output_id = 'gnm-and-ng-removed')

    #various gapfillings
    redo = 0
    while(redo < 5):
        try:
            m2 = Client.probanno_fill(morph, name='gnm-rmv-prob-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m3 = Client.probanno_fill(morph_ng, name='gnm-and-ng-rmv-prob-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m4 = Client.parse_fill(morph, name='gnm-rmv-parse-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m5 = Client.parse_fill(morph_ng, name='gnm-and-ng-rmv-parse-fill')
            redo = 9
        except ServerError:
            redo += 1

    dump(m2, '../data/morph-all1-prob.pkl')
    dump(m3, '../data/morph-all2-prob.pkl')
    dump(m4, '../data/morph-all1-parse.pkl')
    dump(m5, '../data/morph-all2-parse.pkl')

def all_removal():

    #all at once removal
    redo = 0
    morph = make_morph(ws_id=11783)
    while(redo < 5):
        try:
            morph = Client.prepare_supermodel(morph)
            redo = 9
        except ServerError:
            redo += 1
    morph = Client.remove_reactions_by_dict(morph, morph.rxn_labels['gene-no-match'], output_id = 'gnm-removed')
    morph_ng = Client.remove_reactions_by_dict(morph, morph.rxn_labels['no-gene'], output_id = 'gnm-and-ng-removed')

    #various gapfillings
    redo = 0
    while(redo < 5):
        try:
            m2 = Client.probanno_fill(morph, name='gnm-rmv-prob-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m3 = Client.probanno_fill(morph_ng, name='gnm-and-ng-rmv-prob-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m4 = Client.parse_fill(morph, name='gnm-rmv-parse-fill')
            redo = 9
        except ServerError:
            redo += 1
    redo = 0
    while(redo < 5):
        try:
            m5 = Client.parse_fill(morph_ng, name='gnm-and-ng-rmv-parse-fill')
            redo = 9
        except ServerError:
            redo += 1

    dump(m2, '../data/morph-all1-prob.pkl')
    dump(m3, '../data/morph-all2-prob.pkl')
    dump(m4, '../data/morph-all1-parse.pkl')
    dump(m5, '../data/morph-all2-parse.pkl')

def find_kbaliases(filename, genome_id, ws_id):
    """
    given a list of genes, finds their KBase ID aliases and return them as a list
    """
    genes = list()
    with open(filename, "r") as file:
        genes = file.readlines()
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
    for i in range(0, len(genes)):
        # Modify the list in place
        genes[i] = gene_hash[genes[i].split('\n')[0]]
    return genes

def simulate_all_models(ws_id, media_id, phenotype_id):
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
        args['overwrite'] = False
        print args
        Client.fba_client.simulate_phenotypes(args)
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








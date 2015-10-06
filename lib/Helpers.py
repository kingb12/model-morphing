from Morph import Morph
import pickle
import Client
import os
def make_morph():
    args = dict()
    args['genome'] = '3'
    args['src_model'] = '5'
    args['probanno'] = '15'
    args['protcomp'] = '6'
    args['genomews'] = '9145'
    args['src_modelws'] = '9145'
    args['probannows'] = '9145'
    args['protcompws'] = '9145'
    return Morph(args)
def load_morph():
    """
    Loads a morph that has completed
    """
    with open('../data/morph_data.pkl', 'rb') as handle:
          b = pickle.load(handle)
          return b

def modelargs(morph):
    args = dict()
    args['model'] = morph.model
    args['model_workspace'] = morph.ws_id
    return args

def essential_test(morph, rxn_list):
    results = dict()
    for rxn in rxn_list:
        print rxn
        a,b,c = Client._find_alternative(rxn, morph=morph)
        analysis = Client.analyze_alternative(morph, rxn, c)
        results[rxn] = (a, b, c, analysis)
        filename = "essential_results/" + str(rxn) + ".pkl"
        pickle.dump(results, open(filename, "wb"))
        print "saved"
# CURRENT USE SPECIFIC INITIALIZATION CODE
morph = load_morph()

def _unpack_essen_results():
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


import Helpers
from Helpers import *
import Client
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

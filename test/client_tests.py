import unittest
import Client
import Helpers
import Morph
import copy

def refresh_variables(narrative_wsid):
    morph = Helpers.make_morph(ws_id=narrative_wsid)
    supermorph = Client.prepare_supermodel(morph)
    src_filled = copy.deepcopy(morph)
    src_filled.src_modelws = supermorph.src_modelws
    src_filled.src_model = supermorph.src_model
    Helpers.dump(src_filled, '../test/src_filled.pkl')
    Helpers.dump(morph, '../test/morph.pkl')
    Helpers.dump(supermorph, '../test/supermorph.pkl')

test_bank = 12055 #Workspace ID of the test bank
test_space = 12056 #Workspace ID of the test space
fba_client = Client.fba_client
ws_client = Client.ws_client

class ClientTest(unittest.TestCase):

    # Basic Idea, there is a workspace test_bank that stores all the data that
    # will be used in any unit test. All unit tests pull the data they need into
    # the test space via object copying, then perform the functions they need
    # to, and finally clear the test_space workspace. Thus each test gets
    # unmutated data and a 'fresh' workspace to work with

    def setUp(self):
        # pull morphs referencing objects in test bank
        self.src_filled = Helpers.load('../test/src_filled.pkl')
        self.morph = Helpers.load('../test/morph.pkl')
        self.supermorph = Helpers.load('../test/supermorph.pkl')
        # copy objects to the test space (references fixed in helper)
        self.morph = _copy_morph_to_ws(self.morph, test_space)
        self.supermorph =  _copy_morph_to_ws(self.supermorph, test_space)
        self.src_filled = _copy_morph_to_ws(self.src_filled, test_space)
        # get objects
        self.srcobj = Helpers.get_object(self.morph.src_model, self.morph.src_modelws)
        self.superobj = Helpers.get_object(self.supermorph.model, self.supermorph.ws_id)
        self.srcfillobj = Helpers.get_object(self.src_filled.src_model, self.src_filled.src_modelws)

    def tearDown(self):
        _clear_testspace()
        self.src_filled = None
        self.morph = None
        self.supermorph = None
        self.srcobj = None
        self.superobj = None
        self.srcfillobj = None

    def test_translate_features(self):
    # assert trans is a subset of original
        m2 = Client.translate_features(self.morph)
        m2obj = Helpers.get_object(m2.trans_model, m2.ws_id)
        trans_rxns = get_reaction_dict(m2obj)
        src_rxns = get_reaction_dict(self.srcobj)
        genome = Helpers.get_object(m2.genome, m2.genomews)['data']
    # Sanity Checks
        self.assertTrue(m2.ws_id == test_space, msg='test architecture')
        self.assertFalse(m2 is self.morph, msg='aliasing issues')
        self.assertTrue(len(trans_rxns) > 0, msg='there are trans rxns')
    # Translation Definitions and Tests
        trans_set = set(trans_rxns.keys())
        src_set = set(src_rxns.keys())
        self.assertTrue(trans_set.issubset(src_set), msg='translation subset defn')
        # make sure all translation rxns that have genes reference correct
        # genome (no old features)
        # also make sure all trans w/ gene's have gene's in source
        for rxn in trans_rxns:
            prots = trans_rxns[rxn]['modelReactionProteins']
            if len(prots) > 0:
                for p in prots:
                    subunits = p['modelReactionProteinSubunits']
                    for s in subunits:
                        for ftr in s['feature_refs']:
                            self.assertTrue(ftr.split('/')[-1].startswith(genome['id']), msg='incorrect feature ref')
                            self.assertTrue(_has_gene(src_rxns[rxn]),
                                msg='source has no gene for ' + str(rxn) + 'in trans')
        # make sure all removed reactions had genes in source (otherwise, we've
                                # removed a no gene rxn)
        gnm_rxns = src_set - trans_set
        for rxn in gnm_rxns:
            for p in src_rxns[rxn]['modelReactionProteins']:
                self.assertTrue(len(p['modelReactionProteinSubunits']) > 0, msg='removed a rxn w/o a gene in source')

    def test_label_reactions(self):
        m2 = Client.translate_features(self.morph)
        m2obj = Helpers.get_object(m2.trans_model, m2.ws_id)
        trans_rxns = get_reaction_dict(m2obj)
        src_rxns = get_reaction_dict(self.srcobj)
        genome = Helpers.get_object(m2.genome, m2.genomews)['data']
        m4 = Client.reconstruct_genome(m2)
        m3 = Client.label_reactions(m4)
        m3obj = Helpers.get_object(m3.trans_model, m3.ws_id)
        m3obj_recon = Helpers.get_object(m3.recon_model, m3.ws_id)
        trans_rxns = get_reaction_dict(m3obj)
        recon_rxns = get_reaction_dict(m3obj_recon)
        labelnames = {'gene-no-match', 'gene-match', 'no-gene', 'recon'}
        labels = dict()
        for key in m3.rxn_labels:
            labels[key] = set(m3.rxn_labels[key].keys())
        trans_set = set(trans_rxns.keys())
        src_set = set(src_rxns.keys())
        recon_set = set(recon_rxns.keys())
    #Sanity checks
        self.assertTrue(m2.ws_id == test_space, msg='test architecture')
        self.assertFalse(m2 is self.morph, msg='aliasing issues')
        self.assertTrue(len(trans_rxns) > 0, msg='there are trans rxns')
        self.assertTrue(set(labels.keys()) == labelnames)
    #Mechanism Checks
        self.assertTrue(labels['recon'] == recon_set - src_set, msg='incorrect recon mechanism')
        self.assertTrue(labels['gene-no-match'] == src_set - trans_set, msg='incorrect gene-no-match mechanism')
        self.assertTrue(labels['gene-match'] | labels['no-gene'] == src_set - labels['gene-no-match'],
                        msg='incorrect genematch/nogene mechanism')
        self.assertTrue(labels['recon'].isdisjoint(src_set), msg='recon reaction in source')
        items = labels.items()
        for i in range(0, len(items)):
            for j in range(0, len(items)):
                if i != j:
                    self.assertTrue(items[i][1].isdisjoint(items[j][1]), 'a reaction is in two label sets')


    #gnm
    #ng
    #gm
    #

def _has_gene(rxn_struct):
    return len(rxn_struct['modelReactionProteins'][0]['modelReactionProteinSubunits']) > 0

def _copy_morph_to_ws(morph, ws_id):
    orig_ws = morph.ws_id
    morph.ws_id = ws_id
    if morph.model is not None:
        info = Helpers.copy_object(morph.model, orig_ws, ws_id)
        morph.model = info[0]
        morph.ws_id = info[6]
    if morph.src_model is not None:
        info = Helpers.copy_object(morph.src_model, morph.src_modelws, ws_id)
        morph.src_model = info[0]
        morph.src_modelws = info[6]
    if morph.trans_model is not None:
        info = Helpers.copy_object(morph.trans_model, orig_ws, ws_id)
        morph.trans_model = info[0]
        morph.ws_id = info[6]
    return morph

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
def _clear_testspace(ws_id=test_space):
    # get all object ids in {'id': x} form excluding the narrative object
    object_ids = [{'objid': info[0], 'wsid': ws_id} for info in
                 ws_client.list_objects({'ids': [ws_id]}) if not info[2].startswith('KBaseNarrative')]
    if len(object_ids) > 0:
        ws_client.delete_objects(object_ids)

if __name__ == '__main__':
        unittest.main()

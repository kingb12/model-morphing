
import unittest
import lib.Client as Client
import lib.Helpers as Helpers
import lib.Morph as Morph
import copy
from lib.objects import *
import lib.service as service

def refresh_variables(narrative_wsid):
    morph = Helpers.make_morph(ws_id=narrative_wsid)
    supermorph = Client.prepare_supermodel(morph)
    src_filled = copy.deepcopy(morph)
    src_filled.src_modelws = supermorph.src_modelws
    src_filled.src_model = supermorph.src_model
    Helpers.dump(src_filled, '../test/src_filled.pkl')
    Helpers.dump(morph, '../test/morph.pkl')
    Helpers.dump(supermorph, '../test/supermorph.pkl')

test_bank = 12055  # Workspace ID of the test bank
test_space = 12056  # Workspace ID of the test space


class GprTest(unittest.TestCase):

    # Basic Idea, there is a workspace test_bank that stores all the data that
    # will be used in any unit test. All unit tests pull the data they need into
    # the test space via object copying, then perform the functions they need
    # to, and finally clear the test_space workspace. Thus each test gets
    # unmutated data and a 'fresh' workspace to work with

    def __init__(self, *args, **kwargs):
        super(GprTest, self).__init__(*args, **kwargs)
        morph = Helpers.make_morph(12056)
        morph = Client.translate_features(morph)
        morph = Client.reconstruct_genome(morph)
        morph = Client.label_reactions(morph)
        self.recon = morph.objects['recon_model']
        self.trans = morph.objects['trans_model']

    def setUp(self):
        self.recon_gprs = [(Morph.Gpr(r), r) for r in self.recon['modelreactions']]
        self.trans_gprs = [(Morph.Gpr(r), r) for r in self.trans['modelreactions']]

    def tearDown(self):
        _clear_testspace()
        self.src_filled = None
        self.morph = None
        self.supermorph = None
        self.srcobj = None
        self.superobj = None
        self.srcfillobj = None

    def test_gpr_structure(self):
        '''
        verifies that the gpr is structurally sound
        '''
        def test(g):
            structure = False
            for protein in g[0]:
                self.assertEqual(type(protein), frozenset, msg=str(protein) + ' is a protein w/ wrong type')
                for subunit in protein:
                    self.assertEqual(type(subunit), frozenset, msg=str(subunit) + ' is a subunit w/ wrong type')
                    structure = True
                    for feature in subunit:
                        self.assertTrue(feature is not None, msg='feature in ' + str(subunit) + ' is None')
            self.assertTrue(structure, msg='could not iterate over sets to features ' + repr(g))
        for g in self.recon_gprs:
            test(g)
        for g in self.trans_gprs:
            test(g)
    def test_merge(self):
        # the features of each should always be in merged, and mergin should be
        # symmetric
        recon_dict = dict()
        trans_dict = dict()
        for g in self.recon_gprs:
            rxn_id = Helpers.get_rxn_id(g[1])
            recon_dict[rxn_id] = g[0]
        for g in self.trans_gprs:
            rxn_id = Helpers.get_rxn_id(g[1])
            trans_dict[rxn_id] = g[0]
        #TODO: Below code simply ignores rxn00000, a KBase dependency
        for r in recon_dict:
            if r in trans_dict and r is not 'rxn00000':
                merge1 = trans_dict[r].merge(recon_dict[r])
                merge2 = recon_dict[r].merge(trans_dict[r])
                self.assertTrue(merge1 == merge2, msg=str(r) + ' merge ' + repr(merge1) +  'Not symetrical to merge ' + repr(merge2))
                for ftr in recon_dict[r].ftrs:
                    self.assertTrue(ftr in merge1.ftrs and ftr in merge2.ftrs, msg=str(r) + 'merge: ' +  str(ftr) + ' in arguments but not in merge')
                for ftr in trans_dict[r].ftrs:
                    self.assertTrue(ftr in merge1.ftrs and ftr in merge2.ftrs, msg=str(r) + 'merge: ' +  str(ftr) + ' in arguments but not in merge')
                # the proteins of each if they have no features in the other should
                # always be in both
                if len(recon_dict[r].ftrs) > 0 and len(trans_dict[r].ftrs):
                    for protein in recon_dict[r]:
                        if len(recon_dict[r]._unnest_sets(protein) & trans_dict[r].ftrs) == 0:
                            self.assertTrue(merge1.contains_protein(protein), msg=str(protein) + ' has no common features with ' + repr(merge1) + ' but is not in merge')
                    for protein in trans_dict[r]:
                        if len(trans_dict[r]._unnest_sets(protein) & recon_dict[r].ftrs) == 0:
                            self.assertTrue(merge1.contains_protein(protein), msg=str(protein) + ' has no common features with ' + repr(merge1) + ' but is not in merge')
        # prime examples
        a = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb1', 'kb2'])])]))
        b = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb2'])])]))
        c = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb1']), frozenset(['kb2'])])]))
        d = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb1', 'kb3'])])]))
        e = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb4'])])]))
        f = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb4', 'kb2'])])]))
        g = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb1', 'kb2']), frozenset(['kb4'])])]))
        h = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb2'])]), frozenset([frozenset(['kb4'])])]))
        i = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb1', 'kb2'])]), frozenset([frozenset(['kb1']), frozenset(['kb2'])])]))
        j = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb1', 'kb2'])]), frozenset([frozenset(['kb3'])])]))
        k = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb1', 'kb3', 'kb4'])]), frozenset([frozenset(['kb2'])])]))
        l = Morph.Gpr().new_gpr(frozenset([frozenset([frozenset(['kb1', 'kb2', 'kb3', 'kb4'])])]))
        #Homolog cases
        result = a.merge(b)
        self.assertEqual(result, a, msg='failure on adding one homolog T->R')
        result = b.merge(a)
        self.assertEqual(result, a, msg='failure on adding one homolog R->T')
        result = a.merge(d)
        self.assertEqual(result.gpr, frozenset([frozenset([frozenset(['kb1', 'kb2', 'kb3'])])]), msg= repr(result) + " fails merge gpr: " + repr(a) + '\n' + repr(d))
        result = a.merge(a)
        self.assertEqual(result, a, msg='failure on identity case ' + repr(a)  + '\n' + repr(result))
        #Subunit Disagreement cases
        # case you discovered yesterday
        self.assertEqual(j.merge(k), l)
        self.assertEqual(k.merge(j), l)

        # classic case
        self.assertEqual(a.merge(c) , a, msg= repr(a.merge(c)) + '\n' + repr(a))

        #Isofunctional cases
        self.assertEqual(b.merge(e), h)


class StoredObjectTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(StoredObjectTest, self).__init__(*args, **kwargs)
        _clear_testspace(test_space)

    def setUp(self):
        self.test_objects = service.list_objects(test_bank)

    def tearDown(self):
        _clear_testspace()

    def test_overwrite(self):
        obj = StoredObject(self.test_objects[0])
        try:
            obj.object_id = "This should fail"
            self.fail(msg="StoredObject permitted a mutation")
        except MutationError:
            pass

    def test_get_object(self):
        for info in self.test_objects:
            obj = StoredObject(info)
            object = obj.get_object()
            self.assertTrue(object is not None, msg="getting object " + str(info) + " failed")









def _clear_testspace(ws_id=test_space):
    service.clear_workspace(ws_id)

if __name__ == '__main__':
        unittest.main()


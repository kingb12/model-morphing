import unittest
import lib.Client as Client
import lib.Helpers as Helpers
import lib.Morph as Morph
import copy
from lib.objects import *
import lib.service as service


def refresh_variables(narrative_wsid):
    # morph = Helpers.make_morph(ws_id=narrative_wsid)
    # supermorph = Client.prepare_supermodel(morph)
    # src_filled = copy.deepcopy(morph)
    # src_filled.src_modelws = supermorph.src_modelws
    # src_filled.src_model = supermorph.src_model
    # Helpers.dump(src_filled, '../test/src_filled.pkl')
    # Helpers.dump(morph, '../test/morph.pkl')
    # Helpers.dump(supermorph, '../test/supermorph.pkl')
    pass


test_bank = 12055  # Workspace ID of the test bank
test_space = 12056  # Workspace ID of the test space


class StoredObjectTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(StoredObjectTest, self).__init__(*args, **kwargs)
        _clear_testspace(test_space)

    def setUp(self):
        self.test_objects = service.list_objects(test_bank)

    def tearDown(self):
        _clear_testspace()

    def test_overwrite(self):
        t = self.test_objects[0]
        obj = StoredObject(t[0], t[1])
        try:
            obj.object_id = "This should fail"
            self.fail(msg="StoredObject permitted a mutation")
        except MutationError:
            pass

    def test_get_object(self):
        for info in self.test_objects:
            obj = StoredObject(info[0], info[1])
            obj = obj.get_object()
            self.assertTrue(obj is not None, msg="getting object " + str(info) + " failed")

    def test_implicit_get(self):
        for info in self.test_objects:
            obj = StoredObject(info[0], info[1])
            obj = obj.data
            self.assertTrue(obj is not None, msg="getting object " + str(info) + " failed")

    def test_save_object(self):
        for info in self.test_objects[1:4]:
            obj = StoredObject(info[0], info[1])
            obj = obj.get_object()
            a = StoredObject.save(obj, test_space, typestr=info[3])
            self.assertTrue(a.get_object() == obj, msg="improper save")


class FBAModelTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(FBAModelTests, self).__init__(*args, **kwargs)
        _clear_testspace(test_space)

    def setUp(self):
        self.test_objects = service.list_objects(test_bank, typestr='FBAModel')

    def tearDown(self):
        _clear_testspace()

    def test_get_reactions(self):
        for info in self.test_objects:
            model = FBAModel(info[0], info[1])
            reactions = model.get_reactions()
            self.assertTrue(len(reactions) == len(model.data['modelreactions']), msg="got more/less rxns than expected")
            for r in reactions:
                self.assertTrue(isinstance(r, ModelReaction), msg="got a non model reaction")
                self.assertTrue(type(r.rxn_id()) == str or type(r.rxn_id()) == unicode, msg="improper rxn id " +
                                                                                            str(r.rxn_id()))
                self.assertTrue(r.get_equation() is not None, msg="failure on ModelReaction.get_equation()")
                self.assertTrue(isinstance(r.gpr, Gpr), msg="bad gpr type " + str(type(r.gpr)))


class ModelReactionTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(ModelReactionTests, self).__init__(*args, **kwargs)
        _clear_testspace(test_space)

    def setUp(self):
        self.data_normal = Helpers.load('ModelReactionTests/test_data_normal.pkl')
        self.data_no_id = Helpers.load('ModelReactionTests/test_data_no_id.pkl')

    def tearDown(self):
        _clear_testspace()

    def test_normal(self):
        rxn = ModelReaction(self.data_normal)
        self.assertTrue(rxn.rxn_id() == u'rxn02483_c0', msg="Failure on normal rxn_id()")
        self.assertTrue(str(rxn.get_equation()) ==
                        '[-1*(cpd00938_c0), -1*(cpd00067_c0), 1*(cpd00011_c0), 1*(cpd02255_c0)]',
                        msg="failure on get equation of normal ModelReaction:\n" + str(rxn.get_equation()))
        self.assertTrue(rxn.get_removal_id() == u'rxn02483-c0-_c0', msg='incorrect normal removal id ' +
                                                                        str(rxn.get_removal_id()))
        self.assertTrue(rxn.get_biochem_ref() == FBAModel.DEFAULT_BIOCHEM, msg="Got: " +
                                                                               str(rxn.get_biochem_ref().identity) +
                                                                               'Expected: ' +
                                                                               str(FBAModel.DEFAULT_BIOCHEM.identity))
        self.assertTrue(rxn.get_rxn_ref() == 'rxn02483')

    def test_no_id(self):
        rxn = ModelReaction(self.data_no_id)
        self.assertTrue(rxn.rxn_id() == u'rxn07191-c0-_c0', msg="Failure on no id rxn_id() " + rxn.rxn_id())
        self.assertTrue(str(rxn.get_equation() ==
                            '[-1*(cpd00102_c0), 1*(cpd00169_c0), -1*(cpd00001_c0), 3*(cpd00067_c0), -1*(cpd11621_c0), 1*(cpd11620_c0)]'),
                        msg="failure on get equation of no id ModelReaction\n" + str(rxn.get_equation()))
        self.assertTrue(rxn.get_removal_id() == u'rxn07191-c0-_c0', msg='incorrect no_id removal id' +
                                                                        str(rxn.get_removal_id()))
        self.assertTrue(rxn.get_biochem_ref() == FBAModel.DEFAULT_BIOCHEM, msg="Got: " +
                                                                               str(rxn.get_biochem_ref().identity) +
                                                                               'Expected: ' +
                                                                               str(FBAModel.DEFAULT_BIOCHEM.identity))
        self.assertTrue(rxn.get_rxn_ref() == 'rxn00000', msg="bad rxn ref " + str(rxn.get_rxn_ref()))


class CompoundTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(CompoundTests, self).__init__(*args, **kwargs)
        _clear_testspace(test_space)

    def setUp(self):
        self.data_normal = Helpers.load('ModelReactionTests/test_data_normal.pkl')
        self.data_no_id = Helpers.load('ModelReactionTests/test_data_no_id.pkl')

    def tearDown(self):
        _clear_testspace()

    def test_normal(self):
        rxn = ModelReaction(self.data_normal)
        for comp in rxn.get_equation()[0:2]:
            self.assertTrue(isinstance(comp, Compound))
            self.assertTrue(comp.biochem == rxn.get_biochem_ref())
            self.assertTrue(type(comp.coeff) == int)
            self.assertTrue(comp.compound_id.startswith('cpd'))
            self.assertTrue(comp.get_info() is not None, msg=comp.compound_id) # Will Raise an Error if it fails
            self.assertTrue(type(comp.formula()) in [str, unicode], msg=str(type(comp.formula())))


def _clear_testspace(ws_id=test_space):
    service.clear_workspace(ws_id)


if __name__ == '__main__':
    unittest.main()

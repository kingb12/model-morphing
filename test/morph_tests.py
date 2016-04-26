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


# class StructureTests(unittest.TestCase):
# pass


class PrepareSupermodelTests(unittest.TestCase):
    DEBUG = False

    def __init__(self, *args, **kwargs):
        super(PrepareSupermodelTests, self).__init__(*args, **kwargs)
        _clear_testspace(test_space)

    def setUp(self):
        self.mj = Helpers.mari_to_janna(test_space)
        self.ab = Helpers.ace_to_bark(test_space)

    def tearDown(self):
        _clear_testspace()

    def test_fill_source(self):
        """
        WARNING: This Test will take a long time to run. Set module variable to turn off gapfilling tests
        :return:
        """
        if not PrepareSupermodelTests.DEBUG:
            return  # SKIP THIS TEST
        morph = self.ab
        prevSource = morph.src_model
        self._type_check(morph.src_model, FBAModel, 'morph.src_model(pre)')
        morph.fill_src_to_media()
        self.assertTrue(morph.src_model != prevSource, msg='No Change')
        self._type_check(morph.src_model, FBAModel, 'morph.src_model(post)')
        self.assertTrue(morph.runfba(model=morph.src_model).objective > 0, msg='gapfill did not lead to growth')
        self.assertTrue(len(morph.log.actions) == 2, msg='incorrect logging')

    def test_translate_features(self):
        morph = self.mj
        self.assertTrue(morph.trans_model is None, msg='Improper setup, morph.trans_model is not None: \n' +
                                                       repr(morph.trans_model))
        morph.translate_features()
        self._type_check(morph.trans_model, FBAModel, 'morph.trans_model')

        # Testing guarantees of translation
        genome_id = morph.genome.get_genome_id()
        source_rxns = dict([(r.rxn_id(), r) for r in morph.src_model.get_reactions()])
        trans_rxns = dict([(r.rxn_id(), r) for r in morph.trans_model.get_reactions()])
        for r in trans_rxns:
            rxn = trans_rxns[r]
            self._check_feature_reference(genome_id, rxn)
            self.assertTrue(rxn.rxn_id() in source_rxns, msg='reaction: ' + rxn.rxn_id() + ' not in source reactions')
            if len(rxn.gpr.features()) > 0:
                self.assertTrue(len(source_rxns[r].gpr.features()) > 0, msg='a rxn:' + str(r) +
                                                                            ' in trans w/ gene has no-gene in source')
        for r in source_rxns:
            if r not in trans_rxns:
                self.assertTrue(len(source_rxns[r].gpr.features()) > 0, msg='removed a no-gene rxn: ' + str(r))

    def test_reconstruct_genome(self):
        morph = self.mj
        # Set up Check
        self.assertTrue(morph.recon_model is None, msg="improper test setup, morph.recon_model is not None")
        morph.reconstruct_genome()
        # Type Check
        self._type_check(morph.recon_model, FBAModel, 'morph.recon_model')
        recon_rxns = dict([(r.rxn_id(), r) for r in morph.recon_model.get_reactions()])
        genome_id = morph.genome.get_genome_id()
        for r in recon_rxns:
            rxn = recon_rxns[r]
            self._check_feature_reference(genome_id, rxn)

    def test_label_reactions(self):
        morph = self.mj
        # ======= Test specific set up ======== #
        morph.reconstruct_genome()
        morph.translate_features()
        morph.label_reactions()
        # ===================================== #
        source_rxns = dict([(r.rxn_id(), r) for r in morph.src_model.get_reactions()])
        trans_rxns = dict([(r.rxn_id(), r) for r in morph.trans_model.get_reactions()])
        recon_rxns = dict([(r.rxn_id(), r) for r in morph.recon_model.get_reactions()])
        labelnames = {'gene-no-match', 'gene-match', 'no-gene', 'recon', 'common'}

        # Extract the reaction_id's only
        labels = dict()
        for key in morph.rxn_labels:
            if key != 'common':
                labels[key] = set(morph.rxn_labels[key].keys())
        trans_set = set(trans_rxns.keys())
        src_set = set(source_rxns.keys())
        recon_set = set(recon_rxns.keys())
        all_rxns = copy.deepcopy(recon_rxns).update(source_rxns)

        # Definition Checks
        # self.assertTrue(labels['gene-no-match'] == src_set - trans_set - recon_set,
        #                msg='incorrect gene-no-match mechanism')
        items = labels.items()
        for i in range(0, len(items)):
            for j in range(0, len(items)):
                if i != j:
                    self.assertTrue(items[i][1].isdisjoint(items[j][1]), 'a reaction is in two label sets')
        self.assertTrue(labels['gene-no-match'].isdisjoint(trans_set), msg='a gnm reaction is in the trans model')
        basis_rxns = set()
        for rxn in source_rxns:
            if self._has_gene(source_rxns[rxn]) or source_rxns[rxn].gpr.gpr_type in {u'universal', u'spontaneous'}:
                basis_rxns.add(rxn)
        for rxn in recon_rxns:
            basis_rxns.add(rxn)
        self.assertTrue(len(basis_rxns) > 0, msg='NONE of the reactions have a basis')
        self.assertTrue(labels['gene-no-match'].isdisjoint(trans_rxns),
                        msg='a gnm reaction is in the trans model')
        self.assertTrue(labels['gene-no-match'].issubset(basis_rxns),
                        msg='there is a rxn in gnm but has no gene in source')
        self.assertTrue(labels['gene-match'].issubset(basis_rxns),
                        msg='there is a rxn in gene match without a gene in source\n' +
                            str(labels['gene-match'] - basis_rxns) + '\n\n' +
                            self._rxn_info_str(labels['gene-match'] - basis_rxns, lookup=trans_rxns))
        self.assertTrue(labels['gene-match'].isdisjoint(labels['recon']),
                        msg='there is a rxn in gene match without a gene in source\n' +
                            str(labels['gene-match'].intersection(labels['recon'])))
        self.assertTrue(labels['gene-match'].issubset(trans_set.union(recon_set)),
                        msg='there is a rxn in gene match not in tanslation')
        # no gene definition
        self.assertTrue(labels['no-gene'].isdisjoint(basis_rxns), msg='a no gene rxn has a gene in source' +
                                                                      self._rxn_info_str(
                                                                          labels['no-gene'].intersection(basis_rxns),
                                                                          lookup=source_rxns))
        self.assertTrue(labels['no-gene'].issubset(trans_set), msg='a no gene rxn is not in translation')
        self.assertTrue(labels['no-gene'].issubset(src_set), msg='a no gene rxn is not in source')

    def test_build_supermodel(self):
        morph = self.mj
        # ================= SETUP ================== #
        morph.reconstruct_genome()
        morph.translate_features()
        morph.label_reactions()
        morph.build_supermodel()
        # ========================================== #
        source_rxns = dict([(r.rxn_id(), r) for r in morph.src_model.get_reactions()])
        trans_rxns = dict([(r.rxn_id(), r) for r in morph.trans_model.get_reactions()])
        recon_rxns = dict([(r.rxn_id(), r) for r in morph.recon_model.get_reactions()])
        super_rxns = dict([(r.rxn_id(), r) for r in morph.model.get_reactions()])
        labels = dict()
        for key in morph.rxn_labels:
            labels[key] = set(morph.rxn_labels[key].keys())
        trans_set = set(trans_rxns.keys())
        src_set = set(source_rxns.keys())
        recon_set = set(recon_rxns.keys())
        super_set = set(super_rxns.keys())
        # Mechanism Checks
        self.assertTrue(super_set == recon_set | trans_set | labels['gene-no-match'], msg='supermodel mechanism off')
        self.assertTrue(super_set == recon_set | src_set, msg='supermodel mechanism off (definition of trans)')
        # Definition Checks
        for rxn in recon_rxns:
            self.assertTrue(rxn in super_rxns, msg='recon rxn: ' + str(rxn) +
                                                   ' not in super (redundant check, prevent keyerror)')
        for rxn in trans_rxns:
            self.assertTrue(rxn in super_rxns, msg='tran rxn: ' + str(rxn) +
                                                   ' not in super (redundant check, prevent keyerror)')
        for rxn in source_rxns:
            self.assertTrue(rxn in super_rxns, msg='source rxn: ' + str(rxn) +
                                                   ' not in super (redundant check, prevent keyerror)')
        all_rxns = set()
        for key in morph.rxn_labels:
            for r in morph.rxn_labels[key]:
                all_rxns.add(r)
        for rxn in all_rxns:
            self.assertTrue(rxn in super_rxns)
            if rxn in morph.rxn_labels['common']:
                reaction = super_rxns[rxn]
                trans_gpr = morph.rxn_labels['common'][rxn][0].gpr
                recon_gpr = morph.rxn_labels['common'][rxn][1].gpr
                merge_gpr = trans_gpr.merge(recon_gpr)
                super_gpr = reaction.gpr
                self.assertTrue(merge_gpr == super_gpr, msg='rxn: ' + str(rxn) + '\nmerge gpr:\n' + str(merge_gpr) +
                                                            '\nsuper gpr:\n' + str(super_gpr))
            # Failing Because KBase
            elif rxn in trans_rxns:
                self.assertTrue(super_rxns[rxn].gpr == trans_rxns[rxn].gpr, msg=str(rxn) + ' has different GPRS\n' +
                                                                                'super: ' + str(
                    super_rxns[rxn].gpr) + '\ntrans: ' + str(trans_rxns[rxn].gpr))
            elif rxn in recon_rxns:
                self.assertTrue(super_rxns[rxn].gpr == recon_rxns[rxn].gpr, msg=str(rxn) + ' has different GPRS\n' +
                                                                                'super: ' + str(
                    super_rxns[rxn].gpr) + '\nrecon: ' + str(recon_rxns[rxn].gpr))
                # Fin

    # -------------------------------- HELPER FUNCTIONS AND TESTS --------------------------------------------------- #
    def _type_check(self, obj, type_class, obj_name):
        self.assertTrue(isinstance(obj, type_class), msg='Illegal type for ' + obj_name + ' type: ' +
                                                         str(type(obj)) + '\ndata: \n' +
                                                         repr(obj))

    def _check_feature_reference(self, genome_id, model_reaction):
        # Make sure all features in a set of reactions start with a particular genome_id
        for ftr in model_reaction.gpr.features():
            self.assertTrue(ftr.startswith(genome_id), msg='feature: ' + str(ftr) + ' doesnt start with genome: ' +
                                                           str(genome_id) + ' in rxn: ' +
                                                           str(model_reaction.rxn_id()))

    def _has_gene(self, model_reaction):
        return len(model_reaction.gpr.features()) > 0

    def _rxn_info_str(self, reactions, lookup=None):
        result = ""
        if lookup is not None:
            reactions = [lookup[r] for r in reactions]
        for r in reactions:
            result += 'id: ' + str(r.rxn_id()) + ' gpr_type: ' + str(r.gpr.gpr_type)
            result += ' verify gpr_type: \n\t' + str(r.data['pathway']) + '\n'
            result += ' stored id: ' + r.data['id'] + '\n'
        return result


def _clear_testspace(ws_id=test_space):
    service.clear_workspace(ws_id)


if __name__ == '__main__':
    unittest.main()

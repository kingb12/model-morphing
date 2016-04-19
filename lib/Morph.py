# import necesary services
from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import random
import copy
import json
from operator import itemgetter

from lib.log import Log
from lib.objects import *


def _init_clients():
    global ws_client, args
    global fba_client
    global ws_id
    # Get workspace service URL parameter
    with open("./urls/.kbase_workspaceURL", "r") as myfile:
        url = myfile.read().replace('\n', '')
    ws_client = Workspace(url)
    # Get FBA Model Services URL parameter
    with open("./urls/.kbase_fbaModelServicesURL", "r") as myfile:
        url = myfile.read().replace('\n', '')
    fba_client = fbaModelServices(url)
    return ws_client, fba_client


def _init_workspace(ws=None):
    ws_id = ws
    ws_name = 'MMws'
    if ws is None:
        ws_conflict = True
        while ws_conflict:
            create_ws_params = {'workspace': ws_name, 'globalread': 'r', 'description':
                "A workspace for storing the FBA's and meta data of the algorithm"}
            # Try to create a workspace, catch an error if the name is already in use
            try:
                new_ws = ws_client.create_workspace(create_ws_params)
                # new_ws is type workspace_info, a tuple where 0, 1 are id, name
                ws_id = new_ws[0]
                ws_name = new_ws[1]
                ws_conflict = False
            except ServerError:
                ws_name += str(random.randint(1, 9))
    return ws_id, ws_name


class Morph:
    # These are the allowed properties of a morph object. Values not specified
    # by user are set to None
    # This is an Abstract Object representing the Morph of a metabolic model
    # from one species to a close relative genome. It has information related to
    # the source model, the target genome, the reactions in the model as it is
    # morphed from source import to target, and the set of
    properties = {'src_model', 'genome', 'probanno', 'protcomp', 'model', 'rxn_labels', 'ws_id', 'ws_name',
                  'trans_model',
                  'recon_model', 'media', 'probhash', 'log'}

    def __init__(self, *arg_hash, **kwargs):
        """
        Initiializes a Morph! A reasonable Morph needs each of the following set:
            src_model: (FBAModel) the model of the source organism
            genome: (Genome) the genome of the target organism
            probanno: (ReactionProbabilities) the likelihoods for KBase reactions in the target organism
            protcomp: (ProteomeComparison) the proteome comparison between source and target genomes
            media: (Media) the media for the source model and subsequent models
        :param arg_hash: a dictionary with the attributes desired as keys
        :param kwargs: keyword arguments
        :return:
        """
        # TODO: EDIT THIS DOC LINE WITH THE TRUTH
        self.src_model = None
        self.model = None
        self.genome = None
        self.probanno = None
        self.protcomp = None
        self.model = None
        self.rxn_labels = None
        self.media = None
        self.recon_model = None
        self.trans_model = None
        self.probhash = None
        self.log = None
        self.ws_id = None
        self.ws_name = None

        for dictionary in arg_hash:
            for key in dictionary:
                if key in Morph.properties:
                    setattr(self, key, dictionary[key])
        for key in kwargs:
            if key in Morph.properties:
                setattr(self, key, kwargs[key])
        for prop in Morph.properties:
            if not hasattr(self, prop):
                setattr(self, prop, None)
        if self.log is None:
            self.log = Log(self)
        if self.ws_id is None:
            self.ws_id, self.ws_name = _init_workspace()
        self._check_rep()

    def _check_rep(self):
        if self.model is not None:
            assert isinstance(self.model, FBAModel)
        if self.src_model is not None:
            assert isinstance(self.src_model, FBAModel)
        if self.genome is not None:
            assert isinstance(self.genome, Genome)
        if self.media is not None:
            assert isinstance(self.media, Media)
        if self.probanno is not None:
            assert isinstance(self.probanno, ReactionProbabilities)
        if self.protcomp is not None:
            assert isinstance(self.protcomp, ProteomeComparison)
        if self.trans_model is not None:
            assert isinstance(self.trans_model, FBAModel)
        if self.recon_model is not None:
            assert isinstance(self.recon_model, FBAModel)

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True)

    # Overridden Functions to produce unique output
    def __str__(self):
        output = ''
        for key in vars(self):
            attr = getattr(self, key)
            if isinstance(attr, dict):
                attr = attr.keys()
                if len(attr) < 100:
                    output += str(key) + ': ' + str(attr) + '\n'
                else:
                    output += str(key) + ': ' + str(attr[0:100]) + ' ... (more)\n'
            else:
                output += str(key) + ': ' + str(attr) + '\n'
        return output

    def __repr__(self):
        return str(self)

    def __unicode__(self):
        return unicode(str(self))

    def fill_src_to_media(self):
        """
        Gap-fills the Morphs source model to the provided media
        :return: None
        """
        prev = FBAModel(self.src_model.object_id, self.src_model.workspace_id)
        self.src_model = self.src_model.copy(workspace_id=self.ws_id)
        result = service.gapfill_model(self.src_model, self.media,
                                       workspace=self.ws_id, rxn_probs=self.probanno)
        self.src_model = FBAModel(result[0], result[1])
        self.log.add('gapfill', prev, self.src_model, context='fill_src_to_media')

    def runfba(self, model=None, media=None):
        """
        Run FBA on the model in the Morph
        :param model: (optional) FBAModel, default is morph.model
        :param media: (optional) Media, default is morph.media
        :return: FBA
        """
        if model is None:
            model = self.model
        if media is None:
            media = self.media
        objid, wsid = service.runfba(model, media, self.ws_id)
        return FBA(objid, wsid)

    def translate_features(self):
        """
        Translate morph.model using ProteomeComparison (morph.procomp) to a translated model
        :return: None (sets morph.trans_model)
        """
        prev = copy.deepcopy(self.trans_model)
        result = service.translate_model(self.src_model, self.protcomp, workspace=self.ws_id)
        self.trans_model = FBAModel(result[0], result[1])
        self.log.add('model_translation', prev, self.trans_model, context='translate_features')

    def reconstruct_genome(self):
        """
        Reconstruct the genome using automated reconstruction in the service form morph.genome
        :return: None (sets morph.recon_model)
        """
        prev = copy.deepcopy(self.recon_model)
        result = service.reconstruct_genome(self.genome, workspace=self.ws_id)
        self.recon_model = FBAModel(result[0], result[1])
        self.log.add('genome_reconstruction', prev, self.recon_model, context='reconstruct_genome')

    def label_reactions(self):
        """
    Labels morph's reactions from translated model, reconstruction, and source

    Populates the rxn_labels attribute in the Morph object with a Dictionary of four dictionaries of
    reaction_id -> value tuples.
    The first level dicts are named with the keys:
        - gene-match
        - gene-no-match
        - no-gene
        - recon
    Populates morph.probhash with a dictionary of compartment-truncated reaction_ids to their probability relative to
    genome in question (derived from morph.probanno).
    Also populates morph.objects, morph.info with data from KBase objects (advanced use)

    Technical Details (Consult if you are encountering issues):
    For simplicity, an end user can treat the interiors of each dictionary like a set of reaction ids, but it is
    important to note the this is actually a nested dictionary with entries of the form:
    reaction_id -> (model_index, probability)
    Thus, some set behavior works, but some doesn't, the user must account for the fact that these are dictionaries:

        >>>'rxn09073_c0' in morph.rxn_labels['gene-no-match']

            works and returns True or False, indicating whether or not rxn09073_c0 is a gene-no-match reaction

        >>>morph.rxn_labels['gene-match'].add('rxn00456_c0')

        fails and throws an exception. add() is a set method and can't be used on dictionaries. (You could set an
        entry with real or arbitrary value to get around this if you really wished)

    Each inner dictionary is keyed with reaction_ids, hashed to tuples as such: (model_index, probability)
    Where reaction_id is a kbase reaction id with the compartment info appended to the end (e.g. rxn01316_c0),
    model_index is the index of the reaction in the objects['x_model'][modelreactions] list, and the probability
    is the reaction probability associated with each reaction from the probanno object. Reactions not in ProbAnno
    are given an arbitrary probability of -1.0

    Example evaluations:
    >>>rxn_labels['gene-no-match']['rxn01316_c0'][0]
    evaluates to the index of rxn01316_c0 in morph.objects['trans_model']['modelreactions']
    >>>rxn_labels['gene-no-match']['rxn01316_c0'][1]
    evaluates to the reaction probability of rxn01316_c0
    >>>'rxn01316_c0' in rxn_labels['gene-no-match']
    will evaluate True if the reaction is a gene-no-match reaction (an inner dict key)

    Note
    ----
    Function Requirements:
        - morph.probanno is a ReactionProbabilities
        - morph.src_model, morph.recon_model, morph.trans_model are FBAModels
        - morph.ws_id is the ID of a readable/writeable service workspace


    Post-Condition
    -------
    Morph
        a morph in which morph.rxn_labels holds a dictionary with the keys 'gene-match', gene-no-match', 'recon'
        and 'no-gene'. The value of each key holds a dictionary with 0 or more entries of the form:
            reaction_id -> (model_index, probability)
        morph.probhash contains a dictionary of reaction_ids (COMPARTMENT TRUNCATED) mapped to their probabilities
            e.g. rxn09876 -> 0.04545339

    Examples
    --------
    Given a morph of the form (only relevant attributes shown):

    >>>morph = Client.label_reactions(morph)

    would produce something like:

        probannows: 9145
        ws_name: MMws235
        src_modelws: 9145
        src_model: 19
        trans_model: 3
        probhash: ['rxn05653', 'rxn12345', rxn59595', 'rxn45644' ... (more)]
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        recon_model: 4
        probanno: 15
        morph.objects = ['source_model', 'recon_model', 'trans_model', 'probanno']
        morph.info = ['source_model', 'recon_model', 'trans_model', 'probanno']

    These could be examined like so:

    >>>morph.rxn_labels['no-gene'].keys()[1]
    u'rxn10316_c0'

    >>>morph.rxn_labels['no-gene']['rxn10316_c0'][1]
    0.444456666959

    >>>'rxn10316_c0' in morph.rxn_labels['no-gene']
    True
    """
        # get reaction sets
        recon_dict = dict([(r.rxn_id(), r) for r in self.recon_model.get_reactions()])
        trans_dict = dict([(r.rxn_id(), r) for r in self.trans_model.get_reactions()])
        model_dict = dict([(r.rxn_id(), r) for r in self.src_model.get_reactions()])
        # create the rxn_labels dictionary
        self.rxn_labels = {'gene-no-match': dict(),
                           'gene-match': dict(),
                           'no-gene': dict(),
                           'recon': dict(),
                           'common': dict()}
        # Some reference sets
        all_reactions = set(model_dict.keys()).union(recon_dict.keys())  # TODO: runtime will make you cry
        for rxn in all_reactions:
            if rxn in trans_dict and rxn in recon_dict:
                self.rxn_labels['common'][rxn] = (trans_dict[rxn], recon_dict[rxn])
            if rxn in model_dict and rxn not in trans_dict:
                if rxn not in recon_dict:
                    self.rxn_labels['gene-no-match'][rxn] = model_dict[rxn]
                else:
                    self.rxn_labels['recon'][rxn] = recon_dict[rxn]
            if rxn in trans_dict:
                gpr = trans_dict[rxn].gpr
                if gpr.gpr_type == 'no-gene' and rxn not in recon_dict:
                    self.rxn_labels['no-gene'][rxn] = trans_dict[rxn]
                else:
                    if rxn in recon_dict:
                        self.rxn_labels['gene-match'][rxn] = recon_dict[rxn]
                    else:
                        self.rxn_labels['gene-match'][rxn] = trans_dict[rxn]
            if rxn in recon_dict and rxn not in trans_dict and rxn not in model_dict:
                self.rxn_labels['recon'][rxn] = recon_dict[rxn]

    def build_supermodel(self):
        """
    Sets morph.model to a superset of all reactions in morph.rxn_labels

    Note
    ----
    Function Requirements:
        - morph.rxn_labels is a dictionary with the four keys ['gene-match', 'gene-no-match', 'recon', 'no-gene'],
        and it's values are dictionaries with entries of the form: reaction_id -> (model_index, probability)
        - morph.objects contains entries for the 'source_model' and 'recon_model' with data for the models in KBase (this
        is the output of the label_reactions(morph) function)

    Parameters
    ----------
    morph: Morph
        The morph for which you want to build a super_model (initializing morph.model)

    Returns
    -------
    Morph
        a morph object where morph.model, morph.ws_id forms a valid ObjectIdentity for a
        readable/writable KBase model object (the super-model)

    Examples
    --------
    Given a morph like so (only relevant attributes shown):

        ws_name: MMws235
        trans_model: 3
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.objects = ['source_model', 'recon_model', 'trans_model', 'probanno']
        morph.model = None

    >>>morph = Client.build_supermodel(morph)

    would produce something like this:

        ws_name: MMws235
        trans_model: 3
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.objects = ['source_model', 'recon_model', 'trans_model', 'probanno']
        morph.model = 5

    Where morph.model, morph.ws_id forms a valid ObjectIdentity for a model object in KBase (the super model)
    """
        src_rxns = dict([(r.rxn_id(), r) for r in self.src_model.get_reactions()])
        super_rxns = dict()
        specials = list()
        # copy the trans model
        self.model = self.trans_model.copy()
        #  ---->
        # Adding reactions into the translation.
        # First, go through every reaction they have in common and adjust if
        # necessary:
        reactions_to_remove = []
        for rxn_id in self.rxn_labels['common']:
            trans_rxn = self.rxn_labels['common'][rxn_id][0]  # MR
            recon_rxn = self.rxn_labels['common'][rxn_id][1]  # MR
            merge_gpr = trans_rxn.gpr.merge(recon_rxn.gpr)
            direction = _general_direction(trans_rxn, recon_rxn)
            if trans_rxn.gpr != merge_gpr or trans_rxn.get_direction() != direction:
                super_rxns[rxn_id] = (recon_rxn.get_rxn_ref(), recon_rxn.get_comp_ref(), direction, str(merge_gpr))
                removal_id = trans_rxn.get_removal_id()
                reactions_to_remove.append(removal_id)
        # ---->
        # removes the rxns we need to remove in place vs. making a new copy
        service.remove_reactions_in_place(self.model, reactions_to_remove)
        # Next, add all the reactions that aren't already in the translation:
        # Add the GENE_NO_MATCH reactions:
        for rxn_id in self.rxn_labels['gene-no-match']:
            reaction = self.rxn_labels['gene-no-match'][rxn_id]
            if reaction.is_special_ref():
                specials.append(reaction)
            else:
                super_rxns[rxn_id] = (reaction.get_rxn_ref(), reaction.get_comp_ref(), reaction.get_direction())
        # Add the RECON reactions:
        i = 0
        for rxn_id in self.rxn_labels['recon']:
            reaction = self.rxn_labels['recon'][rxn_id]
            if rxn_id in src_rxns:
                direction = _general_direction(reaction, src_rxns[rxn_id])
                i += 1
            else:
                direction = reaction.get_direction()
            if reaction.is_special_ref():
                raise Exception
                specials.append(reaction)
            else:
                if rxn_id not in super_rxns:
                    super_rxns[rxn_id] = (reaction.get_rxn_ref(), reaction.get_comp_ref(), direction, str(reaction.gpr))
        raise Exception(str(i))
        # ---->
        super_rxns = super_rxns.values()
        result = service.add_reactions(self.model, super_rxns, name='super_model')
        self.model = FBAModel(result[0], result[1])
        result = service.add_reactions_manually(self.model, specials, name='super_modelspc')
        self.model = FBAModel(result[0], result[1])


def _general_direction(model_rxn1, model_rxn2):
    """
    picks the more general of the two directions from reactions passed in
    """
    r1d = model_rxn1.get_direction()
    r2d = model_rxn2.get_direction()
    if r1d == r2d:
        return r1d
    elif r1d == '=' or r2d == '=':
        return '='
    else:
        raise Exception('directions are incompatible')


class AbstractGrowthCondition:
    """
    an interface for processing reactions according to some condition

    Subclases must implement evaluate(args) and return true or false. The primary use
    for this class is in the process_reactions method of the Client module, which removes
    reactions iteratiely, and decides to keep or remove a reaction based on the outcome of
    a GrowthCondition.  Here is an example:
        class SimpleCondition(AbstractGrowthCondition):
            def evaluate(args):
                # args must have attribute model
                fba = runfba(args['model'])
                return fba['Objective'] > 0
    This SimpleCondition keeps all reactions that are absolutely necessary for the models growth
    """

    def __init__(self):
        pass

    def evaluate(self, arguments):
        raise NotImplementedError()


class RepresentationError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return str(type(self.value)) + repr(self.value.gpr) + '\n' + repr(self.value.ftrs)


ws_client, fba_client = _init_clients()

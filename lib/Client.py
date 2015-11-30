# import necesary services
from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import random
import argparse
import time
import traceback
import Helpers
import copy
from operator import itemgetter
from Morph import Morph

# This module is the Client for the model morphing service. It performs
# operations necessary to morphing models. It is a module, performing
# functions on models to morph them from a source model to a target genome
# according to the (TODO: Name Algorithm) Algorithm designed by members of
# the Price lab at the Institute for Systems Biology
#
# Abstraction Function: a Client module for using the Morphing service,
# providing a user with functions for morphing models and analyzing the
# results.

# A global value for running requirement and representation checks with each method
debug = True

def _init_clients():
    # Get workspace service URL parameter
    with open ("./urls/.kbase_workspaceURL", "r") as myfile:
        url = myfile.read().replace('\n','')
    ws_client = Workspace(url)
    # Get FBA Model Services URL parameter
    with open ("./urls/.kbase_fbaModelServicesURL", "r") as myfile:
        url = myfile.read().replace('\n','')
    fba_client = fbaModelServices(url)

    return ws_client, fba_client

def prepare_supermodel(morph, fill_src=True):
    """
    Composition of the first several steps in the algorithm

    1) Fill the source model to media using probabilistic gapfilling (can be skipped if fill_src keyword arg is set to False)
    2) Translate the src_model to a pouplate morph.trans_model
    3) Draft reconstruction of target genome fills morph.recon_model field
    4) label reactions in the morph (populates morph.rxn_labels)
    5) Builds a super_model and puts it in the morph.model field. The model is now ready for the process_reactions function

    Note
    ----
    Function Requirements:
        morph.probanno, morph.probannows form a valid ObjectIdentity for a readable RxnProbs object in KBase
        morph.src_model, morph.src_modelws form a valid ObjectIdentity for a readable model object in KBase
        morph.media, morph.mediaws form a valid ObjectIdentity for a readable media object in KBase
        morph.genome, morph.genomews form a valid ObjectIdentity for a readable genome object in KBase
        morph.protcomp, morph.protcompws form a valid ObjectIdentity for a readable protein comparison object in KBase
        morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        An initialized Morph with the following Requirements:
            - morph.ws_id is a valid KBase workspace_id that user has permission
            to write to
            - morph has valid object identities for src_model, genome, probanno,
            protcomp (These 4 are object_ids, their ___ws counterparts are workspace_ids
            of user readable workspaces)
    fill_src: boolean,optional
        a boolean indicating that the src_model should first be filled using probabilistic gapfilling
        Optional, default is true.

    Returns
    -------
    Morph
        A Morph with the following state changes:
            - morph.model has the object_id of the super_model in KBase
            - morph.rxn_labels has the labels for the reactions in the super
            model
            - morph.probhash has the probability hash for reactions related to
            Target genome
            - morph.trans_model has the object_id of the translated model
            - morph.recon model has the object_id of the reconstruction of
            Target Genome

    Examples
    --------
    Suppose you have a morph initialized like so: (this is the output of Helpers.make_morph() )

    >>> morph = make_morph()

        probannows: 9145
        ws_name: MMws235
        src_modelws: 9145
        protcompws: 9145
        src_model: 19
        media: 18
        removed_ids: None
        mediaws: 9145
        trans_model: None
        probhash: None
        genomews: 9145
        essential_ids: None
        genome: 3
        rxn_labels: None
        model: None
        ws_id: 11444
        recon_model: None
        probanno: 15
        protcomp: 6

    A call to Client.prepare_supermodel(morph) will return a morph of this sort of form:

    >>> morph = Client.prepare_supermodel(morph)

        probannows: 9145
        ws_name: MMws235
        src_modelws: 9145
        protcompws: 9145
        src_model: 19
        media: 18
        removed_ids: None
        mediaws: 9145
        trans_model: 3
        probhash: [u'rxn00001' ... u'rxn97854']
        genomews: 9145
        essential_ids: None
        genome: 3
        rxn_labels: ['gene-no-much', 'gene-match', 'recon', 'no-gene']
        model: 5
        ws_id: 11444
        recon_model: 4
        probanno: 15
        protcomp: 6

    See Also
    --------
    fill_src_to_media
    translate_features
    reconstruct_genome
    label_reactions
    build_supermodel

    This functions post condition preps the morph for the Client.process_reactions(morph) function
    """
    morph = copy.deepcopy(morph)
    if(fill_src):
        morph = fill_src_to_media(morph)
    morph = translate_features(morph)
    morph = reconstruct_genome(morph)
    morph = label_reactions(morph)
    morph = build_supermodel(morph)
    return morph

def fill_src_to_media(morph):
    """
    Takes a morph and probanno gap-fills the source model to given media

    Fills a source model to a given media. This is an important step if the source model does not grow initially on the media. NOTE: THIS MODIFIES SRC_MODEL IN THE MORPH
    However, it does not overwrite the actual source model in KBase. It is imported to morph.ws_id, and morph.src_modelws is changed to match ws_id. Therefore, the data of the source
    model is preserved as is, but the reference to its location in KBase is lost from the morph.

    Note
    ----
    Function Requirements:
    morph.src_model, morph.src_modelws form a valid ObjectIdentity for a readable model object in KBase
    morph.media, morph.mediaws form a valid ObjectIdentity for a readable media object in KBase
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        A morph for which you want to gap-fill the source model to it's given media

    Returns
    -------
    Morph
        a morph object where morph.src_model, morph.src_modelws forms a valid ObjectIdentity for a
        readable/writable KBase model object (the gapfilled version of parameter's morph.src_model)

    Raises
    ------
    exception: MediaException
        (TO BE IMPLEMENTED) An exception indicating that the given media could not be gapfilled in such a way that model objective function exceeds 0

    Examples
    --------
    Given a morph of the form (Only relevant attributes shown):

        ws_name: MMws235
        src_modelws: 9145
        src_model: 19
        media: 18
        mediaws: 9145
        ws_id: 11444

    >>> morph = fill_src_to_media(morph)

    Produces something like so:

        ws_name: MMws235
        src_modelws: 11444
        src_model: 2
        media: 18
        mediaws: 9145
        ws_id: 11444

    Note: The reference to the old source model is lost from the model, but it still exists at it's previous KBase location un-modified

    """
    # Copy the source model to the new workspace
    morph = copy.deepcopy(morph)
    srccopy = ws_client.copy_object({'from': {'wsid':morph.src_modelws, 'objid':morph.src_model}, 'to': {'wsid':morph.ws_id, 'name': u'srccopy'}})
    morph.src_model = srccopy[0]
    morph.src_modelws = srccopy[6]

    # Gapfill the model and reset the src_model and src_modelws attributes of
    # the morph
    fba_formulation = {'media': morph.media, 'media_workspace': morph.mediaws}
    gap_formulation = {'probabilisticAnnotation' : morph.probanno, 'probabilisticAnnotation_workspace' : morph.probannows, u'formulation': fba_formulation}
    params = {u'model': morph.src_model, u'model_workspace': morph.src_modelws, u'out_model' : u'source_filled', u'workspace' : morph.ws_id, u'formulation' : gap_formulation, u'integrate_solution' : True, u'gapFill' : u'gf'}
    model_info = fba_client.gapfill_model(params)
    morph.src_model = model_info[0]
    morph.src_modelws = model_info[6]
    return morph
def translate_features(morph):
    """
    Translates the features in the source model to matches in the target genome

    NOTE: The translated model will lose reactions in the source model that did not have matching features
    and is not guaranteed to be useful in FBA simulation.

    Note
    ----
    Function Requirements:
    morph.src_model, morph.src_modelws form a valid ObjectIdentity for a readable model object in KBase
    morph.protcomp, morph.protcompws form a valid ObjectIdentity for a readable protein comparison object in KBase
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        A morph for which you want to produce a translated model to be referenced by morph.trans_model
        (preceeds reaction labelling)

    Returns
    -------
    Morph
        a morph object where morph.trans_model, morph.ws_id forms a valid ObjectIdentity for a
        readable/writable KBase model object (the translation of morph.src_model)
    """
 # import necesary services# import necesary service
    if (debug):
        assert morph.protcomp is not None, "protcomp is None"
        assert morph.protcompws is not None, "protcompws is None"
        assert morph.src_model is not None, "src_model is None"
        assert morph.src_modelws is not None, "src_modelws is None"
    morph = copy.deepcopy(morph)
    trans_params = {'keep_nogene_rxn': 1, 'protcomp': morph.protcomp, 'protcomp_workspace': morph.protcompws, 'model': morph.src_model, 'model_workspace': morph.src_modelws, 'workspace': morph.ws_id}
    morph.trans_model = fba_client.translate_fbamodel(trans_params)[0]
    return morph

def reconstruct_genome(morph):
    """
    Builds a reconstruction from morph.genome and saves it in morph.recon_model

    Note
    ----
    Function Requirements:
    morph.genome, morph.genomews form a valid ObjectIdentity for a genome object in KBase
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        A morph for which you want to build a reconstruction referenced by morph.recon_model
        (preceeds reaction labelling)

    Returns
    -------
    Morph
        a morph object where morph.recon_model, morph.ws_id forms a valid ObjectIdentity for a
        readable/writable KBase model object (the draft reconstruction of morph.genome)

    Examples
    --------
    Given a morph of the form (only relevant attributes shown):

        ws_name: MMws235
        genomews: 9145
        genome: 3
        ws_id: 11444
        recon_model: None

    >>> morph = Client.reconstruct_genome(morph)

    would produce something like this:

        ws_name: MMws235
        genomews: 9145
        genome: 3
        ws_id: 11444
        recon_model: 4

    Where morph.recon_model, morph.ws_id forms a valid ObjectIdentity for a readable/writable KBase model object (the draft reconstruction of morph.genome)

    """
    morph = copy.deepcopy(morph)
    if (debug):
        assert morph.genome is not None, "genome is None"
        assert morph.genomews is not None, "genomews is None"
    recon_params = {'genome': morph.genome, 'genome_workspace': morph.genomews, 'workspace': morph.ws_id}
    morph.recon_model = fba_client.genome_to_fbamodel(recon_params)[0]
    return morph

def _get_objects(morph):
    """
    This is used by the label reactions function and is not intended for end user use, proceed with caution

    Gets the model and probanno objects associated with the morph.

    Populates morph.objects with a dictionary of objects. Three are FBA Models, with keys 'source_model', 'recon_model', 'trans_model'.
    The last is a ProbAnno RxnProbs object, keyed 'probanno'. Corresponding entries are put in a dictionary morph.info, with the object_info
    as it is returned from the workspace service"

    Requires:
        src_model, trans_model, recon_model, and probanno are all attributes of morph and valid IDs for KBase objects. src_modelws and probannows are valid
        workspace ID's corresponding with src_model and probanno attributes of morph to form valid ObjectIdentity's.

   :param morph: the Morph object with the associated objects.
    """
    morph = copy.deepcopy(morph)
    # make the objects and info dictionaries
    if (not isinstance(morph.objects, dict)):
        morph.objects = dict()
    if (not isinstance(morph.info, dict)):
        morph.info = dict()
    # Get the src_model
    obj  = ws_client.get_objects([{'objid': morph.src_model, 'wsid': morph.src_modelws}])[0]
    morph.objects['source_model'] = obj['data']
    morph.info['source_model'] = obj['info']
    # Get the translate model
    obj  = ws_client.get_objects([{'objid': morph.trans_model, 'wsid': morph.ws_id}])[0]
    morph.objects['trans_model'] = obj['data']
    morph.info['trans_model'] = obj['info']
    # Get the recon model
    obj  = ws_client.get_objects([{'objid': morph.recon_model, 'wsid': morph.ws_id}])[0]
    morph.objects['recon_model'] = obj['data']
    morph.info['recon_model'] = obj['info']
    # Get the probanno
    obj  = ws_client.get_objects([{'objid': morph.probanno, 'wsid': morph.probannows}])[0]
    morph.objects['probanno'] = obj['data']
    morph.info['probanno'] = obj['info']
    morph.probhash= dict()
    for rxn in morph.objects['probanno']['reaction_probabilities']:
        morph.probhash[rxn[0]] = rxn[1]
    return morph
def label_reactions(morph):
    """
    Labels morph's reactions from translated model, reconstruction, and source

    Populates the rxn_labels attribute in the Morph object with a Dictionary of four dictionaries of reactrion_id -> value tuples.
    The first level dicts are named with the keys:
        - gene-match
        - gene-no-match
        - no-gene
        - recon
    Populates morph.probhash with a dictionary of compartment-truncated reaction_ids to their probability relative to genome in question
        (derived from morph.probanno).
    Also populates morph.objects, morph.info with data from KBase objects (advanced use)

    Technical Details (Consult if you are encountering issues):
    For simplicity, an end user can treat the interiors of each dictionary like a set of reaction ids, but it is important to note the
    this is actually a nested dictionary with entries of the form: reaction_id -> (model_index, probability)
    Thus, some set like behavior works, but some doesn't, the user must account for the fact that these are dictionaries:

        >>>'rxn09073_c0' in morph.rxn_labels['gene-no'match']

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
    rxn_labels['gene-no-match']['rxn01316_c0'][0] evaluates to the index of rxn01316_c0 in morph.objects['trans_model']['modelreactions']
    rxn_labels['gene-no-match']['rxn01316_c0'][1] evaluates to the reaction probability of rxn01316_c0
    'rxn01316_c0' in rxn_labels['gene-no-match'] will evaluate True if the reaction is a gene-no-match reaction (an inner dict key)

    Note
    ----
    Function Requirements:
    morph.probanno, morph.probannows form a valid ObjectIdentity for a readable RxnProbs object in KBase
    morph.src_model, morph.src_modelws form a valid ObjectIdentity for a readable model object in KBase
    morph.recon_model, morph.ws_id form a valid ObjectIdentity for a readable model object in KBase
    morph.trans_model, morph.ws_id form a valid ObjectIdentity for a readable model object in KBase
    morph.ws_id is the ID of a readable/writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        the morph fo which you want to generate reaction labels

    Returns
    -------
    Morph
        a morph in which morph.rxn_labels holds a dictionary with the keys 'gene-match', gene-no-match', 'recon'
        and 'no-gene'. The value of each key holds a dictionary with 0 or more entries of the form:
            reaction_id -> (model_index, probability)
        morph.probhash contains a dictionary of reaction_ids (COMPARTMENT TRUNCATED) mapped to their probabilities
            e.g. rxn09876 -> 0.04545339

    Examples
    --------
    Given a morph of the form (only relevant attricutes shown):

        probannows: 9145
        ws_name: MMws235
        src_modelws: 9145
        src_model: 19
        trans_model: 3
        probhash: None
        rxn_labels: None
        ws_id: 11444
        recon_model: 4
        probanno: 15
        morph.objects = None
        morph.info = None

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
    morph = _get_objects(morph)
    morph = copy.deepcopy(morph)
    label_time = time.time()
    # create local vars for parts of the morph
    trans_model = morph.objects['trans_model']
    model = morph.objects['source_model']
    recon = morph.objects['recon_model']
    probanno = morph.objects['probanno']
    # create the rxn_labels dictionary
    rxn_labels = {'gene-no-match': dict(), 'gene-match': dict(), 'no-gene': dict(), 'recon': dict()}
    # Build a hash of rxn_ids to probability to save look up time
    prob_hash = dict()
    for rxn in probanno['reaction_probabilities']:
        prob_hash[rxn[0]] = rxn[1]
    # label gene-match and no-gene from translation
    for i in range(len(trans_model['modelreactions'])):
        #anni e.g. rxn_id = rxn01316_c0
        mdlrxn = trans_model['modelreactions'][i]
        rxn_id = mdlrxn['reaction_ref'].split('/')[-1] + '_' + str(mdlrxn['modelcompartment_ref'].split('/')[-1]) # -1 index gets the last in the list
        rxn_prot = mdlrxn['modelReactionProteins']
        if (len(rxn_prot) == 1 and len(rxn_prot[0]['modelReactionProteinSubunits']) == 0 and rxn_prot[0]['note'] != 'spontaneous' and rxn_prot[0]['note'] != 'universal'):
            try:
                rxn_labels['no-gene'][rxn_id] = (i, prob_hash[mdlrxn['reaction_ref'].split('/')[-1]])
            except KeyError:
                rxn_labels['no-gene'][rxn_id] = (i, -1.0)
        else:
            try:
                rxn_labels['gene-match'][rxn_id] = (i, prob_hash[mdlrxn['reaction_ref'].split('/')[-1]])
            except KeyError:
                rxn_labels['gene-match'][rxn_id] = (i, -1.0)
    # Label gene-no-match as those unique to model compared with translation
    for i in range(len(model['modelreactions'])):
        mdlrxn = model['modelreactions'][i]
        rxn_id = mdlrxn['reaction_ref'].split('/')[-1] + '_' + str(mdlrxn['modelcompartment_ref'].split('/')[-1]) # -1 index gets the last in the list
        if (rxn_id not in rxn_labels['gene-match']) and (rxn_id not in rxn_labels['no-gene']):
            try:
                rxn_labels['gene-no-match'][rxn_id] = (i, prob_hash[mdlrxn['reaction_ref'].split('/')[-1]])
            except KeyError:
                rxn_labels['gene-no-match'][rxn_id] = (i, -1.0)
    # label recon as those unique to recon compared with model (model = gene-match + gene-no-match + no-gene)
    for i in range(len(recon['modelreactions'])):
        mdlrxn = recon['modelreactions'][i]
        rxn_id = mdlrxn['reaction_ref'].split('/')[-1] + '_' + str(mdlrxn['modelcompartment_ref'].split('/')[-1]) # -1 index gets the last in the list
        # if the recon reaction is already in the model
        if (rxn_id not in rxn_labels['gene-match']) and (rxn_id not in rxn_labels['no-gene']) and (rxn_id not in rxn_labels['gene-no-match']):
            try:
                rxn_labels['recon'][rxn_id] = (i, prob_hash[rxn_id])
            except KeyError:
                rxn_labels['recon'][rxn_id] = (i, -1.0)
    morph.rxn_labels = rxn_labels
    return morph

def build_supermodel(morph):
    """
    Sets morph.model to a superset of all reactions in morph.rxn_labels

    Note
    ----
    Function Requirements:
    morph.rxn_labels is a dictionary with the four keys ['gene-match', 'gene-no-match', 'recon', 'no-gene'],
        and it's values are dictionaries with entries of the form: reaction_id -> (model_index, probability)
    morph.objects contains entries for the 'source_model' and 'recon_model' with data for the models in KBase (this
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
    morph = copy.deepcopy(morph)
    super_time = time.time()
    model = morph.objects['source_model']
    recon = morph.objects['recon_model']
    super_rxns = list()
    # Add the GENE_NO_MATCH reactions:
    for rxn_id in morph.rxn_labels['gene-no-match']:
        # morph.rxn_labels['gene-no-match'][0] gives the index of the reaction in the model['modelreactions'] list to make this look up O(1) instead of O(n)
        reaction = model['modelreactions'][morph.rxn_labels['gene-no-match'][rxn_id][0]]
        #TODO: See what more you can fix/add (gpr?)
        super_rxns.append((reaction['reaction_ref'].split('/')[-1], str(reaction['modelcompartment_ref'].split('/')[-1][0]), reaction['direction'], 'GENE_NO_MATCH', '', reaction['name']))
    # Add the recon reactions:
    for rxn_id in morph.rxn_labels['recon']:
        # morph.rxn_labels['recon'][0] gives the index of the reaction in the model['modelreactions'] list to make this look up O(1) instead of O(n)
        reaction = recon['modelreactions'][morph.rxn_labels['recon'][rxn_id][0]]
        #TODO: See what more you can fix/add (gpr?)
        super_rxns.append((reaction['reaction_ref'].split('/')[-1], str(reaction['modelcompartment_ref'].split('/')[-1][0]), reaction['direction'], 'recon', '', reaction['name']))
    morph.model = fba_client.add_reactions({'model': morph.trans_model, 'model_workspace': morph.ws_id, 'output_id': 'super_model', 'workspace': morph.ws_id, 'reactions': super_rxns})[0]
    return morph

def removal_list(rxn_dict, list_range=None):
    """
    Generates a removal list from the rxn_dict (e.g. morph.rxn_labels['no-gene'])

    Accepts a dictionary of reactions -> (index, probability) as is given by morph.rxn_labels[label] and returns a sorted
    list of tuples that can be passed as the rxn_list parameter to the process_reactions function. This allows for processing
    of only some reactions from the model at a time. Tuples are in the form (reaction_id, index, probability). Setting the range keyword
    arg returns only a subset of the sorted list given by the index range passed in as a tuple (start, end).

    Note
    ----
    Function Requirements:
    rxn_dict is of a form analagous to rxn_labels[key]

    Parameters
    ----------
    rxn_dict: dict
        A dictionary of the form specified in rxn_labels[key] (see above description)
    list_range: tuple, optional. (start, end)
        A range within the list once it is sorted by probability (start inclusive, end exclusive)

    Returns
    -------
    list
        an sorted list of reactions to remove (sorted from low import probability to high

    Examples
    --------
    Given a rxn_dict with 4 entries (shown below) and list_range of (0, 3)
        >>>rxn_dict = { 'rxn00001_c0': (3, 0.2453)
                     'rxn03035_c0': (78, 0.8744)
                     'rxn00451_c0': (32, 0.0004)
                     'rxn00602_c0': (9, 0.0056) }

        >>>rxn_list = Client.removal_list(rxn_dict, list_range=(0, 3))

    rxn_list would evaluate to: [('rxn00451_c0', (32, 0.0004)), ('rxn00602_c0', (9, 0.0056)), ('rxn00001_c0, (3, 0.2453))]
    """
    # Sort by probanno. items() returns (K, V=(model_index, prob))
    def get_key(item):
        return item[1][1]
    removal_list = sorted(rxn_dict.items(), key=get_key)
    if list_range is not None:
       return removal_list[list_range[0]:list_range[1]]
    return removal_list

def process_reactions(morph, rxn_list=None, name='', process_count=0, get_count=False):
    """
    Attempts removal of rxn_list reactions from morph (i.e. morph.rxn_labels[label])

    Attempts removal of each reaction , keeping it only if it is essential to the objective function.
    Populates morph.essential_ids with the reaction_ids that could not be removed without breaking the model.
    morph.removed_ids is populated with the ids of the reactions that were removed. Both are populated as a dictionary of
    reaction_ids to their model_index and probability (like the entries in rxn_labels).
        e.g. reaction_id -> (model_index, probability)

    if rxn_list is None, the function removes (in low to high probability order) gene-no-match reactions followed by no-gene reactions

    Controlling name and process_count parameters allows the user tocontrol the number of models created in the morph.ws_id

    Note
    ----
    Function Requirements:
    if rxn_list is None, rxn_labels['gene-no-match'] and rxn_labels['no-gene'] must be dictionaries
        with 0 or more entries of the form reaction_id -> (model_index, probability)
    morph.model, morph.ws_id form a valid ObjectIdentity for a readable model object in KBase
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        The morph containing the model (morph.model) from which reactions will be removed
    rxn_list: list, optional
        A sorted list of tuples of the form (reaction_id, (model_index, probability)) which will be processed from morph,model
        (form is the output of removal_list() function)
        (if left out, all gene-no-match follwed by all no-gene reactions will be processed in low to high likelihood order)
    name: String, optional
        Setting name will begin all output model names with name. Default is MM
        output models are named as follows: str(name) + '-' + str(process_count)
    process_count: int, optional
        A number indicating how many reactions have been processed so far, used in nameing output models. Default is 0.
        output models are named as follows: str(name) + '-' + str(process_count)
    get_count: Boolean, optional
        A Boolean flag indicating whether the process_count should be returned with the morph (as a tuple). Used when not processing all
        reactions at once. Deafault is False

    Returns
    -------
    Morph
        A morph with a new model, and essential_ids/removed)ids used to keep track of changes
    int (Only if get_count=True)
        process_count (number of models created, used for name managing)

    Examples
    --------
    Given a morph of the form (only relevant attributes shown):
        ws_name: MMws235
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.model = 5
        morph.essential_ids = None
        morph.removed_ids = None

    >>>morph = Client.process_reactions(morph, rxn_list=Client,removal_list(morph.rxn_labels['gene-no-match'], list_range=(0, 10)))

    would process the 10 least likely gene-no-match reactions, which might produce a morph like this:

        ws_name: MMws235
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.model = 5
        morph.essential_ids = ['rxn10316_c0', 'rxn23434_c0', 'rxn78687_c0']
        morph.removed_ids = ['rxn11123_c0', 'rxn34534_c0', 'rxn90565_c0', 'rxn78987_c0', 'rxn12321_c0', 'rxn89034_c0', 'rxn88888_c0']

    where removal of one of the reactions given by a key in morph.essential_ids would result in a model that has an objective value of 0.000 in FBA simulation
    """
    morph = copy.deepcopy(morph)
    ws = morph.ws_id
    # Sort by probanno. items() returns (K, V=(model_index, prob))
    def get_key(item):
        return item[1][1]
    # label argument behavior
    if (rxn_list is None):
        rxn_dict = morph.rxn_labels['gene-no-match']
        removal_list = sorted(rxn_dict.items(), key=get_key)
        if (debug):
            for i in range(len(removal_list) - 1):
                assert removal_list[i][1][1] <= removal_list[i + 1][1][1]
        rxn_dict = morph.rxn_labels['no-gene']
        removal_list.append(sorted(rxn_dict.items(), key=get_key))
    else:
        removal_list = rxn_list
    # instantiate lists only if needed
    if morph.essential_ids is None:
        morph.essential_ids = dict()
    if morph.removed_ids is None:
        morph.removed_ids = dict()
    #Give objs a general name if none is provided
    if name == '':
        name = 'MM'
    for i in range(len(removal_list)):
        reaction_id = removal_list[i][0]
        print 'Reaction to remove: ' + str(removal_list[i])
        #TODO Find someway to fix the behaivior bug if model_id is not in ws, etc.
        new_model_id = fba_client.remove_reactions({'model': morph.model, 'model_workspace':ws, 'output_id': name + '-' + str(process_count), 'workspace':morph.ws_id, 'reactions':[reaction_id]})[0]
        if iterative_models:
            process_count += 1
        fba_formulation = {'media': morph.media, 'media_workspace': morph.mediaws}
        fba_params = {'fba': str(name + '-FBA-' + str(process_count)), 'workspace': ws, 'model' : new_model_id, 'model_workspace':ws, 'formulation': fba_formulation}
        fbaMeta = fba_client.runfba(fba_params)
        flux = fbaMeta[-1]['Objective']
        print 'FBA Flux: ' + str(flux)
        if (flux > 0.0):
            # removed successfully
            print 'Removed ' + str(reaction_id)
            morph.model = new_model_id
            morph.removed_ids[reaction_id] = removal_list[i][1]
        else:
            # essential
            print str(reaction_id) + ' is Essential'
            morph.essential_ids[reaction_id] = removal_list[i][1]
    if get_count:
        return  morph, process_count
    return morph

def find_alternative(morph, reaction):
    """
    NOT FUNCTIONING AS DESIRED

    attemps to find an alternative to a reaction given by reaction item of the form (reaction_id, (model_index, probanno value))

    An example reaction_item would be morph.essential_ids.items()[0]

    Note
    ----
    Function Requirements:
    reaction corresponds to a Reaction in morph.model
    morph.probanno, morph.probannows form a valid ObjectIdentity for a readable RxnProbs object in KBase
    morph.model, morph.ws_id forms a valid ObjectIdentity for a model object in KBase
    morph.probhash is initialized to a dictionary of compartment truncated reaction_ids -> probabilities
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        contains the morph.model from which will try to find an alternative to a particular reaction
    reaction: String (reaction_id)
        The reaction_id for a reaction to find a more likely alternative for

    Returns
    -------
    Morph
        a morph where morph,model, morph.ws_id is a valid ObjectIdentity for a model object in KBase
        (if a more likely alternative was found for reaction, this will be reflected in the model object)
    list
        new_reactions: a list of the new reactions considered as a potential alternatice to 'reaction'

    Examples
    --------
    Given a morph of the form:
        ws_name: MMws235
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.model = 214
        probannows: 9145
        probhash: ['rxn00001', 'rxn03434', 'rxn67834' ... (more)]
        probanno: 15
        morph.essential_ids = ['rxn10316_c0', 'rxn23434_c0', 'rxn78687_c0']

    Below, morph.essential_ids.keys()[0] evaluates to 'rxn10316_c0' as shown above
    >>>morph, new_rxns = Client.find_alternative(morph, morph.essential_ids.keys()[0])

    would return:

        morph =

        ws_name: MMws235
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.model = 216
        probannows: 9145
        probhash: ['rxn00001', 'rxn03434', 'rxn67834' ... (more)]
        probanno: 15
        morph.essential_ids = ['rxn10316_c0', 'rxn23434_c0', 'rxn78687_c0']

        new_reactions =

        ['rxn34343_c0', 'rxn 78923_c0', 'rxn00002_c0']

    """
    m = copy.deepcopy(morph)
    model_id = m.model
    ws_id = m.ws_id
    fba_formulation = {'media': m.media, 'media_workspace': m.mediaws}
    # blacklist reactions we've already removed and the one we are removing now
    new_model_id = fba_client.remove_reactions({'model': model_id, 'model_workspace': ws_id, 'workspace': ws_id, 'output_id': str('findalt' + str(reaction)), 'reactions': [reaction]})[0]
    m.removed_ids[reaction] = reaction_item[1]
    blacklisted_rxns = [i.split('_')[0] for i in m.removed_ids]
    gap_formulation = {'probabilisticAnnotation' : m.probanno, 'probabilisticAnnotation_workspace' : m.probannows, u'formulation': fba_formulation, u'blacklisted_rxns':blacklisted_rxns}
    for i in gap_formulation['blacklisted_rxns']:
        print i
    print "\n\n" + str(reaction)
    print m.removed_ids[reaction]
    assert reaction.split('_')[0] in gap_formulation['blacklisted_rxns']
    params = {'model': new_model_id, 'model_workspace': ws_id, 'workspace' : ws_id, 'formulation' : gap_formulation, 'integrate_solution' : True, 'gapFill' : 'gf', 'completeGapfill': False}
    model_info = fba_client.gapfill_model(params)
    filled_model = model_info[0]
    object = ws_client.get_objects([{'name': model_info[0] + '.gffba', 'wsid': params['workspace']}])[0]
    gapfill = object['data']['gapfillingSolutions'][0]
    rxn_dicts = gapfill['gapfillingSolutionReactions']
    gap_reactions = [a['reaction_ref'].split('/')[-1] + '_' + a['compartment_ref'].split('/')[-1] + str(0) for a in rxn_dicts]
    modelreactions = set()
    for key in m.rxn_labels:
        if reaction in m.rxn_labels[key].keys():
            del m.rxn_labels[key][reaction]
    for key in m.rxn_labels:
        modelreactions |= set(m.rxn_labels[key].keys())
    new_reactions = [r for r in gap_reactions if r not in modelreactions]
    rmv_probability = m.probhash[reaction.split('_')[0]]
    probsum  = 0.0
    for i in new_reactions:
        try:
            prob = m.probhash[i/split('_')[0]]
            if prob > 0:
                probsum += prob
        except KeyError:
            m.probhash[i.split('_')[0]] = -1.0
    new_probability = probsum / (0.0 + len(new_reactions))
    if new_probability > rmv_probability:
        m.model = filled_model
    a = rmv_probability
    b = new_probability
    return m, new_reactions

def probanno_optimization(morph, rxn_list=None):
    """
    Deprecated, not maintained
    """
    morph = copy.deepcopy(morph)
    if rxn_list is None:
        rxn_list = morph.essential_ids.keys()
    #Build Massive Probanno Hash. Consider doing this elsewhere
    #Its actually not terribly massive, but 2000+ items is still too many to iterate over more than once
    probhash= dict()
    for rxn in morph.objects['probanno']['reaction_probabilities']:
        probhash[rxn[0]] = rxn[1]
    a = 0.0
    b = 0.0
    for reaction in rxn_list:
        print reaction
        m, model_id, new_reactions = find_alternative(reaction, morph=morph)
        rmv_probability = probhash[reaction.split('_')[0]]
        probsum  = 0.0
        for i in new_reactions:
            print i
            print i.split('_')[0]
            try:
                probsum += probhash[i.split('_')[0]]
            except KeyError:
                # THIS NEEDS TO BE CHANGED IF YOU ARE SAVING THE HASH OUTSIDE OF
                # THIS FUNCTION
                probhash[i.split('_')[0]] = 0.0
        new_probability = probsum / (0.0 + len(new_reactions))
        if new_probability > rmv_probability:
            morph = m
        a = rmv_probability
        b = new_probability
    return morph, m, a, b, new_reactions, probhash


# Finishing/Cleanup  Steps
def _finish(morph, save_ws=False):
    if not save_ws:
        ws_client.delete_workspace({'id': morph.ws_id})
    else:
        print 'Saved'

ws_client, fba_client =  _init_clients()

def get_morph_rxns(morph, label=None):
    """
    returns a list of all reactions in the model (by reaction_id)

    Note
    ----
    Function Requirements:
    morph.model, morph.ws_id form a valid ObjectIdentity for a model object in KBase
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        morph which has the model to get reactions from
    label: String, optional
        pass a label to only get reactions of a certain type
        valid labels: essential, removed, gebe-no-match, gene-match, no-gene, recon

    Returns
    -------
    list
        a list of reaction_ids in the morph's model. If a label parameter is set, only that specific type of reaction
        will be added to the list, If none exist, the empty list is returned

    Examples
    --------
    Given a morph of the form:
        ws_name: MMws235
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.model = 216
        probannows: 9145
        probhash: ['rxn00001', 'rxn03434', 'rxn67834' ... (more)]
        probanno: 15
        morph.essential_ids = ['rxn10316_c0', 'rxn23434_c0', 'rxn78687_c0']

    >>>rxn_list = Client.get_morph_rxns(morph)

     would return a list of the form:

    ['rxn00001_c0', 'rxn78786_c0', 'rxn45324_e0' ... (more)]

    >>>rxn_list = Client.get_morph_rxns(morph, label='gene-no-match')

    would return a list of the same form, and would be identical to morph.rxn_labels['gene-no-match'].keys()
    """
    if label is None:
        modelobj = Helpers.get_object(morph.model, morph.ws_id)
        rxnlist = list()
        model = modelobj['data']
        for i in range(len(model['modelreactions'])):
            mdlrxn = model['modelreactions'][i]
            rxn_id = mdlrxn['reaction_ref'].split('/')[-1] + '_' + str(mdlrxn['modelcompartment_ref'].split('/')[-1]) # -1 index gets the last in the list
            rxnlist.append(rxn_id)
    if label is "essential":
        if morph.essential_ids is not None:
            return morph.essential_ids.keys()
        else:
            return []
    if label is "removed":
        if morph.removed_ids is not None:
            return morph.removed_ids.keys()
        else:
            return []
    if label in morph.rxn_labels.keys():
        return morph.rxn_labels[label].keys()
    else:
        return []
def get_reactions(reactions):
    """
    returns reaction info for given reactions in truncated compartment form, including definition and more

    Parmeters
    ---------
    reactions: list,(String) reaction_id
        a list of reaction_ids (compartment truncated) to get definition, equation, and other info on (will update this doc section to be more specific

    Returns
    -------
    list
        a list of dictionaries where the keys are strings describing info about a reaction and the values are that information specific to a reaction

    Examples
    --------
    >>>rxn_info_list = Client.get_reactions(['rxn00001', 'rxn54545',])

    >>>rxn_info_list[0].keys()
    [u'definition', u'direction', u'name', u'enzymes', u'equation', u'deltaGErr', u'reversibility', u'abbrev', u'deltaG', u'id', u'aliases']

    >>>rxn_info_list[0]['definition']
    u'(1) H2O[c] + (1) PPi[c] <=> (2) Phosphate[c] + (2) H+[c]'
    """
    return fba_client.get_reactions({'reactions':reactions})
def get_reaction(reaction):
    """
    returns reaction info for a given reaction in truncated compartment form, including definition and more

    Parmeters
    ---------
    reaction: String, reaction_id
        a reaction_id (compartment truncated) to get definition, equation, and other info on (will update this doc section to be more specific

    Returns
    -------
    dict
        a dictionary where the keys are strings describing info about a reaction and the values are that information specific to a reaction

    Examples
    --------
    >>>rxn_info = Client.get_reaction('rxn00001')

    >>>rxn_info.keys()
    [u'definition', u'direction', u'name', u'enzymes', u'equation', u'deltaGErr', u'reversibility', u'abbrev', u'deltaG', u'id', u'aliases']

    >>>rxn_info['definition']
    u'(1) H2O[c] + (1) PPi[c] <=> (2) Phosphate[c] + (2) H+[c]'
    """
    return fba_client.get_reactions({'reactions':[reaction]})
def get_mdlrxn_info(morph, reaction_ids=None):
    """
    Returns a list of dict() of informtation about reactions in the model and it's relation to the morph

    Returns a dictionary with the following information for each reaction specified. If no reaction_ids key
    word is set, examines all reactions in morph.model

    'label' - the morph classification of the reaction: (gene-match, gene-no-match, no-gene, recon)
    'equation' - the equation for the reaction returned in terms of compound_ids
    'definition' - the equation for the reaction returned in terms of coumpound names
    """
    allreactions = False
    not_in_morph = list()
    if reaction_ids is None:
        all_reactions = True
        reaction_ids = list()
        for key in morph.rxn_labels.keys():
            reaction_ids.extend(morph.rxn_labels[key].items())
    else:
        # Make tuples like rxN_labels[key].items()
        for i in range(len(reaction_ids)):
            rxn = reaction_ids[i]
            for key in morph.rxn_labels.keys():
                if (rxn in morph.rxn_labels[key]):
                    value = morph.rxn_labels[key][rxn]
                    reaction_ids[i] = (rxn, value[0], value[1])
    results = list()
    reactions = fba_client.get_reactions({'reactions': [rxntup[0].split('_')[0] for rxntup in reaction_ids]})
    # WARNING: This assumes lists are identically ordered
    for i in range(len(reaction_ids)):
        try:
            rxn = reactions[i]
        except:
            print i
        mdlrxn = reaction_ids[i]
        model = morph.objects['source_model']
        if (mdlrxn[0] in morph.rxn_labels['recon']):
            model = morph.objects['recon_model']
        ans = dict()
        for key in morph.rxn_labels.keys():
            if (mdlrxn[0] in morph.rxn_labels[key]):
               ans['label'] = str(key)
        ans['id'] = mdlrxn[0]
        ans['definiton'] = rxn['definition']
        ans['equation'] = rxn['equation']
        index = mdlrxn[1]
        results.append(ans)
    return results

def remove_reactions_by_dict(morph, rxn_dict):
    """
    Removes reactions provided in a dictionary of a style like morph.rxn_labels['gene-no-match'] or morph.essential_ids
    """
    morph = copy.deepcopy(morph)
    morph.model  = fba_client.remove_reactions({'model': morph.model, 'model_workspace': morph.ws_id, 'workspace': morph.ws_id, 'output_id': 'removerxnsbydict', 'reactions': rxn_dict.keys()})[0]
    return morph
def probanno_fill(morph, name=None):
    """
    runs probabilistic gapfilling on morph.model

    Note
    ----
    Function Requirements:
    morph.probanno, morph.probannows form a valid ObjectIdentity for a readable RxnProbs object in KBase
    morph.media, morph.mediaws form a valid ObjectIdentity for a readable media object in KBase
    morph.model, morph.ws_id form a valid ObjectIdentity for a readable model object in KBase
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        A morph with model to be filled using probabilistic gapfilling

    Returns
    -------
    Morph
        a morph such that morph.model, morph.ws_id forms an ObjectIdentity for a probalistic gapfilled model in KBase (to media)

    Examples
    --------
    Given a morph of the form:
        ws_name: MMws235
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.model = 214
        probannows: 9145
        probhash: ['rxn00001', 'rxn03434', 'rxn67834' ... (more)]
        probanno: 15
        morph.essential_ids = ['rxn10316_c0', 'rxn23434_c0', 'rxn78687_c0']
        probannows: 9145
        media: 18
        mediaws: 9145

    >>>morph = Client.probanno_fill(morph)

    Might return a morph of form
        ws_name: MMws235
        rxn_labels: ['gene-match', 'gene-no-match', 'no-gene', 'recon']
        ws_id: 11444
        morph.model = 216
        probannows: 9145
        probhash: ['rxn00001', 'rxn03434', 'rxn67834' ... (more)]
        probanno: 15
        morph.essential_ids = ['rxn10316_c0', 'rxn23434_c0', 'rxn78687_c0']
        probannows: 9145
        media: 18
        mediaws: 9145

    Where morph.model now holds the object_id of a model filled to morph.media using porbabilistic gapfilling
    """
    morph = copy.deepcopy(morph)
    # Gapfill the model and reset the model and modelws attributes of
    # the morph
    if name is None:
        name = u'morph-filled'
    fba_formulation = {'media': morph.media, 'media_workspace': morph.mediaws}
    gap_formulation = {'probabilisticAnnotation' : morph.probanno, 'probabilisticAnnotation_workspace' : morph.probannows, u'formulation': fba_formulation}
    params = {u'model': morph.model, u'model_workspace': morph.ws_id, u'out_model' : name, u'workspace' : morph.ws_id, u'formulation' : gap_formulation, u'integrate_solution' : True, u'gapFill' : u'gf'}
    model_info = fba_client.gapfill_model(params)
    morph.model = model_info[0]
    return morph

def remove_reactions(morph, rxn_list):
    """
    remove the reactions in rxn_list from morph.model

    Note
    ----
    Function Requirements:
    morph.model, morph.ws_id form a valid ObjectIdentity for a readable model object in KBase
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        the morph with the model reactions will be removed from
    rxn_list: list (String)
        a list of reaction_ids corresponding to reactions to be removed from the Morph

    Returns
    -------
    Morph
        a morph such that morph,model, morph.ws_id forms an ObjectIdentity for a KBase model object (the previous morph.model minus reactions in rxn_list)

    Examples
    --------
    Given a morph of the form (only relevant attributes shown):
        ws_id: 11444
        model: 214

    >>>morph = Client.remove_reactons(morph, ['rxn00001_c0', 'rxn11898_c0'])

    Might return a morph of the form:
        ws_id: 11444
        model: 216

    where morph.model, morph.ws_id forms an ObjectIdentity for a KBase model object (the previous model minus rxn00001_c0 and rxn11898_c0)
    """
    morph = copy.deepcopy(morph)
    morph.model  = fba_client.remove_reactions({'model': morph.model, 'model_workspace': morph.ws_id, 'workspace': morph.ws_id, 'output_id': 'removerxnsbydict', 'reactions': rxn_list})[0]
    return morph

def remove_reaction(morph, rxn):
    """
    remove the reaction rxn from morph.model

    Note
    ----
    Function Requirements:
    morph.model, morph.ws_id form a valid ObjectIdentity for a readable model object in KBase
    morph.ws_id is the ID of a writeable KBase workspace

    Parameters
    ----------
    morph: Morph
        the morph with the model reactions will be removed from
    rxn: String
        a reaction_id corresponding to a reaction to be removed from the Morph

    Returns
    -------
    Morph
        a morph such that morph,model, morph.ws_id forms an ObjectIdentity for a KBase model object (the previous morph.model minus reaction rxn)

    Examples
    --------
    Given a morph of the form (only relevant attributes shown):
        ws_id: 11444
        model: 214

    >>>morph = Client.remove_reactons(morph, 'rxn11898_c0')

    Might return a morph of the form:
        ws_id: 11444
        model: 216

    where morph.model, morph.ws_id forms an ObjectIdentity for a KBase model object (the previous model minus rxn11898_c0)
    """
    morph = copy.deepcopy(morph)
    morph.model  = fba_client.remove_reactions({'model': morph.model, 'model_workspace': morph.ws_id, 'workspace': morph.ws_id, 'output_id': 'removerxnsbydict', 'reactions': [rxn]})[0]
    return morph

def simple_probanno_morph(morph):
    morph = copy.deepcopy(morph)
    morph = prepare_supermodel(morph, fill_src=False)
    morph = remove_reactions(morph, rxn_list=morph.rxn_labels['gene-no-match'].keys())
    morph = remove_reactions(morph, rxn_list=morph.rxn_labels['no-gene'].keys())
    morph = probanno_fill(morph, name=u'simple_probanno_morph')
    return morph

def compare_morphmodels(list_of_morphs):
    """
    Returns a detailed comparison of the models in a list of models
    """
    args = {'models':[], 'workspaces': []}
    for i in list_of_morphs:
        args['models'].append(i.model)
        args['workspaces'].append(i.ws_id)
    return fba_client.compare_models(args)
def prob_annotate(genome_id, ws_id):
    print "to be implemented"
ws_client, fba_client = _init_clients()
def build_media(filename, ws_id, suppressError=False, objid=None, isMinimal=False):
    """
    Builds a media from a text file in tab delimited format

    Note
    ----
    Function Requirements:
    filename is a valid path to a properly formatted text file representing a media object
    ws_id corresponds to a KBase workspace you are authorized to write to

    Parameters
    ----------
    filename: String
        a String with the path to a txt file representing a media in tab-delimitted format. See MediaDB for an example
    ws_id: Integer
        a workspace_id for the media destination (can also be a string so long as the string represents an integer)
    suppressError: Boolean, optional
        a boolean flag for suppressing improperly formatted compound errors. Best used if one or more compounds does
        not have a kbase compound_id but they are NOT essential to the media. Default is False.
    objid: Integer, optional
        a object_id for where to save the media (optional, for workspace system control)
    isMinimal: Boolean, optional
        a boolean indicating whether the given media is minimal. Default is False

    Returns
    -------
    Tuple: (objid, wsid) (Integer, Integer)
        An object_id and workspace_id in a tuple that correspond to the ObjectIdentity of the media in KBase

    Examples
    --------
    Proper file format example:
        Compound    Concentration   Production_Bound(-)    Uptake_Bound(+)
        cpd00001[e0]    0.01    -1000.000000    1000.000000
        cpd00009[e0]    0.01    -1000.000000    1000.000000
        cpd00011[e0]    0.01    -1000.000000    1000.000000
        cpd00013[e0]    0.01    -1000.000000    1000.000000
        cpd00029[e0]    0.01    -1000.000000    1000.000000
        cpd00030[e0]    0.01    -1000.000000    1000.000000
        cpd00034[e0]    0.01    -1000.000000    1000.000000
    - First Row is headers for the columns below
    - Each row is tab-delimitted ('\\t'), the 4th term ending with a new new
    line

    First, put a media txt file of the above form in a directory sucl as (relative) '../media/mymedia.txt'

    >>>Client.build_media('../media/mymedia.txt', 9145)
    (15, 9145)

    >>>

    """
    compounds = list()
    with open(filename,'r') as file:
        cpds =  file.readlines()[1:]
        for c in cpds:
            item = dict()
            c = c.split('\n')[0]
            values = c.split('\t')
            if values[0][0:3] == "cpd":
                values[0] = values[0][0:8] + '_' + values[0][9:11]
                item[u'compound_ref'] = u'489/6/1/compounds/id/' + values[0]
                item[u'concentration'] = float(values[1])
                item[u'minFlux'] = float(values[2])
                item[u'maxFlux'] = float(values[3])
                compounds.append(item)
            else:
                if suppressError:
                    print values[0] + " not added to media!!"
                else:
                    raise MediaFormatError(values[0] + 'is not a properly formed KBase compound_id')
    name = filename.split('/')[-1].split('.')[0]
    obj = Helpers.load('blankmedia.pkl')
    obj =copy.deepcopy(obj)
    media = obj['data']
    media[u'name'] = name
    media[u'source_id'] = name
    if isMinimal:
        media['isMinimal'] = 1
    media[u'mediacompounds'] = compounds
    info = Helpers.save_object(media, obj['info'][2], ws_id, objid, name=name)[0]
    return (info[0], info[6])


class MediaFormatError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

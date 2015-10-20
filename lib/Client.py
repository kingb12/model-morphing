# import necesary services
from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import random
import argparse
import time
import traceback
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

def morph_model(morph):
    morph = translate_features(morph)
    morph = reconstruct_genome(morph)
    morph = get_objects(morph)
    morph = label_reactions(morph)
    morph = build_supermodel(morph)
    return morph

def fill_to_media(morph):
    """
    Takes a morph and produces a model filled to the provided media

    Fills a model to a given media. This is an important step if the source model does not grow initially on the media. NOTE: THIS OVERWRITES SRC_MODEL IN THE MORPH
    However, it does not overwrite the actual source model in KBase. It is imported to morph.ws_id, and morph.src_modelws is changed to match ws_id. Therefore, the data of the source
    model is preserved as is, but the reference is lost from the morph.

    Requires:
        morph.media, morph.mediaws is a valid ObjectIdentity for a media capable of allowing growth for the organisms in question.
           This is an important distinction that is left to the discretion of the user. If no media is known, complete can likely be used.
        morph.src_model, morph.src_modelws is a valid ObjectIdentity for a source model
    Returns:
        a morph. with morph.src_model set to a filled version of the source organism model. Guarantees that the resulting model can grow on media or
        an Exception is thrown declariing the source model and media incompatible
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

    Requires:
        morph.src_model, morph.src_modelws must form a valid ObjectIdentity. morph.protcomp/morph.protcompws must form a valid ObjectIdentity
    Modifies:
        morph - morph.trans_model
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
    Builds a draft reconstruction from the genome in the morph.

    Builds a draft model reconstruction and saves it in the workspace associated with morph.ws_id.

    Requires:
        morph.genome and morph.genomews must be valid ID's for KBase objects/workspace (respectively), forming a valid
        ObjectIdentity.
    Modifies:
        morph - morph.recon_model
    """
    morph = copy.deepcopy(morph)
    if (debug):
        assert morph.genome is not None, "genome is None"
        assert morph.genomews is not None, "genomews is None"
    recon_params = {'genome': morph.genome, 'genome_workspace': morph.genomews, 'workspace': morph.ws_id}
    morph.recon_model = fba_client.genome_to_fbamodel(recon_params)[0]
    return morph

def get_objects(morph):
    """
    Gets the model and probanno objects associated with the morph.

    Populates morph.objects with a dictionary of objects. Three are FBA Models, with keys 'source_model', 'recon_model', 'trans_model'.
    The last is a ProbAnno RxnProbs object, keyed 'probanno'. Corresponding entries are put in a dictionary morph.info, with the object_info
    as it is returned from the workspace service"

    Modifies:
        morph - morph.objects and morph.info

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
    return morph
def label_reactions(morph):
    """
    Labels the reactions in the translated model and the reconstruction

    Populates the rxn_labels attribute in the Morph object with a Dictionary of four dictionaries of reactrion_id -> value tuples.
    The first level dicts are named with the keys:
        - gene-match
        - gene-no-match
        - no-gene
        - recon

    Each inner dictionary is keyed with reaction_ids, hashed to tuples as such: (model_index, probability)
    Where reaction_id is a kbase reaction id with the compartment info appended to the end (e.g. rxn01316_c0),
    model_index is the index of the reaction in the objects['x_model'][modelreactions] list, and the probability
    is the reaction probability associated with each reaction from the probanno object. Reactions not in ProbAnno
    are given an arbitrary probability of -1.0

    Example uses:
    rxn_labels['gene-no-match']['rxn01316_c0'][0] evaluates to the index of rxn01316_c0 in morph.objects['trans_model']['modelreactions']
    rxn_labels['gene-no-match']['rxn01316_c0'][1] evaluates to the reaction probability of rxn01316_c0
    'rxn01316_c0' in rxn_labels['gene-no-match'] will evaluate True if the reaction is a gene-no-match reaction (an inner dict key)

    Requires: Morph.objects and Morph.info must have valid entries for 'source_model', 'trans_model', 'recon_model',
    and 'probanno' (The post condition of _get_objects()).

    Modifies:
        Morph - Morph.rxn_labels

    :param morph: the Morph object for labeling as described above
    """
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
    Sets morph.model to a superset of all reaction types in morph.rxn_labels

    Builds a model with the reactions unique to the target genome added and sets the morph.model attribute to this 'super-model'

    Requires:
    morph.rxn_labels must be populated to meet the post condition of _label_reactions(). morph.objects and morph,info must meet the post condition of _get_objects()

    Modifies:
    morph.model
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
    Generates a removal list from the import rxn_dict (e.g. morph.rxn_labels['no-gene'])

    Accepts a dictionary of reactions -> (index, probability) as is given by morph.rxn_labels[label] and returns a sorted
    list of tuples that can be passed as the rxn_list parameter to the _process_reactions function. This allows for processing
    of only some reactions from the model at a time. Tuples are in the form (reaction_id, index, probability). Setting the range keyword
    arg returns only a subset of the sorted list given by the index range passed in as a tuple (start, end).
    """
    # Sort by probanno. items() returns (K, V=(model_index, prob))
    def get_key(item):
        return item[1][1]
    removal_list = sorted(rxn_dict.items(), key=get_key)
    if list_range is not None:
        removal_list = removal_list[list_range[0]:list_range[1]]
    return removal_list

def process_reactions(morph, rxn_list=None, model_id=None, name='', ws=None, process_count=0, iterative_models=True):
    """
    Attempts removal of rxn_list reactions from morph (i.e. morph.rxn_labels[label])

    Attempts removal of each reaction , keeping it only if it is essential to the objective function
    of a given FBA formulation. Populates morph.essential_ids with the reaction_ids that could not be removed without breaking FBA
    simulation. morph.removed_ids is populated with the ids of the reactions that were removed. Both are populated as a dictionary of
    reaction_ids to their model_index and probability

    if label is a key to morph.rxn_labels AND rxn_list is None, the function behaves as if morph.rxn_labels[label] was the rxn_list argument

    Requires:
        -You must set keyword label to a valid key of morph.rxn_labels OR keyword rxn_list to a list of tuples
        - if rxn_labels[label] is not None, it must be a dictionary of reaction_ids -> (model_index, probability). Behavior uncertain if model_id
        and rxn_labels are not related to the same model as morph.

    Modifies:
       morph - morph.essential_ids, morph.removed_ids, morph.rxn_labels[label]
    Raises:
        KeyError if (label not in rxn_labels)
    """
    morph = copy.deepcopy(morph)
    if (ws is None):
        ws = morph.ws_id
    if (model_id is None):
        model_id = morph.model
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
        new_model_id = fba_client.remove_reactions({'model': model_id, 'model_workspace':ws, 'output_id': name + '-' + str(process_count), 'workspace':morph.ws_id, 'reactions':[reaction_id]})[0]
        if iterative_models:
            process_count += 1
        fba_params = {'fba': str(name + '-FBA-' + str(process_count)), 'workspace': ws, 'model' : new_model_id, 'model_workspace':ws}
        fbaMeta = fba_client.runfba(fba_params)
        flux = fbaMeta[-1]['Objective']
        print 'FBA Flux: ' + str(flux)
        if (flux > 0.0):
            # removed successfully
            print 'Removed ' + str(reaction_id)
            model_id = new_model_id
            morph.removed_ids[reaction_id] = removal_list[i]
        else:
            # essential
            print str(reaction_id) + ' is Essential'
            morph.essential_ids[reaction_id] = removal_list[i]
    morph.model = model_id
    return  morph, model_id, process_count

def find_alternative(reaction, formulation=None, morph=None, model_id=None, ws_id=None):
    if (morph is not None):
        morph = copy.deepcopy(morph)
        model_id = morph.model
        ws_id = morph.ws_id
    formulation = {'probabilisticAnnotation' : morph.probanno, 'probabilisticAnnotation_workspace' : morph.probannows}
    new_model_id = fba_client.remove_reactions({'model': model_id, 'model_workspace': ws_id, 'workspace': ws_id, 'output_id': 'findalt', 'reactions': [reaction]})[0]
    params = {'model': new_model_id, 'model_workspace': ws_id, 'workspace' : ws_id, 'formulation' : formulation, 'integrate_solution' : True, 'gapFill' : 'gf'}
    model_info = fba_client.gapfill_model(params)
    filled_model = model_info[0]
    object = ws_client.get_objects([{'name': model_info[0] + '.gffba', 'wsid': params['workspace']}])[0]
    gapfill = object['data']['gapfillingSolutions'][0]
    rxn_dicts = gapfill['gapfillingSolutionReactions']
    gap_reactions = [a['reaction_ref'].split('/')[-1] + '_' + a['compartment_ref'].split('/')[-1] + str(0) for a in rxn_dicts]
    modelreactions = set()
    for key in morph.rxn_labels:
        modelreactions |= set(morph.rxn_labels[key].keys())
    new_reactions = [r for r in gap_reactions if r not in modelreactions]
    return morph, params['gapFill'], new_reactions

# Finishing/Cleanup  Steps
def _finish(morph, save_ws=False):
    if not save_ws:
        ws_client.delete_workspace({'id': morph.ws_id})
    else:
        print 'Saved'

ws_client, fba_client =  _init_clients()

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

def analyze_alternative(morph, reaction, new_reactions):
    """
    Returns information related to what was added to the model when finding alternative reactions
    """
    results = dict()
    results['reaction_re-added'] = reaction in new_reactions
    results['removed_reactions'] = [r for r in new_reactions if r in morph.removed_ids]
    results['essential_reactions'] = [r for r in new_reactions if r in morph.essential_ids]
    return results

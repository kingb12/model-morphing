# import necesary services
from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import random
import argparse
import time
import traceback
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
    with open (".kbase_workspaceURL", "r") as myfile:
        url = myfile.read().replace('\n','')
    ws_client = Workspace(url)
    # Get FBA Model Services URL parameter
    with open (".kbase_fbaModelServicesURL", "r") as myfile:
        url = myfile.read().replace('\n','')
    fba_client = fbaModelServices(url)
    return ws_client, fba_client

def morph_model(args):
    morph = Morph()
    _label_reactions(morph)

def _translate_features(morph):
    """
    Translates the features in the source model to matches in the target genome

    NOTE: The translated model will lose reactions in the source model that did not have matching features
    and is not guaranteed to be useful in FBA simulation.

    Requires:
        morph.src_model, morph.src_modelws must form a valid ObjectIdentity. morph.protcomp/morph.protcompws must form a valid ObjectIdentity
    Modifies:
        morph - morph.trans_model
    """
    if (debug):
        assert morph.protcomp is not None, "protcomp is None"
        assert morph.protcompws is not None, "protcompws is None"
        assert morph.src_model is not None, "src_model is None"
        assert morph.src_modelws is not None, "src_modelws is None"
    trans_params = {'keep_nogene_rxn': 1, 'protcomp': morph.protcomp, 'protcomp_workspace': morph.protcompws, 'model': morph.src_model, 'model_workspace': morph.src_modelws, 'workspace': ws_id}
    morph.trans_model = fba_client.translate_fbamodel(trans_params)[0]

def _reconstruct_genome(morph):
    """
    Builds a draft reconstruction from the genome in the morph.

    Builds a draft model reconstruction and saves it in the workspace associated with morph.ws_id.

    Requires:
        morph.genome and morph.genomews must be valid ID's for KBase objects/workspace (respectively), forming a valid
        ObjectIdentity.
    Modifies:
        morph - morph.recon_model
    """
    if (debug):
        assert morph.genome is not None, "genome is None"
        assert morph.genomews is not None, "genomews is None"
    recon_params = {'genome': morph.genome, 'genome_workspace': morph.genomews, 'workspace': ws_id}
    morph.recon_model = fba_client.genome_to_fbamodel(recon_params)[0]

def _get_objects(morph):
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
    # Get the src_model
    obj  = ws_client.get_objects([{'objid': morph.src_model, 'wsid': morph.src_modelws}])[0]
    morph.objects['source_model'] = obj['data']
    morph.info['source_model'] = obj['info']
    # Get the translate model
    obj  = ws_client.get_objects([{'objid': morph.trans_model, 'wsid': ws_id}])[0]
    morph.objects['trans_model'] = obj['data']
    morph.info['trans_model'] = obj['info']
    # Get the recon model
    obj  = ws_client.get_objects([{'objid': morph.recon_model, 'wsid': ws_id}])[0]
    morph.objects['recon_model'] = obj['data']
    morph.info['recon_model'] = obj['info']
    # Get the probanno
    obj  = ws_client.get_objects([{'objid': morph.probanno, 'wsid': morph.probannows}])[0]
    morph.objects['probanno'] = obj['data']
    morph.info['probanno'] = obj['info']

def _label_reactions(morph):
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
    label_time = time.time()
    # create local vars for parts of the morph
    trans_model = morph.objects['trans_model']
    model = morph.objects['source_model']
    recon = morph.objects['recon_model']
    probanno = morph.objects['probanno']
    # create the rxn_labels dictionary 
    rxn_labels = {'gene-no-match': dict(), 'gene-match': dict(), 'no-gene': dict(), 'recon': dict()}
    # Build a hash of rxn_ids to probability to save look up time
    probanno_hash = dict()
    for rxn in probanno['reaction_probabilities']:
        prob_hash[rxn[0]] = rxn[1]
    # label gene-match and no-gene from translation
    for i in range(len(trans_model['modelreactions'])):
        # e.g. rxn_id = rxn01316_c0
        mdlrxn = trans_model['modelreactions'][i]
        rxn_id = mdlrxn['reaction_ref'].split('/')[-1] + '_' + str(mdlrxn['modelcompartment_ref'].split('/')[-1]) # -1 index gets the last in the list
        rxn_prot = mdlrxn['modelReactionProteins']
        if (len(rxn_prot) == 1 and len(rxn_prot[0]['modelReactionProteinSubunits']) == 0 and rxn_prot[0]['note'] != 'spontaneous' and rxn_prot[0]['note'] != 'universal'):
            try:
                rxn_labels['no-gene'][rxn_id] = (i, prob_hash[rxn_id])
            except KeyError:
                rxn_labels['no-gene'][rxn_id] = (i, -1.0)
        else:
            try:
                rxn_labels['gene-match'][rxn_id] = (i, prob_hash[rxn_id])
            except KeyError:
                rxn_labels['gene-match'][rxn_id] = (i, -1.0)
    # Label gene-no-match as those unique to model compared with translation
    for i in range(len(model['modelreactions'])):
        mdlrxn = model['modelreactions'][i]
        rxn_id = mdlrxn['reaction_ref'].split('/')[-1] + '_' + str(mdlrxn['modelcompartment_ref'].split('/')[-1]) # -1 index gets the last in the list
        if (rxn_id not in rxn_labels['gene-match']) and (rxn_id not in rxn_labels['no-gene']):
            try:
                rxn_labels['gene-no-match'][rxn_id] = (i, prob_hash[rxn_id])
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

def _build_supermodel(morph):
    """
    Sets morph.model to a superset of all reaction types in rxn_labels

    Builds a model with the reactions unique to the target genome added and sets the morph.model attribute to this 'super-model'

    Requires:
    morph.rxn_labels must be populated to meet the post condition of _label_reactions(). morph.objects and morph,info must meet the post condition of _get_objects()

    Modifies:
    morph.model
    """
    super_time = time.time()
    model = morph.objects['source_model']
    super_rxns = list()
    # Add the GENE_NO_MATCH reactions:
    for rxn_id in rxn_labels['gene-no-match']:
        # rxn_labels['gene-no-match'][0] gives the index of the reaction in the model['modelreactions'] list to make this look up O(1) instead of O(n)
        reaction = model['modelreactions'][rxn_id[0]]
        #TODO: See what more you can fix/add (gpr?)
        super_rxns.append((reaction['reaction_ref'].split('/')[-1], str(reaction['modelcompartment_ref'].split('/')[-1][0]), reaction['direction'], 'GENE_NO_MATCH', '', reaction['name']))
    # Add the recon reactions:
    for rxn_id in rxn_labels['recon']:
        # rxn_labels['recon'][0] gives the index of the reaction in the model['modelreactions'] list to make this look up O(1) instead of O(n)
        reaction = model['modelreactions'][rxn_id[0]]
        #TODO: See what more you can fix/add (gpr?)
        super_rxns.append((reaction['reaction_ref'].split('/')[-1], str(reaction['modelcompartment_ref'].split('/')[-1][0]), reaction['direction'], 'recon', '', reaction['name']))
    morph.model = fba_client.add_reactions({'model': trans_model_id, 'model_workspace': ws_id, 'output_id': 'super_model', 'workspace': ws_id, 'reactions': super_rxns})[0]
ws_client, fba_client  _init_clients()

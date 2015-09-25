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

def _label_reactions(morph):
    label_time = time.time()
    # translate the model using Kbase translation (This translates the features
    # for the model as well)
    # Get required objects from workspaces
    object = ws_client.get_objects([{'objid': morph.model, 'wsid': morph.modelws}])[0]
    morph.objects['source_model'] = object['data']
    morph.info['source_model'] = object['info']
    translation_params = {'keep_nogene_rxn': 1, 'protcomp': morph.protcomp, 'protcomp_workspace': morph.protcompws, 'model': morph.model, 'model_workspace': morph.modelws, 'workspace': ws_id}
    morph.trans_model =  fba_client.translate_fbamodel(translation_params)[0]
ws_client, fba_client = _init_clients()

# import necesary services
from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import random
import argparse
import time
import traceback
from operator import itemgetter

class Morph():

    def __init__(self, model, genome, protcomp=None, probanno=None, modelws=None, genomews=None, protcompws=None, probannows=None):
        self.args = _parse_arguments(model, genome, protcomp, probanno, modelws, genomews, protcompws, probannows)
        # rxn_labels, essentials, removed, _id_hash. Make properties?
        self.ws_id = _init_workspace()
        self.save_ws = False
        rxn_labels = {'gene-match' :


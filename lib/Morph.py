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

    # These are the allowed properties of a morph object. VAlues not specified
    # by user are set to None
    # This is an Abstract Object representing the Morph of a metabolic model
    # from one species to a close relative genome. It has information related to
    # the source model, the target genome, the reactions in the model as it is
    # morphed from source import to target, and the set of
    ws_client, fba_client = _init_clients()
    properties = set(['model', 'modelws', 'genome', 'genomews', 'probanno', 'probannows', 'protcomp', 'protcompws', 'rxn_labels', 'objects', 'info', 'essential_ids', 'removed_ids', 'ws_id', 'ws_name'])

    def __init__(self, *arg_hash, **kwargs):
          for dictionary in arg_hash:
              for key in dictionary:
                  if (key in Morph.properties):
                      setattr(self, key, dictionary[key])
          for key in kwargs:
              if (key in Morph.properties):
                  setattr(self, key, kwargs[key])
          for prop in Morph.properties:
              if (not hasattr(self, prop)):
                  setattr(self, prop, None)
          if (self.ws_id is None):
              self.ws_id, self.ws_name = _init_workspace()

    def _init_workspace(ws = None):
        global ws_id
        global ws_name
        ws_id = ws
        ws_name = 'MMws'
        if (ws is None):
            ws_conflict = True
            while (ws_conflict):
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
                     ws_name += str(random.randint(1,9))
        return ws_id, ws_name

    # Remove if unused
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

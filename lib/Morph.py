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
    properties = set('model', 'modelws', 'genome', 'genomews', 'probanno', 'probannows', 'protcomp', 'protcompws', 'rxn_labels')
    def __init__(self, *arg_hash, **kwargs):
          for dictionary in arg_hash:
              for key in dictionary:
                  if (key in Morph.properties):
                      setattr(self, key, dictionary[key])
          for key in kwargs:
              if (key in Morph.properties):
                  setattr(self, key, kwargs[key])
          for prop in Morph.properties:
              if (!hasattr(self, prop)):
                  print prop
                  setattr(self, prop, None)

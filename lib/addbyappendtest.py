from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.GenomeComparison.Client import GenomeComparison
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.userandjobstate.client import UserAndJobState

import copy
import random
import sys
import argparse
import time
from operator import itemgetter

with open (".kbase_fbaModelServicesURL", "r") as myfile:
    url = myfile.read().replace('\n','')
fba_client = fbaModelServices(url)
ws_client = Workspace()
modelx = ws_client.get_objects([{'wsid' : '9145', 'objid' : '9'}])[0]
model = modelx['data']
modtype = modelx['info']
model['modelreactions'].append({'direction': '<', 'name': 'N-Acetylneuraminate pyruvate-lyase (pyruvate-phosphorylating)_c0', 'probability': 0, 'modelcompartment_ref': '~/modelcompartments/id/c0', 'protons': 0, 'modelReactionReagents': [{'coefficient': -1, 'modelcompound_ref': '~/modelcompounds/id/cpd00067_c0'}, {'coefficient': 1, 'modelcompound_ref': '~/modelcompounds/id/cpd00001_c0'}, {'coefficient': -1, 'modelcompound_ref': '~/modelcompounds/id/cpd00009_c0'}, {'coefficient': 1, 'modelcompound_ref': '~/modelcompounds/id/cpd00492_c0'}, {'coefficient': -1, 'modelcompound_ref': '~/modelcompounds/id/cpd00232_c0'}, {'coefficient': 1, 'modelcompound_ref': '~/modelcompounds/id/cpd00061_c0'}], 'reaction_ref': '489/6/1/reactions/id/rxn01316', 'aliases': [], 'id': 'rxn01316_c0', 'modelReactionProteins': [{'note': '', 'complex_ref': '489/5/2/complexes/id/mscpx.741', 'modelReactionProteinSubunits': [{'note': '', 'feature_refs': ['2150/2468/1/features/id/kb|g.186.peg.2209'], 'role': 'N-acetylneuraminate synthase (EC 2.5.1.56)', 'triggering': 1, 'optionalSubunit': 0}]}]})
ws_client.save_objects({'workspace' : '9145', 'id' : '42', 'objects': [{'type': modtype, 'data' : model}]})
fba_client.runfba({'model' : '42', 'model_workspace' : '9145', 'workspace' : '9145', 'fba' : 'testFBAaddReactionAPPENDED'})

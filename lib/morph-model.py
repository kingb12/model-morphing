# This is the script that initiates the model morphing function. Interfacing to be addressed soon. TODO: Fix this description as you determine what this file will actioally do

#import necesary services
from biokbase.workspace.client import Workspace 
from biokbase.workspace.client import ServerError 
from clients.GenomeComparisonClient import GenomeComparison

# TODO remove these imports
import list_objs_by_ID
# from biokbase.workspace.ScriptHelpers import user_workspace

import random
import sys
import argparse
desc1 = '''
NAME
	 mm-morphmodel - 'Morph' a model an existing model to fit a related genome

SYNOPSIS
'''
desc2 = '''
DESCRIPTION
	Takes a Model/Genome for a 'source' organism and 'morphs' it to the 'target' organism. It keeps reactions from the 'source' model for which there is at least one matching gene feature, and attempts to remove those in the source model that do not have BLAST hits in the second genome. It also adds reactions annotated within the target genome that are unique to the target organism. The model is never broken in the process, such that the Biomass reaction and any other specified reaction [insert how you flag this option] must always carry flux. If a reaction can't be removed without breaking the model, it is flagged again and a few options are available for handling these reactions [insert said options] 
'''
	#TODO: Fill in EXAMPLES section with actual model/genomek
desc3 = '''
EXAMPLES
      Generate 'morphed' model for : 
      > mm-morph-model kb|m.1.model kb|g.0.genome 
AUTHORS
	Brendan King, Vangelis Simeonidis, Matt Richards
'''
#Parse incoming arguments
	#TODO: replace sys.argv with appropriate replacement from bash script interface
parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, prog='mm-morph-model', epilog=desc3)	
parser.add_argument('genome', help='ID of the Genome object', action='store', default=None)
parser.add_argument('model', help='ID of the Model object', action='store', default=None)
parser.add_argument('--genomews', help='Workspace of the Genome object', action='store', default=None, dest='genomews')
parser.add_argument('--modelws', help='Workspace of the Model object', action='store', default=None, dest='modelws')
	#TODO: ADD OTHER OPTION ARGUMENTS
usage = parser.format_usage()
parser.description = desc1 + '	' + usage + desc2
parser.usage = argparse.SUPPRESS
input_args = parser.parse_args()

#Prepare the argument dictionary
args = dict()
args['genome'] = input_args.genome
args['model'] = input_args.model
if input_args.genomews is None:
#FIXME: fix this functionality
# 	args['genomews'] = user_workspace()
	args['genomews'] = 8730
else:
	args['genomews'] = input_args.genomews
if input_args.modelws is None:
#	args['modelws'] = user_workspace()
	args['modelws'] = 8730
else:
	args['modelws'] = input_args.modelws

# some global vars
ws_id = None
ws_name = 'MMws'

# initiate a new workspace
ws_client = Workspace()
ws_conflict = True
while (ws_conflict):
	create_ws_params = {'workspace' : ws_name}
	# Try to create a workspace, catch an error if the name is already in use
	try:
		new_ws = ws_client.create_workspace(create_ws_params)
		# new_ws is type workspace_info, a tuple where 0, 1 are id, name
		ws_id = new_ws[0]
		ws_name = new_ws[1]
		ws_conflict = False
	except ServerError:
		ws_name += str(random.randint(1,9))

# Copy objs from current ws to temp one (for argument simplicity to later method calls)
	#TODO: verify this does not significantly impact performance
	ws_client.copy_object({'ws_id' : args['genomews'], 'obj_id' : args['genome']}, {'ws_id' : ws_id, 'obj_id' : args['genome']}) #ObjectIdentity old => ObjectIdentity new
	ws_client.copy_object({'ws_id' : args['modelws'], 'obj_id' : args['model']}, {'ws_id' : ws_id, 'obj_id' : args['model']}) #ObjectIdentity old => ObjectIdentity new
	list_objs_by_ID()
#Prepare Proteome Comparison
gencomp_client = GenomeComparison()
	# Set up parameters
blast_proteomes_params = dict()
blast_proteomes_params['genome1ws'] = ws_id
blast_proteomes_params['genome2ws'] = ws_id
blast_proteomes_params['genome1id'] = args['genome'] #genome 1 = the input genome


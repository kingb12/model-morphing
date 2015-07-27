# This is the script that initiates the model morphing function. Interfacing to be addressed soon. TODO: Fix this description as you determine what this file will actioally do

#import necesary services
from biokbase.workspace.client import Workspace 
from biokbase.workspace.client import ServerError 
from clients.GenomeComparisonClient import GenomeComparison
from biokbase.fbaModelServices.Client import fbaModelServices
# TODO remove these imports
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

# =================================================================================================
# 						Functions
# 
# functions used in the algorithm script TODO: make these private
# =================================================================================================

#Prepare Proteome Comparison
def blast_proteomes(args):
		# Set up parameters
	blast_proteomes_params = dict()
	blast_proteomes_params['genome1ws'] = args['genomews']
	blast_proteomes_params['genome2ws'] = args['modelws']
	blast_proteomes_params['genome1id'] = args['genome'] #genome 1 = the input genome
	get_models_params = dict()
	get_models_params['models'] = list() 
	get_models_params['models'].append(args['model']) 
	get_models_params['workspaces'] = list()
	get_models_params['workspaces'].append(args['modelws']) 
	print args
	print ws_client.get_workspace_info({'id' : ws_id})
	print ws_client.list_objects({'ids' : [ws_id]})
	model = fba_client.get_models(get_models_params)
	# model = fba_client.get_models({'models' : [args['model']], 'workspaces' : args['modelws']})
	# model.
	
# Parses Command Line arguments and TODO: assigns all values to ids for ease of use
def parse_arguments():
	#TODO: replace sys.argv with appropriate replacement from bash script interface
	#FIXME: make it so arguments can be passed as names, then find a way to convert interior to IDs
	parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, prog='mm-morph-model', epilog=desc3)	
	parser.add_argument('genome', type=int,  help='ID of the Genome object', action='store', default=None)
	parser.add_argument('model', type=int, help='ID of the Model object', action='store', default=None)
	parser.add_argument('--genomews', type=int, help='Workspace of the Genome object', action='store', default=None, dest='genomews')
	parser.add_argument('--modelws', type=int, help='Workspace of the Model object', action='store', default=None, dest='modelws')
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
	return args

#Initiate Clients Objects for Function
def init_clients():
	clients = dict()
	clients['ws'] = Workspace()
		# Get FBA Model Services URL parameter
	with open (".kbase_fbaModelServicesURL", "r") as myfile:
		url = myfile.read().replace('\n','')
	clients['fba'] = fbaModelServices(url)
		# Get Genome Comparison URL parameter
	with open (".kbase_genomecomparisonURL", "r") as myfile:
		url = myfile.read().replace('\n','')
	clients['gencomp'] = GenomeComparison(url)
	return clients

# initiate MMws workspace and clone in all needed objects
def init_workspace():
	ws_conflict = True
	global ws_name
	global ws_id
	ws_name = 'MMws'
	while (ws_conflict):
		clone_ws_params = {'workspace' : ws_name, 'wsi' : {'id' : args['genomews']}, 'globalread' : 'r'}
		# Try to create a workspace, catch an error if the name is already in use
		try:
			new_ws = ws_client.clone_workspace(clone_ws_params)
			# new_ws is type workspace_info, a tuple where 0, 1 are id, name
			ws_id = new_ws[0]
			ws_name = new_ws[1]
			ws_conflict = False
		except ServerError:
			 ws_name += str(random.randint(1,9))

# Finishing/Cleanup  Steps 
def finish():
	ws_client.delete_workspace({'id' : ws_id})

# =================================================================================================
# 						Script
#
# the scripted algorithm steps in order
# =================================================================================================

#parse command args
args = parse_arguments()

#initiate clients
clients = init_clients()
fba_client = clients['fba']
ws_client = clients['ws']
gencomp_client = clients['gencomp']

#initiate Model Morphing workspace
init_workspace() # creates global vars ws_name and ws_id

#Blast Proteomes
blast_proteomes(args)
# Clean up/Finish
finish()

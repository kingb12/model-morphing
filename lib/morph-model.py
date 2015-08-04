# This is the script that initiates the model morphing function. Interfacing to be addressed soon. TODO: Fix this description as you determine what this file will actioally do

#import necesary services
from biokbase.workspace.client import Workspace 
from biokbase.workspace.client import ServerError 
from biokbase.GenomeComparison.Client import GenomeComparison
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.userandjobstate.client import UserAndJobState
# TODO remove these imports
# from biokbase.workspace.ScriptHelpers import user_workspace

import random
import sys
import argparse
import time
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

# Build a model composed of ALL reactions (gene matching, non matching, no-gne, and recon rxns)
def build_supermodel(): #model, recon, trans_model, rxn_labels, id_hash):
	tuplist = list()
	print id_hash.keys()
	for rxn_id in rxn_labels['gene-no-match']:
		print rxn_id
		tuplist.append( get_add_rxn_tuple(model['reactions'][id_hash['model'][rxn_id]]))
	print model.keys()
	for i in model.keys():
		if type(model[i]) is list and len(model[i]) > 0:
			print model[i][0]
		else:
			print model[i]
	rxn_keys = model['reactions'][1].keys()
	print rxn_keys
	for i in rxn_keys:
		print model['reactions'][1][i]
# 	fba_client.add_reactions({'model' : trans_model_id, 'model_workspace' : ws_id, 'output_id' : '7', 'workspace' : ws_id, 'reactions' : tuplist})
	print 'added no match'
#	print fba_client.compare_models({'models' : [args['model'], trans_model_id], 'workspaces' : [args['modelws'], ws_id]})
# label reactions in each model
def label_reactions(model, recon):
	# Dictionaries for ids => model.reactions list indices
	model_rxn_ids = dict()
	trans_model_rxn_ids = dict()
	recon_rxn_ids = dict()
	# Dictionary for all rxn ids and their 'labels'
	rxn_labels = dict()
	rxn_labels['recon'] = set()
	rxn_labels['common'] = set() #Check: Is this a subset of gene-match (should be)
	rxn_labels['gene-match'] = set() 
	rxn_labels['gene-no-match'] = set() 
	rxn_labels['no-gene'] = set() 
	#Build a dictionary of rxn_ids to their index in the list so future look ups can be run in constant-time instead of O(n)
	trans_model_id = fba_client.translate_fbamodel({'keep_nogene_rxn' : 1, 'protcomp' : args['protcomp'], 'protcomp_workspace' : args['protcompws'], 'model' : args['model'], 'model_workspace' : args['modelws'], 'workspace' : ws_id})[0]
	trans_model = fba_client.get_models({'models' : [trans_model_id], 'workspaces' : [ws_id]})[0]
	for i in range(len(trans_model['reactions'])):
		trans_model_rxn_ids[trans_model['reactions'][i]['reaction']] = i
		mdlrxn = trans_model['reactions'][i]
		if len(mdlrxn['features']) == 0:
			rxn_labels['no-gene'].add(mdlrxn['reaction'])
		else: 
			rxn_labels['gene-match'].add(mdlrxn['reaction'])
	for j in range(len(model['reactions'])):
		model_rxn_ids[model['reactions'][j]['reaction']] = j
		mdlrxn = model['reactions'][j]
		if mdlrxn['reaction'] not in trans_model_rxn_ids:
			rxn_labels['gene-no-match'].add(mdlrxn['reaction'])
	for k in range(len(recon['reactions'])):
		recon_rxn_ids[recon['reactions'][k]['reaction']] = k
		mdlrxn = recon['reactions'][k]
		# if the recon reaction is already in the model
		if mdlrxn['reaction'] in model_rxn_ids:
			rxn_labels['common'].add(mdlrxn['reaction'])
		else:
			rxn_labels['recon'].add(mdlrxn['reaction'])
	with open("mmlog.txt", "a") as log:
		log.write("MODEL LABELING STATISTICS: ")
		log.write(str(len(model_rxn_ids)) + ' model reactions')
		log.write(str(len(trans_model_rxn_ids)) + ' translated model reactions')
		log.write(str(len(rxn_labels['recon'])) + ' reconstructed reactions')
		log.write(str(len(rxn_labels['gene-no-match'])) + ' no-match reactions')
		log.write(str(len(rxn_labels['gene-match'])) + ' gene-match reactions')
		log.write(str(len(rxn_labels['no-gene'])) + ' no-gene reactions')
	ids = {'model' : model_rxn_ids, 'trans' : trans_model_rxn_ids, 'recon' : recon_rxn_ids}
	return trans_model, rxn_labels, ids
def get_add_rxn_tuple(rxn):
	return (rxn['reaction'], rxn['compartment'], rxn['direction']) 
# Parses Command Line arguments and TODO: assigns all values to ids for ease of use
def parse_arguments():
	#TODO: replace sys.argv with appropriate replacement from bash script interface
	#FIXME: make it so arguments can be passed as names, then find a way to convert interior to IDs
	parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, prog='mm-morph-model', epilog=desc3)	
	parser.add_argument('model', type=int, help='ID of the Model object', action='store', default=None)
	parser.add_argument('genome', type=int,  help='ID of the Genome object', action='store', default=None)
	parser.add_argument('protcomp', type=int,  help='ID of the Proteome Comparison object', action='store', default=None)
	parser.add_argument('--genomews', type=int, help='Workspace of the Genome object', action='store', default=None, dest='genomews')
	parser.add_argument('--modelws', type=int, help='Workspace of the Model object', action='store', default=None, dest='modelws')
	parser.add_argument('--protcompws', type=int, help='Workspace of the Proteome Comparison object', action='store', default=None, dest='protcompws')
		#TODO: ADD OTHER OPTION ARGUMENTS
	usage = parser.format_usage()
	parser.description = desc1 + '	' + usage + desc2
	parser.usage = argparse.SUPPRESS
	input_args = parser.parse_args()
	
	#Prepare the argument dictionary
	args = dict()
	args['genome'] = input_args.genome
	args['model'] = input_args.model
	args['protcomp'] = input_args.protcomp
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
	if input_args.protcompws is None:
	#	args['protcompws'] = user_workspace()
		args['protcompws'] = 8730
	else:
		args['protcompws'] = input_args.protcompws
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
	clients['ujs'] = UserAndJobState()
	return clients

# initiate MMws workspace and clone in all needed objects
def init_workspace():
	ws_conflict = True
	global ws_name
	global ws_id
	ws_name = 'MMws'
	while (ws_conflict):
		create_ws_params = {'workspace' : ws_name, 'globalread' : 'r', 'description' : "A workspace for storing the FBA's and meta data of the algorithm"}
		# Try to create a workspace, catch an error if the name is already in use
		try:
			new_ws = ws_client.create_workspace(create_ws_params)
			# new_ws is type workspace_info, a tuple where 0, 1 are id, name
			ws_id = new_ws[0]
			ws_name = new_ws[1]
			ws_conflict = False
		except ServerError:
			 ws_name += str(random.randint(1,9))

# Get the reactions for the comparison 'recon' model in Genome B
def build_models():
	get_models_params = {'models' : [args['model']], 'workspaces' : [args['modelws']]}
	model = fba_client.get_models(get_models_params)[0]
	recon_params = {'genome' : args['genome'], 'genome_workspace' : args['genomews'], 'workspace' : ws_id}
	recon_id = fba_client.genome_to_fbamodel(recon_params)[0]
	get_models_params = {'models' : [recon_id], 'workspaces' : [ws_id]}
	recon = fba_client.get_models(get_models_params)[0]
	trans_model_id = fba_client.translate_fbamodel({'keep_nogene_rxn' : 1, 'protcomp' : args['protcomp'], 'protcomp_workspace' : args['protcompws'], 'model' : args['model'], 'model_workspace' : args['modelws'], 'workspace' : ws_id})[0]
	trans_model = fba_client.get_models({'models' : [trans_model_id], 'workspaces' : [ws_id]})[0]
	return [model, recon, trans_model, trans_model_id]

# label reactions in each model
def label_reactions(): #model, recon, trans_model):
	# Dictionaries for ids => model.reactions list indices
	model_rxn_ids = dict()
	trans_model_rxn_ids = dict()
	recon_rxn_ids = dict()
	# Dictionary for all rxn ids and their 'labels'
	rxn_labels = dict()
	rxn_labels['recon'] = set()
	rxn_labels['common'] = set() #Check: Is this a subset of gene-match (should be)
	rxn_labels['gene-match'] = set() 
	rxn_labels['gene-no-match'] = set() 
	rxn_labels['no-gene'] = set() 
	#Build a dictionary of rxn_ids to their index in the list so future look ups can be run in constant-time instead of O(n)
	for i in range(len(trans_model['reactions'])):
		trans_model_rxn_ids[trans_model['reactions'][i]['reaction']] = i
		mdlrxn = trans_model['reactions'][i]
		if len(mdlrxn['features']) == 0:
			rxn_labels['no-gene'].add(mdlrxn['reaction'])
		else: 
			rxn_labels['gene-match'].add(mdlrxn['reaction'])
	for j in range(len(model['reactions'])):
		model_rxn_ids[model['reactions'][j]['reaction']] = j
		mdlrxn = model['reactions'][j]
		if mdlrxn['reaction'] not in trans_model_rxn_ids:
			rxn_labels['gene-no-match'].add(mdlrxn['reaction'])
	for k in range(len(recon['reactions'])):
		recon_rxn_ids[recon['reactions'][k]['reaction']] = k
		mdlrxn = recon['reactions'][k]
		# if the recon reaction is already in the model
		if mdlrxn['reaction'] in model_rxn_ids:
			rxn_labels['common'].add(mdlrxn['reaction'])
		else:
			rxn_labels['recon'].add(mdlrxn['reaction'])
	with open(".mmlog.txt", "a") as log:
		log.write("MODEL LABELING STATISTICS: ")
		log.write(str(len(model_rxn_ids)) + ' model reactions')
		log.write(str(len(trans_model_rxn_ids)) + ' translated model reactions')
		log.write(str(len(rxn_labels['recon'])) + ' reconstructed reactions')
		log.write(str(len(rxn_labels['gene-no-match'])) + ' no-match reactions')
		log.write(str(len(rxn_labels['gene-match'])) + ' gene-match reactions')
		log.write(str(len(rxn_labels['no-gene'])) + ' no-gene reactions')
	ids = {'model' : model_rxn_ids, 'trans' : trans_model_rxn_ids, 'recon' : recon_rxn_ids}
	return rxn_labels, ids

# Finishing/Cleanup  Steps 
def finish():
	if ws_id is not None:
		ws_client.delete_workspace({'id' : ws_id})
	else:
		print "ws_id" is None

# =================================================================================================
# 						Script
#
# the scripted algorithm steps in order
# =================================================================================================
# init variables
model = None
recon = None
ws_id = None
try:
	#parse command args
	print 'parsing args'
	args = parse_arguments()
	#initiate clientslk
	print 'init clients...'
	clients = init_clients()
	fba_client = clients['fba']
	ws_client = clients['ws']
	gencomp_client = clients['gencomp']
	ujs_client = clients['ujs']
	#initiate Model Morphing workspace
	print 'init ws...'
	init_workspace() # creates global vars ws_name and ws_id
	# [args['protcomp'], args['protcompws']] = blast_proteomes()
	print 'build models'
	[model, recon, trans_model, trans_model_id] = build_models()
	print 'label rxns...'
	[rxn_labels, id_hash] = label_reactions() #model, recon)
	print 'build supermodel...'
	morphed_model = build_supermodel()#model, recon, trans_model)
finally:
	finish()
# Clean up/Finish

# This is the script that initiates the model morphing function. Interfacing to be addressed soon. TODO: Fix this description as you determine what this file will actioally do

#import necesary services
from biokbase.workspace.client import Workspace 
from biokbase.workspace.client import ServerError 
from biokbase.GenomeComparison.Client import GenomeComparison
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.userandjobstate.client import UserAndJobState

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

#Save model to geven workspace
def save_model(model, workspace, name, model_type):
	save_data = dict()
	save_data['type'] = model_type
	save_data['data'] = model
	save_data['name'] = name
#	ws_client.save_objects({'id' : workspace, 'objects' : [save_data]})
	print model['biomasses']
	print model['biomasses'][0]

# Parses Command Line arguments and TODO: assigns all values to ids for ease of use
def parse_arguments():
	#TODO: replace sys.argv with appropriate replacement from bash script interface
	#FIXME: make it so arguments can be passed as names, then find a way to convert interior to IDs
	parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, prog='mm-morph-model', epilog=desc3)	
	parser.add_argument('model', type=int, help='ID of the Model object', action='store', default=None)
	parser.add_argument('genome', type=int,  help='ID of the Genome object', action='store', default=None)
	parser.add_argument('protcomp', type=int,  help='ID of the Proteome Comparison object', action='store', default=None)
	parser.add_argument('probanno', type=int,  help='ID of the ProbAnno object', action='store', default=None)
	parser.add_argument('--genomews', type=int, help='Workspace of the Genome object', action='store', default=None, dest='genomews')
	parser.add_argument('--modelws', type=int, help='Workspace of the Model object', action='store', default=None, dest='modelws')
	parser.add_argument('--protcompws', type=int, help='Workspace of the Proteome Comparison object', action='store', default=None, dest='protcompws')
	parser.add_argument('--probannows', type=int, help='Workspace of the ProbAnno object', action='store', default=None, dest='probannows')
	parser.add_argument('--outputws', type=int, help='Workspace for the morphed Model object', action='store', default=None, dest='outputws')
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
	args['probanno'] = input_args.probanno
	if input_args.genomews is None:
	# FIXME: fix this functionality
	# 	args['genomews'] = user_workspace()
		args['genomews'] = 9145
	else:
		args['genomews'] = input_args.genomews
	if input_args.modelws is None:
	#	args['modelws'] = user_workspace()
		args['modelws'] = 9145
	else:
		args['modelws'] = input_args.modelws
	if input_args.protcompws is None:
	#	args['protcompws'] = user_workspace()
		args['protcompws'] = 9145
	else:
		args['protcompws'] = input_args.protcompws
	if input_args.probannows is None:
	#	args['probannows'] = user_workspace()
		args['probannows'] = 9145
	else:
		args['probannows'] = input_args.probannows
	if input_args.outputws is None:
	#	args['outputws'] = user_workspace()
		args['outputws'] = 9145
	else:
		args['outputws'] = input_args.outputws
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
	type_string = ws_client.get_object_info_new({'objects' : [{'objid' : args['model'], 'wsid' : args['modelws']}]})[0][2]
	print type_string
	return [model, recon, trans_model, trans_model_id, type_string]

# label reactions in each model
def label_reactions(): #model, recon, trans_model
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
		log.write("MODEL LABELING STATISTICS: \n")
		log.write(str(len(model_rxn_ids)) + ' model reactions\n')
		log.write(str(len(trans_model_rxn_ids)) + ' translated model reactions\n')
		log.write(str(len(rxn_labels['recon'])) + ' reconstructed reactions\n')
		log.write(str(len(rxn_labels['gene-no-match'])) + ' no-match reactions\n')
		log.write(str(len(rxn_labels['gene-match'])) + ' gene-match reactions\n')
		log.write(str(len(rxn_labels['no-gene'])) + ' no-gene reactions\n')
	ids = {'model' : model_rxn_ids, 'trans' : trans_model_rxn_ids, 'recon' : recon_rxn_ids}
	return rxn_labels, ids

# Build a model composed of ALL reactions (gene matching, non matching, no-gne, and recon rxns)
def build_supermodel(): #model, recon, trans_model, rxn_labels, id_hash
#Add the GENE_NO_MATCH reactions:
	i = 0 #an index for how many reactions are addded so provide new indices to the id hash
	orig_size = len(trans_model['reactions'])
	for rxn_id in rxn_labels['gene-no-match']:
		# id_hash['model][rxn key] gives the index of the reaction in the model['reactions'] list to make this look up O(1) instead of O(n)
		reaction = model['reactions'][id_hash['model'][rxn_id]]
		trans_model['reactions'].append(reaction)
		id_hash['trans'][rxn_id] = orig_size + i
		# assert trans_model['reactions'][id_hash['trans'][rxn_id]]['reaction'] is rxn_id #assert that value got added to end of the list and no changes occured 
		i += 1
	with open(".mmlog.txt", "a") as log:
		log.write('MODEL REACTION ADDITION STATISTICS: \n')
		log.write('Added ' +str(i) + '  gene-no-match reactions to translated model: ' + str(rxn_labels['gene-no-match']) + '\n')
#Add the RECON reactions:
	i = 0 #an index for how many reactions are addded so provide new indices to the id hash
	orig_size = len(trans_model['reactions'])
	for rxn_id in rxn_labels['recon']:
		# id_hash['model][rxn key] gives the index of the reaction in the model['reactions'] list to make this look up O(1) instead of O(n)
		reaction = recon['reactions'][id_hash['recon'][rxn_id]]
		trans_model['reactions'].append(reaction)
		id_hash['trans'][rxn_id] = orig_size + i
		assert trans_model['reactions'][id_hash['trans'][rxn_id]]['reaction'] is rxn_id #assert that value got added to end of the list and no changes occured 
		i += 1
	with open(".mmlog.txt", "a") as log:
		log.write('Added ' + str(i) + ' recon reactions to translated model: ' + str(rxn_labels['recon']) + '\n')
		log.write('SUPERMODEL STATE: \n')
		log.write('NAME: ' + trans_model['name'] + '\n')
		numftrs = 0
		for rxn in trans_model['reactions']:
			numftrs += len(rxn['features'])	
		log.write('REACTIONS: ' + str(len(trans_model['reactions'])))
		log.write('FEATURES: ' + str(numftrs))
	return trans_model

# Finishing/Cleanup  Steps 
def finish():
	with open('.mmlog.txt', 'r') as log:
		print log.read()
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
	[model, recon, trans_model, trans_model_id, model_type] = build_models()
	print 'label rxns...'
	[rxn_labels, id_hash] = label_reactions() #model, recon)
	print 'build supermodel...'
	morphed_model = build_supermodel()
	save_model(morphed_model, args['outputws'], 'MM-' + str(args['genome']), model_type)
finally:
	finish()
# Clean up/Finish

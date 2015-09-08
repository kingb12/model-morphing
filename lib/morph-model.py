# This is the script that initiates the model morphing function. Interfacing to be addressed soon. TODO: Fix this description as you determine what this file will actioally do

#IDEA: When findling more likely alternatives to gene-no-match and no-gene reactions, the new more likely pathway must be more likely than it's replacement AND and reactions that are now dead ends as a result of it's use

#import necesary services
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

desc1 = '''
NAME
     mm-morphmodel - 'Morph' a model an existing model to fit a related genome

SYNOPSIS
'''
desc2 = '''
DESCRIPTION
    Takes a Model/Genome for a 'source' organism and 'morphs' it to the 'target' organism. It keeps reactions from the 'source' model for which there is at least one matching gene feature, and attempts to remove those in the source model that do not have BLAST hits in the second genome. It also adds reactions annotated within the target genome that are unique to the target organism. The model is never broken in the process, such that the Biomass reaction and any other specified reaction [insert how you flag this option] must always carry flux. If a reaction can't be removed without breaking the model, it is flagged again and a few options are available for handling these reactions [insert said options]
'''
    #TODO: Fill in EXAMPLES section with actual model/genome
desc3 = '''
EXAMPLES
      Generate 'morphed' model for :
      > mm-morph-model kb|m.1.model kb|g.0.genome
AUTHORS
    Brendan King, Vangelis Simeonidis, Matt Richards
'''

# =================================================================================================
#                         Functions
#
# functions used in the algorithm script TODO: make these private
# =================================================================================================

#Create Reaction removal lists (these are ordered)
def removal_tuples(label_list):
    tup_list = list()
    for rxn in label_list:
        rxnprob = mm_ids[rxn] #THIS IS WHERE YOU ADD PROBANNO SUPPORT
        tup_list.append((rxn, mm_ids[rxn], rxnprob))
    removal_tuples = sorted(tup_list, key=itemgetter(2))
    return removal_tuples

# Process the reactions THIS METHOD IS IN DEVELOPMENT
def process_reactions(model_id, rxn_list, name = '', process_count=0, ws=None):
    # Call generate model stats if more info is wanted
    removed_ids = list()
    if ws is None:
        ws = ws_id
    for i in range(len(rxn_list)):
        print 'TO BE REMOVED: ' + str(rxn_list[i][0])
        new_model_id = fba_client.remove_reactions({'model' : model_id, 'model_workspace' : ws, 'output_id' : name + '-' + str(process_count), 'workspace' : ws_id, 'reactions' : [rxn_list[i][0]]})[0]
        fba_params = dict()
        fba_params['fba'] = name + '-FBA-' + str(process_count)
        fba_params['workspace'] = ws_id
        fba_params['model'] = new_model_id
        fba_params['model_workspace'] = ws
        print fba_params
        fbaMeta = fba_client.runfba(fba_params)
        flux = fbaMeta[-1]['Objective']
        print flux
        if (flux > 0):
                print "Removed Successfully"
                model_id = new_model_id
                # rxn12345_c0 in rxn_list => rxn12345 in removed_ids
                removed_ids.append(rxn_list[i][0].split('_')[0])
        else:
            print "Reaction is Essential. Not Removed"
        process_count += 1
    return model_id, process_count, removed_ids


# Parses Command Line arguments into an argument dictionary for further use.
# REQURED ARGUMENTS: model, genome, protcomp, probanno
def parse_arguments():
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
    parser.description = desc1 + '    ' + usage + desc2
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
    #     args['genomews'] = user_workspace()
        args['genomews'] = 9145
    else:
        args['genomews'] = input_args.genomews
    if input_args.modelws is None:
    #    args['modelws'] = user_workspace()
        args['modelws'] = 9145
    else:
        args['modelws'] = input_args.modelws
    if input_args.protcompws is None:
    #    args['protcompws'] = user_workspace()
        args['protcompws'] = 9145
    else:
        args['protcompws'] = input_args.protcompws
    if input_args.probannows is None:
    #    args['probannows'] = user_workspace()
        args['probannows'] = 9145
    else:
        args['probannows'] = input_args.probannows
    if input_args.outputws is None:
    #    args['outputws'] = user_workspace()
        args['outputws'] = 9145
    else:
        args['outputws'] = input_args.outputws
    return args

#Initiate Clients Objects for KBase API services, including workspace, KBaseFBAModeling, Probabilistic Annotation, and UserandJobState. All are configured to the production URL, https://kbase.us
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

# initiate MMws workspace for implementation side objects and pieces for analysis
def init_workspace(ws = None):
    global ws_id
    global ws_name
    ws_id = ws
    ws_name = 'MMws'
    if (ws is None):
        ws_conflict = True
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
    return ws_id, ws_name

# builds model data structures from the KBase API for each model object: The
# model to be morph-translated, a general translation template, and a
# reconstruction for the target genome.
def build_models():
    model = ws_client.get_objects([{'objid' : args['model'], 'wsid' : args['modelws']}])[0]
    trans_params = {'keep_nogene_rxn' : 1, 'protcomp' : args['protcomp'], 'protcomp_workspace' : args['protcompws'], 'model' : args['model'], 'model_workspace' : args['modelws'], 'workspace' : ws_id}
    trans_model_id = fba_client.translate_fbamodel(trans_params)[0]
    trans_model = ws_client.get_objects([{'objid' : trans_model_id, 'wsid' : ws_id}])[0]
    recon_params = {'genome' : 'Methanosarcina_barkeri_str._fusaro', 'genome_workspace' : args['genomews'], 'workspace' : ws_id}
    recon_id = fba_client.genome_to_fbamodel(recon_params)[0]
    recon = ws_client.get_objects([{'objid' : recon_id, 'wsid' : ws_id}])[0]
    return [model['data'], recon['data'], trans_model['data'], model['info'], recon['info'], trans_model['info'], trans_model_id]

# label reactions in each model triplet. Must be passed a Source model, a
# translated model, and a reconstruction (the output of build_models)
def label_reactions(): #model, recon, trans_model
    # Dictionaries for ids => model.reactions list indices
    model = ws_client.get_objects([{'objid' : args['model'], 'wsid' : args['modelws']}])[0]
    trans_params = {'keep_nogene_rxn' : 1, 'protcomp' : args['protcomp'], 'protcomp_workspace' : args['protcompws'], 'model' : args['model'], 'model_workspace' : args['modelws'], 'workspace' : ws_id}
    trans_model_id = fba_client.translate_fbamodel(trans_params)[0]
    trans_model = ws_client.get_objects([{'objid' : trans_model_id, 'wsid' : ws_id}])[0]
    recon_params = {'genome' : 'Methanosarcina_barkeri_str._fusaro', 'genome_workspace' : args['genomews'], 'workspace' : ws_id}
    recon_id = fba_client.genome_to_fbamodel(recon_params)[0]
    recon = ws_client.get_objects([{'objid' : recon_id, 'wsid' : ws_id}])[0]
    model_info = model['info']
    recon_info = recon['info']
    trans_info = trans_model['info']
    model = model['data']
    recon = recon['data']
    trans_model = trans_model['data']
    model_rxn_ids = dict()
    recon_rxn_ids = dict()
    trans_model_rxn_ids = dict()
    # Dictionary for all rxn ids and their 'labels'
    rxn_labels = dict()
    common_rxns = set() #Check: Is this a subset of gene-match (should be)
    rxn_labels['recon'] = set()
    rxn_labels['gene-match'] = set()
    rxn_labels['gene-no-match'] = set()
    rxn_labels['no-gene'] = set()
    #Build a dictionary of rxn_ids to their index in the list so future look ups can be run in constant-time instead of O(n)
    for i in range(len(trans_model['modelreactions'])):
        rxn_id = trans_model['modelreactions'][i]['reaction_ref'].split('/')[-1] + '_' + str(trans_model['modelreactions'][i]['modelcompartment_ref'].split('/')[-1]) #-1 index gets the last in the list
        trans_model_rxn_ids[rxn_id] = i
        mdlrxn = trans_model['modelreactions'][i]
        rxn_prot = mdlrxn['modelReactionProteins']
        if len(rxn_prot) == 1 and len(rxn_prot[0]['modelReactionProteinSubunits']) == 0 and rxn_prot[0]['note'] != 'spontaneous' and rxn_prot[0]['note'] != 'universal':
            rxn_labels['no-gene'].add(rxn_id)
        else:
            rxn_labels['gene-match'].add(rxn_id)
    for i in range(len(model['modelreactions'])):
        rxn_id = model['modelreactions'][i]['reaction_ref'].split('/')[-1] + '_' + str(model['modelreactions'][i]['modelcompartment_ref'].split('/')[-1]) #-1 index gets the last in the list
        model_rxn_ids[rxn_id] = i
        mdlrxn = model['modelreactions'][i]
        if rxn_id not in trans_model_rxn_ids:
                rxn_labels['gene-no-match'].add(rxn_id)
    for i in range(len(recon['modelreactions'])):
        rxn_id = recon['modelreactions'][i]['reaction_ref'].split('/')[-1] + '_' + str(recon['modelreactions'][i]['modelcompartment_ref'].split('/')[-1]) #-1 index gets the last in the list
        recon_rxn_ids[rxn_id] = i
        mdlrxn = recon['modelreactions'][i]
        # if the recon reaction is already in the model
        if rxn_id in model_rxn_ids:
            common_rxns.add(rxn_id)
        else:
            rxn_labels['recon'].add(rxn_id)
    # Log Results
    with open(".mmlog.txt", "a") as log:
        log.write("MODEL LABELING STATISTICS: \n")
        log.write(str(len(model_rxn_ids)) + ' model reactions\n')
        log.write(str(len(trans_model_rxn_ids)) + ' translated model reactions\n')
        log.write(str(len(rxn_labels['recon'])) + ' reconstructed reactions out of ' + str(len(recon['modelreactions'])) + '. ' + str(len(common_rxns)) + ' in common\n')
        log.write(str(len(rxn_labels['gene-no-match'])) + ' no-match reactions\n')
        log.write(str(len(rxn_labels['gene-match'])) + ' gene-match reactions\n')
        log.write(str(len(rxn_labels['no-gene'])) + ' no-gene reactions\n')
    #Return Results, include look-up ID dict for future use
    id_hash = {'model' : model_rxn_ids, 'trans' : trans_model_rxn_ids, 'recon' : recon_rxn_ids}
    return model, recon, trans_model, model_info, recon_info, trans_info, trans_model_id, rxn_labels, id_hash

# Build a model composed of ALL reactions (gene matching, non matching, no-gne, and recon rxns)
def build_supermodel(): #model, recon, trans_model, rxn_labels, id_hash
#Add the GENE_NO_MATCH reactions:
    i = 0 #an index for how many reactions are addded so provide new indices to the id hash
    orig_size = len(trans_model['modelreactions'])
    super_rxns = list()
    for rxn_id in rxn_labels['gene-no-match']:
        # id_hash['model][rxn key] gives the index of the reaction in the model['modelreactions'] list to make this look up O(1) instead of O(n)
        reaction = model['modelreactions'][id_hash['model'][rxn_id]]
        # aliasing destroys features in 'model'. This might be ok, but consider copying
        # eliminate features from reaction
        reaction['modelReactionProteins'] = [{'note':'Manually Specified GPR', 'complex_ref' : '', 'modelReactionProteinSubunits' : []}]
        super_rxns.append((reaction['reaction_ref'].split('/')[-1], str(reaction['modelcompartment_ref'].split('/')[-1][0]), reaction['direction'], 'GENE_NO_MATCH', '', reaction['name']))
        trans_model['modelreactions'].append(reaction)
        id_hash['trans'][rxn_id] = orig_size + i
        assert trans_model['modelreactions'][id_hash['trans'][rxn_id]]['reaction_ref'].split('/')[-1]+ '_' + trans_model['modelreactions'][id_hash['trans'][rxn_id]]['modelcompartment_ref'].split('/')[-1] == rxn_id #assert that value got added to end of the list and no changes occured
        i += 1
    with open(".mmlog.txt", "a") as log:
        log.write('MODEL REACTION ADDITION STATISTICS: \n')
        log.write('Added ' +str(i) + '  gene-no-match reactions to translated model: ' + str(rxn_labels['gene-no-match']) + '\n')
#Add the RECON reactions:
    i = 0 #an index for how many reactions are addded so provide new indices to the id hash
    orig_size = len(trans_model['modelreactions'])
    for rxn_id in rxn_labels['recon']:
        # id_hash['model][rxn key] gives the index of the reaction in the model['modelreactions'] list to make this look up O(1) instead of O(n)
        reaction = recon['modelreactions'][id_hash['recon'][rxn_id]]
        super_rxns.append((reaction['reaction_ref'].split('/')[-1], str(reaction['modelcompartment_ref'].split('/')[-1][0]), reaction['direction'], 'RECON', '', reaction['name']))
        trans_model['modelreactions'].append(reaction)
        id_hash['trans'][rxn_id] = orig_size + i
        assert trans_model['modelreactions'][id_hash['trans'][rxn_id]]['reaction_ref'].split('/')[-1]+ '_' + trans_model['modelreactions'][id_hash['trans'][rxn_id]]['modelcompartment_ref'].split('/')[-1] == rxn_id #assert that value got added to end of the list and no changes occured
        i += 1
    super_model_id = fba_client.add_reactions({'model': trans_model_id, 'model_workspace': ws_id, 'output_id': 'super_model', 'workspace': ws_id, 'reactions': super_rxns})[0]
    with open(".mmlog.txt", "a") as log:
        log.write('Added ' + str(i) + ' recon reactions to translated model: ' + str(rxn_labels['recon']) + '\n')
        log.write('SUPERMODEL STATE: \n')
        log.write('    NAME: ' + trans_model['name'] + '\n')
        numprots = 0
        for rxn in trans_model['modelreactions']:
            numprots += len(rxn['modelReactionProteins'])
        log.write('    REACTIONS: ' + str(len(trans_model['modelreactions'])))
        log.write('. PROTEINS: ' + str(numprots) + '\n')
    return trans_model, id_hash['trans'], super_model_id

#Save model to geven workspace
def save_model(model, workspace, name, model_type):
    save_data = {'type' : model_type, 'name' : name,  'data' : model}
    ws_client.save_objects({'id' : workspace, 'objects' : [save_data]})
    return {'model' : name, 'model_workspace' : workspace}

# Finishing/Cleanup  Steps
def finish(save_ws=False):
    with open('.mmlog.txt', 'r') as log:
        print 'Finished'
    if not save_ws:
        ws_client.delete_workspace({'id' : ws_id})
    else:
        print 'ERROR:  workspace is None'

# =================================================================================================
#                         Script
#
# the scripted algorithm steps in order
# =================================================================================================
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
    ws_id, ws_name = init_workspace(ws='10232')
    # [args['protcomp'], args['protcompws']] = blast_proteomes()
    print 'label reactions...'
    [model, recon, trans_model, model_info, recon_info, trans_info, trans_model_id, rxn_labels, id_hash] = label_reactions()
    print 'build supermodel...'
    morphed_model, mm_ids, super_model_id = build_supermodel()
    print 'get reaction removal lists...'''
    gene_no_match_tuples = removal_tuples(rxn_labels['gene-no-match'])
    no_gene_tuples = removal_tuples(rxn_labels['no-gene'])
    # info[2] is 'type'
    model_id = save_model(model, ws_id, 'MM-0', model_info[2])
    print 'process reactions...'
    removed_ids = list()
    super_model_id, gnm, removed = process_reactions(super_model_id, gene_no_match_tuples, name = 'MM')
    removed_ids.append(removed)
    super_model_id, total, removed = process_reactions(super_model_id, no_gene_tuples, name = 'MM', process_count=gnm)
    removed_ids.append(removed)
    removed_reactions = fba_client.get_reactions({'reactions' : removed_ids})
    with open('.mmlog.txt', 'a') as log:
        log.write('\n\n Removed ' + str(total) + ' Reactions: ' + str(gnm) + ' Gene-no-match, ' + str(total - gnm) + ' No-Gene')
        for rxn in removed_reactions:
            log.write('id: ' + str(rxn['id']) + ' name: ' + str(rxn['name']) + ' Equation: ' + str(rxn['equation']) )
    print 'output model...'
    ws_client.copy_objects({'from' : {'objid' : super_model_id, 'wsid' : ws_id}, 'to' : {'objid' : super_model_id, 'wsid' : args['modelws']}})
    print 'further analysis...'
    removed_ids = list()
    i=0
    super_model_id, i, removed = process_reactions(args['model'], gene_no_match_tuples, name = 'Aonly')
    removed_ids.append(removed)
    super_model_id, i, removed = process_reactions(args['model'], no_gene_tuples, name = 'Aonly', process_count=i)
    removed_ids.append(removed)
finally:
    finish(save_ws=True)

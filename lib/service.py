from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import LocalFbaModelServices
import random

ws_url = 'https://kbase.us/services/ws'
fba_url = 'https://kbase.us/services/KBaseFBAModeling/'


def _init_clients():
    ws_c = Workspace(ws_url)
    # The KBase web-client instantiated below is being deprecated and no longer works for all functions.
    # A local work-around has been added using kb_sdk and fba_tools
    fba_c = LocalFbaModelServices
    return ws_c, fba_c


ws_client, fba_client = _init_clients()


# =====================================================================================================================
# Type Strings From KBase

def types():
    return {'FBAModel': 'KBaseFBA.FBAModel',
            'Biochemistry': 'KBaseBiochem.Biochemistry',
            'Genome': 'KBaseGenomes.Genome',
            'FBA': 'KBaseFBA.FBA',
            'ReactionProbabilities': 'ProbabilisticAnnotation.RxnProbs',
            'ProteomeComparison': 'GenomeComparison.ProteomeComparison',
            'Media': 'KBaseBiochem.Media'
            }


# TODO: More pragmatic approach, have these stored in some sort of XML or plain text file that can be read in
# =====================================================================================================================


def get_object(objid, wsid, name=None):
    """
    Returns an object and it's associated KBase information

    Returns an ObjectData (dictionary) like what is returned in the workspace service 'get_objects' function:

	/* The data and supplemental info for an object.

		UnspecifiedObject data - the object's data or subset data.
		object_info info - information about the object.
		list<ProvenanceAction> provenance - the object's provenance.
		username creator - the user that first saved the object to the
			workspace.
		timestamp created - the date the object was first saved to the
			workspace.
		list<obj_ref> - the references contained within the object.
		obj_ref copied - the reference of the source object if this object is
			a copy and the copy source exists and is accessible.
			null otherwise.
		boolean copy_source_inaccessible - true if the object was copied from
			another object, but that object is no longer accessible to the
			user. False otherwise.
		mapping<id_type, list<extracted_id>> extracted_ids - any ids extracted
			from the object.
		string handle_error - if an error occurs while setting ACLs on
			embedded handle IDs, it will be reported here.
		string handle_stacktrace - the stacktrace for handle_error.

	*/
	typedef structure {
		UnspecifiedObject data;
		object_info info;
		list<ProvenanceAction> provenance;
		username creator;
		timestamp created;
		list<obj_ref> refs;
		obj_ref copied;
		boolean copy_source_inaccessible;
		mapping<id_type, list<extracted_id>> extracted_ids;
		string handle_error;
		string handle_stacktrace;
	} ObjectData;

    :param name: (optional) the name for the object to be retrieved. if included, favored over ID
    :param wsid: the workspace to retrieve the object from
    :param objid: the id of the object to be retrieved

    """
    result = None
    if name is None:
        result = ws_client.get_objects([{'objid': objid, 'wsid': wsid}])[0]
    else:
        result = ws_client.get_objects([{'name': name, 'wsid': wsid}])[0]
    return result['data'], result['info']


def get_info(objid, wsid, name=None):
    if name is None:
        return ws_client.get_object_info_new({'objects': [{'objid': objid, 'wsid': wsid}]})[0]
    else:
        return ws_client.get_object_info_new({'objects': [{'name': name, 'wsid': wsid}]})[0]


def save_object(data, type, wsid, objid=None, name=None):
    """
    Saves an object in KBase

    :param data: data representing the object to be saved
    :param type: a string representing the KBase type of the object
    :param wsid: destination workspace
    :param objid: (optional) ID for location of object to be saved (use with care, overwriting/failures are at KBase's
        discretion).
    :param name: (optional) string name for the pbject to be saved
    :return: a list of information about the object as it is stored in KBase
    """
    sv = {u'data': data, u'type': type, u'name': name}
    if objid is not None:
        sv[u'objid'] = objid
    if name is not None:
        sv[u'name'] = name
    info = ws_client.save_objects({u'id': wsid, u'objects': [sv]})[0]
    return info[0], info[6]


def list_objects(workspace_id, typestr=None):
    """
    returns a list of all the objects within a workspace in tuples (obj_id, ws_id, object_name)

    :rtype: list
    :param typestr: (optional) if set, lists only objects of this type (filter over default case)
    :param workspace_id: the workspace to list the objects from
    :return: a list of tuples of objects
    """
    objects = ws_client.list_objects({'ids': [workspace_id]})
    result = list()
    for obj in objects:
        object_type = obj[2]
        if typestr is None or typestr in object_type or types()[typestr] in object_type:  # type filtering of our list
            result.append((obj[0], obj[6], obj[1], obj[2]))
    return result


def clear_workspace(workspace_id):
    """
    clear all objects in a workspace (except for a Narrative object if applicable)
    :param workspace_id: workspace to clear
    :return: None
    """
    object_ids = [{'objid': info[0], 'wsid': workspace_id} for info in
                  ws_client.list_objects({'ids': [workspace_id]}) if not info[2].startswith('KBaseNarrative')]
    if len(object_ids) > 0:
        ws_client.delete_objects(object_ids)


def delete_objects(object_tuples):
    """
    delete objects
    :param object_tuples: list of tuples representing objects to delete of the form (obj_id, ws_id)
    :return: None
    """
    object_ids = [{'objid': info[0], 'wsid': info[1]} for info in object_tuples]
    if len(object_ids) > 0:
        ws_client.delete_objects(object_ids)


def copy_object(from_tuple, to_tuple):
    """
    Copies an object in the service to another location in the service

    :param from_tuple: (objid, wsid) of the object to be copied
    :param to_tuple: (name, wsid) of the destination. workspace may differ. NOTE NAME IS A STRING
    :return: a tuple with information on the new objectmodel
    """
    info = ws_client.copy_object({'from': {'wsid': from_tuple[1],
                                           'objid': from_tuple[0]},
                                  'to': {'wsid': to_tuple[1], 'name': to_tuple[0]}})
    return info[0], info[6]


def gapfill_model(model, media, workspace=None, rxn_probs=None):
    """

    :param model: FBAModel to gapfill
    :param media: Media to gapfill the model to
    :param workspace: destination workspace for new model and gapfill object
    :param rxn_probs: (optional) ReactionProbabilities for probabilisitc gapfilling
    :param name: (optional) name for new model. KBase will overwrite original if left unspecified.
    :return: the information for a new gap-filled model
    """
    gap_formulation = {u'formulation': media.fba_formulation()}
    if workspace is None:
        workspace = model.workspace_id
    if rxn_probs is not None:
        gap_formulation.update({u'probabilisticAnnotation': rxn_probs.object_id,

                                u'probabilisticAnnotation_workspace': rxn_probs.workspace_id})

    params = {u'fbamodel_id': str(model.object_id), u'fbamodel_workspace': str(model.workspace_id), u'fbamodel_output_id': str(model.name),
              u'workspace': workspace,
              u'media_id': media.object_id, u'media_workspace': media.workspace_id,
              u'comprehensive_gapfill': False}
    fba_client.gapfill_model(params)
    return model.object_id, model.workspace_id


def _gapfill_solution(fba):
        """
        If this FBA was created as a gapfilling solution, then this returns a list of reactions to be added/adjusted
        :return: list(tuple) (rxn_id, direction, etc.)
        """
        # For now, naively assume first = best = only gap-filling solution
        solutions = fba['gapfillingSolutions']
        if len(solutions) < 1:
            raise ValueError("This is not a gapfilling solution")
        gsol = solutions[0]['gapfillingSolutionReactions']
        result = []
        for r in gsol:
            reaction_id = r['reaction_ref'].split('/')[-1] + '_' + \
                          r['compartment_ref'].split('/')[-1] + str(r['compartmentIndex'])
            direction = r['direction']
            result.append((reaction_id, direction))
        return result


def fba_formulation(media):
    return {u'media_id': str(media.object_id), u'media_workspace': str(media.workspace_id)}


def runfba(model, media, workspace=None):
    """
    runs Flux Balance Analysis on an FBAModel in the fba modeling service

    :param model: FBAModel to run flux balance analysis on
    :param media: Media to run FBA with
    :param workspace: (optional) workspace for the FBA object to be left in, default is model workspace
    :return: tuple identity of the FBA stored in the service
    """
    if workspace is None:
        workspace = model.workspace_id
    fba_params = fba_formulation(media).update({u'workspace': workspace,
                                                u'fbamodel_id': model.object_id,
                                                u'fbamodel_workspace': model.workspace_id})
    info = fba_client.runfba(fba_params)
    ws, name = info[3].split('/')
    info = get_info(None, ws, name)
    return info[0], info[6]


def runfva(model, media, workspace=None):
    """
    runs Flux Balance Analysis on an FBAModel in the fba modeling service

    :param model: FBAModel to run flux balance analysis on
    :param media: Media to run FBA with
    :param workspace: (optional) workspace for the FBA object to be left in, default is model workspace
    :return: tuple identity of the FBA stored in the service
    """
    if workspace is None:
        workspace = model.workspace_id
    fba_params = fba_formulation(media).update({u'workspace': workspace,
                                                u'fbamodel_id': model.object_id,
                                                u'fbamodel_workspace': model.workspace_id,
                                                u'fva': True})
    info = fba_client.runfba(fba_params)
    ws, name = info[3].split('/')
    info = get_info(None, ws, name)
    return info[0], info[6]


def translate_model(src_model, protcomp, workspace=None, translation_name=None):
    """
    Uses the service to translate an FBAModel to a close genome relative
    :param protcomp: ProteomeComparison with source and target Genome
    :param src_model: FBAModel of source
    return: tuple identity of the translated model stored in the service
    """
    if translation_name is None:
        translation_name = src_model.name + '_translation'
    if workspace is None:
        workspace = src_model.workspace_id
    trans_params = {u'keep_nogene_rxn': 1,
                    u'proteincomparison_id': protcomp.object_id,
                    u'proteincomparison_workspace': protcomp.workspace_id,
                    u'fbamodel_id': src_model.object_id,
                    u'fbamodel_workspace': src_model.workspace_id,
                    u'fbamodel_output_id': translation_name,
                    u'workspace': workspace}
    info = fba_client.translate_model(trans_params)
    ws, name = info['new_fbamodel_ref'].split('/')
    info = get_info(None, ws, name)
    return info[0], info[6]


def reconstruct_genome(genome, workspace=None):
    """
    Reconstructs a genome and returns the identity of a stored draft recon model (FBAModel)
    :param workspace: (optional) destination workspace. Default is genome.workspace_id
    :param genome: Genome to draft a reconstruction for
    :return: tuple identity of the draft model stored in the service (FBAModel)
    """
    if workspace is None:
        workspace = genome.workspace_id
    recon_params = {u'genome_id': genome.object_id, u'genome_workspace': genome.workspace_id, u'workspace': workspace}
    info = fba_client.build_metabolic_model(recon_params)
    ws, name = info['new_fbamodel_ref'].split('/')
    info = get_info(None, ws, name)
    return info[0], info[6]


def remove_reactions_in_place(model, reactions_to_remove):
    """
    Removes reactions from an FBAModel IN PLACE (changes object as it is stored)

    Recommended to make a copy first

    :param model: FBAModel to remove reactions form
    :param reactions_to_remove: reactions to remove (removal_id's)
    :return:
    """
    removal_args = {'fbamodel_id': model.object_id,
                    'fbamodel_workspace': model.workspace_id,
                    'data': {'reactions_to_remove': reactions_to_remove}}
    fba_client.remove_reactions(removal_args)


def remove_reaction(model, reaction, workspace=None, output_id=None, in_place=False):
    """

    :param workspace: (optional) (str) the workspace to put the new model in, disregarded if in_place=True, default is
        model.workspace_id
    :param model: FBAModel to remove the reaction from
    :param reaction: removal_id (str) of the reaction to remove
    :param output_id: (optional) (str) of the new name for the output model
    :param in_place: (optional) set to true if you want to remove the reaction from the model in place instead of making
        a new model. Will disregard output_id argument if set to true
    :return: info tuple for the new FBAModel in the stored environment
    """

    if in_place:
        remove_reactions_in_place(model, [reaction])
    if output_id is None:
        i = 0
        output_id = model.name + '-' + str(i)
        names = set([info[3] for info in list_objects(model.workspace_id)])
        while output_id in names:
            i += 1
            output_id = model.name + '-' + str(i)

    info = fba_client.remove_reactions({'fbamodel_id': model.object_id,
                                        'fbamodel_workspace':model.workspace_id,
                                        'fbamodel_output_id': output_id,
                                        'workspace': model.workspace_id,
                                        'data': {'reactions_to_remove': [reaction]}})
    ws, name = info['new_fbamodel_ref'].split('/')
    info = get_info(None, ws, name)
    return info[0], info[6]


def add_reactions(model, new_reactions, workspace=None, name=None):
    """
    adds reactions to an FBAModel, in place or with a copy (set name to a new name)
    :param model: FBAModel to add reactions to
    :param new_reactions: list of tuples of the form (rxn_id, rxn_comp, direction, gpr) (gpr is optional)
    :param workspace: (optional) destination workspace, default is model.workspace_id
    :param name: output name for the new model. use to make a new one or modify in place
    :return: tuple identity of the model stored in the service (FBAModel)
    """
    if workspace is None:
        workspace = model.workspace_id
    args = {'fbamodel_id': model.object_id,
            'fbamodel_workspace': model.workspace_id,
            'workspace': workspace,
            'data': {'reactions_to_add': new_reactions}}
    if name is not None:
        args['output_id'] = name
    info = fba_client.add_reactions(args)
    ws, name = info['new_fbamodel_ref'].split('/')
    info = get_info(None, ws, name)
    return info[0], info[6]


def add_reactions_manually(model, reactions, workspace=None, name=None):
    """
    Manually fix special reactions within the the object itself (use with caution)
    :param name: what to name the model when it is saved
    :param workspace: workspace to save the new FBAModel in
    :param reactions: (list<ModelReaction>) list of reactions to add manually
    :param model: FBAModel to add the reactions to
    """
    model.get_object()
    if workspace is None:
        workspace = model.workspace_id
    obj = model.data
    cpds = dict([(c['id'], c) for c in obj['modelcompounds']])
    for r in reactions:
        obj['modelreactions'].append(r.data)
        for cpd in r.data['modelReactionReagents']:
            c = cpd['modelcompound_ref'].split('/')[-1]
            if c not in cpds:
                compound = {'id': c,
                            'name': c,
                            'aliases': [u'mdlid:' + c.split('_')[0]],
                            'charge': 0,
                            'compound_ref': '489/6/6/compounds/id/cpd00000',
                            'modelcompartment_ref': '~/modelcompartments/id/' + c.split('_')[-1],
                            'formula': ''
                            }
                obj['modelcompounds'].append(compound)
                cpds = dict([(c['id'], c) for c in obj['modelcompounds']])
    if name is not None:
        return save_object(obj, types()['FBAModel'], workspace, name=name)
    return save_object(obj, types()['FBAModel'], workspace, objid=model.object_id)


def adjust_gprs(model, adjustments):
    adjust_args = {'fbamodel_id': model.object_id,
                   'fbamodel_workspace': model.workspace_id,
                   'workspace': model.workspace_id,
                   'data': {'reactions_to_modify': [(r[0], r[0].split('_')[-1][0], r[1], r[2]) for r in adjustments]},
                   }
    fba_client.edit_reactions(adjust_args)


def adjust_directions(model, adjustments):
    """
    adjusts directions for reactions in an FBAModel
    :param model: FBAModel to adjust directions for
    :param adjustments: list<tuple> (rxn_id, direction). if rxn_id is not already in the model, it may be added
    :return: None
    """
    adjust_args = {'fbamodel_id': model.object_id,
                    'fbamodel_workspace': model.workspace_id,
                   'workspace': model.workspace_id,
                   'data': {'reactions_to_modify': [(r[0], r[0].split('_')[-1][0], r[1]) for r in adjustments]}
                    }
    fba_client.edit_reactions(adjust_args)


def _integrate_gapfill(model, solution_fba, workspace=None):
    changes = _gapfill_solution(solution_fba)
    reactions = dict([(r.rxn_id(), r) for r in model.get_reactions()])
    dirs = []
    additions = []
    for r in changes:
        if r[0] in reactions:
            dirs.append((reactions[r[0]].get_removal_id(), r[1]))
        else:
            temp = r[0].split('_')
            rxn_id = temp[0]
            rxn_comp = temp[1]
            additions.append((rxn_id, rxn_comp, r[1]))
    adjust_directions(model, dirs)
    info = add_reactions(model, additions, workspace=workspace)
    return info


def model_info(model):
    comp = fba_client.compare_models({'models': [model.object_id], 'workspaces': [model.workspace_id]})
    return (comp['model_comparisons'], dict([(r['reaction'], r) for r in comp['reaction_comparisons']]))

def init_workspace(ws=None, name=None):
    ws_id = ws
    ws_name = name
    if ws_name is None:
        ws_name = 'MMws'
    if ws is None:
        ws_conflict = True
        while ws_conflict:
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
                ws_name += str(random.randint(1, 9))
    return ws_id, ws_name
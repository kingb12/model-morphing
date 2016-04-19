from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import sys
import objects


def _init_clients():
    # Get workspace service URL parameter
    with open("./urls/.kbase_workspaceURL", "r") as myfile:
        url = myfile.read().replace('\n', '')
    ws_c = Workspace(url)
    # Get FBA Model Services URL parameter
    with open("./urls/.kbase_fbaModelServicesURL", "r") as myfile:
        url = myfile.read().replace('\n', '')
    fba_c = fbaModelServices(url)
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
        return ws_client.get_objects_info_new({'objects': [{'name': name, 'wsid': wsid}]})[0]


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
    :param type: (optional) the type of object to return (a filter of default case)
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


def gapfill_model(model, media, workspace=None, rxn_probs=None, name=None):
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
    if name is None:
        name = model.name
    params = {u'model': str(model.object_id), u'model_workspace': str(model.workspace_id), u'out_model': name,
              u'workspace': workspace, u'formulation': gap_formulation,
              u'gapFill': name + u'-gf'}
    for key in params:
        print str(key) + u': ' + str(params[key])
    print u'\n\n\n'
    model_info = fba_client.gapfill_model(params)
    return model_info[0], model_info[6]


def fba_formulation(media):
    return {u'media': str(media.object_id), u'media_workspace': str(media.workspace_id)}


def runfba(model, media, workspace=None):
    """
    Runs Flux Balance Analysis on an FBAModel in the fba modeling service

    :param model: FBAModel to run flux balance analysis on
    :param media: Media to run FBA with
    :param workspace: (optional) workspace for the FBA object to be left in, default is model workspace
    :return: tuple identity of the FBA stored in the service
    """
    if workspace is None:
        workspace = model.workspace_id
    fba_params = {u'workspace': workspace, u'model': model.object_id, u'model_workspace': model.workspace_id,
                  u'formulation': fba_formulation(media)}
    info = fba_client.runfba(fba_params)
    return info[0], info[6]


def translate_model(src_model, protcomp, workspace=None):
    """
    Uses the service to translate an FBAModel to a close genome relative
    :param protcomp: ProteomeComparison with source and target Genome
    :param src_model: FBAModel of source
    return: tuple identity of the translated model stored in the service
    """
    if workspace is None:
        workspace = src_model.workspace_id
    trans_params = {u'keep_nogene_rxn': 1,
                    u'protcomp': protcomp.object_id, u'protcomp_workspace': protcomp.workspace_id,
                    u'model': src_model.object_id, u'model_workspace': src_model.workspace_id,
                    u'workspace': workspace}
    info = fba_client.translate_fbamodel(trans_params)
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
    recon_params = {u'genome': genome.object_id, u'genome_workspace': genome.workspace_id, u'workspace': workspace}
    info = fba_client.genome_to_fbamodel(recon_params)
    return info[0], info[6]


def remove_reactions_in_place(model, reactions_to_remove):
    """
    Removes reactions from an FBAModel IN PLACE (changes object as it is stored)

    Recommended to make a copy first

    :param model: FBAModel to remove reactions form
    :param reactions_to_remove: reactions to remove (removal_id's)
    :return:
    """
    removal_args = {'model': model.object_id,
                    'workspace': model.workspace_id,
                    'reaction': reactions_to_remove,
                    'removeReaction': True}
    fba_client.adjust_model_reaction(removal_args)


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
    args = {'model': model.object_id,
            'model_workspace': model.workspace_id,
            'workspace': workspace,
            'reactions': new_reactions}
    if name is not None:
        args['output_id'] = name
    info = fba_client.add_reactions(args)
    return info[0], info[6]

def add_reactions_manually(model, specials, workspace=None, name=None):
    """
    Manually fix special reactions within the the object itself
    """
    model.get_object()
    if workspace is None:
        workspace = model.workspace_id
    obj = model.data
    cpds = dict([(c['id'], c) for c in obj['modelcompounds']])
    for r in specials:
        obj['modelreactions'].append(r.data)
        for cpd in r.data['modelReactionReagents']:
            c = cpd['modelcompound_ref'].split('/')[-1]
            if c not in cpds:
                compound = {'id': c,
                            'name': c,
                            'aliases': [u'mdlid:' + c.split('_')[0]],
                            'charge': 0,
                            'compound_ref': '489/6/6/compounds/id/cpd00000',
                            'modelcompartment_ref': '~/modelcompartments/id/' + c.split('_')[-1]
                            }
                obj['modelcompounds'].append(compound)
                cpds = dict([(c['id'], c) for c in obj['modelcompounds']])
    if name is not None:
        return save_object(obj, types()['FBAModel'], workspace, name=name)
    return save_object(obj, types()['FBAModel'], workspace, objid=model.object_id)

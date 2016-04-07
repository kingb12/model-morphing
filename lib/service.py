from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import sys


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
    if name is None:
        return ws_client.get_objects([{'objid': objid, 'wsid': wsid}])[0]
    else:
        return ws_client.get_objects([{'name': name, 'wsid': wsid}])[0]


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
    return ws_client.save_objects({u'id': wsid, u'objects': [sv]})[0]


def list_objects(workspace_id, typestr=None):
    """
    returns a list of all the objects within a workspace in tuples (obj_id, ws_id, object_name)

    :param type: (optional) the type of object to return (a filter of default case)
    :param workspace_id: the workspace to list the objects from
    :return: a list of tuples of objects
    """
    objects = ws_client.list_objects({'ids': [workspace_id]})
    for obj in objects:
        object_type = obj[2]
        result = list()
        if typestr is None or typestr in object_type:  # type filtering of our list
            result.append((obj[0], obj[6], obj[1]))
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
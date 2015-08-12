
from biokbase.workspace.client import Workspace 
from biokbase.workspace.client import ServerError 
from biokbase.GenomeComparison.Client import GenomeComparison
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.userandjobstate.client import UserAndJobState
import sys
ws = Workspace()
ws_next = Workspace(url='https://next.kbase.us/services/ws/')
old_id = sys.argv[1]
new_id = sys.argv[2]
objlist = ws.list_objects({'ids' : [old_id]})
object_ids = list()
for obj in objlist:
	object_ids.append({'wsid' : old_id, 'objid' : obj[0]})
objects = ws.get_objects(object_ids)
save_data = list()
for obj in objects:
	if obj['info'][0] != 4:
		save_data.append({'data' : obj['data'], 'type': obj['info'][2]})
ws_next.save_objects({'id' : new_id, 'objects' : save_data})

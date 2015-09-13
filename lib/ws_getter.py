from biokbase.workspace.client import Workspace
ws_client = Workspace()
ws_next_client = Workspace(url='https://next.kbase.us/services/ws')
a, b = ws_next_client.get_objects([{'objid' : '4', 'wsid' : '68'}, {'objid' : '5', 'wsid' : '68'}])[0:2]
a_params = {'type' : a['info'][2], 'data': a['data']}
b_params = {'type' : b['info'][2], 'data': b['data']}
ws_client.save_objects({'id': '9145', 'objects': [a_params, b_params]})


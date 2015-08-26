from biokbase.workspace.client import Workspace 
params = dict()
# params['workspace'] = 'eureka'
ws = Workspace(url='https://kbase.us/services/ws')
# b = ws.create_workspace(params)
a = ws.list_workspace_info(params)
for i in a:
	print i


from biokbase.workspace.client import Workspace 
params = dict()
# params['workspace'] = 'eureka'
ws = Workspace()
# b = ws.create_workspace(params)
a = ws.list_workspace_info(params)
for i in a:
	print i[0:2]


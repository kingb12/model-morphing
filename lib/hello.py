from biokbase.workspace.client import Workspace 
params = dict()
# params['workspace'] = 'eureka'
ws = Workspace()
# ws.create_workspace(params)
a = ws.list_workspace_info(params)
for i in range(len(a)):
	print a[i][1]

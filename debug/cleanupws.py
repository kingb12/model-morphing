
from biokbase.workspace.client import Workspace 
ws = Workspace()
for i in range(9160,9167):
	params = {'id' : i}
	ws.delete_workspace(params) 


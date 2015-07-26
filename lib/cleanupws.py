
from biokbase.workspace.client import Workspace 
ws = Workspace()
for i in range(9107,9111):
	params = {'id' : i}
	ws.delete_workspace(params) 


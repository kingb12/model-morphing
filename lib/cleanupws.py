
from biokbase.workspace.client import Workspace 
ws = Workspace()
for i in range(9056,9057):
	params = {'id' : i}
	ws.delete_workspace(params) 



from biokbase.workspace.client import Workspace 
ws = Workspace()

with open('ws.txt' , 'r') as f:
    for line in f:
        int_list = [int(i) for i in line.split()]
        print int_list
for i in int_list:
	ws.delete_workspace({'id' : str(i)})


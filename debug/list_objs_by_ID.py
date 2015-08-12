from biokbase.workspace.client import Workspace 
import sys
ws_client = Workspace()
ws_idsa = list()
ws_idsa.append(sys.argv[1])
obj_list = ws_client.list_objects({'ids' : ws_idsa}) 
for obj in obj_list:
	print obj

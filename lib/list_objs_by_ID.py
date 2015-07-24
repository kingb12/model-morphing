from biokbase.workspace.client import Workspace 
import sys
ws_client = Workspace()
ws_idsa = sys.argv[1:len(sys.argv)]
obj_list = ws_client.list_objects({'ids' : ws_idsa}) 
for obj in obj_list:
	print obj[0:4]

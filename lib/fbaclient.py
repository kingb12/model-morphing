# A class for wrapping the functionality of the KBase Clients, primarily FBA and
# WS
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError

class FBAClient():
    def __init__(self, url=None):
        if url is None:
            with open ("./urls/.kbase_workspaceURL", "r") as myfile:
                url = myfile.read().replace('\n','')
            with open ("./urls/.kbase_fbaModelServicesURL", "r") as myfile:
                url = myfile.read().replace('\n','')


from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import Model
import Reaction
import random
import Helpers
from Morph import AbstractGrowthCondition
import Client
from operator import itemgetter


class SimpleCondition(AbstractGrowthCondition):
    '''
    a growth conditon for absolute growth (objective > 0)

    Required attributes of args:
        - morph
        - model
        - fba_name
    '''
    def evaluate(self, args):
        morph = args['morph']
        ws = morph.ws_id
        fba_formulation = {'media': morph.media, 'media_workspace': morph.mediaws}
        fba_params = {'fba': args['fba_name'], 'workspace': ws, 'model' : args['model'], 'model_workspace':ws, 'formulation': fba_formulation}
        fbaMeta = Client.fba_client.runfba(fba_params)
        flux = fbaMeta[-1]['Objective']
        print "Flux is " + str(flux)
        return flux > 0.0
class BarkeriCondition(AbstractGrowthCondition):
    '''
    a growth condition for barkeri (3 media)
    '''
    def evaluate(self, args):
        morph = args['morph']
        ws = morph.ws_id
        fba_formulation = {'media': morph.media, 'media_workspace': morph.mediaws}
        fba_params = {'fba': args['fba_name'], 'workspace': ws, 'model' : args['model'], 'model_workspace':ws, 'formulation': fba_formulation}
        fbaMeta = Client.fba_client.runfba(fba_params)
        flux = fbaMeta[-1]['Objective']
        print "Flux is " + str(flux)
        raise NotImplementedError()
class AllMedia(AbstractGrowthCondition):
    def __init__(self, media):
        self.media = media
    def evaluate(self, args):
        morph = args['morph']
        ws = morph.ws_id
        for med in self.media:
            fba_formulation = {'media': med[0], 'media_workspace': med[1]}
            fba_params = {'fba': args['fba_name'], 'workspace': ws, 'model' : args['model'], 'model_workspace':ws, 'formulation': fba_formulation}
            fbaMeta = Client.fba_client.runfba(fba_params)
            flux = fbaMeta[-1]['Objective']
            print "Flux is " + str(flux)
            if flux <= 0.0:
                return False
        return True




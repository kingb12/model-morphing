import service

from lib import objects


class AbstractGrowthCondition:
    """
    an interface for processing reactions according to some condition

    Subclases must implement evaluate(args) and return true or false. The primary use
    for this class is in the process_reactions method of the Client module, which removes
    reactions iteratiely, and decides to keep or remove a reaction based on the outcome of
    a GrowthCondition.  Here is an example:
        class SimpleCondition(AbstractGrowthCondition):
            def evaluate(args):
                # args must have attribute model
                fba = runfba(args['model'])
                return fba['Objective'] > 0
    This SimpleCondition keeps all reactions that are absolutely necessary for the models growth
    """

    def __init__(self):
        self.fba = None

    def evaluate(self, arguments):
        raise NotImplementedError()

class SimpleCondition(AbstractGrowthCondition):
    """
    a growth conditon for absolute growth (objective > 0)

    Required attributes of args:
        - morph
        - model
        - fba_name
    """

    def evaluate(self, arguments):
        morph = arguments['morph']
        model = arguments['model'] if 'model' in arguments else morph.model
        info = service.runfba(model, morph.media, workspace=morph.ws_id)
        self.fba = objects.FBA(info[0], info[1])
        return self.fba.objective > 0.0


class BarkeriCondition(AbstractGrowthCondition):
    """
    a growth condition for barkeri (3 media)
    """
    def evaluate(self, arguments):
        raise NotImplementedError()
        morph = arguments['morph']
        model = arguments['model'] if 'model' in arguments else morph.model
        info = service.runfba(model, morph.media, workspace=morph.ws_id)
        fba = objects.FBA(info[0], info[1])


class AllMedia(AbstractGrowthCondition):

    def __init__(self, media):
        AbstractGrowthCondition.__init__(self)
        self.media = media

    def evaluate(self, args):
        morph = args['morph']
        ws = morph.ws_id
        for med in self.media:
            morph = args['morph']
            model = args['model'] if 'model' in args else morph.model
            info = service.runfba(model, med, workspace=morph.ws_id)
            fba = objects.FBA(info[0], info[1])
            self.fba = fba
            if not fba.objective > 0.0:
                return False
        return True





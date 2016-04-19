import copy


class Log:
    """
    A simple class for logging what occurs to a mutatable object. Object is responsible for composing a log, and making
    the required edits to it.

    A log is a simple ADT, it is a sequence of Actions that occurred to an object. The first action is 'initialize' and
    stores the initial state of the object as encoded by the object. Default behavior for this might be to just perform
    a deepcopy on the object and store a reference to it, etc. This is at the discretion of the composing object.
    """
    # TODO: A LOG SHOULD BE ABLE TO WRITE ITSELF IN MARKDOWN
    def __init__(self, object):
        """

        :param object: the object to be logged. Remember, a log does not update itself
        :return: Log with initial action
        """
        try:
            a = object.logInitialize()
        except AttributeError:
            a = copy.deepcopy(object)
        self.actions = []
        self.actions.append(Action('initialize', [a], context=str(Log.__init__)))

    def add(self, action_type, inputs, outputs, context=None, notes=None):
        a = Action(action_type, {'in': inputs, 'out': outputs}, context=context, notes=notes)
        self.actions.append(a)


class Action:
    """
    a class representing an action to be documented in a log. How it is used is up to the discretion of the object

    an action with type attribute 'initialize' is special, it's members field is the singleton list of whatever it was
    encoded as in log.__init__
    """
    def __init__(self, action_type, members, context=None, notes=None):
        """

        :param action_type: String describing the broad class/type of action. e.g. log uses 'initialize' for __init__
        :param members: the involved parties, in a dict 'in', 'out' as keys
        :param context: (optional) at object discretion. context in which action was performed
        :param notes: (optional) at object discretion. You can put ANYTHING here
        :return: Action as initialized with parameters
        """
        self.type = action_type
        self.members = members
        self.context = context
        self.notes = notes

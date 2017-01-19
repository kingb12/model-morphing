import copy
import Plotter
import objects
import json

#String constants
TYPE_STR = 'type'
MEMBERS = 'members'
CONTEXT = 'context'
NOTES = 'notes'
class Log:
    """
    A simple class for logging what occurs to a mutatable object. Object is responsible for composing a log, and making
    the required edits to it.

    A log is a simple ADT, it is a sequence of Actions that occurred to an object. The first action is 'initialize' and
    stores the initial state of the object as encoded by the object. Default behavior for this might be to just perform
    a deepcopy on the object and store a reference to it, etc. This is at the discretion of the composing object.
    """
    def __init__(self, object=None):
        """

        :param object: the object to be logged. Remember, a log does not update itself
        :return: Log with initial action
        """
        if object is not None:
            try:
                a = object.logInitialize()
            except AttributeError:
                a = copy.deepcopy(object)
            self.actions = []
            self.actions.append(Action('initialize', [a], context=str(Log.__init__)))
        else:
            self.actions = []

    def add(self, action_type, inputs, outputs, context=None, notes=None):
        self._make_references(inputs)
        self._make_references(outputs)
        a = Action(action_type, {'in': inputs, 'out': outputs}, context=context, notes=notes)
        self.actions.append(a)

    def markdown(self):
        header = ('Action', 'Members', 'Context', 'Notes')
        actions = Plotter.SimpleTable(header)
        for a in self.actions:
            actions.add((a.md_tuple()))
        return actions.markdown()

    def _make_references(self, item):
        # Mutates list in place replacing StoredObjects with their references
        if type(item) == list:
            for i in range(len(item)):
                if isinstance(item[i], objects.StoredObject):
                    item[i] = item[i].reference()
        else:
            if isinstance(item, objects.StoredObject):
                    item = item.reference()

    def to_json(self):
        j_actions = [a.to_json() for a in self.actions]
        return json.dumps(j_actions)

    @staticmethod
    def from_json(json_str):
        actions = [Action.from_json(a) for a in json.loads(json_str)]
        log = Log()
        log.actions = [a for a in actions if a is not None]

        return log




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

    def __str__(self):
        result = 'Action: ' + str(self.type) + '\nMembers: \n'
        if type(self.members) == dict:
            for i in self.members:
                result += str(i) + ': \n' + str(self.members[i])
        else:
            result += str(self.members)
        result += '\nContext: ' + str(self.context)
        return result

    def md_tuple(self):
        return str(self.type), str(self.members), str(self.context), str(self.notes)

    def to_json(self):
        optional = {}
        if self.context is not None:
            optional[CONTEXT] = self.context
        if self.notes is not None:
            optional[NOTES] = self.notes
        return json.dumps({TYPE_STR: self.type, MEMBERS: self.members}.update(optional))

    @staticmethod
    def from_json(json_str):
        data = json.loads(json_str)
        if data is None:
            return None
        return Action(data[TYPE_STR], data[MEMBERS], context=data[CONTEXT] if CONTEXT in data else None, notes=data[NOTES] if NOTES in data else None)
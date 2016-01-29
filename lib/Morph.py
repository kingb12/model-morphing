# import necesary services
from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError
from biokbase.fbaModelServices.Client import fbaModelServices
import random
import json
import argparse
import time
import traceback
from operator import itemgetter

def _init_clients():
    global ws_client, args
    global fba_client
    global ws_id
    # Get workspace service URL parameter
    with open ("./urls/.kbase_workspaceURL", "r") as myfile:
        url = myfile.read().replace('\n','')
    ws_client = Workspace(url)
    # Get FBA Model Services URL parameter
    with open ("./urls/.kbase_fbaModelServicesURL", "r") as myfile:
        url = myfile.read().replace('\n','')
    fba_client = fbaModelServices(url)
    return ws_client, fba_client

def _init_workspace(ws = None):
    ws_id = ws
    ws_name = 'MMws'
    if (ws is None):
        ws_conflict = True
        while (ws_conflict):
            create_ws_params = {'workspace': ws_name, 'globalread': 'r', 'description':
                                "A workspace for storing the FBA's and meta data of the algorithm"}
            # Try to create a workspace, catch an error if the name is already in use
            try:
                new_ws = ws_client.create_workspace(create_ws_params)
                # new_ws is type workspace_info, a tuple where 0, 1 are id, name
                ws_id = new_ws[0]
                ws_name = new_ws[1]
                ws_conflict = False
            except ServerError:
                 ws_name += str(random.randint(1,9))
    return ws_id, ws_name

class Morph():

    # These are the allowed properties of a morph object. VAlues not specified
    # by user are set to None
    # This is an Abstract Object representing the Morph of a metabolic model
    # from one species to a close relative genome. It has information related to
    # the source model, the target genome, the reactions in the model as it is
    # morphed from source import to target, and the set of
    properties = set(['src_model', 'src_modelws', 'genome', 'genomews', 'probanno', 'probannows', 'protcomp', 'protcompws',
                      'model', 'rxn_labels', 'objects', 'info', 'essential_ids', 'removed_ids', 'ws_id', 'ws_name', 'trans_model', 'recon_model', 'media', 'mediaws', 'probhash'])


    def __init__(self, *arg_hash, **kwargs):
          for dictionary in arg_hash:
              for key in dictionary:
                  if (key in Morph.properties):
                      setattr(self, key, dictionary[key])
          for key in kwargs:
              if (key in Morph.properties):
                  setattr(self, key, kwargs[key])
          for prop in Morph.properties:
              if (not hasattr(self, prop)):
                  setattr(self, prop, None)
          if (self.ws_id is None):
              self.ws_id, self.ws_name = _init_workspace()

    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True)
    # Overridden Functions to produce unique output
    def __str__(self):
        output = ''
        for key in vars(self):
            attr = getattr(self, key)
            if (isinstance(attr, dict)):
                attr = attr.keys()
                if len(attr) < 100:
                    output += str(key) + ': ' + str(attr) + '\n'
                else:
                    output += str(key) + ': ' + str(attr[0:100]) + ' ... (more)\n'
            else:
                output += str(key) + ': ' + str(attr) + '\n'
        return output
    def __repr__(self):
        return str(self)
    def __unicode__(self):
        return unicode(str(self))

class Gpr():
    def __init__(self, reaction=None):
        '''
        creates a GPR object t represent the gene-protein-reaction relationship for a rxn/model
        '''
        if reaction is not None:
            self.gpr = self._gpr_set(reaction)
            self.ftrs = self._feature_set()
        else:
            self.gpr = None
        self._check_rep()

    def __str__(self):
        '''
        returns the gpr represented as a string, both Human and KBase readable
        '''
        gpr_set = self.gpr
        proteins = list(gpr_set)
        proteins_str = ""
        for i in range(0, len(proteins)):
            units_str = ""
            subunits = list(proteins[i])
            for j in range(0, len(subunits)):
                unit = list(subunits[j])
                features = ""
                for k in range(0, len(unit)):
                    feature = unit[k]
                    if (k > 0):
                        feature = " or " + feature
                    features += feature
                unit_str = "(" + features + ")"
                if j > 0:
                    unit_str = " and " + unit_str
                units_str += unit_str
            protein = "(" + units_str + ")"
            if (i > 0):
                protein = " or " + protein
            proteins_str += protein
        gpr = "(" + proteins_str + ")"
        return gpr
    def __repr__(self):
        return str(self.gpr)
    def __unicode__(self):
        return unicode(str(self))
    def _gpr_set(self, rxn_object):
        reaction = rxn_object
        rxn_proteins = reaction['modelReactionProteins']
        prots = set()
        for i in range(0, len(rxn_proteins)):
            prot = set()
            subunits = rxn_proteins[i]['modelReactionProteinSubunits']
            for j in range(0, len(subunits)):
                unit = subunits[j]
                ftrs = [f.split('/')[-1] for f in unit['feature_refs']]
                if len(ftrs) > 0:
                    prot.add(frozenset(ftrs))
            if len(prot) > 0:
                prots.add(frozenset(prot))
        if len(prots) > 0:
            return frozenset(prots)
        return None

    def _feature_set(self):
        features = set()
        for protein in self.gpr:
            for sub in protein:
                for f in sub:
                    features.add(f)
        return frozenset(features)
    def features(self):
        '''
        returns a set of the features in this gpr
        '''
        # ok to return attibute, it's a frozen set
        self._check_rep()
        return self.ftrs

    def contains_feature(self, feature):
        '''
        returns true if the feature is somewhere in the gpr (features of the string form 'kb|g.587.peg.1234')
        '''
        return feature in self.ftrs
    def contains_protein(self, protein):
        '''
        returns true if the protein is in the gpr (proteins are frozen sets of subunits, which in turn are frozen sets of features)
        '''
        return protein in self.gpr
    def contains_subunit(self, subunit):
        '''
        returns true if the subunit is in a protein in the gpr (subunits are frozensets of features)
        '''
        for protein in self.gpr:
            if subunit in protein:
                return True
        return False

    def merge(self, other_gpr):
        # if at least one is None, attempt to return a possibly non-None one
        gpr_set1 = self.gpr
        gpr_set2 = other_gpr.gpr
        if gpr_set1 is None or gpr_set2 is None:
            if gpr_set1 is None:
                return gpr_set2
            return gpr_set1
        g1 = gpr_set1
        g2 = set(gpr_set2)
        # enclosing set is the set of proteins
        for protein in g1:
            if protein not in g2:
                matchedProtein = False
                proteinsToRemove = set()
                proteinsToAdd = set()
                # Look for a nearly matching protein
                for g2_protein in g2:
                    #if they share a subunit or any feature (catches homolog and subunit cases)
                    if  len(self._unnest_sets(protein) & (self._unnest_sets(g2_protein))) != 0:
                        proteinsToRemove.add(g2_protein)
                        prot = set(g2_protein)
                        matchedProtein = True
                        # Look for a matching subunit
                        matchedSub = False
                        for subunit in protein:
                            if subunit not in g2_protein:
                                for other in g2_protein:
                                    if len(subunit & other) != 0:
                                        matchedSub = True
                                        prot.remove(other)
                                        prot.add(frozenset(subunit.union(other)))
                                        if frozenset(prot) not in proteinsToAdd:
                                            proteinsToAdd.add(frozenset(prot))
                        if not matchedSub:
                            proteinsToRemove.remove(g2_protein)
                            # do nothing, but better other solutions should be
                            # implememted
                                # recon  - do nothing
                                # trans - always push
                                # stronger/weaker
                assert (len(proteinsToRemove) > 0 or not matchedProtein or (matchedProtein and not matchedSub))
                g2 = g2 - set(proteinsToRemove)
                g2 |= set(proteinsToAdd)
                #Simple Case, proteins don't conflict
                if not matchedProtein or not matchedSub:
                    g2.add(protein)
        return_gpr = self.new_gpr(frozenset(g2))
        return_gpr._check_rep()
        return return_gpr

    def new_gpr(self, gpr_set):
        '''
        returns a new gpr with the given gpr_set
        '''
        newgpr = Gpr()
        newgpr.gpr = gpr_set
        newgpr.ftrs = newgpr._feature_set()
        newgpr._check_rep()
        return newgpr

    def _check_rep(self):
        '''
        check rep invariant
        '''
        if self.gpr is not None and self.ftrs is None:
            raise RepresentationError(self)
    def _unnest_sets(self, nested_set):
        '''
        Takes anything in a nested set form and returns a single set of it's non-set elements

        i.e. {{{a}}{{b, c}}{{d}{e}}} -->  {a, b, c, d, e,}
        '''
        single_set = set()
        for item in nested_set:
            if type(item) is set or type(item) is frozenset:
                single_set |= self._unnest_sets(item)
            else:
                single_set.add(item)
        return single_set



class RepresentationError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return str(type(self.value)) + repr(self.value.gpr) + '\n' + repr(self.value.ftrs)




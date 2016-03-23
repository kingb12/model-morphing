import service
import copy


class StoredObject:
    """
    A class representing any object stored in our environment (D.o.E. KBase), having an object_id and a workspace_id
    """

    # The current environment is KBase, where all calls to the API refer to their objects by ids and workspace_ids
    def __init__(self, object_id, workspace_id):
        self.object_id = object_id
        self.workspace_id = workspace_id
        self.initialized = True

    def __setattr__(self, key, value):
        if self.initialized:
            raise NotImplementedError
        self.__dict__[key] = value

    def get_object(self):
        """
        Returns the representation of our object from it's stored environment
        :return dict: the data representing this object

        Use this at your own peril. You can't overwrite anything without editing it in the saved environment but
        the interior elements of this dictionary are NOT the same as the classes representing them. If you want data
        from the interior of the object, it is best to use the abstractions provided by it's more specific type.
        """
        return service.get_object(self.object_id, self.workspace_id)['data']
        # TODO: Implement clone or copy


class Biochemistry(StoredObject):
    """
    a class representing a Biochemistry Object
    """
    pass


class FBAModel(StoredObject):
    """
    A class representing an FBA Model in the stored environment

    Parent Classes: StoredObject -> FBAModel
    """

    DEFAULT_BIOCHEM = Biochemistry(6, 489)

    def get_reactions(self):
        """
        Returns a list of ModelReaction objects representing this model's reactions
        """
        model_obj = self.get_object()
        return [ModelReaction(r) for r in model_obj['modelreactions']]  # Potentially Wasteful but maintains D.R.Y.


class ModelReaction:
    """
    a class representing a reaction in an FBAModel. Has a GPR, Compounds, etc.
    """

    def __init__(self, reaction_obj):
        """

        :param reaction_obj: the dictionary representing the reaction as it is in the FBAModel
        :return:
        """
        self.gpr = Gpr(reaction_obj)
        self.data = reaction_obj

    # TODO add features from Reactions Module

    def rxn_id(self):
        """
        returns the reaction id of a model reaction object (i.e. rxn34565_c0)
        """
        rxn_id = str(self.data['reaction_ref'].split('/')[-1]) + '_' + self.data['modelcompartment_ref'].split('/')[-1]
        if rxn_id.startswith('rxn00000'):
            return self.data['id']
        return rxn_id

    def get_equation(self):
        """

        :return: equation of the reaction
        """
        return [Compound(cpd) for cpd in self.data['modelReactionReagents']]

    def get_rxn_ref(self):
        """
        gets the reaction object reference if it is associated with a KBase Biochemistry. 'rxn00000' if not in a biochem
        """
        return self.data['reaction_ref'].split('/')[-1]

    def get_removal_id(self):
        """
        returns the ID which can be used by the service to remove the reaction from a model
        :return: the ID which can be used to address and remove the reaction in the model by the service
        """
        return self.data['id']

    def get_biochem_ref(self):
        """
        returns a tuple with the reference to a
        """
        ref = self.data['reaction_ref'].split('/')
        return Biochemistry(ref[1], ref[0])


class Compound(object):
    """
    small wrapper class for compound coefficient and ID.
    """

    def __init__(self, compound_obj, biochem=FBAModel.DEFAULT_BIOCHEM):
        """
        initializes a Compound Object

        :param compound_obj:
        :param biochem: (optional, default is KBase Default Biochem) the biochemistry in which a compound can be
            referenced. If it is not found in the biochemistry
        :return:
        """
        self.coeff = compound_obj['coefficient']
        self.compound_id = compound_obj['modelcompound_ref'].split('/')[-1]
        self.biochem = biochem

    def get_info(self):
        """
        returns the information about this compound frmo its referenced biochemistry

        :return: (type=dict) the information on the biochemistry for this compound, with the following keys:
            [u'cues',
            u'name',
            u'pkas',
            u'deltaGErr',
            u'pkbs',
            u'abbreviation',
            u'mass',
            u'isCofactor',
            u'deltaG',
            u'formula',
            u'id',
            u'defaultCharge',
            u'unchargedFormula']
        """
        biochem = self.biochem.get_object()
        for cpd in biochem['compounds']:
            if cpd['id'] == self.compound_id:
                return cpd
        raise BiochemistryError("Information not found in provided biochemistry")

    def __str__(self):
        """
        returns a str readable format of compound, e.g: -3*(cpd01111)
        :return:
        """
        return str(self.coeff) + '*(' + str(self.compound_id) + ')'

    def formula(self):
        """
        returns the formula for the compound

        :return: the formula for the compound
        """
        return self.get_info()['formula']


class Gpr:
    """
    a class representing the Gene -> Protein -> Reaction relationship for a ModelReaction in a model
    """
    def __init__(self, reaction=None):
        """
        creates a GPR object that represent the gene-protein-reaction relationship for a ModelReaction
        """
        if reaction is not None:
            self.gpr, self.gpr_type = self._gpr_set(reaction)
            self.ftrs = self._feature_set()
            if self.gpr_type is None:
                if len(self.ftrs) > 0:
                    self.gpr_type = 'genes'
                else:
                    self.gpr_type = 'no-gene'
        else:
            self.gpr = None
        self._check_rep()

    def __str__(self):
        """
        returns the gpr represented as a string, both Human and service KBase readable
        """
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
                    if k > 0:
                        feature = " or " + feature
                    features += feature
                unit_str = "(" + features + ")"
                if j > 0:
                    unit_str = " and " + unit_str
                units_str += unit_str
            protein = "(" + units_str + ")"
            if i > 0:
                protein = " or " + protein
            proteins_str += protein
        gpr = "(" + proteins_str + ")"
        return gpr

    def __repr__(self):
        return repr(self.gpr)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.gpr == other.gpr
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __unicode__(self):
        return unicode(str(self))

    @staticmethod
    def _gpr_set(rxn_object):
        """
        creates the gpr set for the self.gpr field. Intended to only be called once

        A gpr with no features is a frozenset(frozenset(frozenset()))
        """
        reaction = rxn_object
        rxn_proteins = reaction['modelReactionProteins']
        prots = set()
        for i in range(0, len(rxn_proteins)):
            prot = set()
            if rxn_proteins[i]['note'] == u'spontaneous' or rxn_proteins[i]['note'] == u'universal':
                return frozenset([frozenset([frozenset([])])]), rxn_proteins[i]['note']
            subunits = rxn_proteins[i]['modelReactionProteinSubunits']
            for j in range(0, len(subunits)):
                unit = subunits[j]
                ftrs = [f.split('/')[-1] for f in unit['feature_refs']]
                if len(ftrs) > 0:
                    prot.add(frozenset(ftrs))
            if len(prot) > 0:
                prots.add(frozenset(prot))
        if len(prots) > 0:
            return frozenset(prots), None
        return frozenset([frozenset([frozenset([])])]), None

    def __iter__(self):
        """
        returns an interator over the gpr attribute
        """
        if self.gpr is not None:
            return self.gpr.__iter__()

    def _feature_set(self):
        features = set()
        for protein in self.gpr:
            for sub in protein:
                for f in sub:
                    features.add(f)
        return frozenset(features)

    def features(self):
        """
        returns a set of the features in this gpr
        """
        # ok to return attribute, it's a frozen set
        self._check_rep()
        return self.ftrs

    def contains_feature(self, feature):
        """
        returns true if the feature is somewhere in the gpr (features of the string form 'kb|g.587.peg.1234')
        :param feature: the feature in question. e.g: 'kb|g.587.peg.123'
        """
        return feature in self.ftrs

    def contains_protein(self, protein):
        """
        returns true if the protein is in the gpr (proteins are frozen sets of subunits, which in turn are frozen sets
        of features)
        :param protein: the protein in question e.g. frozenset(frozenset('kb|g.587.peg.1234'))
        """
        return protein in self.gpr

    def contains_subunit(self, subunit):
        """
        returns true if the subunit is in a protein in the gpr (subunits are frozen sets of features)
        :param subunit:
        """
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
        examine_gpr = False
        # enclosing set is the set of proteins
        for protein in g1:
            if protein not in g2:
                matched_protein = False
                proteins_to_remove = set()
                proteins_to_add = set()
                # Look for a nearly matching protein
                for g2_protein in g2:
                    # if they share a subunit or any feature (catches homolog and subunit cases)
                    # AND they are equal in number of subunits
                    if len(protein) == len(g2_protein) and len(self._unnest_sets(protein) &
                                                               self._unnest_sets(g2_protein)) != 0:
                        proteins_to_remove.add(g2_protein)
                        prot = set(g2_protein)
                        matched_protein = True
                        # Look for a matching subunit
                        matched_sub = False
                        for subunit in protein:
                            if subunit not in g2_protein:
                                for other in g2_protein:
                                    if len(subunit & other) != 0:
                                        matched_sub = True
                                        prot.remove(other)
                                        new_sub = subunit.union(other)
                                        if new_sub in prot:
                                            examine_gpr = True
                                        else:
                                            prot.add(frozenset(subunit.union(other)))
                                        if frozenset(prot) not in proteins_to_add:
                                            proteins_to_add.add(frozenset(prot))
                        if not matched_sub:
                            proteins_to_remove.remove(g2_protein)
                            # do nothing, but better other solutions should be
                            # implemented
                            # recon  - do nothing
                            # trans - always push
                            # stronger/weaker

                assert (len(proteins_to_remove) > 0 or not matched_protein or (matched_protein and not matched_sub))
                g2 = g2 - set(proteins_to_remove)
                g2 |= set(proteins_to_add)
                # Simple Case, proteins don't conflict
                if not matched_protein or not matched_sub:
                    g2.add(protein)
        return_gpr = self.new_gpr(frozenset(g2))
        return_gpr.remove_redundancy()
        if examine_gpr:
            return_gpr.gpr_type = u'potential merge conflict'
        else:
            return_gpr.gpr_type = u'merge'
        return_gpr.parents = (self, other_gpr)
        return_gpr._check_rep()
        return return_gpr

    def remove_redundancy(self):
        """
        finds proteins that are subsets of each other and removes the smaller one

        e.g. ((a or b)) or ((a or b or c or d)) ==> ((a or b or c or d))
        performing this check prevents redundancy and helps ensure symmetry in T.merge(R) == R.merge(T)
        """
        matches = dict()
        for protein in self.gpr:
            for protein2 in self.gpr:
                if protein != protein2 and self._unnest_sets(protein).issubset(self._unnest_sets(protein2)):
                    matches[protein] = protein2
        for protein in matches:
            # we know that all features in prot are in protein2
            protein2 = copy.deepcopy(matches[protein])
            matched_subs = set()
            for sub in protein:
                # sub should be a subset of a sub in p2
                for sub2 in protein2:
                    if sub.issubset(sub2):
                        matched_subs.add(sub2)
                        protein2 = protein2 - sub2
            if len(matched_subs) == len(matches[protein]):
                assert matched_subs == matches[protein]
                i = len(self.gpr)
                self.gpr = self.gpr - frozenset([protein])
                assert len(self.gpr) == i - 1

    def is_empty(self):
        self._check_rep()
        return len(self.ftrs) == 0

    @staticmethod
    def new_gpr(gpr_set):
        """
        returns a new gpr with the given gpr_set
        :param gpr_set: the gpr_set to make a gpr from
        """
        newgpr = Gpr()
        newgpr.gpr = gpr_set
        newgpr.ftrs = newgpr._feature_set()
        newgpr._check_rep()
        return newgpr

    def _check_rep(self):
        # check rep invariant

        if self.gpr is not None and self.ftrs is None:
            raise RepresentationError(self)
        if hasattr(self, 'ftrs'):
            feature_set = self._unnest_sets(self.gpr)
            if len(self.ftrs) == 0 and self.gpr != frozenset([frozenset([frozenset([])])]):
                raise RepresentationError(self)
            for f in self.ftrs:
                if f not in feature_set:
                    raise RepresentationError(self)
            for f in feature_set:
                if not self.contains_feature(f):
                    raise RepresentationError(self)
        # verify structure
        if hasattr(self, 'gpr') and self.gpr is not None:
            assert type(self.gpr) is frozenset
            for protein in self.gpr:
                assert type(protein) is frozenset
                for sub in protein:
                    assert type(sub) is frozenset
                    for ftr in sub:
                        assert type(ftr) is str or type(ftr) is unicode

    def _unnest_sets(self, nested_set):
        """
        Takes anything in a nested set form and returns a single set of it's non-set elements

        i.e. {{{a}}{{b, c}}{{d}{e}}} -->  {a, b, c, d, e,}
        """
        single_set = set()
        for item in nested_set:
            if type(item) is set or type(item) is frozenset:
                single_set |= self._unnest_sets(item)
            else:
                single_set.add(item)
        return single_set


class Genome(StoredObject):
    """
    a class representing a genome in the stored environment
    """
    pass


class Media(StoredObject):
    """
    a class representing a media in the stored environment
    """
    pass


class FBA(StoredObject):
    """
    a class representing an FBA result in the stored environment
    """
    pass


class ProteomeComparison(StoredObject):
    """
    a class representing a Porteome Comparison in the stored environment
    """
    pass


class ReactionProbabilities(StoredObject):
    """
    a class representing a ReactionProbabilities in the stored environment
    """
    pass


class BiochemistryError(Exception):
    """
    an error for when a look-up in a Biochemistry Object fails
    """
    pass


class RepresentationError(Exception):
    # an error type for internal use to label when a rep invariant has been broken

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return str(type(self.value)) + repr(self.value.gpr) + '\n' + repr(self.value.ftrs)

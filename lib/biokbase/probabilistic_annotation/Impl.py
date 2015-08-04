#BEGIN_HEADER
import sys
import tempfile
import os
import traceback
import time
import re
from biokbase.probabilistic_annotation.DataParser import DataParser, NotReadyError
from biokbase.probabilistic_annotation.Helpers import timestamp, make_object_identity, make_job_directory, ProbAnnoType, RxnProbsType, ServiceVersion
from biokbase.workspace.client import Workspace
from biokbase.fbaModelServices.Client import *
from biokbase.cdmi.client import CDMI_EntityAPI
from biokbase.userandjobstate.client import UserAndJobState
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase import log

# Exception thrown when static database file is missing from Shock.
class MissingFileError(Exception):
    pass

# Exception thrown when role not found in roleToTotalProb dictionary
class RoleNotFoundEror(Exception):
    pass

# Exception thrown when object version is not valid
class WrongVersionError(Exception):
    pass
#END_HEADER


class ProbabilisticAnnotation:
    '''
    Module Name:
    ProbabilisticAnnotation

    Module Description:
    The purpose of the Probabilistic Annotation service is to provide users with
alternative annotations for genes, each attached to a likelihood score, and to
translate these likelihood scores into likelihood scores for the existence of
reactions in metabolic models.  With the Probabilistic Annotation service:

- Users can quickly assess the quality of an annotation.

- Reaction likelihood computations allow users to estimate the quality of
  metabolic networks generated using the automated reconstruction tools in
  other services.

- Combining reaction likelihoods with gapfilling both directly incorporates
  available genetic evidence into the gapfilling process and provides putative
  gene annotations automatically, reducing the effort needed to search for
  evidence for gapfilled reactions.
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    def _checkInputArguments(self, ctx, input, requiredArgs, defaultArgDict):
        ''' Check that required input arguments are present and set defaults for non-required arguments.
        
            If a key in defaultArgDict is not found in the input it is added with the
            specified default value.

            @param ctx Current context object
            @param input Dictionary keyed by argument name of argument values
            @param requiredArgs List of required arguments in the input dictionary
            @param defaultArgDict Dictionary keyed by argument name of default argument values
            @return Dictionary keyed by argument name of argument values (including default arguments)
            @raise ValueError when required argument is not found
        '''
        if requiredArgs is not None:
            for arg in requiredArgs:
                if arg not in input:
                    message = "Required argument %s not found" %(arg)
                    ctx.log_err(message)
                    raise ValueError(message)
        if defaultArgDict is not None:
            for arg in defaultArgDict:
                if arg in input:
                    continue
                else:
                    input[arg] = defaultArgDict[arg]

        return input

    def _rolesetProbabilitiesToRoleProbabilities(self, ctx, input, genome, queryToTuplist, workFolder):
        ''' Compute probability of each role from the rolesets for each query protein.

            At the moment the strategy is to take any set of rolestrings containing
            the same roles and add their probabilities.  So if we have hits to both
            a bifunctional enzyme with R1 and R2, and hits to a monofunctional enzyme
            with only R1, R1 ends up with a greater probability than R2.

            I had tried to normalize to the previous sum but I need to be more careful
            than that (I'll put it on my TODO list) because if you have e.g. one hit
            to R1R2 and one hit to R3 then the probability of R1 and R2 will be unfairly
            brought down due to the normalization scheme.

            @param ctx: Current context object
            @param input: Dictionary of input parameters to calculate() function
            @param genome: Genome ID string
            @param queryToTuplist: Dictionary keyed by query gene of list of tuples with roleset and likelihood
            @param workFolder: Path to directory in which to store temporary files
            @return List of tuples with query gene, role, and likelihood
        '''
    
        ctx.log_debug('Started computing role probabilities from roleset probabilities for '+genome)

        # Start with an empty list.
        roleProbs = list()

        # Iterate over all of the query genes in the dictionary.
        # querygene -> [ (roleset1, likelihood_1), (roleset2, likelihood_2), ...]
        for query in queryToTuplist:
            # This section actually does the conversion of likelihoods.
            # See equation 3 in the paper ("Calculating reaction likelihoods" section).
            queryRolesToProbs = dict()
            for tup in queryToTuplist[query]:
                rolelist = tup[0].split(self.config["separator"])
                # Add up all the instances of each particular role on the list.
                for role in rolelist:
                    if role in queryRolesToProbs:
                        queryRolesToProbs[role] += tup[1]
                    else:
                        queryRolesToProbs[role] = tup[1]
    
            # Add them to the array.
            for role in queryRolesToProbs:
                roleProbs.append( (query, role, queryRolesToProbs[role]) )
    
        # Save the generated data when debug is turned on.
        if ctx.get_log_level() >= log.DEBUG2:
            role_probability_file = os.path.join(workFolder, "%s.roleprobs" %(genome))
            fid = open(role_probability_file, "w")
            for tuple in roleProbs:
                fid.write("%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2]))
            fid.close()

        ctx.log_debug('Finished computing role probabilities from roleset probabilities for '+genome)
            
        return roleProbs
    
    def _totalRoleProbabilities(self, ctx, input, genome, roleProbs, workFolder):
        ''' Given the likelihood that each gene has each role, estimate the likelihood
            that the entire ORGANISM has that role.

            To avoid exploding the likelihoods with noise, I just take the maximum
            likelihood of any query gene having a function and use that as the
            likelihood that the function exists in the cell.

            A gene is assigned to a role if it is within DILUTION_PERCENT of the maximum
            probability. DILUTION_PERCENT can be adjusted in the config file. For each
            role the maximum likelihood and the estimated set of genes that perform that
            role are linked with an OR relationship to form a Boolean Gene-Function
            relationship.

            @param ctx Current context object
            @param input Dictionary of input parameters to calculate() function
            @param genome Genome ID string
            @param roleProbs List of tuples with query gene, role, and likelihood
            @param workFolder Path to directory in which to store temporary files
            @return List of tuples with role, likelihood, and estimated set of genes that perform the role
            @raise RoleNotFoundError when role is not placed properly in roleToTotalProb dictionary
        '''
    
        ctx.log_debug('Started generating whole-cell role probability file for '+genome)
    
        # Find maximum likelihood among all query genes for each role.
        # This is assumed to be the likelihood of that role occurring in the organism as a whole.
        roleToTotalProb = dict()
        for tuple in roleProbs:
            if tuple[1] in roleToTotalProb:
                if float(tuple[2]) > roleToTotalProb[tuple[1]]:
                    roleToTotalProb[tuple[1]] = float(tuple[2])
            else:
                roleToTotalProb[tuple[1]] = float(tuple[2])
    
        # Get the genes within DILUTION_PERCENT percent of the maximum
        # likelihood and assert that these are the most likely genes responsible for that role.
        # (note - DILUTION_PERCENT is defined in the config file)
        # This produces a dictionary from role to a list of genes
        # See equation 4 in the paper ("Calculating reaction likelihoods" section).
        roleToGeneList = dict()
        for tuple in roleProbs:
            if tuple[1] not in roleToTotalProb:
                message = "Role %s not placed properly in roleToTotalProb dictionary?" %(tuple[1])
                ctx.log_err(message)
                raise RoleNotFoundError(message)
            if float(tuple[2]) >= float(self.config["dilution_percent"])/100.0 * roleToTotalProb[tuple[1]]:
                if tuple[1] in roleToGeneList:
                    roleToGeneList[tuple[1]].append(tuple[0])
                else:
                    roleToGeneList[tuple[1]] = [ tuple[0] ]
        
        # Build the array of total role probabilities.     
        totalRoleProbs = list()
        for role in roleToTotalProb:
            gpr = " or ".join(list(set(roleToGeneList[role])))
            # We only need to group these if there is more than one of them (avoids extra parenthesis when computing complexes)
            if len(list(set(roleToGeneList[role]))) > 1:
                gpr = "(" + gpr + ")"
            totalRoleProbs.append( (role, roleToTotalProb[role], gpr ) )   
    
        # Save the generated data when debug is turned on.
        if ctx.get_log_level() >= log.DEBUG2:
            total_role_probability_file = os.path.join(workFolder, "%s.cellroleprob" %(genome))
            fid = open(total_role_probability_file, "w")
            for tuple in totalRoleProbs:
                fid.write("%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2]))
            fid.close()
        
        ctx.log_debug('Finished generating whole-cell role probability file for '+genome)
            
        return totalRoleProbs
    
    def _complexProbabilities(self, ctx, input, genome, totalRoleProbs, workFolder, complexesToRequiredRoles = None):
        ''' Compute the likelihood of each protein complex from the likelihood of each role.

            A protein complex represents a set functional roles that must all be present
            for a complex to exist.  The likelihood of the existence of a complex is
            computed as the minimum likelihood of the roles within that complex (ignoring
            roles not represented in the subsystems).

            For each protein complex, the likelihood, type, list of roles not in the
            organism, and list of roles not in subsystems is returned.  The type is a
            string with one of the following values:
        
            CPLX_FULL - All roles found in organism and utilized in the complex
            CPLX_PARTIAL - Only some roles found in organism and only those roles that
                were found were utilized. Note this does not distinguish between not
                there and not represented for roles that were not found
            CPLX_NOTTHERE - Likelihood is 0 because the genes aren't there for any of
                the subunits
            CPLX_NOREPS - Likelihood is 0 because there are no representative genes in
                the subsystems for any of the subunits
            CPLX_NOREPS_AND_NOTTHERE - Likelihood is 0 because some genes aren't there
                for any of the subunits and some genes have no representatives

            @param ctx: Current context object
            @param input: Dictionary of input parameters to calculate() function
            @param genome: Genome ID string
            @param totalRoleProbs: List of tuples with role, likelihood, and estimated set
                of genes that perform the role
            @param workFolder: Path to directory in which to store temporary files
            @param complexesToRequiredRoles: Dictionary keyed by complex ID to the roles
                involved in forming that complex. If it is None we read it from the CDMI
                files we downloaded, otherwise we use the provided dictionary. This is
                included for template model support in the future
            @return List of tuples with complex ID, likelihood, type, list of roles not in
                organism, list of roles not in subsystems, and boolean Gene-Protein
                relationship
        '''
    
        ctx.log_debug('Started computing complex probabilities for '+genome)
    
        # Get the mapping from complexes to roles if it isn't already provided.
        if complexesToRequiredRoles is None:
            complexesToRequiredRoles = self.dataParser.readComplexRoles()
        
        # Get the subsystem roles (used to distinguish between NOTTHERE and NOREPS).
        otu_fidsToRoles, otu_rolesToFids = self.dataParser.readFilteredOtuRoles()
        allroles = set()
        for fid in otu_fidsToRoles:
            for role in otu_fidsToRoles[fid]:
                allroles.add(role)
    
        # Build two dictionaries, both keyed by role, one mapping the role to its
        # likelihood and one mapping to the gene list.
        rolesToProbabilities = dict()
        rolesToGeneList = dict()
        for tuple in totalRoleProbs:
            rolesToProbabilities[tuple[0]] = float(tuple[1]) # can skip the float()?
            rolesToGeneList[tuple[0]] = tuple[2]
    
        # Iterate over complexes and compute complex probabilities from role probabilities.
        # Separate out cases where no genes seem to exist in the organism for the reaction
        # from cases where there is a database deficiency.
        # See equation 5 in the paper ("Calculating reaction likelihoods" section).
        SEPARATOR = self.config["separator"]
        complexProbs = list()
        for cplx in complexesToRequiredRoles:
            allCplxRoles = complexesToRequiredRoles[cplx]
            availRoles = list() # Roles that may have representatives in the query organism
            unavailRoles = list() # Roles that have representatives but that are not apparently in the query organism
            noexistRoles = list() # Roles with no representatives in the subsystems
            for role in complexesToRequiredRoles[cplx]:
                if role not in allroles:
                    noexistRoles.append(role)
                elif role not in rolesToProbabilities:
                    unavailRoles.append(role)
                else:
                    availRoles.append(role)
            TYPE = ""
            GPR = ""
            if len(noexistRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOREPS"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            if len(unavailRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOTTHERE"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            # Some had no representatives and the rest were not found in the cell
            if len(unavailRoles) + len(noexistRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOREPS_AND_NOTTHERE"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            # Otherwise at least one of them is available
            if len(availRoles) == len(allCplxRoles):
                TYPE = "CPLX_FULL"
            elif len(availRoles) < len(allCplxRoles):
                TYPE = "CPLX_PARTIAL_%d_of_%d" %(len(availRoles), len(allCplxRoles))

            # Link individual functions in complex with an AND relationship to form a
            # Boolean Gene-Protein relationship.
#            partialGprList = [ "(" + s + ")" for s in [ rolesToGeneList[f] for f in availRoles ] ]
            partialGprList = [ rolesToGeneList[f] for f in availRoles ]
            GPR = " and ".join( list(set(partialGprList)) )

            if GPR != "" and len(list(set(partialGprList))) > 1:
                GPR = "(" + GPR + ")"

            # Find the minimum probability of the different available roles (ignoring ones
            # that are apparently missing) and call that the complex likelihood.
            minp = 1000
            for role in availRoles:
                if rolesToProbabilities[role] < minp:
                    minp = rolesToProbabilities[role]
            complexProbs.append( (cplx, minp, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )

        # Save the generated data when debug is turned on.
        if ctx.get_log_level() >= log.DEBUG2:
            complex_probability_file = os.path.join(workFolder, "%s.complexprob" %(genome))
            fid = open(complex_probability_file, "w")
            for tuple in complexProbs:
                fid.write("%s\t%1.4f\t%s\t%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4], tuple[5]))
            fid.close()
        
        ctx.log_debug('Finished computing complex probabilities for '+genome)
        return complexProbs
    
    def _reactionProbabilities(self, ctx, input, genome, complexProbs, workFolder, rxnsToComplexes = None):
        ''' Estimate the likelihood of reactions from the likelihood of complexes.

            The reaction likelihood is computed as the maximum likelihood of complexes
            that perform that reaction.

            If the reaction has no complexes it won't even be in this file because of the way
            I set up the call... I could probably change this so that I get a list of ALL reactions
            and make it easier to catch issues with reaction --> complex links in the database.
            Some of the infrastructure is already there (with the TYPE).

            @param ctx: Current context object
            @param input: Dictionary of input parameters to calculate() function
            @param genome: Genome ID string
            @param complexProbs: List of tuples with complex ID, likelihood, type, list of
                roles not in organism, list of roles not in subsystems, and boolean
                Gene-Protein relationship
            @param workFolder: Path to directory in which to store temporary files
            @param rxnsToComplexes: Dictionary keyed by reaction ID to a list of catalyzing
                complexes. If it is None we read it from the CDMI files we downloaded,
                otherwise we use the provided dictionary. This is included for template
                model support in the future
            @return List of tuples with reaction ID, likelihood, reaction type, complex info,
                and gene-protein-reaction relationship
        '''
    
        ctx.log_debug('Started computing reaction probabilities for '+genome)
        
        # Build a dictionary keyed by complex ID of tuples with likelihood, type, and GPR.
        # Note we don't need to use the list of roles not in organism and list of roles
        # not in subsystems.
        # cplx --> {likelihood, type, GPR}
        cplxToTuple = dict()
        for tuple in complexProbs:
            cplxToTuple[tuple[0]] = ( tuple[1], tuple[2], tuple[5] )
        
        # Get the mapping from reactions to complexes if it isn't already provided.
        if rxnsToComplexes is None:
            rxnsToComplexes = self.dataParser.readReactionComplex()

        # Take the MAXIMUM likelihood of complexes catalyzing a particular reaction
        # and call that the reaction likelihood.
        # See equation 6 in the paper ("Calculating reaction likelihoods" section).
        reactionProbs = list()
        for rxn in rxnsToComplexes:
            TYPE = "NOCOMPLEXES"
            rxnComplexes = rxnsToComplexes[rxn]
            maxProb = 0
            GPR = ""
            complexList = list()
            for cplx in rxnComplexes:
                if cplx in cplxToTuple:
                    # Complex1 (P1; TYPE1) ///Complex2 (P2; TYPE2) ...
                    complexList.append( [ cplx, cplxToTuple[cplx][0], cplxToTuple[cplx][1] ])
                    TYPE = 'HASCOMPLEXES'
            complexString = ''
            if len(complexList) > 0:
                complexList.sort(key=lambda tup: tup[1], reverse=True)
                maxProb = complexList[0][1]
                for complex in complexList:
                    complexString += '%s (%1.4f; %s)%s' %(complex[0], complex[1], complex[2], self.config['separator'])
                complexString = complexString[:-len(self.config['separator'])] # Remove the final separator

            # Iterate separately to get a GPR. We want to apply a cutoff here too to avoid
            # a complex with 80% probability being linked by OR to another with a 5%
            # probability.  For now I've implemented using the same cutoff as we used for
            # which genes go with a role.
            cplxGprs = []
            for cplx in rxnComplexes:
                if cplx in cplxToTuple:
                    if cplxToTuple[cplx][0] < maxProb * float(self.config["dilution_percent"])/100.0:
                        continue
                    cplxGprs.append(cplxToTuple[cplx][2])
            if len(cplxGprs) > 0:
                GPR = " or ".join( list(set(cplxGprs)) )

            # Use a list so that we can modify the reaction IDs if needed to translate to ModelSEED IDs
            reactionProbs.append( [rxn, maxProb, TYPE, complexString, GPR] )
    
        # Save the generated data when debug is turned on.
        if ctx.get_log_level() >= log.DEBUG2:
            reaction_probability_file = os.path.join(workFolder, "%s.rxnprobs" %(genome))
            fid = open(reaction_probability_file, "w")
            for tuple in reactionProbs:
                fid.write("%s\t%1.4f\t%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4]))
            fid.close()
    
        ctx.log_debug('Finished computing reaction probabilities for '+genome)
        return reactionProbs
    
    def _metaboliteWeights(input, model):
        '''Given a model object, computes an S-matrix.
     
        model: Model object

        This function returns three things:
        - A sparse matrix object (coo_matrix) from scipy with indexes matching the lists
        - A dictionary from metabolite UUIDs to their index in the matrix
        - A dictionary from reaction UUIDs to their index in the matrix
     
        The lists should have the same order as the input model.
     
        If absval = True, returns the absolute value of the S-matrix (absolute value of every
        term in the S matrix) rather than S itself.
         
        TODO - put in biomass equation too.
        '''
     
        print(input)
        metIdToIdx = {}
        rxnIdToIdx = {}
     
        idx = 0
        for compound in model["modelcompounds"]:
            metIdToIdx[compound["uuid"]] = idx
            idx += 1
         
        i = []
        j = []
        data = []
        costs = []
        idx = 0
        numZeroReagents = 0
        for reaction in model["modelreactions"]:
            if "modelReactionReagents" not in reaction:
                numZeroReagents += 1
                continue
            for reagent in reaction["modelReactionReagents"]:
                met_uuid = reagent["modelcompound_uuid"]
                coefficient = reagent["coefficient"]
                met_idx = metIdToIdx[met_uuid]
                i.append(met_idx)
                j.append(idx)
                if input["absval"]:
                    coefficient = abs(coefficient)
                data.append(coefficient)
            costs.append(1.0 - reaction["probability"])
            rxnIdToIdx[reaction["uuid"]] = idx
            idx += 1
     
        S_matrix = sparse.coo_matrix( ( data, ( i, j ) ) )
        
        sys.stderr.write("%d model reactions had zero reagents\n" %(numZeroReagents))
        sys.stderr.write("len i=%d len j=%d len data=%d\n" %(len(i), len(j), len(data)))
        sys.stderr.write("S_matrix rows=%d, cols=%d\n" %(S_matrix.shape[0], S_matrix.shape[1]))
    
        rxncosts = numpy.array(costs)    
        '''
        Given an S matrix (sparse) S, and a vector of reaction 
        probabilities rxnprobs, calls the scipy least-squares solver
        to obtain our best estimate for the metabolite weights
        (gamma)
     
        The reaction weights (probabilities) and metabolite weights
        are related by the equation
        R = |S|^T * gamma
     
        where R is the vector of reaction weights, |S| means the absolute value
        of S and gamma is the vector of metabolite weights. Solving in the least-squares
        sense is the best we can do since S is not a square matrix.
     
        Returns a list of metabolite weights with the same indexing as the rows of S.
        '''
         
        S_prime = S_matrix.transpose(copy=True)
        res = linalg.lsqr(S_prime, rxncosts)
        
        for index in range(len(res[0])):
            cpd = model["modelcompounds"][index]
            cpd["weight"] = res[0][index]
            
        return True
    
    def _checkDatabaseFiles(self, ctx):
        ''' Check the status of the static database files.

            @param ctx Current context object
            @return Nothing
            @raise NotReadyError if the database has not been loaded correctly.
        '''
        try:
            status = self.dataParser.readStatusFile()
            if status != "ready":
                message = "Static database files are not ready.  Current status is '%s'." %(status)
                ctx.log_err(message)
                raise NotReadyError(message)
        except IOError:
            message = "Static database files are not ready.  Failed to open status file '%s'." %(self.dataParser.StatusFiles['status_file'])
            ctx.log_err(message)
            raise NotReadyError(message)
        return

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        ''' Constructor for ProbabilisticAnnotation object.

            @param config Dictionary of configuration variables
            @return Nothing
            @raise ValueError when a valid configuration was not provided
        '''
        if config == None:
            # There needs to be a config for the server to work.
            raise ValueError("__init__: A valid configuration was not provided.  Check KB_DEPLOYMENT_CONFIG and KB_SERVICE_NAME environment variables.")
        else:
            self.config = config
        
        submod = os.environ.get('KB_SERVICE_NAME', 'probabilistic_annotation')
        self.mylog = log.log(submod, ip_address=True, authuser=True, module=True, method=True,
            call_id=True, config=os.getenv('KB_DEPLOYMENT_CONFIG'))
        self.mylog.log_message(log.NOTICE, 'Server started, version is '+ServiceVersion)
        configValues = 'shock_url='+self.config['shock_url']
        configValues += ', userandjobstate_url='+self.config['userandjobstate_url']
        configValues += ', workspace_url='+self.config['workspace_url']
        configValues += ', cdmi_url='+self.config['cdmi_url']
        configValues += ', work_folder_path='+self.config['work_folder_path']
        configValues += ', data_folder_path='+self.config['data_folder_path']
        configValues += ', load_data_option='+self.config['load_data_option']
        configValues += ', separator='+self.config['separator']
        configValues += ', dilution_percent='+self.config['dilution_percent']
        configValues += ', pseudo_count='+self.config['pseudo_count']
        configValues += ', job_queue='+self.config['job_queue']
        configValues += ', search_program='+self.config['search_program']
        configValues += ', search_program_path='+self.config['search_program_path']
        configValues += ', blast_threads='+self.config['blast_threads']
        configValues += ', usearch_accel='+self.config['usearch_accel']
        self.mylog.log_message(log.NOTICE, configValues)

        # Create a DataParser object for working with the static database files (the
        # data folder is created if it does not exist).
        self.dataParser = DataParser(self.config)

        # Get the static database files.  If the files do not exist and they are downloaded
        # from Shock, it can take a few minutes before the server is ready.
        testDataPath = os.path.join(os.environ['KB_SERVICE_DIR'], 'testdata')
        self.config['load_data_option'] = self.dataParser.getDatabaseFiles(self.mylog, testDataPath)

        # Validate the value of the job_queue variable.  Currently the only supported value is 'local'.
        # Force it to a valid value to avoid an error trying to submit a job later.
        if self.config['job_queue'] != 'local':
            self.mylog.log_message(log.NOTICE, 'Configuration variable job_queue='+self.config['job_queue']+' switched to local')
            self.config['job_queue'] = 'local'
        #END_CONSTRUCTOR
        pass

    def version(self, ctx):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN version
        ''' Return the name and version number of the service.

            @param ctx Current context object
            @return List with service name string and version number string
        '''

        returnVal = [ os.environ.get('KB_SERVICE_NAME'), ServiceVersion ]
        #END version

        # At some point might do deeper type checking...
        if not isinstance(returnVal, list):
            raise ValueError('Method version return value ' +
                             'returnVal is not type list as required.')
        # return the results
        return [returnVal]

    def annotate(self, ctx, input):
        # ctx is the context object
        # return variables are: jobid
        #BEGIN annotate
        ''' Compute probabilistic annotations from the specified genome object.

            The input dictionary must contain the following keys:
            genome: Name of genome object
            genome_workspace: Workspace from which to grab the Genome object
            probanno: Name of probanno object to output
            probanno_workspace: Workspace to which to save the ProbAnno object

            The following keys are optional:
            verbose: Print lots of messages on the progress of the algorithm

            @param ctx Current context object
            @param input Dictionary with input parameters for function
            @return Job ID of job started to compute annotation likelihoods
        '''

        input = self._checkInputArguments(ctx, input, 
                                          [ "genome", "genome_workspace", "probanno", "probanno_workspace"],
                                          { "verbose" : False }
                                          )
        
        # Make sure the static database files are ready.
        self._checkDatabaseFiles(ctx)
        
        # Set log level to INFO when verbose parameter is enabled.
        if input['verbose']:
            ctx.set_log_level(log.DEBUG)

        # Make sure the Genome object is available.
        wsClient = Workspace(self.config["workspace_url"], token=ctx['token'])
        genomeIdentity = make_object_identity(input['genome_workspace'], input['genome'])
        wsClient.get_object_info( [ genomeIdentity ], 0 )

        # Create a user and job state client and authenticate as the user.
        ujsClient = UserAndJobState(self.config['userandjobstate_url'], token=ctx['token'])

        # Create a job to track running probabilistic annotation.
        description = 'pa-annotate for genome %s to probanno %s for user %s' %(input['genome'], input['probanno'], ctx['user_id'])
        progress = { 'ptype': 'task', 'max': 5 }
        jobid = ujsClient.create_and_start_job(ctx['token'], 'initializing', description, progress, timestamp(3600))
        ctx.log_info('Job '+jobid+' started for genome '+input['genome']+' to probanno '+input['probanno'])

        # Run the job on the local machine.
        if self.config["job_queue"] == "local":
            # Create working directory for job and build file names.
            jobDirectory = make_job_directory(self.config['work_folder_path'], jobid)
            jobDataFilename = os.path.join(jobDirectory, 'jobdata.json')
            outputFilename = os.path.join(jobDirectory, 'stdout.log')
            errorFilename = os.path.join(jobDirectory, 'stderr.log')
    
            # Save data required for running the job.
            jobData = { 'id': jobid, 'input': input, 'context': ctx, 'config': self.config }
            json.dump(jobData, open(jobDataFilename, "w"), indent=4)
    
            # Start worker to run the job.
            jobScript = os.path.join(os.environ['KB_TOP'], 'bin/pa-runjob')
            cmdline = "nohup %s %s >%s 2>%s &" %(jobScript, jobDirectory, outputFilename, errorFilename)
            status = os.system(cmdline)
            ctx.log_info('Job %s is running on local host, status %d' %(jobid, status))

        #END annotate

        # At some point might do deeper type checking...
        if not isinstance(jobid, basestring):
            raise ValueError('Method annotate return value ' +
                             'jobid is not type basestring as required.')
        # return the results
        return [jobid]

    def calculate(self, ctx, input):
        # ctx is the context object
        # return variables are: output
        #BEGIN calculate
        ''' Compute reaction probabilities from a probabilistic annotation.

            The input dictionary must contain the following keys:
            probanno: Name of ProbAnno object to input
            probanno_workspace: Workspace from which to grab the ProbAnno object
            rxnprobs: Name of RxnProbs object
            rxnprobs_workspace: Workspace to which to save the RxnProbs object

            The following keys are optional:
            verbose: Print lots of messages on the progress of the algorithm
            template_model: Name of TemplateModel object
            template_workspace: Workspace from which to grab TemplateModel object

            @param ctx Current context object
            @param input Dictionary with input parameters for function
            @return Object info for RxnProbs object
            @raise WrongVersionError when ProbAnno object version number is invalid
            @raise ValueError when template_workspace input argument is not specified
        '''

        # Sanity check on input arguments
        input = self._checkInputArguments(ctx, input, 
                                          ["probanno", "probanno_workspace", "rxnprobs", "rxnprobs_workspace"], 
                                          { "verbose" : False ,
                                            "template_model" : None,
                                            "template_workspace" : None
                                          }
                                         )

        # Make sure the static database files are ready.
        self._checkDatabaseFiles(ctx)

        # Set log level to INFO when verbose parameter is enabled.
        if input['verbose']:
            ctx.set_log_level(log.DEBUG)
        
        # Create a workspace client.
        wsClient = Workspace(self.config["workspace_url"], token=ctx['token'])
        
        # Get the ProbAnno object from the specified workspace.
        probannoObjectId = make_object_identity(input["probanno_workspace"], input["probanno"])
        objectList = wsClient.get_objects( [ probannoObjectId ] )
        probannoObject = objectList[0]
        if probannoObject['info'][2] != ProbAnnoType:
            message = "ProbAnno object type %s is not %s for object %s" %(probannoObject['info'][2], ProbAnnoType, probannoObject['info'][1])
            ctx.log_err(message)
            raise WrongVersionError(message)
        genome = probannoObject["data"]["genome"]
        
        # Create a temporary directory for storing intermediate files when debug is turned on.
        if ctx.get_log_level() >= log.DEBUG2:
            workFolder = tempfile.mkdtemp("", "calculate-%s-" %(genome), self.config["work_folder_path"])
            ctx.log_debug('Intermediate files saved in '+workFolder)
        else:
            workFolder = None

        # When a template model is specified, use it to build dictionaries for roles,
        # complexes, and reactions instead of retrieving from static database files.
        complexesToRoles = None
        reactionsToComplexes = None
        if input["template_model"] is not None or input["template_workspace"] is not None:
            if not(input["template_model"] is not None and input["template_workspace"] is not None) :
                message = "Template model workspace is required if template model ID is provided"
                ctx.log_err(message)
                raise ValueError(message)

            # Create a dictionary to map a complex to a list of roles and a dictionary
            # to map a reaction to a list of complexes.  The dictionaries are specific to
            # the specified template model instead of covering everything in the central
            # data model.
            complexesToRoles = dict()
            reactionsToComplexes = dict()

            # Get the list of RoleComplexReactions for the template model from the
            # fba modeling service.  The RoleComplexReactions structure has a list
            # of ComplexReactions structures for the given role.  And each ComplexReactions
            # structure has a list of reactions for the given complex.
            fbaClient = fbaModelServices(self.config['fbamodeling_url'], token=ctx['token'])
            roleComplexReactionsList = fbaClient.role_to_reactions( { 'templateModel': input['template_model'], 'workspace': input['template_workspace'] } )

            # Build the two dictionaries from the returned list.
            for rcr in roleComplexReactionsList:
                for complex in rcr['complexes']:
                    complexId = re.sub(r'cpx0*(\d+)', r'kb|cpx.\1', complex['name']) # Convert ModelSEED format to KBase format
                    if complexId in complexesToRoles:
                        complexesToRoles[complexId].append(rcr['name'])
                    else:
                        complexesToRoles[complexId] = [ rcr['name'] ]
                    for reaction in complex['reactions']:
                        reactionId = reaction['reaction']
                        if reactionId in reactionsToComplexes:
                            reactionsToComplexes[reactionId].append(complexId)
                        else:
                            reactionsToComplexes[reactionId] = [ complexId ]

        # Calculate per-gene role probabilities.
        roleProbs = self._rolesetProbabilitiesToRoleProbabilities(ctx, input, genome, probannoObject["data"]["roleset_probabilities"], workFolder)

        # Calculate whole cell role probabilities.
        # Note - eventually workFolder will be replaced with a rolesToReactions call
        totalRoleProbs = self._totalRoleProbabilities(ctx, input, genome, roleProbs, workFolder)

        # Calculate complex probabilities.
        complexProbs = self._complexProbabilities(ctx, input, genome, totalRoleProbs, workFolder, complexesToRequiredRoles = complexesToRoles)

        # Calculate reaction probabilities.
        reactionProbs = self._reactionProbabilities(ctx, input, genome, complexProbs, workFolder, rxnsToComplexes = reactionsToComplexes)

        # If the reaction probabilities were not calculated using the data from the fba modeling service
        # via the template model, we need to convert from the KBase ID format to the ModelSEED format.
        if input["template_model"] is None:
            reactionList = list()
            for index in range(len(reactionProbs)):
                reactionList.append(reactionProbs[index][0])
            EntityAPI = CDMI_EntityAPI(self.config["cdmi_url"])
            numAttempts = 4
            while numAttempts > 0:
                try:
                    numAttempts -= 1
                    reactionData = EntityAPI.get_entity_Reaction( reactionList, [ "source_id" ] )
                    if len(reactionList) == len(reactionData):
                        numAttempts = 0
                except HTTPError as e:
                    pass
            for index in range(len(reactionProbs)):
                rxnId = reactionProbs[index][0]
                reactionProbs[index][0] = reactionData[rxnId]['source_id']
 
        # Create a reaction probability object
        objectData = dict()
        objectData["genome"] = probannoObject["data"]["genome"]
        objectData['genome_workspace'] = probannoObject['data']['genome_workspace']
        if input["template_model"] is None:
            objectData['template_model'] = 'None'
        else:
            objectData["template_model"] = input["template_model"]
        if input["template_workspace"] is None:
            objectData['template_workspace'] = 'None'
        else:
            objectData["template_workspace"] = input["template_workspace"]
        objectData["probanno"] = input['probanno']
        objectData['probanno_workspace'] = input['probanno_workspace']
        objectData["id"] = input["rxnprobs"]
        objectData["reaction_probabilities"] = reactionProbs

        objectMetaData = { "num_reaction_probs": len(objectData["reaction_probabilities"]) }
        objectProvData = dict()
        objectProvData['time'] = timestamp(0)
        objectProvData['service'] = os.environ['KB_SERVICE_NAME']
        objectProvData['service_ver'] = ServiceVersion
        objectProvData['method'] = 'calculate'
        objectProvData['method_params'] = input.items()
        objectProvData['input_ws_objects'] = [ '%s/%s/%d' %(probannoObject['info'][7], probannoObject['info'][1], probannoObject['info'][4]) ]
        objectSaveData = dict();
        objectSaveData['type'] = RxnProbsType
        objectSaveData['name'] = input["rxnprobs"]
        objectSaveData['data'] = objectData
        objectSaveData['meta'] = objectMetaData
        objectSaveData['provenance'] = [ objectProvData ]
        objectInfo = wsClient.save_objects( { 'workspace': input["rxnprobs_workspace"], 'objects': [ objectSaveData ] } )
        output = objectInfo[0]
        
        #END calculate

        # At some point might do deeper type checking...
        if not isinstance(output, list):
            raise ValueError('Method calculate return value ' +
                             'output is not type list as required.')
        # return the results
        return [output]

    def get_rxnprobs(self, ctx, input):
        # ctx is the context object
        # return variables are: output
        #BEGIN get_rxnprobs
        ''' Convert a reaction probability object into a human-readable table.

            @param ctx Current context object
            @param input Dictionary with input parameters for function
            @return List of reaction_probability tuples
            @raise WrongVersionError when RxnProbs object version number is invalid
        '''

        # Sanity check on input arguments
        input = self._checkInputArguments(ctx, input, 
                                          [ "rxnprobs", "rxnprobs_workspace" ], 
                                          { 'rxnprobs_version': None, 'sort_field': 'rxnid' }
                                          )

        wsClient = Workspace(self.config["workspace_url"], token=ctx['token'])
        rxnProbsObjectId = make_object_identity(input["rxnprobs_workspace"], input["rxnprobs"], input['rxnprobs_version'])
        objectList = wsClient.get_objects( [ rxnProbsObjectId ] )
        rxnProbsObject = objectList[0]
        if rxnProbsObject['info'][2] != RxnProbsType:
            message = 'RxnProbs object type %s is not %s for object %s' %(rxnProbsObject['info'][2], RxnProbsType, rxnProbsObject['info'][1])
            ctx.log_err(message)
            raise WrongVersionError(message)
        output = rxnProbsObject["data"]["reaction_probabilities"]
        if input['sort_field'] == 'rxnid':
            output.sort(key=lambda tup: tup[0])
        elif input['sort_field'] == 'probability':
            output.sort(key=lambda tup: tup[1], reverse=True)
        #END get_rxnprobs

        # At some point might do deeper type checking...
        if not isinstance(output, list):
            raise ValueError('Method get_rxnprobs return value ' +
                             'output is not type list as required.')
        # return the results
        return [output]

    def get_probanno(self, ctx, input):
        # ctx is the context object
        # return variables are: output
        #BEGIN get_probanno
        ''' Convert a probabilistic annotation object into a human-readbable table.

            @param ctx Current context object
            @param input Dictionary with input parameters for function
            @return Dictionary keyed by gene to a list of tuples with roleset and likelihood
            @raise WrongVersionError when ProbAnno object version number is invalid
        '''

        input = self._checkInputArguments(ctx, input,
                                          ['probanno', 'probanno_workspace'],
                                          { 'probanno_version': None }
                                          )

        wsClient = Workspace(self.config["workspace_url"], token=ctx['token'])
        probAnnoObjectId = make_object_identity(input["probanno_workspace"], input["probanno"], input['probanno_version'])
        objectList = wsClient.get_objects( [ probAnnoObjectId ] )
        probAnnoObject = objectList[0]
        if probAnnoObject['info'][2] != ProbAnnoType:
            message = 'ProbAnno object type %s is not %s for object %s' %(probAnnoObject['info'][2], ProbAnnoType, probAnnoObject['info'][1])
            ctx.log_err(message)
            raise WrongVersionError(message)
        output = probAnnoObject["data"]["roleset_probabilities"]

        #END get_probanno

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method get_probanno return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

#!/usr/bin/python

# Read and write data files
import os
import sys
import math
import subprocess
import json
import traceback
import time
from shock import Client as ShockClient
from biokbase import log
from biokbase.probabilistic_annotation.Helpers import now

# E values of less than 1E-200 are treated as 1E-200 to avoid log of 0 issues.
MIN_EVALUE = 1E-200

# Exception thrown when makeblastdb command failed
class MakeblastdbError(Exception):
    pass

# Exception thrown when static database file is missing from Shock.
class MissingFileError(Exception):
    pass

# Exception thrown when static database files are not ready
class NotReadyError(Exception):
    pass

''' Read and write data files. '''

class DataParser:
    
    def __init__(self, config):
        ''' Initialize the object.
        
            @param config Dictionary of configuration variables
        '''

        # Save the configuration variables related to data files.
        self.dataFolderPath = config['data_folder_path']
        self.separator = config['separator']
        self.searchProgram = config['search_program']
        self.searchProgramPath = config['search_program_path']
        self.shockURL = config['shock_url']
        self.loadDataOption = config['load_data_option']
       
        # Paths to files for tracking status of static database files.
        self.StatusFiles = dict()
        self.StatusFiles['status_file'] = os.path.join(self.dataFolderPath, 'staticdata.status')
        self.StatusFiles['cache_file'] = os.path.join(self.dataFolderPath, 'staticdata.cache')

        # Paths to files with source data.
        self.DataFiles = dict()
        self.DataFiles['otu_id_file'] = os.path.join(self.dataFolderPath, 'OTU_GENOME_IDS')
        self.DataFiles['subsystem_fid_file'] = os.path.join(self.dataFolderPath, 'SUBSYSTEM_FIDS')
        self.DataFiles['dlit_fid_file'] = os.path.join(self.dataFolderPath, 'DLIT_FIDS')
        self.DataFiles['concatenated_fid_file'] = os.path.join(self.dataFolderPath, 'ALL_FIDS')
        self.DataFiles['concatenated_fid_role_file'] = os.path.join(self.dataFolderPath, 'ALL_FID_ROLES')
        self.DataFiles['subsystem_otu_fid_roles_file'] = os.path.join(self.dataFolderPath, 'SUBSYSTEM_OTU_FID_ROLES')
        self.DataFiles['subsystem_otu_fasta_file'] = os.path.join(self.dataFolderPath, 'SUBSYSTEM_FASTA')
        self.DataFiles['complexes_roles_file'] = os.path.join(self.dataFolderPath, 'COMPLEXES_ROLES')
        self.DataFiles['reaction_complexes_file'] = os.path.join(self.dataFolderPath, 'REACTIONS_COMPLEXES')
        
        # Paths to files for searching for proteins.
        self.SearchFiles = dict()
        if self.searchProgram == 'usearch':
            self.SearchFiles['subsystem_udb_file'] = os.path.join(self.dataFolderPath, 'SUBSYSTEM.udb')
        else:
            self.SearchFiles['subsystem_otu_index_file'] = os.path.join(self.dataFolderPath, 'SUBSYSTEM_FASTA.pin')
            self.SearchFiles['subsystem_otu_sequence_file'] = os.path.join(self.dataFolderPath, 'SUBSYSTEM_FASTA.psq')
            self.SearchFiles['subsystem_otu_header_file'] = os.path.join(self.dataFolderPath, 'SUBSYSTEM_FASTA.phr')

        # Create the data folder if it does not exist.
        if not os.path.exists(config["data_folder_path"]):
            os.makedirs(config["data_folder_path"], 0775)

        return

    # The OTU ID file is a list of representative OTU genome IDs.  Each line has these fields:
    #   1. Genome ID in KBase format (e.g. kb|g.0)
    #   2. Flag indicating if the genome is a prokaryote (1 means yes, 0 means no)
    
    def readOtuData(self):
        ''' Read data from the representative OTU genome ID file.
        
            @return List of all OTU genome IDs, list of prokaryote OTU genome IDs
        '''
    
        fid = open(self.DataFiles['otu_id_file'], 'r')
        otus = list()
        prokotus = list()
        for line in fid:
            spl = line.strip("\r\n").split("\t")
            otus.append(spl[0])
            if int(spl[1]) == 1:
                prokotus.append(spl[0])
        fid.close()
        return otus, prokotus
    
    def writeOtuData(self, otus, prokotus):
        ''' Write data to the representative OTU genome ID file.
        
            @param otus List of all OTU genome IDs
            @param prokouts List of prokaryote OTU genome IDs
            @return Nothing
        '''
    
        fid = open(self.DataFiles['otu_id_file'], 'w')
        for otu in otus:
            if otu in prokotus:
                fid.write("%s\t%d\n" %(otu, 1))
            else:
                fid.write("%s\t%d\n" %(otu, 0))
        fid.close()
        return
    
    # The subsystem feature ID file is a list of feature IDs from SEED subsystems.  Each line has
    # one field that is the feature ID in KBase format (e.g. kb|g.3.peg.541).
    
    def readSubsystemFids(self):
        ''' Read data from the subsystem feature ID file.
        
            @return List of feature IDs from SEED subystems
        '''
    
        fid = open(self.DataFiles['subsystem_fid_file'], 'r')
        sub_fids = list()
        for line in fid:
            spl = line.strip("\r\n")
            sub_fids.append(spl)
        fid.close()
        return sub_fids
    
    def writeSubsystemFids(self, sub_fids):
        ''' Write data to the subsystem feature ID file.
        
            @param sub_fids List of feature IDs from SEED subsystems
            @return Nothing
        '''
    
        fid = open(self.DataFiles['subsystem_fid_file'], 'w')
        for f in sub_fids:
            fid.write("%s\n" %(f))
        fid.close()
        return
    
    # The direct literature-supported feature ID file is a list of feature IDs identified in the
    # literature.  Each line has one field that is the feature ID in KBase format (e.g. kb|g.428.peg.6254).
    
    def readDlitFids(self):
        ''' Read data from the direct literature-supported feature ID file.
        
            @return List of feature IDs from literature
        '''
    
        fid = open(self.DataFiles['dlit_fid_file'], 'r')
        otu_fids = list()
        for line in fid:
            spl = line.strip("\r\n")
            otu_fids.append(spl)
        fid.close()
        return otu_fids
    
    def writeDlitFids(self, otu_fids):
        ''' Write data to the direct literature-supported feature ID file.
        
            @param otu_fids List of feature IDs from literature
            @return Nothing
        '''
    
        fid = open(self.DataFiles['dlit_fid_file'], 'w')
        for f in otu_fids:
            fid.write("%s\n" %(f))
        fid.close()
        return
    
    # The concatenated feature ID file is a list of all feature IDs from both
    # SEED subsystems and literature.  Each line has one field that is the feature
    # ID in KBase format (e.g. kb|g.428.peg.6254).

    def readAllFids(self):
        ''' Read data from the concatenated feature ID file.
        
            @return List of all feature IDs
        '''
        
        fid = open(self.DataFiles['concatenated_fid_file'], 'r')
        all_fids = list()
        for line in fid:
            spl = line.strip('\r\n')
            all_fids.append(spl)
        fid.close()
        return all_fids
        
    def writeAllFids(self, all_fids):
        ''' Write data to the concatenated feature ID file.
        
            @param all_fids List of all feature IDs
            @return Nothing
        '''
        
        fid = open(self.DataFiles['concatenated_fid_file'], 'w')
        for f in all_fids:
            fid.write("%s\n" %(f))
        fid.close()
        return
    
    # The concatenated feature ID to role file is a mapping of feature IDs to functional roles.
    # Each line has these fields:
    #   1. Feature ID in KBase format (e.g. kb|g.0.peg.2094)
    #   2. List of names of functional roles (e.g. Conserved ATP-binding protein YghS)
    #
    # Note that functional roles must be separated by a string that does not occur in any role.
    
    def readAllFidRoles(self):
        ''' Read data from the concatenated feature ID to role file.
        
            @return Dictionary mapping a feature ID to list of names of roles, dictionary mapping a role to feature ID
        '''
    
        fid = open(self.DataFiles['concatenated_fid_role_file'], 'r')
        all_fidsToRoles = dict()
        for line in fid:
            spl = line.strip("\r\n").split("\t")
            roles = spl[1].split(self.separator)
            if spl[0] in all_fidsToRoles:
                all_fidsToRoles[spl[0]] += roles
            else:
                all_fidsToRoles[spl[0]] = roles
        fid.close()
    
        all_rolesToFids = dict()
        for fid in all_fidsToRoles:
            roles = all_fidsToRoles[fid]
            for role in roles:
                if role in all_rolesToFids:
                    all_rolesToFids[role].append(fid)
                else:
                    all_rolesToFids[role] = [ fid ]
    
        return all_fidsToRoles, all_rolesToFids
    
    def writeAllFidRoles(self, otu_fidsToRoles):
        ''' Write data to the concatenated feature ID to role file.
        
            @param otu_fidsToRoles Dictionary mapping a feature ID to list of names of roles
            @return Nothing
        '''
    
        fid = open(self.DataFiles['concatenated_fid_role_file'], 'w')
        for f in otu_fidsToRoles:
            fid.write("%s\t%s\n" %(f, self.separator.join(otu_fidsToRoles[f])))
        fid.close()
        return
    
    # The filtered feature ID to roles file is a mapping of feature IDs to functional roles where
    # there is one protein from each OTU for each functional role.  Each line has these fields:
    #   1. Feature ID in KBase format (e.g. kb|g.0.peg.2094)
    #   2. List of names of functional roles (e.g. Conserved ATP-binding protein YghS)
    #
    # Note that functional roles must be separated by a string that does not occur in any role.
    
    def readFilteredOtuRoles(self):
        ''' Read data from the filtered feature ID to roles file.
        
            @return Dictionary mapping a feature ID to list of names of roles, dictionary mapping a role to feature ID
        '''
    
        fid = open(self.DataFiles['subsystem_otu_fid_roles_file'], 'r')
        otu_fidsToRoles = dict()
        for line in fid:
            spl = line.strip("\r\n").split("\t")
            roles = spl[1].split(self.separator)
            if spl[0] in otu_fidsToRoles:
                otu_fidsToRoles[spl[0]] += roles
            else:
                otu_fidsToRoles[spl[0]] = roles
        fid.close()
    
        otu_rolesToFids = dict()
        for fid in otu_fidsToRoles:
            roles = otu_fidsToRoles[fid]
            for role in roles:
                if role in otu_rolesToFids:
                    otu_rolesToFids[role].append(fid)
                else:
                    otu_rolesToFids[role] = [ fid ]
    
        return otu_fidsToRoles, otu_rolesToFids
    
    def writeFilteredOtuRoles(self, otu_fidsToRoles):
        ''' Write data to the filtered feature ID to roles file.
        
            @param otu_fidsToRoles Dictionary mapping a feature ID to list of names of roles
            @return Nothing
        '''
    
        fid = open(self.DataFiles['subsystem_otu_fid_roles_file'], 'w')
        for f in otu_fidsToRoles:
            fid.write("%s\t%s\n" %(f, self.separator.join(otu_fidsToRoles[f])))
        fid.close()
        return
    
    # The subsystem FASTA file contains the amino acid sequences for a set of feature IDs.  When the file
    # is written a search database is automatically generated.
    
    def readSubsystemFasta(self):
        ''' Read data from the subsystem FASTA file.
        
            @note This function does not return data (not sure it is really needed except for symmetry).
            @return Nothing
        '''
    
        fid = open(self.DataFiles['subsystem_otu_fasta_file'], 'r')
        fid.close()
        return
    
    def writeSubsystemFasta(self, fidsToSeqs):
        ''' Write data to the subsystem FASTA file.
        
            @param fidsToSeqs Dictionary mapping a feature ID to amino acid sequence
            @return Nothing
        '''
    
        fid = open(self.DataFiles['subsystem_otu_fasta_file'], 'w')
        # Sort the fids so that fasta files containing the same proteins hash to the same MD5 (for
        # data provenance purposes)
        for fids in sorted(fidsToSeqs.keys()):
            fid.write(">%s\n%s\n" %(fids, fidsToSeqs[fids]))
        fid.close()
        return
    
    def buildSearchDatabase(self):
        ''' Build a search database for the configured search program.

            @note Make sure the subsystem FASTA file is available.
            @return Nothing
        '''

        # Build the command based on the configured search program.
        if self.searchProgram == 'usearch':
            args = [ self.searchProgramPath, '-makeudb_ublast', self.DataFiles['subsystem_otu_fasta_file'], '-output', self.SearchFiles['subsystem_udb_file'] ]
        else:
            args = [ "/usr/bin/makeblastdb", "-in", self.DataFiles['subsystem_otu_fasta_file'], "-dbtype", "prot" ]
    
        # Run the command to compile the database from the subsystem fasta file.
        try:
            proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            (stdout, stderr) = proc.communicate()
            if proc.returncode < 0:
                cmd = ' '.join(args)
                raise MakeblastdbError("'%s' was terminated by signal %d" %(cmd, -proc.returncode))
            else:
                if proc.returncode > 0:
                    cmd = ' '.join(args)
                    details = "'%s' failed with return code %d:\nCommand: '%s'\nStdout: '%s'\nStderr: '%s'" \
                        %(args[0], proc.returncode, cmd, stdout, stderr)
                    raise MakeblastdbError(details)
        except OSError as e:
            cmd = ' '.join(args)
            raise MakeblastdbError("Failed to run '%s': %s" %(cmd, e.strerror))
        return
    
    def parseBlastOutput(self, blastResultsPath):
        ''' Read BLAST results file and store in a convenient structure.

            The results file is in BLAST output format 6 where each line describes an alignment
            found by the search program.  A line has 12 tab delimited fields: (1) query label,
            (2) target label, (3) percent identity, (4) alignment length, (5) number of
            mismatches, (6) number of gap opens, (7) 1-based position of start in query,
            (8) 1-based position of end in query, (9) 1-based position of start in target,
            (10) 1-based position of end in target, (11) e-value, and (12) bit score.
 
            @note Score is the negative log E-value
            @param blastResultsPath Path to BLAST results file
            @return Dictionary mapping query ID to tuple of target ID and score
        '''
    
        idToTargetList = dict()
        for line in open(blastResultsPath, 'r'):
            fields = line.strip('\r\n').split('\t')
            queryid = fields[0]
            targetid = fields[1]
            if float(fields[11]) < 0.0: # Throw out alignments with a negative bit score
                print 'throwing out %s' %(line)
                continue
            logeval = -1.0 * math.log10(float(fields[10]) + MIN_EVALUE)
            tup = ( targetid, logeval )
            if queryid in idToTargetList:
                idToTargetList[queryid].append( tup )
            else:
                idToTargetList[queryid] = [ tup ]
        return idToTargetList
    
    # The complexes to roles file contains a mapping of complex IDs to functional roles.
    # Each line has these fields:
    #   1. Complex ID in KBase format (e.g. kb|cpx.1048)
    #   2. List of functional roles for the complex
    #
    # Note that functional role names must be separated by a string that does not occur in any name.

    def readComplexRoles(self):
        ''' Read data from the complex to roles file.
        
            @return Dictionary mapping a complex ID to list of names of functional roles
        '''
    
        fid = open(self.DataFiles['complexes_roles_file'], 'r')
        complexToRequiredRoles = dict()
        for line in fid:
            spl = line.strip("\r\n").split("\t")
            complexes = spl[0]
            roles = spl[1].split(self.separator)
            # This shouldn't be necessary but just to be safe...
            if complexes in complexToRequiredRoles:
                complexToRequiredRoles[complexes] += roles
            else:
                complexToRequiredRoles[complexes]  = roles
        fid.close()
        return complexToRequiredRoles
    
    def writeComplexRoles(self, complexToRequiredRoles):
        ''' Write data to the complex to roles file.
        
            @param complexToRequiredRoles Dictionary mapping a complex ID to list of names of functional roles
            @return Nothing
        '''
    
        fid = open(self.DataFiles['complexes_roles_file'], 'w')
        for complexes in complexToRequiredRoles:
            fid.write("%s\t%s\n" %(complexes, self.separator.join(complexToRequiredRoles[complexes])))
        fid.close()
        return
    
    # The reaction to complexes file contains a mapping of reaction IDs to complex IDs.
    # Each line has these fields:
    #   1. Reaction ID in KBase format (e.g. kb|rxn.5682)
    #   2. List of complex IDs in KBase format (e.g. kb|cpx.1507///kb|cpx.1813)
    #
    # Note that complex IDs must be separated by a string that does not occur in any complex ID.
    
    def readReactionComplex(self):
        ''' Read data from the reaction to complexes file.
        
            @return Dictionary mapping a reaction ID to list of complex IDs
        '''
    
        fid = open(self.DataFiles['reaction_complexes_file'], 'r')
        rxnToComplexes = dict()
        for line in fid:
            spl = line.strip("\r\n").split("\t")
            rxn = spl[0]
            cplxlist = spl[1].split(self.separator)
            # This shouldn't be necessary but just to be safe...
            if rxn in rxnToComplexes:
                rxnToComplexes[rxn] += cplxlist
            else:
                rxnToComplexes[rxn] = cplxlist
        fid.close()
        return rxnToComplexes
    
    def writeReactionComplex(self, rxnToComplexes):
        ''' Write data to reaction to complexes file.
        
            @param rxnToComplexes Dictionary mapping a reaction ID to list of complex IDs
            @return Nothing
        '''
    
        fid = open(self.DataFiles['reaction_complexes_file'], 'w')
        for rxn in rxnToComplexes:
            fid.write("%s\t%s\n" %(rxn, self.separator.join(rxnToComplexes[rxn])))
        fid.close()
        return
    
    def readRolesetProbabilityFile(self, roleset_probability_file):
        ''' Read the roleset probability file.
        
            @param roleset_probability_file Path to roleset probability file
            @return Dictionary mapping query ID to list of tuples with list of roles and probability
        '''
    
        queryToTuplist = dict()
        for line in open(roleset_probability_file, 'r'):
            spl = line.strip("\r\n").split("\t")
            if spl[0] in queryToTuplist:
                queryToTuplist[spl[0]].append( (spl[1], float(spl[2])) )
            else:
                queryToTuplist[spl[0]] = [ (spl[1], float(spl[2])) ]
        return queryToTuplist
    
    # The status file is used to track the status of setting up the static database files when
    # the server starts.  The first line of the file contains the status which is one of
    # these values:
    #   1. 'building' when the pa-gendata command is building the files
    #   2. 'running' when the server initialization is in progress
    #   3. 'ready' when the server initialization is complete or a build is complete
    #   4. 'failed' when there was an error building or loading the files
    #
    # The second line has the timestamp of when the status was last changed.
    
    def readStatusFile(self):
        ''' Read the current status value from the status file.
        
            @return Current status string
        '''
    
        fid = open(self.StatusFiles['status_file'], 'r')
        statusLine = fid.readline()
        fid.close()
        return statusLine.strip("\r\n")
    
    def writeStatusFile(self, status):
        ''' Write new status value to the status file.
        
            @param status New status value
            @return Nothing
        '''
    
        fid = open(self.StatusFiles['status_file'], 'w')
        fid.write("%s\nupdated at %s\n" %(status, now()))
        fid.close()
        return
    
    def checkIfDatabaseFilesExist(self):
        ''' Check for existence of all of the database files.
        
            @raise NotReadyError: A database file does not exist
            @return Nothing
        '''
    
        for path in self.DataFiles.values():
            if not os.path.exists(path):
                raise NotReadyError("Static database file '%s' does not exist" %(path))
        for path in self.SearchFiles.values():
            if not os.path.exists(path):
                raise NotReadyError("Static database file '%s' does not exist" %(path))
        return

    def loadDatabaseFiles(self, mylog):
        ''' Load the static database files from Shock.

            The static database files are stored in the directory specified by the
            data_folder_path configuration variable.  A file is only downloaded if
            the file is not available on this system or the file has been updated
            in Shock.

            @param mylog Log object for messages
            @return Nothing
            @raise MissingFileError when database file is not found in Shock
        '''
        
        # Get the current info about the static database files from the cache file.
        cacheFilename = self.StatusFiles['cache_file']
        if os.path.exists(cacheFilename):
            fileCache = json.load(open(cacheFilename, "r"))
        else:
            fileCache = dict()
        
        # Create a shock client.
        shockClient = ShockClient(self.shockURL)

        # See if the static database files on this system are up-to-date with files stored in Shock.
        shockFiles = dict(self.DataFiles.items() + self.SearchFiles.items())
        for key in shockFiles:
            # Get info about the file stored in Shock.
            localPath = shockFiles[key]
            name = os.path.basename(localPath)
            nodelist = shockClient.query_node( { 'lookupname': 'ProbAnnoData/'+name } )
            if len(nodelist) == 0:
                message = "Database file %s is not available from %s\n" %(name, self.shockURL)
                mylog.log_message(log.ERR, message) # MBM
                raise MissingFileError(message)
            node = nodelist[0]
            
            # Download the file if the checksum does not match or the file is not available on this system.
            download = False
            if key in fileCache:
                if node['file']['checksum']['md5'] != fileCache[key]['file']['checksum']['md5']:
                    download = True
            else:
                download = True
            if os.path.exists(localPath) == False:
                download = True
            if download:
                sys.stderr.write("Downloading %s to %s\n" %(key, localPath))
                shockClient.download_to_path(node["id"], localPath)
                fileCache[key] = node
                mylog.log_message(log.INFO, 'Downloaded %s to %s' %(key, localPath))
                
        # Save the updated cache file.
        json.dump(fileCache, open(cacheFilename, "w"), indent=4)
        return
     
    def storeDatabaseFiles(self, token):
        ''' Store the static database files to Shock.

            @param token: Authorization token for authenticating to shock
            @return Nothing
        '''
        
        # Create a shock client.
        shockClient = ShockClient(self.shockURL, token=token)
        
        # Upload all of the static database files to shock.
        fileCache = dict()
        shockFiles = dict(self.DataFiles.items() + self.SearchFiles.items())
        for key in shockFiles:
            localPath = shockFiles[key]
            name = os.path.basename(localPath)
            if os.path.exists(localPath):
                sys.stderr.write("Saving '%s'..." %(localPath))
                
                # See if the file already exists in Shock.
                query = { 'lookupname': 'ProbAnnoData/'+name }
                nodelist = shockClient.query_node(query)
                
                # Remove all instances of the file in Shock.
                if nodelist != None:
                    for node in nodelist:
                        shockClient.delete_node(node["id"])
     
                # Build the attributes for this file and store as json in a separate file.
                moddate = time.ctime(os.path.getmtime(localPath))           
                attr = { "lookupname": "ProbAnnoData/"+name, "moddate": moddate }
                attrFilename = os.path.join(self.dataFolderPath, name+".attr")
                attrFid = open(attrFilename, "w")
                json.dump(attr, attrFid, indent=4)
                attrFid.close()
                
                # Upload the file to Shock.
                metadata = shockClient.create_node(localPath, attrFilename)
                fileCache[key] = metadata
                os.remove(attrFilename)
                
                # Remove the list of users from the read ACL to give the file public read permission.
                # Note this needs to change for Shock version 0.9.5 but not sure how to set public ACLs.
                readacl = shockClient.get_acl(metadata["id"])
                shockClient.delete_acl(metadata['id'], 'read', readacl['read'][0])
                sys.stderr.write("done\n")
                
            else:
                sys.stderr.write("Could not find '%s' so it was not saved\n" %(localPath))
                
        # Save the metadata on all of the database files.
        cacheFilename = os.path.join(self.dataFolderPath, self.StatusFiles["cache_file"])
        json.dump(fileCache, open(cacheFilename, "w"), indent=4)

        return

    def getDatabaseFiles(self, mylog, testDataPath):
        ''' Get the static database files.

            The static database files come from one of three places: (1) Shock,
            (2) a local data directory, (3) a test data directory.  If the files
            are not found in Shock, the local data directory is searched.  If the
            file are not found in the local data directory, the test data directory
            is used.

            @param mylog: Log object for messages
            @param testDataPath: Path to directory with test database files
            @return Current value of load data option which indicates which of the
                three places is being used for the static database files
        '''

        # Update the status file to indicate that the static database files are being updated.
        self.writeStatusFile('running')
        status = 'failed'

        # Get the static database files from Shock (only missing or changed files are downloaded).
        if self.loadDataOption == 'shock':
            try:
                self.loadDatabaseFiles(mylog)
                status = 'ready'
                sys.stderr.write('All static database files loaded from Shock to %s.\n' %(self.dataFolderPath))
                mylog.log_message(log.INFO, 'All static database files loaded from Shock to %s' %(self.dataFolderPath))
            except:
                traceback.print_exc(file=sys.stderr)
                sys.stderr.write('WARNING: Failed to load static database files from Shock. Checking current files but they might not be the latest!\n')
                mylog.log_message(log.NOTICE, 'Failed to load static database files from Shock. Checking current files...')
                self.loadDataOption = 'preload'

        # Get the static database files from the data directory specified in the configuration.
        if self.loadDataOption == 'preload':
            try:
                self.checkIfDatabaseFilesExist()
                status = 'ready'
                sys.stderr.write('All static database files are available in %s.\n' %(self.dataFolderPath))
                mylog.log_message(log.INFO, 'All static database files are available in %s' %(self.dataFolderPath))
            except:
                # There is a problem with at least one of the static database files so switch
                # to the test data.
                status = 'ready'
                self.loadDataOption = 'test'
                self.dataFolderPath = testDataPath
                traceback.print_exc(file=sys.stderr)
                sys.stderr.write('WARNING: Static database files are missing. Switched to test database files in %s.\n' %(testDataPath))
                mylog.log_message(log.NOTICE, 'Static database files are missing. Switched to test database files in %s' %(testDataPath))

        # Update the status file to indicate that the static database files updating is done.
        self.writeStatusFile(status)
        return self.loadDataOption

    #####################
    # OTU neighborhoods #
    #####################
    
    #def readOtuNeighborhoods(folder):
    #    fid = open(os.path.join(folder, OTU_NEIGHBORHOOD_FILE), "r")##
    
    #    tuplist = []
    #    fidToRoles = {}
    #    for line in fid:
    #        spl = line.strip("\r\n").split("\t")
    #        roles = spl[4].split(SEPARATOR)
    #        if spl[1] in fidToRoles:
    #            fidToRoles[spl[1]] += roles
    #        else:
    #            fidToRoles[spl[1]] = roles
    #        tuplist.append( (spl[0], spl[1], spl[2], spl[3],) )
    #    fid.close()
    #    return tuplist, fidToRoles
    
    #def writeOtuNeighborhoods(tuplist, fidToRoles, verbose, fname):
    #    fid = open(fname, "w")
    #    for f in tuplist:
    #        if f[1] in fidToRoles:
    #            roles = fidToRoles[f[1]]
    #        else:
    #            if verbose:
    #                sys.stderr.write("WARNING: Fid %s has no role despite being a neighbor of an OTU gene!\n" %(f[1]) )
    #            roles = ""
    #        try:
    #            fid.write("%s\t%s\t%s\t%s\t%s\n" %(f[0], f[1], f[2], f[3], SEPARATOR.join(roles)))
    #        except UnicodeEncodeError:
    #            sys.stderr.write("ERROR: encountered roles that contain non-ASCII characters?\n")
    #            sys.stderr.write("In gene ID %s\n" %(f[1]))
    #            sys.stderr.write("Skipping...\n")
    #            continue
    #    fid.close()
    #    return
    
    #####################
    # Complex --> roles #
    #####################
    
    # The complex to roles file contains a mapping of complex IDs to functional roles.
    # Each line has these fields:
    #   1. Complex ID in KBase format (e.g. kb|cpx.51240)
    #   2. List of names of functional roles (e.g. V-type ATP synthase subunit H (EC 3.6.3.14))
    #
    # Note that functional roles must be separated by a string that does not occur in any role.
    


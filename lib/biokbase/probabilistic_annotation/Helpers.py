#! /usr/bin/python

import os
import sys
import time
from biokbase.auth import kb_config
from ConfigParser import ConfigParser

# Default URL for production server
DefaultURL = 'https://kbase.us/services/probabilistic_annotation/'

# Current version number of ProbAnno object
ProbAnnoType = 'ProbabilisticAnnotation.ProbAnno-1.0'

# Current version number of RxnProbs object
RxnProbsType = 'ProbabilisticAnnotation.RxnProbs-1.0'

# Current version of service.
ServiceVersion = '1.1.0'

def read_config(filename=None):
    ''' Read a configuration file.

        If filename is None, the default configuration file is read.  If the
        probabilistic_annotation section is not present in the configuration
        file, it is automatically added.

        @param filename Path to configuration file.
        @return Config object
    '''

    # Use default config file if one is not specified.
    if filename == None:
        filename = os.path.join(os.environ['KB_TOP'], 'deployment.cfg')

    # Read the config file.
    config = ConfigParser()
    try:
        config.read(filename)
    except Exception as e:
        print "Error while reading config file %s: %s" % (filename, e)

    # Make sure there is a probabilistic_annotation section in the config file.
    if not config.has_section('probabilistic_annotation'):
        config.add_section('probabilistic_annotation')
        with open(filename, 'w') as configfile:
            config.write(configfile)

    return config

def get_config(filename=None):
    ''' Get the probabilistic_annotation section from a configuration file.

        If filename is None, the default configuration file is read.

        @param filename Path to configuration file.
        @return Dictionary mapping configuration variables to values
    '''

    # Read the config file.
    config = read_config(filename)

    # Extract the probabilistic annotation section.
    sectionConfig = dict()
    for nameval in config.items("probabilistic_annotation"):
        sectionConfig[nameval[0]] = nameval[1]
    return sectionConfig

def get_url():
    ''' Get the current URL for the service.

        @return Current URL string
    '''

    # Just return the default URL when running in IRIS.
    if 'KB_RUNNING_IN_IRIS' in os.environ:
        return DefaultURL

    # Get the URL from the config file or use the default if it is not set.
    config = get_config(kb_config)
    if 'url' in config:
        currentURL = config['url']
    else:
        currentURL = DefaultURL;
    return currentURL

def set_url(newURL):
    ''' Set the current URL for the service and store in config file.

        @param newURL New value for URL
        @return New URL string
    '''

    # Check for special value for the default URL.
    if newURL ==  'default':
        newURL = DefaultURL

    # Just return the URL when running in IRIS. There is no place to save it.
    if 'KB_RUNNING_IN_IRIS' in os.environ:
        return newURL

    # Save the new URL to the config file.
    config = read_config(kb_config)
    config.set('probabilistic_annotation', 'url', newURL)
    with open(kb_config, 'w') as configfile:
        config.write(configfile)
    return newURL

def timestamp(deltaSeconds):
    ''' Get a timestamp in the format required by user and job state service.

        @param deltaSeconds Seconds added to the current time to get a time in the future
        @return Formatted timestamp string
    '''

    # Use UTC timestamps to avoid timezone issues.
    now = time.time() + deltaSeconds
    ts = time.gmtime(time.time() + deltaSeconds)
    return time.strftime('%Y-%m-%dT%H:%M:%S+0000', ts)

def now():
    ''' Get the current time as a printable string.

        @return Formatted timestamp string in local timezone
    '''

    return '%s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

def job_info_dict(infoTuple):
    ''' Convert a job info tuple into a dictionary.

        @param infoTuple Job info tuple returned by user and job state service functions
        @returns Dictionary of job info data from tuple
    '''

    info = dict()
    info['id'] = infoTuple[0]
    info['service'] = infoTuple[1]
    info['stage'] = infoTuple[2]
    info['started'] = infoTuple[3]
    info['status'] = infoTuple[4]
    info['last_update'] = infoTuple[5]
    info['total_progress'] = infoTuple[6]
    info['max_progress'] = infoTuple[7]
    info['progress_type'] = infoTuple[8]
    info['est_complete'] = infoTuple[9]
    info['complete'] = infoTuple[10]
    info['error'] = infoTuple[11]
    info['description'] = infoTuple[12]
    info['results'] = infoTuple[13]
    return info

def make_object_identity(workspace, object, ver=None):
    ''' Make an object identity structure.

        @param workspace Name or number of workspace containing object
        @param object Name or number of object
        @param ver Optional version number of object
        @returns ObjectIdentity structure for workspace APIs
    '''

    objectIdentity = dict()
    if workspace.isdigit():
        objectIdentity['wsid'] = workspace
    else:
        objectIdentity['workspace'] = workspace
    if object.isdigit():
        objectIdentity['objid'] = object
    else:
        objectIdentity['name'] = object
    if ver is not None:
        objectIdentity['ver'] = ver
    return objectIdentity

def make_job_directory(workDirectory, jobID):
    ''' Make working directory for a job.

        @param workDirectory Path to base working directory
        @param jobID Job identifier
        @returns Path to job directory
    '''

    jobDirectory = os.path.join(workDirectory, jobID)
    if not os.path.exists(jobDirectory):
        os.makedirs(jobDirectory, 0775)
    return jobDirectory

# Our Project And KBase (v2)

Our project currently depends on KBase, a set of software tools for metabolic modeling, genome annotation, and other
biology computing tools. Their website is here: [https://kbase.us/](https://kbase.us/)

## How KBase Services Work
KBase organizes their tools into services, which are included into their overall project. A service is composed of
an API using standard KBase data types, definitions for new types (ReactionProbabilities in the case of ProbAnno), 
and other code required for working with their infrastructure.

### KBase Services we use
- **KBaseFBAModeling**: The service we use to to run Flux Balance Analysis (FBA), add and remove reactions to a model, model translation, etc.
- **GenomeComparison**: The service we use to create a ProteomeComparison as part of our process (see graphic)
- **ProbabilisticAnnotation**: The service for assigning likelihoods to reactions in a database for inclusion in a model

##### Process Graphic:
<img src="https://github.com/kingb12/model-morphing/blob/master/KBaseDep.png" alt="Dependencies" style="width: 200px;"/>



#### Using a KBase Service
KBase services can be accessed in a number of ways:
- **Narrative Interface** The web interface for the services. [Narrative Log-in](https://narrative.kbase.us/#login)
- **Command Line Interface** We have the KBase CLI Installed on our AWS Instance
- **Language Clients** Programming API's are supported for Python, Java, and Perl for most services. Our code is written as a client to the Python API's for each service we use

The Python API's are used with instances of clients for each service, initialized with a URL for the service in question. Here's an example for initializing a client for
the *KBaseFBAModeling* service. The URLs for each service are found at: [https://kbase.us/services/](https://kbase.us/services/). There are also development and
Continuous Integration versions of each service, found at [https://next.kbase.us/services/](https://next.kbase.us/services/) and [https://ci.kbase.us/services/](https://ci.kbase.us/services/) respectively.

Here's an example: 

```python
from biokbase.fbaModelServices.Client import fbaModelServices
fba_client = fbaModelServices('https://kbase.us/services/KBaseFBAModeling/')
```

Arguments to each function for a service are almost always passed as dictionaries, with the keys specified in the .spec file for the client (in the repository).

Here's another example with the KBaseFBAModeling Service. The KBaseFBAModeling service is on [GitHub](https://github.com/kbase/KBaseFBAModeling),
and can be cloned. Here's a code sample that will create an FBA Model from a genome. The
workspace in question can be seen in the narrative interface here: 
[https://narrative.kbase.us/narrative/notebooks/ws.14581.obj.1](https://narrative.kbase.us/narrative/notebooks/ws.14581.obj.1)


Here are the headers for the methods we'll call:


```
/* Input parameters for the "genome_to_fbamodel" function.
	
		genome_id genome - ID of the genome for which a model is to be built (a required argument)
		workspace_id genome_workspace - ID of the workspace containing the target genome (an optional argument; default is the workspace argument)
		template_id templatemodel - 
		workspace_id templatemodel_workspace - 
		bool probannoOnly - a boolean indicating if only the probabilistic annotation should be used in building the model (an optional argument; default is '0')
		fbamodel_id model - ID that should be used for the newly constructed model (an optional argument; default is 'undef')
		bool coremodel - indicates that a core model should be constructed instead of a genome scale model (an optional argument; default is '0')
		workspace_id workspace - ID of the workspace where the newly developed model will be stored; also the default assumed workspace for input objects (a required argument)
		string auth - the authentication token of the KBase account changing workspace permissions; must have 'admin' privelages to workspace (an optional argument; user is "public" if auth is not provided)
		
	*/
    typedef structure {
		genome_id genome;
		workspace_id genome_workspace;
		template_id templatemodel;
		workspace_id templatemodel_workspace;
		fbamodel_id model;
		bool coremodel;
		workspace_id workspace;
		string auth;
		bool fulldb;
    } genome_to_fbamodel_params;
    /*
        Build a genome-scale metabolic model based on annotations in an input genome typed object
    */
    authentication required;
    funcdef genome_to_fbamodel(genome_to_fbamodel_params input) returns (object_metadata modelMeta);
 ```
 
 
 In the workspace linked above, we have a genome 'Methanococcus_maripaludis_S2'. 
 The arguments we need to include are the `genome`, `genome_workspace`,`model`, `workspace`. The `genome` and `genome_workspace`
  arguments identify the genome object as input, and the `model` and `workspace` arguments identify the FBAModel that will result from the call.
 
 ```python
from biokbase.fbaModelServices.Client import fbaModelServices
fba_client = fbaModelServices('https://kbase.us/services/KBaseFBAModeling/')
args = {'genome': 2,
        'genome_workspace': 14581,
        'model': 'MyNewMaripaludisModel',
        'workspace': 14581
       }
meta_info = fba_client.genome_to_fbamodel(args)
```

Now there wil be an FBA Model called 'MyNewMaripaludisModel' in the workspace 14581, which can be seen form the narrative interface
as linked above. `meta_info` will be a list of information about the model, including it's object_id, workspace, type, creator, etc.
 
 
## Probabilistic Annotation and KBase
Probabalistic Annotation (or ProbAnno) is a tool written by Matt Benedict and Mike Mundy that uses
genome annotations to assign likelihoods for metabolic reactions in a database. For this reason, 
ProbAnno is currently dependent on KBase as the service that runs the algorithm and the database for
the reactions and genomes. The only problem is that it isn't available on the Narrative service and 
is incompatible with updates to KBase: It crashes when run using the KBase CLI.


[Probabilistic Annotation Repo](https://github.com/kbase/probabilistic_annotation)


[Paper describing Algorithm and its Applications](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003882)









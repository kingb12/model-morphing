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
![KBase Dependencies](https://github.com/kingb12/model-morphing/blob/master/KBaseDep.png)


## Probabilistic Annotation and KBase
Probabalistic Annotation (or ProbAnno) is a tool written by Matt Benedict and Mike Mundy that uses
genome annotations to assign likelihoods for metabolic reactions in a database. For this reason, 
ProbAnno is currently dependent on KBase as the service that runs the algorithm and the database for
the reactions and genomes. The only problem is that it isn't available on the Narrative service and 
is incompatible with updates to KBase: It crashes when run using the KBase CLI.








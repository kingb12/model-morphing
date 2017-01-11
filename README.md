A tool for morphing a metabolic model from one genome to a closely related one, without ever disrupting growth

# Morphing a Model for Organism 'A' to Organism 'B'

Takes a Model/Genome for a 'source' organism and 'morphs' it to the 'target' organism. It keeps reactions from the 'source' model for which there is at least one matching gene feature, and attempts to remove those in the source model that do not have BLAST hits in the second genome. It also adds reactions annotated within the target genome that are unique to the target organism. The model is never broken in the process, such that the Biomass reaction and any other specified reaction [insert how you flag this option] must always carry flux. If a reaction can't be removed without breaking the model, it is removed and gap-filled using Prob-Anno gapfilling [other options?]

## Requirements

This tool runs within a Docker container, so it a working installation of Docker is a pre-requisite. Information on how to install Docker for your system can be found at [docker.com](https://www.docker.com/) and elsewhere online.

## Installation
1. Clone this repository: `git clone https://github.com/kingb12/model-morphing.git`
2. Navigate to `model-morphing/`
3. Make sure the docker daemon is running (i.e. you can run `docker ...` commands)
4. Build the image: `docker build -t kingb12/mm .`

## Running
5. Run the image:  `docker run -t -i kingb12/mm .`
6. With no other arguments, this will bring you into an **iPython console** with appropriate modules and classes from `model_morphing` imported. To run the container with some other command (e.g. bash), use: `docker run -t -i kingb12/mm /bin/bash`. This can be customized as needed.

## Tutorial
A walk-through of the steps for completing a morph, and touching on possible customizations

### Instantiating a Morph
For easy management of the necessary data for morphing, a Morph class exists in the Morph module, which has instance data relevant to your particular morph, and methods for each step in the process.

These are the minimal required arguments to instantiate a morph. Many of these are subclass instances of StoredObject's, a wrapper abstraction for objects in KBase. If an object is stored in KBase, then it can be uniquely referenced by its object_id and workspace_id. 
 - `genome`: a target genome (Genome object). A Genome is initialized with an object_id, workspace_id for it storage in KBase (e.g. `my_genome = Genome(object_id, workspace_id)`).
 - `src_model`: a source model (FBAModel object). A Model is initialized with an object_id, workspace_id for it storage in KBase (e.g. `my_model = FBAModel(object_id, workspace_id)`).
 - `probanno`: reaction probabilities for the target genome (ReactionProbabailities object). This holds likelihoods for the presence in the target organism of all KBase reactions. A ReactionProbabilities is initialized with an object_id, workspace_id for it storage in KBase (e.g. `my_rxn_probs = ReactionProbabilities(object_id, workspace_id)`).
 - `protcomp`: a proteome comparison between the genomes of the source and target organism (ProteomeComparison object). A ProteomeComparison is initialized with an object_id, workspace_id for it storage in KBase (e.g. `my_prot_comp = ProteomeComparison(object_id, workspace_id)`).
 - `media`: media that the source model grows on (Media object). If the source model does not grow on any media, you can specify one you wish the model to grow on when completed for the target organism, and gapfill using `fill_src_to_media` method or the `fill_src=True` flag in `prepare_supermodel`.
 
 Each of these arguments are passed as keyword arguments, in a dictionary of arguments, or a combination. Keyword arguments can be used to override a dictionary argument. e.g:
 
 ```
 # As keyword arguments:
 morph = Morph(src_model=FBAModel(58, 9145), genome=Genome(26, 9145), protcomp=ProteomeComparison(50, 9145), probanno=ReactionProbabilities(30, 9145), media=Media(24, 9145)) 
 
 # As a dictionary
 morph = Morph({'src_model': FBAModel(58, 9145), 'genome': Genome(26, 9145), 'media': Media(24, 9145), 'probanno': ReactionProbabilities(30, 9145), 'protcomp': ProteomeComparison(50, 9145)})
 
 # Combination (with overwriting)
 morph = Morph({'src_model': FBAModel(58, 9145), 'genome': Genome(26, 9145), 'media': Media(24, 9145), 'probanno': ReactionProbabilities(30, 9145), 'protcomp': ProteomeComparison(50, 9145)}, media=Media(12, 14623))
 
 # morph.media will be Media(12, 14623), allowing customizations via keyword arguments on top of a basic morph specified as an argument dictionary. 
 ```
 
##### ws_id attribute
 useful optional attribute to set is `ws_id`. Setting `ws_id` allows you to control what workspace the morph bases itself in.
Default behavior is to create a new workspace for the user upon instantiation, but this doesn't work with KBase narrative. To view morph data in a KBase narrative, the simplest way is to create a new narrative in KBase and supply it's workspace_id to the Morph constructor via `ws_id`. e.g. `morph = Morph(arg_dictionary, ws_id=12345)`. If you have a narrative whose worksace is 12345, you can now see all objects created in the process of the morph.

### Building a SuperModel

The SuperModel is the model composed of all reactions in your source model, with all the unique ones we can recover from a draft reconstruction of the target genome. It is all the reactions that will be considered for the final model. Additionally, in the process of building a supermodel, a morph is able to categorize reactions as having genetic basis in one,both or no organisms. In the simplest case, you can build a "super model" with a single method call:
```
# morph is a valid Morph instance as specified above
morph.prepare_supermodel()
```

The `fill_src` flag can be set to true to add a gap-filling step prior to translating features, setting the src_model of the morph to a gapfilled version of the one originally supplied.
This function composes a few others into a single command:
```
        self.translate_features()
        self.reconstruct_genome()
        self.label_reactions()
        self.build_supermodel()
```
These steps can be run individually as well. The translation and reconstruction steps set `morph.trans_model` and `morph.recon_model` respectively, and can be run in any order. Alternatively, one could choose to assign these values instead of running the functions, but the use-cases for this are left to the discretion of the user.

The labeling step and building of the super model must be performed after the translation and reconstruction steps.

### Translating Media (optional):

A common use case occurs when the source organism and target organism don't grow on identical media. In this case, one can specify the source media upon Morph construction as described above, 
and prior to processing reactions, use the `translate_media` method. The `translate_media` method takes a Media object as an argument and if necessary, gapfills the model so it can grow on this new media, changing `morph.media` to the supplied argument.


### Processing Reactions (Iterative FBA Based Removal)

A next step would be to try iterative FBA Based reaction removal. This is done using the process_reactions method. In the simplest case, this can be done using a call with no arguments:

```
    # morph has a model attribute as set by build_supermodel/prepare_supermodel
    morph.process_reactions()
```

In the default case, all 'gene-no-match' reactions and 'no-gene' reactions are tested for removability, in order of increasing reaction likelihood. A reaction is considered removable if
an FBA solution for a model in which the reaction has been removed exhibits growth (f > 0.0). Any reactions that are essential to growth on morph.media are kept, regardless of source or genetic basis.

This process van be customized with a few keyword arguments:
self, rxn_list=None, name='', process_count=0, get_count=False, iterative_models=True, growth_condition=GrowthConditions.SimpleCondition()
- `rxn_list`: the reactions to test for removal in order (a list of reaction IDs). Default case is gene-no-match sorted by increasing likelihood, then no-gene by increasing likelihood.
- `growth_condition`: a GrowthCondition object, which has an evaluate function that returns true or false, indicating when to keep or remove a reaction. Default is SimpleCondition, which takes a morph/model, runs FBA on it, and returns True if flux > 0.0, False otherwise. Other options are available, including one for ensuring growth on multiple media. Subclasses can be used for desired process effects, see the GrowthConditions module for details.
- `name`: (non-algorithmic) The prefix name for all intermediate models. Default is 'MM'. Can be changed if you plan on using/inspecting intermiediate models.
- `process_count`: (non-algorithmic) an increasing count of models through the removal process, default is 0. e.g, the first intermediate model will be 'MM-0', then 'MM-1', etc. allowing inspection of intermediates. If you've already removed n reactions, it can be useful to set process_count to n.
- `get_count`: (non-algorithmic) A Boolean flag indicating whether the process_count should be returned with the morph (as a tuple). Used when not processing all reactions at once. Default is False.
- `iterative_models`: (non-algorithmic) A Boolean flag indicating whether or not to use a new model for each removal (e.g. count MM-0, MM-1, etc. or just use a single model overwritten each time). 


This function will run for a while, and will give updates on progress to stdout as it is able to remove/keep reactions. When it complete, the morphing process is finished! The final model is the one referenced by `morph.model`

## Questions, Concerns, Feature Requests, etc.

Feel free to contact us via GitHub here, via issues, pull requests, messaging, etc.




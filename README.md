A tool for morphing a metabolic model from one genome to a closely related one, without ever disrupting growth

# Morphing a Model for Organism 'A' to Organism 'B'

Takes a Model/Genome for a 'source' organism and 'morphs' it to the 'target' organism. It keeps reactions from the 'source' model for which there is at least one matching gene feature, and attempts to remove those in the source model that do not have BLAST hits in the second genome. It also adds reactions annotated within the target genome that are unique to the target organism. The model is never broken in the process, such that the Biomass reaction and any other specified reaction [insert how you flag this option] must always carry flux. If a reaction can't be removed without breaking the model, it is removed and gap-filled using Prob-Anno gapfilling [other options?]

## Installation
1. Clone this repository: `git clone https://github.com/kingb12/model-morphing.git`
2. Add `path-to-the-cloned-repo/model-morphing/lib/` to your `PYTHONPATH` environment variable
3. From the model-morphing directory, run in shell: `python lib/httplib2/setup.py install`
4. run in shell `sudo pip install python-firebase`

## Running from the Interpretter
1. Navigate to `model-morphing/lib/`
2. run in shell `ipython`
3. `from Helpers import *`

## Tutorial
A Walkthrough of some of the main processes used in the Client module and how to use it's functions
### Instantiating a Morph
All functions in the Client module that participate in the morphing process take a morph as a parameter. A Morph manages the state of the model as it is translated from one organism's genome to another, as well as other useful information for use in our algorithm. [More about Morphs can be read here](Morph.md)
These are the minimal required arguments to instantiate a morph:
- genome - an object_id for a KBase genome object (the target genome)
- genomews - an workspace_id for a KBase genome object (the target genome)
- src_model - an object_id for a KBase model object (the source model)
- src_modelws - an workspace_id for a KBase model object (the source model)
- media - an object_id for a KBase media object 
- mediaws - an workspace_id for a KBase media object
- protcomp - an object_id for a KBase protein comparison object (source and target genomes)
- protcompws - an workspace_id for a KBase protein comparison object (source and target genomes)
- probanno - an object_id for a KBase RxnProbs object 
- probannows - an workspace_id for a KBase RxnProbs object
##### ws_id attribute
 useful optional attribute to set is `ws_id`. Setting `ws_id` allows you to control what workspace the morph bases itself in.
Default behavior is to create a new workspace for the user upon instantiation, but this doesn't work with KBase narrative. To view morph data in a KBase narrative, the simplest way is to create a new narrative in KBase and supply it's workspace_id to the Morph constructor.
First, a morph can be instantiated with either a dictionary of arguments, a set of keyword arguments, or some combination of the two. Here's an example function in the Helpers module which constructs and returns a morph:
```python
def make_morph(ws_id=None):
    args = dict()
    args['genome'] = '3'    # genome object_id
    args['src_model'] = '19'    # source model object_id
    args['probanno'] = '15'    # probanno object_id
    args['protcomp'] = '6'    # protcomp object_id
    args['genomews'] = '9145'    # genome workspace_id
    args['src_modelws'] = '9145'    # source model workspace_id
    args['probannows'] = '9145'    # probanno workspace_id
    args['protcompws'] = '9145'    # protcomp workspace_id
    args['mediaws'] = '9145'    # media workspace_id
    args['media'] = '18'    # media object_id
    args['ws_id'] = ws_id    # morph workspace_id (used for all new created objects)
    return Morph(args)
```
This function instantiates a Morph with the arguments included in the args dictionary, but keywords could also be used:
```python 
morph = Morph(genome=3, src_model=19, probanno=15, protcomp=6, media=18
    genomews=9145, src_modelws=9145, probannows=9145, protcompws=9145, mediaws=9145)
```

### Getting a SuperModel
Getting a super model for a givenn argument source model and target genome can be done with just one Client function. To start from scratch, first instantiate a Morph with the minimal required parameters. This can then be follwed by a call to "prepare_supermodel":
```python
# For simplicity, 'args' is understood to be the dictionary as it is shown above in the make_morph() function
>>>from Helpers import *
>>>morph = Morph(args)
>>>morph = Client.prepare_supermodel(morph)
>>>morph
    probhash: [u'rxn10122', u'rxn10123', u'rxn10120', u'rxn10121', u'rxn10126', u'rxn10127', u'rxn10124', u'rxn10125', u'rxn14133', u'rxn10824', u'rxn04659', u'rxn04658', u'rxn28302', u'rxn09265', u'rxn01506', u'rxn04657', u'rxn04656', u'rxn10324', u'rxn01951', u'rxn01406', u'rxn01025', u'rxn01021', u'rxn01022', u'rxn01029', u'rxn01028', u'rxn01073', u'rxn00377', u'rxn00375', u'rxn00371', u'rxn01205', u'rxn01204', u'rxn01753', u'rxn01200', u'rxn01751', u'rxn00379', u'rxn05040', u'rxn20068', u'rxn05687', u'rxn05681', u'rxn05680', u'rxn05683', u'rxn05682', u'rxn01111', u'rxn01115', u'rxn01114', u'rxn01117', u'rxn01116', u'rxn00578', u'rxn00575', u'rxn00577', u'rxn00571', u'rxn10258', u'rxn10259', u'rxn11650', u'rxn10253', u'rxn10254', u'rxn10255', u'rxn10256', u'rxn10257', u'rxn23469', u'rxn01682', u'rxn01685', u'rxn02937', u'rxn09240', u'rxn00790', u'rxn00791', u'rxn09113', u'rxn08199', u'rxn09114', u'rxn00799', u'rxn00400', u'rxn00405', u'rxn00159', u'rxn00408', u'rxn00157', u'rxn00152', u'rxn00151', u'rxn01977', u'rxn01975', u'rxn01974', u'rxn01973', u'rxn01972', u'rxn23799', u'rxn23798', u'rxn23796', u'rxn15173', u'rxn24928', u'rxn00088', u'rxn11938', u'rxn00085', u'rxn00086', u'rxn00082', u'rxn00083', u'rxn01132', u'rxn12822', u'rxn06802', u'rxn22688', u'rxn22689', u'rxn03838', u'rxn03839'] ... (more)
    protcompws: 9145
    removed_ids: None
    trans_model: 4
    objects: ['source_model', 'recon_model', 'probanno', 'trans_model']
    ws_id: 11498
    probanno: 15
    info: ['source_model', 'recon_model', 'probanno', 'trans_model']
    probannows: 9145
    src_modelws: 11498
    src_model: 3
    media: 18
    mediaws: 9145
    genomews: 9145
    essential_ids: None
    genome: 3 
    rxn_labels: ['gene-no-match', 'gene-match', 'recon', 'no-gene']
    model: 6
    ws_name: MMws264
    recon_model: 5
    protcomp: 6
>>>
```
Here, `morph.model` refers to the `object_id` of the super model in workspace `morph.ws_id`.
####Skipping gapfilling in prepare_supermodel
Adding the keyword argument `fill_src` and setting it to `False` skips a preliminary gap-filling step that runs probanno based gap-filling on the source model prior to reaction labelling. This can be done like so:

```python
>>>morph = Client.prepare_supermodel(morph, fill_src=False)
>>>
```
The prepare_supermodel function wraps several functions that can be called individually if necessary:
```python
def prepare_supermodel(morph, fill_src=True):
    morph = copy.deepcopy(morph)
    if fill_src:
        morph = fill_src_to_media(morph)
    morph = translate_features(morph)
    morph = reconstruct_genome(morph)
    morph = label_reactions(morph)
    morph = build_supermodel(morph)
    return morph
```
### Processing Reactions (Iterative FBA Based Removal)
A next step would be to try iterative FBA Based reaction removal. This is done using the process_reactions module in the CLient. A few scenarios to consider:
####Standard Removal Procedure
Without passing an extra argument using key word `rxn_list`, the `process_reactions` method by default removes all reactions labelled 'gene-no-match' and all reactions labelled 'no-gene'. The morph models can then be compared The commands to do so are shown below. 
NOTE: The `process_reactions` function requires in this scenario that `morph.rxn_labels` have valid states for key's `gene-no-match` and `no-gene`, as well as a valid `morph.model`. Details can be found in the docs for the Client module. The short answer is, the functions `label_reactions` and `build_supermodel` must proceed `process_reactions` (these two are a part of `prepare_supermodel`).
```python
>>>morph_processed = Client.process_reactions(morph)
# Note the [] around the arguments, means the morphs are passed as a list
>>>comparison_data = Client.compare_morphmodels([morph, morph_processed]) 
>>>
```
`comparison_meta` now holds comparison information between the models in the two morphs as well as comparison data for each reaction in each model. `morph_processed` is now a model such that every reaction that is in the model and doesn't have a gene association for the new genome is essential.
#### Processing reactions by list
If you only want to process a specific set of reactions, you can pass them as a list to the key word argument `rxn_list` in `process_reactions`. Processing reactions requires the list in a slightly complex form:
`reaction_id -> (model_index, probability)` *will be fixed in the future to be simpler*. Thanksfully, if you are trying to remove reactions of a particular label, a helper function exists to make this process easier:
##### Helpful Function `removal_list`
`removal_list` takes a dictionary of the form: `(K) reaction_id -> (V) (model_index, probability)` and turns it into a list of reactions to remove of the form required by process reactions, sorted by probability low to high. So for example, of you want to process only the `gene-no-match` reactions, this can be done like so:
```python
>>>morph_gnm_processed = Client.process_reactions(morph, rxn_list=removal_list(morph.rxn_labels['gene-no-match']))
```
`removal_list` even has a key word parameter `list_range` for getting only a range within the removal_list. `list_range accepts a tuple with the first index to get and one beyond the last index you want (matches range() behavior). For example, if you wanted to remove the first 10 'gene-no-match' reactions, then the next 5, you could do so like this (removal list put on it's own line for readability):
```python
>>>rmv1to10 = removal_list(morph.rxn_labels['gene-no-match'], list_range=(0,10)
>>>morph_gnm1to10_processed = Client.process_reactions(morph, rxn_list=rmv1to10)
>>>rmv10to15 = removal_list(morph.rxn_labels['gene-no-match'], list_range=(10,15)
>>>morph_gnm1to15_processed = Client.process_reactions(morph_gnm1to10_processed, rxn_list=rmv10to15)
>>>
```
#### Other parameters to `process_reactions`
Other key word parameters can be set in process_reacitons, and these are generally used together to control output of the process_reactions in a KBase workspace so they can be more easily analyzed in the narrative. A guide for this will be written in the future, but this is generally an advanced topic and can be explored in the Client module docs.
### Simple Probanno Morph
A simple function can be used for creating a morph with a shorter algorithm. `simple_probanno_morph` takes an initialized morph, prepares a supermodel, then removes all 'gene-no-match' and 'no-gene' reactions, and finally fills in the model with probanno based gapfilling. It's very straight forward, here is an example of it in use, follwed by a `compare_morphmodels` call comparing it with the processed morph and the supermodel:
```python
>>>m = make_morph()
>>>prob_morph = Client.simple_probanno_morph(m)
>>>comparison_data = Client.compare_morphmodels([morph, prob_morph, processed_morph])
```
### Other useful functions:
Here are some other useful functions that I will tutorial more deeply soon. As best as possible, I formatted the input and output forms of these to be similar or identical to those above. Some are fairly self-explanatory, and most are well documented in the Client module:
- `find_alternative` *not currently functioning as desired*
- `probanno_fill` accepts a morph and returns a morph where the model has been filled using probanno based gapfilling
- `remove_reaction` accepts a `reaction_id` and a morph and removes the given reaction from the morph
- `remove_reactions` accepts a list of `reaction_id` and a morph and removes the all of the listed reactions from the morph
- `get_reaction` accepts a `reaction_id` and returns a dictionary of information about the reaction, including it's equation and more
- `get_reactions` like `get_reaction`, bit accepts a list of reactions and returns a list of dictionaries with reaction info.
- `get_morph_rxns` accepts a morph (and optionally a label) and returns all the reactions in the morph's model ([] if no model)


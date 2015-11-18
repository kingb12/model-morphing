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
A useful optional attribute to set is ws_id. Setting ws_id allows you to control what workspace the morph bases itself in.
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
```
Here, `morph.model` refers to the object_id of the super model in workspace `morph.ws_id`.
####Skipping gapfilling in prepare_supermodel
Adding the keyword argument `fill_src` and setting it to `False` skips a preliminary gap-filling step that runs probanno based gap-filling on the source model prior to reaction labelling. This can be done like so:

```python
>>>morph = Client.prepare_supermodel(morph, fill_src=False)
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
A next step

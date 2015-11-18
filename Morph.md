# Understanding the Morph Class

## General Idea
Morph is an instantiable class used as a state manager for the morphing algorithm. It's primary component is a reference
to a KBase FBA Model, but it includes attributes that describe the current state of the morphing process. The goal of the 
Morph class is to make managing the state of the process simple. There is no need to keep track of reaction labels, or
removal/essentiality data in other variables or contexts. This helps prevent inconsistent states, such as when a KBase
object_id anf workspace_id don't correctly pair to an ObjectIdentity. It also makes function calls simple: every step in the 
model morphing process accepts a morph as a parameter and returns a new morph as a result (allowing the user to manage references
to stages in the process).

## Usage
### Instantiating a Morph
A Morph is instantiated by declaring it's properties. The constructor is flexible, so this can be done with a dictionary of
attribute names to their values, keyword arguments, or a suitable combination of the two. These are the minimal required arguments
to create a morph:
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
Default behavior is to create a new workspace upon instantiation.

Here's an example function in the Helpers module which makes and constructs a morph:
```python
def make_morph(ws_id=None):
    args = dict()
    args['genome'] = '3'
    args['src_model'] = '19'
    args['probanno'] = '15'
    args['protcomp'] = '6'
    args['genomews'] = '9145'
    args['src_modelws'] = '9145'
    args['probannows'] = '9145'
    args['protcompws'] = '9145'
    args['mediaws'] = '9145'
    args['media'] = '18'
    args['ws_id'] = ws_id
    return Morph(args)
```
This function instantiates a Morph with the arguments included in the args dictionary, but keywords could also be used:
```python 
morph = Morph(genome=3, src_model=19, probanno=15, protcomp=6, media=18
    genomews=9145, src_modelws=9145, probannows=9145, protcompws=9145, mediaws=9145)
```
Helper functions in the Helpers module can simplify this step by returning new morphs with some default attribute values
or morphs representing a state somewhere aling the process

### Using a Morph in Client module
Using a Morph in the client module was designed to be as simple as possible. Most Client functions accept a morph and return a morph.
These are reference safe so that if the return value is assigned to a different variable name than the morph argument, they will have
distinct state values. In some cases, a morph is passed to a function that might intuitively accept a model. These functions make KBase
calls that modify morph.model, among other attributes. morph.model can be thought of as the primary data component of the morph; it 
references the KBase model as we change it and is the primary output of the overall algorithm. The Morph class has additional attributes
to make state management easy, to provide information for analysis, and to reduce the number of required arguments to each Client function
to a reasonable number.
Here are a few example calls to Client methods:
```python
# A function in Helpers module that makes a default morph for acetivorans -> barkeri
morph = Helpers,make_morph()
morph = Client.prepare_supermodel(morph)
fully_processed_morph = Client.process_reactions(morph)
gnm_processed_morph = Client.process_reactions(morph, 
    rxn_list=removal_list(morph.rxn_labels['gene-no-match']))
```
The above code takes advantage of the fact that argument morphs keep their original state, and a function returns a new Morph with a resulting state. First, a morph is instantiated using a Helper function. Then, the Client function prepare_supermodel runs through the first steps in the algorith, creating a morph with labelled reactions and a super_model regerence in morph.model. After this, we take advantage of state preservation to diverge into two new morphs from just the one argument, one called 'fully_processed_morph' which has run processing on all gene-no-match and no-gene reactions, and one called 'gnm_processed_morph' which only processed the 'gene-no-match' reactions. In the end, we can reference three distinct morph states: morph (pre-processing), gnm_processed_morph (mid-processing), and fully_processed_morph (post-processing).

### Modifying and Reading Attributes
####Modifying
Morph Attributes are modified through direct manipulation: just treat each attribute like you would any other variable
in the namespace. *Overriding of attribute setting methods in the Morph class will prevent logcally inconsistent states
for each attribute (i.e. an object_id like morph.genome can't be set to a Double or Boolean)

*NOT YET IMPLEMENTED

Here's an example showing how to change the model attribute (done frequently in Client module):
```python
morph.model = 78
```
#### Reading
Attributes of a morph are also read directly. If an attribute is a structure of some sort, like a dictionary, it's components
can also be read directly, and reader methods can be used to look at the data in specific ways. For instance, here is how
the 5th 'no-gene' reaction can be accessed:
```python
#Note the use of the keys() method to get an ordered list, getting element at index 4 (fifth in list)
fifth_nogene_reaction = morph.rxn_labels['no-gene'].keys()[4]
```
#### Aliasing Warning!
Note: Morph is mutable data structure, and the namespace handles only references to each one. This can create some confusion for
those inexperienced with programming. Here's an example of where aliasing causes an issue that might not be expected:
```python
morph.model = 18
print morph.model
# creates a possibly inadvertent variable name alias
m2 = morph
m2.model = 25
print morph.model
print m2.model
```
One not familiar with aliasing might expect this code to output:
```
18
18
25
```
When it really outputs:
```
18
25
25
```
This is because m2 is just an alias for morph; both variable names refer to the same instance in memory. To correct for this in
your own code, use the copy module (This is taken care of in all Client module functions):
```python
import copy

morph.model = 18
print morph.model
# creates a new instance in memory
m2 = copy.deepcopy(morph)
m2.model = 25
print morph.model
print m2.model
```


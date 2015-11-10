A tool for morphing a metabolic model from one genome to a closely related one, without ever disrupting growth

## Morphing a Model for Organism 'A' to Organism 'B'

Takes a Model/Genome for a 'source' organism and 'morphs' it to the 'target' organism. It keeps reactions from the 'source' model for which there is at least one matching gene feature, and attempts to remove those in the source model that do not have BLAST hits in the second genome. It also adds reactions annotated within the target genome that are unique to the target organism. The model is never broken in the process, such that the Biomass reaction and any other specified reaction [insert how you flag this option] must always carry flux. If a reaction can't be removed without breaking the model, it is removed and gap-filled using Prob-Anno gapfilling [other options?]

### Installation
1. Clone this repository: `git clone https://github.com/kingb12/model-morphing.git`
2. Add path-to-the-cloned-repo/model-morphing/lib/ to your PYTHONPATH environment variable
3. From the model-morphing directory, run in shell: `python lib/httplib2/setup.py install`
4. run in shell `sudo pip install python-firebase`

### Running from the Interpretter
1. Navigate to model-morphing/lib/
2. run in shell `python`
3. `from Helpers import *`


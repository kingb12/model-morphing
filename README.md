A tool for morphing a metabolic model from one genome to a closely related one, without ever disrupting growth

## Morphing a Model for Organism 'A' to Organism 'B'

Takes a Model/Genome for a 'source' organism and 'morphs' it to the 'target' organism. It keeps reactions from the 'source' model for which there is at least one matching gene feature, and attempts to remove those in the source model that do not have BLAST hits in the second genome. It also adds reactions annotated within the target genome that are unique to the target organism. The model is never broken in the process, such that the Biomass reaction and any other specified reaction [insert how you flag this option] must always carry flux. If a reaction can't be removed without breaking the model, it is removed and gap-filled using Prob-Anno gapfilling [other options?]

### Inputs:
	- A functioning Model for Organism A
	- A Genome for Organism A
	- A Genome for Organism B

### Output:
	- A functioning Model for Organism 


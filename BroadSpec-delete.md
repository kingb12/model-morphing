#Broad Specification:

##Parameters and Outputs:

####Inputs: 
- Model
  - Can be gap filled, or just a reconstruction, but the more complete input, the better the output
- A Proteome Comparison
  - Probably use the method and structure already in KBase, or else re-implement a similar idea elsewhere
- ProbAnno for New Genome
  - The new Genome

####Output:
- Model for New Genome

####Requires:
- Proteome Comparison compares Genome on the Model to a closely related genome

##Broad Algorithm:
1. Accept a Model, Proteome Comparison, and ProbAnno Object for the ‘other’ genome
  1. (Try to incorporate these into a ‘toolbox’ but not the actual model-morphing function)
2. Accept these Options: (?)
  1. ‘Label’ Reactions (Proteins?) in functional model as ‘gene-matched’, ‘gene-unmatched’, and ‘no-gene'
  2. ‘gene-matched’: The reaction has an associated gene annotation in both the model AND the other genome (determined from proteome comparison)
  3. ‘gene-unmatched’: The reaction has an associated gene annotation in the model but NOT in the other genome
  4. ‘no-gene’: the reaction has no gene association in the model 
  5. Should I also check here for an annotation in the other genome?
3. Build a reconstruction from the ‘other’ genome
  1. Save these reactions for comparison with ones already in the model
  2. ‘reconstruction’: reactions that come from the annotations directly from the ‘other’ genome
4. Build a ‘super-model'
  1. ‘super-model’: a model composed of all (now labelled) reactions in the old model and the addition of all reactions in the new model’s reconstruction
5. Settle conflicts between reconstruction reactions and those form the old model
  1. i.e. ‘reconstruction’ == ‘gene-matched’ > ‘no-gene’ >> ‘gene-unmatched'
6.  Try to remove all ‘gene-unmatched’ reactions
  1. One at a time
  2. Don’t let the model break
  3. If it breaks, try to replace it with a high likelihood reaction from ProbAnno
  4. This will be a long cycle of remove a reaction, run fba, and replace the reaction with a gap-filling algorithm if the model breaks, keeping the ‘gene-unmatched’ reaction only if no satisfying gap-fill reaction can be found
7. Try removing ‘no-gene’ reactions
  1. This should be simpler, ProbAnno rxn for ‘other’ genome > gap-fill for old model
8. Output the model, and a stat sheet explaining the state of the model’s reactions, and maybe a calculated confidence score for each reaction

##Environment Needs
- Access to these functions:
  - Part of the function:
    - running FBA
  - Part of the workflow:
    - reconstructing a model from a genome
    - run ProbAnno
    - run proteome comparison
 - Ability to examine the model at each step
  - For debugging AND curation
  - If you never break the model, you should ALWAYS be able to export it and look at it

##Documentation:
- Specification:
  - A description of all functions of the tool, their requirements, inputs, outputs, critical dependencies, etc.
  - This is for users
- Implementation Guide:
  - A guide to understanding the final implementation of the tool when I’m done working on it
  - This is for future developers (and myself as I go along)
  - This helps keep the tool up to date and makes it finish-able in the case it is not finished at the end of August. Helps explain where new features can be added, etc.
- Usage Guide
  - A detailed explanation of how to use the tool, applications of the results, etc. 
  - Include a few Sample Workflows

###Goals, Ideals, Intuitions, Avoids
- Ideal: Function can be ‘used’ with as little duplication as possible on KBase,  Matlab, and CobraPy
- Avoid: Don’t write the code in Matlab
- Ideal: written in Python so that Matt or others can easily make changes to the code as it aged
- Maybe talk to Buruk about other ways of making sure this works
- Goal: No one should need to understand the KBase architecture, ProbAnno, Proteome Comparison, etc. to understand the function code
- Goal: Describe all parts of the KBase code that the implementer needs to know in ‘Implementation’ Docs

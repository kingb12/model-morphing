# FVA Analysis
 automated analysis of the sparsity of our network, comparing our morph to the supermodel. Further analysis will included the Source Model, other morphs, etc.
### Blockage Data for each model

| Model | Blocked Reactions | Total Reactions | % Blocked | % blocked with genes | 
| --- | --- | --- | --- | --- | 
| Super Model | 329 | 905 | 36.35 | 74.47 | 
| Morphed Model | 340 | 632 | 53.80 | 90.00 | 

### Change in reaction blocking

| Blocked In:  | Morph Only | Both | Supermodel Only | 
| --- | --- | --- | --- | 
| M. jannschii | 78 | 262 | 67 | 

the sets that inform these reactions are available, e.g. what are the qualities of the reactions blocked in only the Morphed Model?
##### Blocked Only in Morphed Model
| Number of Reactions | Percent Gene Based | Percent Gene-Match | 
| --- | --- | --- | 
| 78 | 87.18 | 32.05 | 

 Full Table: 

| Reaction | genes in morph? | # of genes | morph labels | class | sublclass | subsystem | 
| --- | --- | --- | --- | --- | --- | --- | 
| rxn10481_c0 | False | 0 | gene-no-match | None, Clustering-based subsystems,  | , None,  | None, CBSS-196620.1.peg.2477,  | 
| rxn04046_c0 | True | 2 | gene-match, common,  | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn05651_c0 | False | 0 | recon | None, Sulfur Metabolism, Amino Acids and Derivatives,  | Lysine, threonine, methionine, and cysteine, Inorganic sulfur assimilation, None,  | None, Cysteine_Biosynthesis, Inorganic_Sulfur_Assimilation,  | 
| rxn04048_c0 | True | 2 | gene-match, common,  | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn03086_c0 | False | 0 | gene-no-match | Amino Acids and Derivatives, None,  | Lysine, threonine, methionine, and cysteine, None,  | None, Lysine_Biosynthesis_DAP_Pathway,  | 
| rxn03437_c0 | True | 1 | recon | None, Experimental Subsystems,  | , None,  | None, Iron-sulfur_experimental,  | 
| rxn00646_c0 | True | 1 | gene-match | None, Stress Response,  | None, Oxidative stress,  | Glutathione_synthesis, None,  | 
| rxn00881_c0 | True | 1 | recon | Carbohydrates, None,  | None, Sugar alcohols,  | None, Di-Inositol-Phosphate_biosynthesis,  | 
| ADTHOR_c0 | True | 1 | gene-match | None | None | None | 
| rxn01678_c0 | True | 2 | gene-match, common,  | None, Nucleosides and Nucleotides,  | None, Purines,  | None, Purine_conversions,  | 
| rxn01485_c0 | True | 1 | recon | None, Cell Wall and Capsule,  | , None,  | UDP-N-acetylmuramate_from_Fructose-6-phosphate_Biosynthesis, None,  | 
| rxn03075_c0 | False | 0 | gene-no-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | Thiamine and thiamine pyrophosphate, None,  | Thiamin_biosynthesis, None,  | 
| rxn04047_c0 | True | 1 | gene-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn04050_c0 | True | 1 | gene-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn03491_c0 | False | 0 | gene-no-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Coenzyme_B12_biosynthesis,  | 
| rxn12845_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| 2OPs_c0 | True | 1 | gene-match | None | None | None | 
| rxn05287_c0 | True | 1 | gene-match | None, Fatty Acids, Lipids, and Isoprenoids,  | Isoprenoids, None,  | Polyprenyl_Diphosphate_Biosynthesis, None,  | 
| rxn07588_c0 | False | 0 | gene-no-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Coenzyme_B12_biosynthesis,  | 
| rxn05232_c0 | True | 1 | recon | None, Nucleosides and Nucleotides,  | , None,  | Ribonucleotide_reduction, None,  | 
| rxn06874-c0-_c0 | True | 4 | gene-match | None | None | None | 
| rxn01353_c0 | True | 2 | gene-match, common,  | None, Nucleosides and Nucleotides,  | None, Purines,  | None, Purine_conversions,  | 
| ADTHs_c0 | True | 1 | gene-match | None | None | None | 
| rxn12846_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| rxn00707_c0 | True | 2 | gene-match, common,  | None, Experimental Subsystems,  | , None,  | None, 271-Bsub,  | 
| rxn01316_c0 | True | 1 | recon | None, Virulence,  | , None,  | None, Streptococcus_agalactiae_virulome,  | 
| rxn05109-c0-_c0 | True | 1 | gene-match | None | None | None | 
| rxn03892_c0 | True | 1 | gene-match | None, Fatty Acids, Lipids, and Isoprenoids,  | Isoprenoids, None,  | Polyprenyl_Diphosphate_Biosynthesis, None,  | 
| rxn05201_c0 | True | 3 | gene-match, common,  | None, Nucleosides and Nucleotides,  | , None,  | Purine_Utilization, None,  | 
| rxn07589_c0 | False | 0 | gene-no-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Coenzyme_B12_biosynthesis,  | 
| rxn12844_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| rxn12646_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| rxn04042_c0 | True | 1 | gene-match | ,  | ,  | ,  | 
| rxn05054_c0 | False | 0 | gene-no-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn08309_c0 | True | 1 | recon | None, Fatty Acids, Lipids, and Isoprenoids,  | None, Phospholipids,  | Glycerolipid_and_Glycerophospholipid_Metabolism_in_Bacteria, None,  | 
| DKFPs1_c0 | True | 1 | gene-match | None | None | None | 
| DKFPs2_c0 | True | 1 | gene-match | None | None | None | 
| rxn07587_c0 | True | 1 | gene-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn00692_c0 | True | 1 | recon | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | Folate and pterines, None,  | None, 5-FCL-like_protein,  | 
| rxn02775_c0 | True | 2 | gene-match, common,  | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn03512_c0 | True | 2 | gene-match, common,  | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn00274_c0 | True | 1 | gene-match | Amino Acids and Derivatives, None,  | Alanine, serine, and glycine, None,  | None, Glycine_and_Serine_Utilization,  | 
| rxn11544_c0 | True | 3 | gene-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn03492_c0 | True | 1 | gene-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn00835_c0 | True | 2 | gene-match, common,  | None, Nucleosides and Nucleotides,  | None, Purines,  | None, Purine_conversions,  | 
| rxn12635_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| rxn00293_c0 | True | 1 | recon | None, Experimental Subsystems, Cell Wall and Capsule,  | , None,  | None, UDP-N-acetylmuramate_from_Fructose-6-phosphate_Biosynthesis, Peptidoglycan_Biosynthesis_experimental,  | 
| rxn00988_c0 | True | 1 | gene-match | None, Amino Acids and Derivatives, Fatty Acids, Lipids, and Isoprenoids,  | , None, Fatty acids, Branched-chain amino acids,  | HMG_CoA_Synthesis, None, Polyhydroxybutyrate_metabolism, Fatty_acid_degradation_regulons,  | 
| rxn01257_c0 | True | 1 | recon | None, Cofactors, Vitamins, Prosthetic Groups, Pigments, Amino Acids and Derivatives,  | Folate and pterines, None, Aromatic amino acids and derivatives,  | None, Folate_Biosynthesis, Chorismate:_Intermediate_for_synthesis_of_PAPA_antibiotics,_PABA,_anthranilate,_3-hydroxyanthranilate_and_more.,  | 
| rxn12641_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| rxn12847_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| rxn07586_c0 | False | 0 | gene-no-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn00832_c0 | True | 2 | gene-match, common,  | None, Nucleosides and Nucleotides, Experimental Subsystems,  | , None, Purines,  | ZZ_gjo_scratch, None, De_Novo_Purine_Biosynthesis,  | 
| rxn12639_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| F6PG3Pl_c0 | True | 1 | gene-match | None | None | None | 
| rxn00669_c0 | True | 1 | gene-match | Protein Metabolism, Carbohydrates, None,  | Protein processing and modification, None, Central carbohydrate metabolism, Organic acids,  | None, Pyruvate_metabolism_II:_acetyl-CoA,_acetogenesis_from_pyruvate, Propionate-CoA_to_Succinate_Module, Protein_Acetylation_and_Deacetylation_in_Bacteria,  | 
| rxn03127_c0 | True | 8 | gene-match, common,  | Carbohydrates, None,  | None, One-carbon Metabolism,  | None, Methanogenesis,  | 
| rxn00650_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| rxn00515_c0 | True | 2 | gene-match, common,  | None, Nucleosides and Nucleotides,  | None, Purines,  | None, Purine_conversions,  | 
| rxn03513_c0 | True | 1 | gene-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn12642_c0 | True | 3 | gene-match, common,  | Protein Metabolism, None, Carbohydrates,  | Protein degradation, Central carbohydrate metabolism, None, Protein biosynthesis,  | Dehydrogenase_complexes, None, Translation_factors_bacterial, Protein_degradation,  | 
| rxn03638_c0 | True | 1 | recon | None, Experimental Subsystems, Cell Wall and Capsule,  | , None,  | UDP-N-acetylmuramate_from_Fructose-6-phosphate_Biosynthesis, None, Peptidoglycan_Biosynthesis_experimental,  | 
| rxn01069_c0 | True | 1 | recon | Amino Acids and Derivatives, None,  | Lysine, threonine, methionine, and cysteine, None,  | Threonine_and_Homoserine_Biosynthesis, None,  | 
| rxn04052_c0 | True | 2 | gene-match, common,  | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn05179_c0 | True | 4 | recon | None, Virulence, Experimental Subsystems,  | , None,  | Transporters_In_Models, None, C_jejuni_colonization_of_chick_caeca,  | 
| rxn03514_c0 | True | 2 | gene-match, common,  | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn03030_c0 | True | 1 | gene-match | None, Experimental Subsystems,  | , None,  | None, 271-Bsub,  | 
| rxn02774_c0 | True | 4 | gene-match, common,  | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | Experimental_tye, None,  | 
| rxn03136_c0 | True | 1 | recon | None, Nucleosides and Nucleotides,  | None, Purines,  | None, Purine_conversions,  | 
| rxn04045-c0-_c0 | True | 1 | gene-match | None | None | None | 
| rxn02113_c0 | True | 1 | gene-match | Carbohydrates, None,  | None, Organic acids,  | None, Alpha-acetolactate_operon,  | 
| rxn00883_c0 | True | 1 | recon | Carbohydrates, None,  | None, Sugar alcohols,  | None, Di-Inositol-Phosphate_biosynthesis,  | 
| rxn00351_c0 | False | 0 | gene-no-match | None, Stress Response,  | None, Oxidative stress,  | Glutathione_synthesis, None,  | 
| rxn03147_c0 | True | 2 | gene-match, common,  | None, Clustering-based subsystems,  | , None,  | CBSS-251221.1.peg.1863, None,  | 
| rxn06979_c0 | True | 3 | gene-match | None, Cofactors, Vitamins, Prosthetic Groups, Pigments,  | None, Tetrapyrroles,  | None, Cobalamin_synthesis,  | 
| rxn00914_c0 | True | 2 | gene-match, common,  | None, Nucleosides and Nucleotides,  | None, Purines,  | None, Purine_conversions,  | 
| rxn01300_c0 | True | 2 | gene-match, common,  | None, Experimental Subsystems,  | , None,  | None, Methionine_Biosynthesis_short,  | 
| rxn01575_c0 | True | 2 | gene-match, common,  | Carbohydrates, None,  | None, Central carbohydrate metabolism,  | Dehydrogenase_complexes, None,  | 


### Blockage Data for each model

| Model | Blocked Reactions | Total Reactions | % Blocked | % blocked with genes | 
| --- | --- | --- | --- | --- | 
| Super Model | 329 | 905 | 36.35 | 74.47 | 
| Morphed Model | 340 | 632 | 53.80 | 90.00 | 

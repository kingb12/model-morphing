# Model-Morphing: Analyses Thus Far

### Morphing Organisms Covered:
* M. maripaludis to M. jannaschii
* M. maripaludis to M. stadtmanae
* M. maripaludis to M. barkeri
* M. maripaludis to M. maripaludis

### Analyses for All Morphs:
3-way Venn comparison of Reaction Sources


### M. maripaludis to M. jannaschii Data

| Model  | Genes | Metabolites  | Reactions |
| ------------- | ------------- | ------------- | ------------- |
| M. maripaludis Source          | 340     | 710           |         636 | 
| M. jannaschii Morph            | 567     | 841           | 631         | 
| M. jannaschii Reconstruction   | 285     | 708           | 665         | 

#### Reaction Overlap Analysis
![alt text](https://raw.githubusercontent.com/kingb12/model-morphing/master/analysisimages/marijannarxnsvenn.png "marijannarxnsvenn.png")

#### Gene Overlap Analysis
![alt text](https://raw.githubusercontent.com/kingb12/model-morphing/master/analysisimages/marijannagenesvenn.png "marijannagenesvenn.png")

### M. maripaludis to M. stadtmanae Data

| Model  | Genes | Metabolites  | Reactions |
| ------------- | ------------- | ------------- | ------------- |
| M. maripaludis Source          | 222 | 710 | 636 |
| M. stadtmanae Morph            | 464 | 930 | 696 |
| M. stadtmanae Reconstruction   | 296 | 766 | 691 |

#### Reaction Overlap Analysis
![alt text](https://raw.githubusercontent.com/kingb12/model-morphing/master/analysisimages/maristadtrxnsvenn.png "maristadtrxnsvenn.png")

#### Gene Overlap Analysis
![alt text](https://raw.githubusercontent.com/kingb12/model-morphing/master/analysisimages/maristadtgenesvenn.png "maristadtgenesvenn.png")

### M. maripaludis to M. barkeri Data

| Model  | Genes | Metabolites  | Reactions |
| ------------- | ------------- | ------------- | ------------- |
| M. maripaludis Source       | 290 | 710 | 636 |
| M. barkeri Morph            | 701 | 892 | 741 |
| M. barkeri Reconstruction   | 452 | 797 | 728 |

#### Reaction Overlap Analysis
![alt text](https://raw.githubusercontent.com/kingb12/model-morphing/master/analysisimages/maribarkrxnsvenn.png "maribarkrxnsvenn.png")

#### Gene Overlap Analysis
![alt text](https://raw.githubusercontent.com/kingb12/model-morphing/master/analysisimages/maribarkgenesvenn.png "maribarkgenesvenn.png")

### M. maripaludis to M. maripaludis Data

| Model  | Genes | Metabolites  | Reactions |
| ------------- | ------------- | ------------- | ------------- |
| M. maripaludis Source       | 528 | 710 | 636 |
| M. maripaludis Morph            | 767 | 742 | 597 |
| M. maripaludis Reconstruction   | 322 | 707 | 530 |

#### Reaction Overlap Analysis
![alt text](https://raw.githubusercontent.com/kingb12/model-morphing/master/analysisimages/marimarirxnsvenn.png "marimarirxnsvenn.png")

#### Gene Overlap Analysis
![alt text](https://raw.githubusercontent.com/kingb12/model-morphing/master/analysisimages/marimarigenesvenn.png "marimarigenesvenn.png")



###Pairwise Comparisons

####*E. coli* vs. Our Models

| Model | Model Only | Shared With *E. Coli* Model | Percent Shared |
| --- | --- | --- | --- |
| MaripaludisModel | 239 | 397 | 62.4213836478% |
| mari_to_mari_morph | 193 | 403 | 67.6174496644% |
| E.coli_model | 1413 | 1413 | 50.0% |
| MaripaludisReconstruction | 133 | 497 | 78.8888888889% |
| StadtmanaeReconstruction | 131 | 555 | 80.9037900875% |
| mari_to_bark_morph_3media | 195 | 554 | 73.9652870494% |
| JannashciiReconstruction | 125 | 509 | 80.2839116719% |
| Mari_to_Janna_Morph | 193 | 438 | 69.4136291601% |
| Mari_to_Stadt_Morph | 191 | 504 | 72.5179856115% |
| BarkeriReconstruction | 139 | 584 | 80.7745504841% |

####*M. maripaludis* Source vs. The Morphs and Reconstructions

| Model | Model Only | Shared With MaripaludisModel | Percent Shared |
| --- | --- | --- | --- |
| Mari_to_Stadt_Morph | 156 | 539 | 77.5539568345% |
| Mari_to_Janna_Morph | 85 | 546 | 86.529318542% |
| mari_to_mari_morph | 37 | 559 | 93.7919463087% |
| mari_to_bark_morph_3media | 209 | 540 | 72.0961281709% |
| StadtmanaeReconstruction | 312 | 374 | 54.5189504373% |
| MaripaludisReconstruction | 217 | 413 | 65.5555555556% |
| JannashciiReconstruction | 258 | 376 | 59.3059936909% |
| MaripaludisModel | 636 | 636 | 50.0% |
| BarkeriReconstruction | 346 | 377 | 52.1438450899% |


## Qualitative Analyses

#### HdrABC_c0 and Hdr-formate_c0
HdrABC and Hdr-formate are *M. maripaludis* specific reactions, must other organisms use a veryt similar alternative pathway with rxn03126_c0 

```python
In [12]: for m in models:
    model = get_object(m[0], m[1])['data']
    print Model.has_reaction(model, 'Hdr-formate_c0')
   ....:     print m[2]
   ....:
False
JannashciiReconstruction
True
Mari_to_Janna_Morph
True
MaripaludisModel
True
Mari_to_Stadt_Morph
False
StadtmanaeReconstruction
True
mari_to_bark_morph_3media
False
BarkeriReconstruction
True
mari_to_mari_morph
False
MaripaludisReconstruction
```
rxn03126_c0 was found in all of the morphs as well, as it was introduced in each reconstruction (M. mariapludis included). Our source model was the only model to exclude this reaction.

What we found is that Hdr-formate_c0 was carried over to the morphs erronerously. In the course of our algorithm, this is unavoidable without drastic change
The ProteomeComparison informs us that each of the methanogens should be able to produce the enzymes for HdrABC and Hdr-formate despite it not belonging in the
target model. To detect this, we would have to distrust the ProteomeComparison. Is this normal or an edge case? 

Looking at our reaction venn comparisons, we see that on average ~ 185 reactions are brought to the morphed model that are not found in an automatic reconstruction. DEpending on phylogeny 50-80% of these had a gene match in our target organism. Hdr-formate would be one of these, but we could go through these lists and see how much information was correctly vs. erronerously gained.

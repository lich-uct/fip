# Feature Interrelation Profiling library (fip)

fip is a Python library for feature interrelation profiling. The library is intended predominantly for use by the cheminformatic community (i.e. for analyzing structural keys, fingerprints, and other chemically relevant binary feature vectors), though it should be able to process any given binary feature vectors due to the generic information theory origin of the [underlying principles](https://en.wikipedia.org/wiki/Pointwise_mutual_information). fip contains a collection of scripts and tools used to create feature interrelation profiles of given chemical structure sets, and to compare said profiles between each other as well as against specific chemical structures. For added convenience, fip also includes several pre-computed interrelation profiles of the DrugBank, ChEMBL, PubChem and ZINC databases.

The current form of fip is still a proof-of-concept, "alpha" version. The extent and direction of further development will be determined by its future use, if any. Therefore, the contained objects, methods, dependencies etc. may still be subject to change. If you have any feature requests, recommendations or other input, please let us know.

## Feature Interrelation Profiling

Feature Interrelation Profiling aims to quantify how much do the binary features observed in a given set of feature vectors affect each other, i.e. how much does the presence of each individual feature affect the presence of all others in the same feature vector. These interrelations can range from negative, where a pair of features appears together less than what could be expected from their individual occurrence probabilities, through neutral, where the feature pair appears together as could be expected from their individual occurrence probabilities (indicating no mutual relation), up to positive interrelations, where a feature pair occurs together more than could be expected if they were independent. The type and strength of such interrelations can be quantified using [pointwise mutual information](https://en.wikipedia.org/wiki/Pointwise_mutual_information), as well as other derivative measurements.

To give a more specific example, many contemporary cheminformatic methodologies use binary features to predict properties of chemical structures. For instance, methods such as [SAScore](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8) and [SYBA](https://github.com/lich-uct/syba) monitor the presence of specific structural patterns to estimate synthetic accessibility of given chemical structures. While this approach is tried and true, there are still the edge cases of quite common structural patterns that seldom appear together in known chemical structures, and, *vice versa*, unusual structural patterns that can be often observed together. Possible causes for negative interrelations could involve structural patterns that are likely to induce self-reaction when present within the same structure. Conversely, positive interrelation may indicate the presence of some overarching structural motifs within the characterized feature set.

To borrow a different interrelation example from the field of [computational linguistics](https://www.aclweb.org/anthology/W13-1504/), some regularly used words also appear together in the same sentence less that could be expected if their use was mutually independent. For example, words like 'sleep' and 'furious' are both relatively common, but seldom seen together in a meaningful sentence. Texts that contain many such unusual combinations might be harder to understand, or possibly even completely nonsensical. Conversely, some words are used together often, such as 'car' and 'driver'.

fip provides tools to quantify such interrelations in given sets of binary feature vectors. It currently implements:
* Ways to create interrelation profiles for given sets of feature vectors:
  * A co-occurrence relational matrix (CORM) that contains the raw co-occurrence counts
  * A co-occurrence probability relational matrix (COPRM) that contains the probability of features co-occurring within a feature vector randomly selected form the characterized set
  * A pointwise mutual information relational matrix (PMIRM) containing [PMI values](https://en.wikipedia.org/wiki/Pointwise_mutual_information) for each feature pair within the characterized set
  * A Z-scored pointwise mutual information relational matrix (ZPMIRM) with PMI values normalized by [Z-score](https://en.wikipedia.org/wiki/Standard_score), in which each feature pair is assigned a value indicating how much does its PMI value compare to all others within the characterized set in terms of standard deviation
* Ways to compare interrelation profiles:
  * Relative feature tightness (RFT) quantifying how well do the feature co-occurrences in a measured interrelation profile or a single feature vector match the interrelations within a given reference interrelation profile
  * Z-scored relative feature tightness (ZRFT) is similar to RFT, but instead of the co-occurrences being matched to absolute PMI values within the reference profile, they are instead matched to ZPMI values that reflect how favorable are the co-occurrences relative to all others within the reference profile.
* Pre-computed data, since creating interrelation profiles for very large datasets can get computationally expensive:
  * A CORM for [DrugBank](https://www.drugbank.ca/), [ChEMBL](https://www.ebi.ac.uk/chembl/), [PubChem](https://pubchem.ncbi.nlm.nih.gov/), [ZINC](http://zinc15.docking.org/), and all four combined, for the MACCS Keys, PubChem Substructure Keys, ECFP4-1024 and ECFP6-1024 feature vectors; 80 CORMs in total
  * Direct conversion of the above CORMs to COPRMs, PMIRMs and ZPMIRMs

## Installation

fip is written in Python3 only, so just cloning this repository and adding it to pythonpath should work fine.
The main dependencies are currently:
* [Pandas](https://pandas.pydata.org/) for data analysis
* [RDKit](https://www.rdkit.org/) for chemistry-related functionality, especially structural fingerprints

As optional tools to complement the above, we recommend:
* [Jupyter](https://jupyter.org/) for easy, interactive use
* [Seaborn](https://seaborn.pydata.org/) or just [matplotlib](https://matplotlib.org/) for visualization
* [ChemFP](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0398-8) for structure fingerprinting

## Getting started

Basic use and interpretation of the implemented interrelation matrices is covered in an included [Jupyter notebook](https://github.com/lich-uct/fip/blob/master/docs/notebooks/fip_basic_usage.ipynb).

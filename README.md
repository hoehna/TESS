# TESS

TESS provides a flexible framework for specifying diversification models---where diversification rates are constant, vary continuously, or change episodically through time (including explicit models of mass-extinction events)---and implements numerical methods to estimate parameters of these models from estimated phylogenies.
A major feature of TESS is the ability to include various methods of incomplete taxon sampling.
Additionally, we provide robust Bayesian methods for assessing the relative fit of these models of lineage diversification to a given study tree---e.g., where stepping-stone simulation is used to estimate the marginal likelihoods of competing models, which can then be compared using Bayes factors.
We also provide Bayesian methods for evaluating the absolute fit of these branching-process models to a given study tree---i.e., where posterior-predictive simulation is used to assess the ability of a candidate model to generate the observed phylogenetic data.

# Installation

At present, TESS is downloadable from https://github.com/hoehna/TESS. There are multiple ways to download it. The easiest is to use devtools and install from GitHub.

### Installing from GitHub using devtools
Run the following code from your R console:

```{r eval=FALSE}
install.packages("devtools")
library(devtools)
install_github("hoehna/TESS")
library(TESS)
```

# Manual and Tutorials

All functions in TESS are documented according to the CRAN standard and the manual can be downloaded from https://raw.githubusercontent.com/hoehna/TESS/master/docs/manual.pdf

Additionally, we provide an elaborate vignette containing several tutorials on how to use TESS available from https://raw.githubusercontent.com/hoehna/TESS/master/docs/Bayesian_Diversification_Rate_Analysis.pdf


# Citations

The idea and methods behind TESS has been published in several papers, listed below. The main paper about TESS is the application note in Bioinformatics. Please cite the software and papers according to your type of analysis.

* Höhna, S. Fast simulation of reconstructed phylogenies under global time-dependent birth-death processes Bioinformatics, Oxford Univ Press, 2013, 29, 1367-1374
* Höhna, S. Likelihood Inference of Non-Constant Diversification Rates with Incomplete Taxon Sampling PLoS One, Public Library of Science, 2014, 9, e84184
* Höhna, S.; Stadler, T.; Ronquist, F. & Britton, T. Inferring speciation and extinction rates under different species sampling schemes Molecular Biology and Evolution, 2011, 28, 2577-2589
* Höhna, S. The time-dependent reconstructed evolutionary process with a key-role for mass-extinction events Journal of Theoretical Biology, 2015, 380, 321-331
* Höhna, S.; May, M. R. & Moore, B. R. TESS: an R package for efficiently simulating phylogenetic trees and performing Bayesian inference of lineage diversification rates Bioinformatics, 2016, 32, 789-791
* May, M. R.; Höhna, S. & Moore, B. R. A Bayesian Approach for Detecting the Impact of Mass-Extinction Events on Molecular Phylogenies When Rates of Lineage Diversification May Vary Methods in Ecology and Evolution, 2016, 7, 947-959

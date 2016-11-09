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

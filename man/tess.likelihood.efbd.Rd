\name{tess.likelihood.efbd}
\alias{tess.likelihood.efbd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.likelihood.efbd: Probability density of a tree under a tree-wide time-dependent fossilized birth-death process}
\description{tess.likelihood computes the probability of a reconstructed phylogenetic tree under time-dependent diversification rates. The rates may be any positive function of time or a constant. Additionally, mass-extinction event can be provided and a uniform taxon sampling probability. You have several options for the start of the process (origin vs MRCA) and the condition of the process (time, survival or taxa; note that survival and taxa implicitly condition on the time too!). Note that constant rates lead to much faster computations.
The likelihood can be computed for incompletely sampled trees if you give a sampling probability != 1.0. You have two options for the sampling strategy: uniform|diversified. The detailed description of these can be found in the references. More information can be obtained in the vignette about how to apply this likelihood function.}
\usage{
tess.likelihood.efbd(nodes,
                     lambda,
                     mu,
                     phi,
                     rateChangeTimesLambda=c(),
                     rateChangeTimesMu=c(),
                     rateChangeTimesPhi=c(),
                     massExtinctionTimes=c(),
                     massExtinctionSurvivalProbabilities=c(),
                     samplingStrategy="uniform",
                     samplingProbability=1.0,
                     MRCA=TRUE,
                     CONDITION="survival",
                     log=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodes}{Information for each node, if it's a sampled ancestor, a fossil tip, or an extant tip. Also node age and parent ages Calculated using tess.branching.times(phy).}
  \item{lambda}{The speciation rate constant or vector.}
  \item{mu}{The extinction rate constant or vector.}
  \item{phi}{The fossilization rate constant or vector.}
  \item{rateChangeTimesLambda}{The speciation change times}
  \item{rateChangeTimesMu}{The extinction change times}
  \item{rateChangeTimesPhi}{The fossilization change times}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{samplingStrategy}{The strategy how samples were obtained. Options are: uniform|diversified|age.}
  \item{MRCA}{Does the process start with the most recent common ancestor? If not, the tree must have a root edge!}
  \item{CONDITION}{do we condition the process on time|survival|taxa?}
  \item{log}{Should we log-transform the likelihood?}
}
\value{Returns the (log) probability of the tree, i.e. the likelihood of the parameters given the tree.}
\references{
S. Hoehna: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374

S. Hoehna: Likelihood Inference of Non-Constant Diversification Rates with Incomplete Taxon Sampling. 2014, PLoS one, Public Library of Science, 9, e84184.

S. Hoehna: The time-dependent reconstructed evolutionary process with a key-role for mass-extinction events. 2015, Journal of Theoretical Biology, 380, 321-331.
}
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# load a test data set
data(cettiidae)

# convert the phylogeny into the branching times
times <- tess.branching.times(cettiidae)

## Set up some rate values
lambda <- c(1, 2, 1)
mu <- c(0.95, 0.5, 0.95)
phi <- c(0.1, 0.2, 0.1)

## Assign the times at which the rates change
changetimes <- c(0.5, 0.3)

## Calculate likelihood
tess.likelihood.efbd(times, lambda, mu, phi,
                     rateChangeTimesLambda = changetimes,
                     rateChangeTimesMu = changetimes,
                     rateChangeTimesPhi = changetimes)

}
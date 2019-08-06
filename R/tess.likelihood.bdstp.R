################################################################################
#
# tess.likelihood.rateshift.R
#
# Copyright (c) 2012- Sebastian Hoehna
#
# This file is part of TESS.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# TESS is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  TESS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with TESS; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA
#
################################################################################



################################################################################
#
# @brief Computation of the likelihood for a given tree under an episodic fossilized-birth-death model (i.e. piecewise constant rates).
#
# @date Last modified: 2018-05-25
# @author Sebastian Hoehna
# @version 5.0
# @since 2018-05-25, version 3.0
#
# @param    nodes                                         list          all node times in the tree (from tess.branching.times)
# @param    lambda                                        scalar        speciation (birth) rate
# @param    mu                                            scalar        extinction (death) rate
# @param    phi                                           scalar        serial sampling rate
# @param    r                                             scalar        conditional probability that a lineage dies upon (serial) sampling
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: uniform|diversified|age
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|sampleAtLeastOneLineage?
# @param    log                                           boolean       likelhood in log-scale?

# @return                                                 scalar        probability of the speciation times
#
################################################################################

tess.likelihood.bdstp <- function( nodes,
                                   lambda,
                                   mu,
                                   phi,
                                   r,
                                   samplingProbability,
                                   samplingStrategy = "uniform",
                                   MRCA=TRUE,
                                   CONDITION="survival",
                                   log=TRUE) {
  
    if ( samplingStrategy != "uniform" ) {
      stop("Only uniform sampling currently supported for BDSTP.")
    }
    
    # If most recent sampling time > 0, shift all nodes upwards
    if (min(nodes$age) > 0) {
      offset <- min(nodes$age)
      nodes$age <- nodes$age - offset
      nodes$age_parent <- nodes$age_parent - offset
    }
    
    # recover()
    
    # Important information
    if (samplingProbability < .Machine$double.eps) {
      # If there is no sampling at the present, all tips are fossil tips
      fossil_times <- nodes$age[nodes$fossil_tip | nodes$tip]
      tip_times <- c()
    } else {
      fossil_times <- nodes$age[nodes$fossil_tip]
      tip_times <- nodes$age[nodes$tip]
    }
    branching_times <- nodes$age[!(nodes$tip | nodes$fossil_tip | nodes$sampled_ancestor)]
    
    n_tips <- length(tip_times)
    n_fossils <- length(fossil_times)
    n_SA <- sum(nodes$sampled_ancestor)
    
    # Necessary functions
    A <- abs(sqrt((lambda - mu - phi)^2 + 4*lambda*phi))
    if (abs(1 - samplingProbability) > .Machine$double.eps) {
      C <- 1 - samplingProbability
    } else {
      C <- 1
    }
    B <- ((1 - 2*C) * lambda + mu + phi)/A
    
    E <- function(t) {
      plus_minus <- (1 + B - exp(-A*t)*(1-B)) / (1 + B + exp(-A*t)*(1-B))
      E_t <- (lambda + mu + phi - A*plus_minus) / (2 * lambda)
      return(E_t)
    }
    
    if ( abs(samplingProbability -1) > .Machine$double.eps ) {
      D <- function(t) {
        (1 - samplingProbability) * (4 * (exp(-A*t))) / 
          ((1 + B + exp(-A*t)*(1-B))^2)
      }
    } else {
      D <- function(t) {
        (4 * (exp(-A*t))) / 
          ((1 + B + exp(-A*t)*(1-B))^2)
      }
    }
    
    # Components of log-likelihood
    if ( samplingProbability > 0 && samplingProbability < 1) {
      lnPr_tip_samples <- n_tips * log(samplingProbability) # log-probability of lineages sampled at the present (if samplingProbability == 0.0 || samplingProbability == 1.0, nothing to compute)
    } else {
      lnPr_tip_samples <- 0
    }
    if (n_fossils > 0) {
      lnPr_fossils <- sum( log(phi * r + phi*(1 - r)*E(fossil_times)) ) # log-probability density of all fossils
    } else {
      lnPr_fossils <- 0
    }
    if (n_SA > 0) {
      lnPr_SA <- n_SA * log(phi * (1 - r)) # log-probability density of all sampled ancestors
    } else {
      lnPr_SA <- 0
    }
    lnPr_births <- (n_tips+n_fossils-1) * log(lambda) # log-probability density of all birth events
    lnPr_branch_segments <- log(D(max(nodes$age))) + sum(log(D(branching_times))) - sum(log(D(c(tip_times,fossil_times)))) # log-probability density of all branch segments
    
    lnl <- lnPr_fossils + lnPr_SA + lnPr_births + lnPr_branch_segments
    
    # Conditioning
    if ( CONDITION == "survival" ) {
      root_age <- max(nodes$age)
      lnl <- lnl - (1 + MRCA) * log(1 - E(root_age))
    } else if ( CONDITION == "sampleAtLeastOneLineage" ) {
      root_age <- max(nodes$age)
      lnl <- lnl - log(1 - E(root_age))
    }
    
    
    if ( log == FALSE ) {
      lnl <- exp(lnl)
    }
    
    return(lnl)
}

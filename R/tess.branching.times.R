################################################################################
#
# tess.branching.times.R
#
# Copyright (c) 2018- Sebastian Hoehna
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
# Recursively compute branching times for phylogenetic tree.
#  Allows for non-ultrametric (fossil) trees.
#
# This code is adapted from bammtools/BAMMtools/R/NU.branching.times.R
#
################################################################################

#' Title
#'
#' @param phy the phylogeny
#' @param tip.age.threshold threshold to be considered an extant tip (not fossil)
#'
#' @return a list with items c("sampled_ancestor", "fossil_tip", "age_parent", "age", "tip")
#' @export
#'
#' @examples
#' data(conifers)
#'
#' nodes <- tess.branching.times(conifers)
tess.branching.times <- function(phy, tip.age.threshold=1e-5) {

  # Do it recursively
  fx <- function(phy, node, cur_time, nodes) {

    children <- phy$edge[,2][phy$edge[,1] == node]

    for (i in 1:length(children)) {
      child <- children[i]
      child_edge <- which(phy$edge[,2] == child)
      branch_time <- cur_time + phy$edge.length[child_edge]

      nodes$age[child]              <- branch_time
      nodes$age_parent[child]       <- nodes$age[node]
      nodes$sampled_ancestor[child] <- ( length(phy$edge[,2][phy$edge[,1] == child]) == 1 )
      nodes$tip[child]              <- ( length(phy$edge[,2][phy$edge[,1] == child]) == 0 )

      if (child > length(phy$tip.label)) {
        nodes <- fx(phy, node = child, cur_time = branch_time, nodes)
      }

    }

    return (nodes)
  }

  nodes <- list(sampled_ancestor=c(),fossil_tip=c(),age_parent=c(),age=c(),tip=c())

  # do the recursive call starting with the root
  nodes$age[length(phy$tip.label) + 1]              <- 0.0
  nodes$age_parent[length(phy$tip.label) + 1]       <- Inf
  nodes$tip[length(phy$tip.label) + 1]              <- FALSE
  nodes$sampled_ancestor[length(phy$tip.label) + 1] <- FALSE
  nodes$fossil_tip[length(phy$tip.label) + 1]       <- FALSE
  nodes <- fx(phy, node = length(phy$tip.label) + 1, cur_time = 0, nodes)

  max_bt <- max(nodes$age)
  nodes$age        <- max_bt - nodes$age
  nodes$age_parent <- max_bt - nodes$age_parent

  nodes$age_parent[length(phy$tip.label) + 1]       <- Inf

  nodes$age[ nodes$age < tip.age.threshold ] <- 0.0
  nodes$fossil_tip <- nodes$age > tip.age.threshold & nodes$tip
  nodes$tip <- nodes$fossil_tip == FALSE & nodes$tip

  return( nodes )

}


################################################################################
# Compute how many lineages are active at time t in the reconstructed phylogeny.
#  Allows for non-ultrametric (fossil) trees.
#
# Uses nodes object from tess.branching.times.
#
################################################################################

tess.num.active.lineages <- function(nodes,t) {

  epsilon <- ifelse( t==0, 1E-3, 1E-5 )
  if ((t-epsilon) > max(nodes$age)) {
    return(0)
  } else {
    n_born_before <- sum(nodes$age[!(nodes$tip | nodes$fossil_tip | nodes$sampled_ancestor)] > (t-epsilon))
    n_died_before <- sum(nodes$age[nodes$fossil_tip] > (t+epsilon))
    return(1 + n_born_before - n_died_before)
  }
}


################################################################################
# Compute probability of any arbitrary topology under any non-lineage-specific BDP.
#  Allows for non-ultrametric (fossil) trees.
#
# Uses nodes object from tess.branching.times.
#
################################################################################

tess.topology.probability <- function(nodes, log=TRUE) {
  lnl <- (sum(nodes$tip) + sum(nodes$fossil_tip) - 1) * log(2) -
    lfactorial(sum(nodes$sampled_ancestor) + sum(nodes$tip) + sum(nodes$fossil_tip))

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }

  return(lnl)
}

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
#' tess.likelihood.ebdstp
#'
#' @description Computation of the likelihood for a given tree under an episodic fossilized-birth-death model (i.e. piecewise constant rates).
#'
#' @param nodes node times from tess.branching.times
#' @param lambda birth (speciation or infection) rates
#' @param mu death (extinction or becoming non-infectious without treatment) rates
#' @param phi serial sampling (fossilization) rates
#' @param r treatment probability, Pr(death | sample) (does not apply to samples take at time 0)
#' @param rateChangeTimesLambda times at which birth rates change
#' @param rateChangeTimesMu times at which death rates change
#' @param rateChangeTimesPhi times at which serial sampling rates change
#' @param rateChangeTimesR times at which treatment probabilities change
#' @param massDeathTimes time at which mass-deaths (mass-extinctions) happen
#' @param massDeathProbabilities probability of a lineage dying in a mass-death event
#' @param burstBirthTimes
#' @param burstBirthProbabilities
#' @param eventSamplingTimes time at which every lineage in the tree may be sampled
#' @param eventSamplingProbabilities probability of a lineage being sampled at an event sampling time
#' @param samplingStrategyAtPresent Which strategy was used to obtain the samples (taxa). Options are: uniform|diversified|age
#' @param samplingProbabilityAtPresent probability of uniform sampling at present
#' @param MRCA does the tree start at the mrca?
#' @param CONDITION do we condition the process on nothing|survival|taxa?
#' @param log likelhood in log-scale?
#'
#' @return probability of the speciation times
#' @export
#'
#' @examples
#' data(conifers)
#'
#' nodes <- tess.branching.times(conifers)
#'
#' lambda <- c(0.2, 0.1, 0.3)
#' mu <- c(0.1, 0.05, 0.25)
#' phi <- c(0.1, 0.2, 0.05)
#'
#' changetimes <- c(100, 200)
#'
#' tess.likelihood.ebdstp(nodes,
#'                        lambda = lambda,
#'                        mu = mu,
#'                        phi = phi,
#'                        r = 0.0,
#'                        rateChangeTimesLambda = changetimes,
#'                        rateChangeTimesMu = changetimes,
#'                        rateChangeTimesPhi = changetimes,
#'                        samplingProbability = 1.0
#' )
tess.likelihood.ebdstp <- function( nodes,
                                    lambda,
                                    mu,
                                    phi,
                                    r,
                                    samplingProbabilityAtPresent,
                                    rateChangeTimesLambda = c(),
                                    rateChangeTimesMu = c(),
                                    rateChangeTimesPhi = c(),
                                    rateChangeTimesR = c(),
                                    massDeathTimes = c(),
                                    massDeathProbabilities = c(),
                                    burstBirthTimes = c(),
                                    burstBirthProbabilities = c(),
                                    eventSamplingTimes = c(),
                                    eventSamplingProbabilities = c(),
                                    samplingStrategyAtPresent = "uniform",
                                    MRCA=TRUE,
                                    CONDITION="survival",
                                    log=TRUE) {

  if ( length(lambda) != (length(rateChangeTimesLambda)+1) || length(mu) != (length(rateChangeTimesMu)+1) || length(phi) != (length(rateChangeTimesPhi)+1) || length(r) != (length(rateChangeTimesR)+1) ) {
    stop("Number of rate-change times needs to be one less than the number of rates!")
  }

  if ( length(massDeathTimes) != length(massDeathProbabilities) || length(burstBirthTimes) != length(burstBirthProbabilities) || length(eventSamplingTimes) != length(eventSamplingProbabilities) ) {
    stop("Number of mass-extinction times needs to equal the number of mass-extinction survival probabilities!")
  }

  # if ( CONDITION != "time" && CONDITION != "survival" && CONDITION != "taxa" ) {
  #    stop("Wrong choice of argument for \"CONDITION\". Possible option are time|survival|taxa.")
  # }
  #
  # if ( samplingStrategyAtPresent != "uniform" && samplingStrategyAtPresent != "diversified") {
  #    stop("Wrong choice of argument for \"samplingStrategyAtPresent\". Possible option are uniform|diversified.")
  # }

  if ( CONDITION != "time" && CONDITION != "survival" && CONDITION != "sampleAtLeastOneLineage" ) {
    stop("Wrong choice of argument for \"CONDITION\". Possible option are time|survival|sampleAtLeastOneLineage")
  }

  if ( samplingStrategyAtPresent != "uniform") {
    stop("Wrong choice of argument for \"samplingStrategyAtPresent\". Possible options are uniform")
  }

  # make sure the times and values are sorted
  if ( length(rateChangeTimesLambda) > 0 ) {
    sortedRateChangeTimesLambda <- sort( rateChangeTimesLambda )
    if ( !(all(sortedRateChangeTimesLambda == rateChangeTimesLambda)) ) {
      stop("Rate change times must be sorted in increasing order")
    }
  }
  if ( length(rateChangeTimesMu) > 0 ) {
    sortedRateChangeTimesMu <- sort( rateChangeTimesMu )
    if ( !(all(sortedRateChangeTimesMu == rateChangeTimesMu)) ) {
      stop("Rate change times must be sorted in increasing order")
    }
  }
  if ( length(rateChangeTimesPhi) > 0 ) {
    sortedRateChangeTimesPhi <- sort( rateChangeTimesPhi )
    if ( !(all(sortedRateChangeTimesPhi == rateChangeTimesPhi)) ) {
      stop("Rate change times must be sorted in increasing order")
    }
  }
  if ( length(rateChangeTimesR) > 0 ) {
    sortedRateChangeTimesR <- sort( rateChangeTimesR )
    if ( !(all(sortedRateChangeTimesR == rateChangeTimesR)) ) {
      stop("Rate change times must be sorted in increasing order")
    }
  }
  if ( length(burstBirthTimes) > 0 ) {
    sortedBurstBirthTimes <- sort( burstBirthTimes )
    if ( !(all(sortedBurstBirthTimes == burstBirthTimes)) ) {
      stop("Event times must be sorted in increasing order")
    }
  }
  if ( length(massDeathTimes) > 0 ) {
    sortedMassDeathTimes <- sort( massDeathTimes )
    if ( !(all(sortedMassDeathTimes == massDeathTimes)) ) {
      stop("Event times must be sorted in increasing order")
    }
  }
  if ( length(eventSamplingTimes) > 0 ) {
    sortedEventSamplingTimes <- sort( eventSamplingTimes )
    if ( !(all(sortedEventSamplingTimes == eventSamplingTimes)) ) {
      stop("Event times must be sorted in increasing order")
    }
  }

  # If most recent sampling time > 0, shift all nodes upwards
  if (min(nodes$age) > 0) {
    offset <- min(nodes$age)
    nodes$age <- nodes$age - offset
    nodes$age_parent <- nodes$age_parent - offset
  }

  # recover()

  # join the times of the rate changes and the birth/death/sampling events
  if ( length( rateChangeTimesLambda ) > 0 ||  length( rateChangeTimesMu ) > 0 ||  length( rateChangeTimesPhi ) > 0 || length( rateChangeTimesR ) > 0 || length( massDeathTimes ) > 0 || length( burstBirthTimes ) > 0 || length( eventSamplingTimes ) > 0 ) {
    changeTimes <- sort( unique( c( rateChangeTimesLambda, rateChangeTimesMu, rateChangeTimesPhi, rateChangeTimesR, massDeathTimes, burstBirthTimes, eventSamplingTimes ) ) )
  } else {
    changeTimes <- c()
  }
  birth_rate            <- rep(NaN,length(changeTimes)+1)
  death_rate            <- rep(NaN,length(changeTimes)+1)
  sampling_rate         <- rep(NaN,length(changeTimes)+1)
  treatment_probability <- rep(NaN,length(changeTimes)+1)
  bbp <- rep(NaN,length(changeTimes))
  mdp <- rep(NaN,length(changeTimes))
  esp <- rep(NaN,length(changeTimes))
  birth_rate[1] <- lambda[1]
  if ( length(lambda) > 1 ) {
    birth_rate[ match(rateChangeTimesLambda,changeTimes)+1 ] <- lambda[ 2:length(lambda) ]
  }
  death_rate[1] <- mu[1]
  if ( length(mu) > 1 ) {
    death_rate[ match(rateChangeTimesMu,changeTimes)+1 ] <- mu[ 2:length(mu) ]
  }
  sampling_rate[1] <- phi[1]
  if ( length(phi) > 1 ) {
    sampling_rate[ match(rateChangeTimesPhi,changeTimes)+1 ] <- phi[ 2:length(phi) ]
  }
  treatment_probability[1] <- r[1]
  if ( length(r) > 1 ) {
    treatment_probability[ match(rateChangeTimesR,changeTimes)+1 ] <- r[ 2:length(r) ]
  }
  if ( length( burstBirthTimes ) > 0 ) {
    bbp[ match(burstBirthTimes,changeTimes) ] <- burstBirthProbabilities[ 1:length(burstBirthProbabilities) ]
  }
  if ( length( massDeathTimes ) > 0 ) {
    mdp[ match(massDeathTimes,changeTimes) ] <- massDeathProbabilities[ 1:length(massDeathProbabilities) ]
  }
  if ( length( eventSamplingTimes ) > 0 ) {
    esp[ match(eventSamplingTimes,changeTimes) ] <- eventSamplingProbabilities[ 1:length(eventSamplingProbabilities) ]
  }
  for ( i in seq_len(length(changeTimes)) ) {
    if ( is.null(birth_rate[i+1]) || !is.finite(birth_rate[i+1]) ) {
       birth_rate[i+1] <- birth_rate[i]
    }
    if ( is.null(death_rate[i+1]) || !is.finite(death_rate[i+1]) ) {
       death_rate[i+1] <- death_rate[i]
    }
    if ( is.null(sampling_rate[i+1]) || !is.finite(sampling_rate[i+1]) ) {
       sampling_rate[i+1] <- sampling_rate[i]
    }
    if ( is.null(treatment_probability[i+1]) || !is.finite(treatment_probability[i+1]) ) {
      treatment_probability[i+1] <- sampling_rate[i]
    }
    if ( is.null(bbp[i]) || !is.finite(bbp[i]) ) {
      bbp[i] <- 0.0
    }
    if ( is.null(mdp[i]) || !is.finite(mdp[i]) ) {
      mdp[i] <- 0.0
    }
    if ( is.null(esp[i]) || !is.finite(esp[i]) ) {
      esp[i] <- 0.0
    }
  }

  # set the uniform taxon sampling probability
  if (samplingStrategyAtPresent == "uniform") {
    rho <- samplingProbabilityAtPresent
  } else {
    rho <- 1.0
  }

  lambda <- birth_rate
  mu     <- death_rate
  phi    <- sampling_rate
  r      <- treatment_probability
  burstBirthProbabilities <- c(0.0,bbp)
  massDeathProbabilities <- c(0.0,mdp)
  eventSamplingProbabilities <- c(rho,esp)

  # recover()

  # Precompute vectors
  ABCDE <- precomputeVectors(changeTimes,lambda,mu,phi,r,burstBirthProbabilities,massDeathProbabilities,eventSamplingProbabilities)

  # Classify tips and bifurcations, do they belong to an event or are they between events?
  event_tips_and_fossils <- vector("list",length(changeTimes)+1)
  event_sampled_ancestors <- vector("list",length(changeTimes)+1)
  serial_tips_and_fossils <- rep(NA,length(nodes$age))
  serial_sampled_ancestors <- rep(NA,length(nodes$age))
  event_bifurcations <- vector("list",length(changeTimes)+1)
  serial_bifurcations <- rep(NA,length(nodes$age))
  for (i in 1:length(nodes$age)) {
    idx <- findInterval(nodes$age[i],changeTimes,left.open=TRUE)+1
    ti  <- ifelse( idx <= 1, 0.0, changeTimes[idx-1] )
#    at_time <- abs(nodes$age[i] - ti) < .Machine$double.eps
    # there might be rounding errors in branch lengths, so we need to be a bit flexible
    at_time <- abs(nodes$age[i] - ti) < ifelse( idx <= 1, 1E-3, 1E-5 )
    if ( length(at_time) != 1 ) {
      cat("ti =",ti,"\n")
      cat("idx =",idx,"\n")
      cat("nodes$age[i] =",nodes$age[i],"\n")
    }

result = tryCatch({
    if ( at_time == FALSE ) {
      test <- at_time == FALSE
    }
}, warning = function(w) {

}, error = function(e) {
  cat("ti =",ti,"\n")
  cat("idx =",idx,"\n")
  cat("nodes$age[i] =",nodes$age[i],"\n")
  cat("at_time =",at_time,"\n")
  cat("changeTimes =",changeTimes,"\n")
}, finally = {
})

    if ( at_time == FALSE ) {
      if (nodes$tip[i] | nodes$fossil_tip[i]) {
        serial_tips_and_fossils[i] <- nodes$age[i]
      } else if (nodes$sampled_ancestor[i]) {
        serial_sampled_ancestors[i] <- nodes$age[i]
      } else {
        serial_bifurcations[i] <- nodes$age[i]
      }

    } else {

      if (nodes$tip[i] | nodes$fossil_tip[i]) {
        if (eventSamplingProbabilities[idx] >= .Machine$double.eps) {
          event_tips_and_fossils[[idx]] <- c(event_tips_and_fossils[[idx]],nodes$age[i])
        } else {
          serial_tips_and_fossils[i] <- nodes$age[i]
        }
      } else if (nodes$sampled_ancestor[i]) {
        if (eventSamplingProbabilities[idx] >= .Machine$double.eps) {
          event_sampled_ancestors[[idx]] <- c(event_sampled_ancestors[[idx]],nodes$age[i])
        } else {
          serial_sampled_ancestors[i] <- nodes$age[i]
        }
      } else {
        if (burstBirthProbabilities[idx] >= .Machine$double.eps) {
          event_bifurcations[[idx]] <- c(event_bifurcations[[idx]],nodes$age[i])
        } else {
          serial_bifurcations[i] <- nodes$age[i]
        }
      }
    }
  }
  serial_tips_and_fossils <- serial_tips_and_fossils[!is.na(serial_tips_and_fossils)]
  serial_sampled_ancestors <- serial_sampled_ancestors[!is.na(serial_sampled_ancestors)]
  serial_bifurcations <- serial_bifurcations[!is.na(serial_bifurcations)]

  # recover()

  # initialize the log likelihood
  lnl <- 0

  # Sampling event probabilities
  A0 <- tess.num.active.lineages(nodes,0)
  N0 <- length(event_sampled_ancestors[[1]]) + length(event_tips_and_fossils[[1]])
  # Sampling events with Phi = 0 or Phi = 1 are special cases
  if ( samplingProbabilityAtPresent >= .Machine$double.eps && abs(samplingProbabilityAtPresent - 1.0) >= .Machine$double.eps ) {
    if ( A0 == N0 ) {
      lnPr_event_samples_at_present <- samplingProbabilityAtPresent^N0
    } else {
      lnPr_event_samples_at_present <- -Inf
    }
  } else if ( abs(samplingProbabilityAtPresent - 1.0) < .Machine$double.eps ) {
    if ( A0 == N0 ) {
      lnPr_event_samples_at_present <- 0.0
    } else {
      lnPr_event_samples_at_present <- -Inf
    }
  } else {
    # There are no tips to be event-sampled here
    lnPr_event_samples_at_present <- 0.0
  }
  if ( length(changeTimes) > 0 ) {
    lnPr_event_samples_not_at_present <- sum(sapply(2:(length(changeTimes)+1),function(i){
      Ai <- tess.num.active.lineages(nodes,changeTimes[i-1])
      Ni <- length(event_sampled_ancestors[[i]]) + length(event_tips_and_fossils[[i]])
      Ri <- length(event_sampled_ancestors[[i]])
      if (Ni > Ai) {
        return(-Inf)
      } else {
        Ei <- E.ebdstp(i,changeTimes[i], lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE )
        # If sampling probability here is 0 or 1, need to handle separately
        if ( eventSamplingProbabilities[i] >= .Machine$double.eps && abs(eventSamplingProbabilities[i] - 1 ) >= .Machine$double.eps ) {
          pr_sampling_event <- (1 - eventSamplingProbabilities[i])^(Ai-Ni) * (eventSamplingProbabilities[i]^Ni) * ((1 - r[i])^Ri) * ((r[i]) + (((1 - r[i])*Ei)^(Ni-Ri)))
          return(log(pr_sampling_event))
        } else if ( eventSamplingProbabilities[i] < .Machine$double.eps ) {
          # if there are samples and there is no sampling probability, this is impossible
          if ( Ni > 0) {
            return(-Inf)
          } else {
            return(0.0)
          }
        } else if ( abs(eventSamplingProbabilities[i] - 1 ) < .Machine$double.eps ) {
          # If not all active lineages are sampled, this is impossible
          if (Ni != Ai) {
            return(-Inf)
          } else {
            pr_sampling_event <- ((1 - r[i])^Ri) * ((r[i]) + (((1 - r[i])*Ei)^(Ni-Ri)))
            return(log(pr_sampling_event))
          }
        }
      }
    }))
    lnPr_event_samples <- lnPr_event_samples_at_present + lnPr_event_samples_not_at_present
  } else {
    lnPr_event_samples <- lnPr_event_samples_at_present
  }

  # Serially sampled tip/fossil probabilities
  if ( length(serial_tips_and_fossils) > 0 ) {
    lnPr_serial_tips <- sum(sapply(serial_tips_and_fossils,function(t){
      idx <- findInterval(t,changeTimes,left.open=TRUE)+1
      log(phi[idx] * (r[idx] + (1 - r[idx])*E.ebdstp(idx,t,lambda,mu,phi,r,burstBirthProbabilities,massDeathProbabilities,eventSamplingProbabilities,changeTimes,ABCDE)))
    }))
  } else {
    lnPr_serial_tips <- 0.0
  }

  # Serially sampled ancestor probabilities
  if ( length(serial_sampled_ancestors) > 0 ) {
    lnPr_serial_sampled_ancestors <- sum(sapply(serial_sampled_ancestors,function(t){
      idx <- findInterval(t,changeTimes,left.open=TRUE)+1
      log(phi[idx] * (1 - r[idx]))
    }))
  } else {
    lnPr_serial_sampled_ancestors <- 0.0
  }

  # Burst birth probabilities
  if ( length(changeTimes) > 0 ) {
    lnPr_burst_births <- sum(sapply(2:(length(changeTimes)+1),function(i){
      if (burstBirthProbabilities[i] > .Machine$double.eps) {
        Ki <- length(event_bifurcations[[i]])
        Ai <- tess.num.active.lineages(nodes,changeTimes[i])
        Ei <- E.ebdstp(i,changeTimes[i], lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE )
        pr_sampling_event <- (burstBirthProbabilities[i]^Ki) * ((burstBirthProbabilities[i]^(Ai - Ki))*Ei + ((1 - burstBirthProbabilities[i])^(Ai - Ki)))
        return(log(pr_sampling_event))
      } else {
        return(0.0)
      }
    }))
  } else{
    lnPr_burst_births <- 0.0
  }

  # Serial birth probabilities
  lnPr_births <- sum(sapply(serial_bifurcations,function(t){
    idx <- findInterval(t,changeTimes,left.open=TRUE)+1
    log(lambda[idx])
  }))

  # Branch segment probabilities
  root_age <- max(nodes$age)
  root_idx <- findInterval(root_age,changeTimes,left.open=TRUE)+1
  lnPr_branch_segments <- sum(sapply(1:length(nodes$age),function(i){
    if ( nodes$tip[i] || nodes$fossil_tip[i] ) {
      t <- nodes$age[i]
      idx <- findInterval(t,changeTimes,left.open=TRUE)+1
      return(-log(D.ebdstp(idx, t, lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE )))
    } else if ( !(nodes$tip[i] || nodes$fossil_tip[i] || nodes$sampled_ancestor[i]) ) {
      t <- nodes$age[i]
      idx <- findInterval(t,changeTimes,left.open=TRUE)+1
      return(log(D.ebdstp(idx, t, lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE )))
    } else {
      return(0)
    }
  })) + log(D.ebdstp(root_idx, root_age, lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE ))

  if (is.nan(lnl)) lnl <- -Inf

  lnl <- lnl + lnPr_event_samples
  lnl <- lnl + lnPr_serial_tips
  lnl <- lnl + lnPr_serial_sampled_ancestors
  lnl <- lnl + lnPr_burst_births
  lnl <- lnl + lnPr_births
  lnl <- lnl + lnPr_branch_segments

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }

  return (lnl)
}



E.ebdstp <- function( idx, t, lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE ) {

   # idx <- findInterval(t,rateChangeTimes,left.open=TRUE)+1

   # get the parameters
   birth       <- lambda[idx]
   death       <- mu[idx]
   sampling    <- phi[idx]
   treatment   <- r[idx]

   # BIRTH       <- burstBirthProbabilities[idx]
   # DEATH       <- massDeathProbabilities[idx]
   # SAMPLING    <- eventSamplingProbabilities[idx]

   A <- ABCDE$A[idx]
   B <- ABCDE$B[idx]

   ti  <- ifelse( idx <= 1, 0.0, changeTimes[idx-1] )

   diff <- birth - death - sampling
   dt   <- t - ti

   e <- exp(-A*dt)
   tmp <- birth + death + sampling - A *
     ((1.0+B)-e*(1.0-B))/((1.0+B)+e*(1.0-B))

   return ( tmp / (2.0*birth) )
}



D.ebdstp <- function( idx, t, lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE ) {

  # idx <- findInterval(t,rateChangeTimes,left.open=TRUE)+1

  # get the parameters
  birth       <- lambda[idx]
  death       <- mu[idx]
  sampling    <- phi[idx]
  treatment   <- r[idx]

  BIRTH       <- burstBirthProbabilities[idx]
  DEATH       <- massDeathProbabilities[idx]
  SAMPLING    <- eventSamplingProbabilities[idx]

  A <- ABCDE$A[idx]
  B <- ABCDE$B[idx]

  D_minus <- ABCDE$D_minus[idx]
  E_minus <- ABCDE$E_minus[idx]

  ti  <- ifelse( idx <= 1, 0.0, changeTimes[idx-1] )
  dt <- t - ti

  tmp_event <- D_minus * ifelse(abs(SAMPLING - 1) > .Machine$double.eps, 1 - SAMPLING, 1.0) * (1 - DEATH) * (1 - BIRTH + 2*BIRTH*E_minus)
  e   <- exp(A*dt)
  tmp_rate <- (1.0+B) + e*(1.0-B)


  return ( tmp_event * 4.0*e / (tmp_rate*tmp_rate) )
}

precomputeVectors <- function(changeTimes,
                              lambda,
                              mu,
                              phi,
                              r,
                              burstBirthProbabilities,
                              massDeathProbabilities,
                              eventSamplingProbabilities) {

  ABCDE <- list()

  ABCDE$A  <- numeric(length(changeTimes)+1)
  ABCDE$B  <- numeric(length(changeTimes)+1)
  ABCDE$C  <- numeric(length(changeTimes)+1)
  ABCDE$D_minus <- numeric(length(changeTimes)+1)
  ABCDE$E_minus <- numeric(length(changeTimes)+1)

  ABCDE$A[1] <- sqrt((lambda[1] - mu[1] - phi[1])^2 + 4 * lambda[1] * phi[1])
  if ( abs(eventSamplingProbabilities[1] - 1) > .Machine$double.eps ) {
    ABCDE$C[1] <- 1.0 - eventSamplingProbabilities[1]
  } else {
    ABCDE$C[1] <- 1.0
  }
  ABCDE$B[1] <- ((1 - 2 * ABCDE$C[1]) * lambda[1] + mu[1] + phi[1]) / ABCDE$A[1]
  ABCDE$D_minus[1] <- 1.0
  ABCDE$E_minus[1] <- 1.0

  if ( length(changeTimes) > 1 ) {
    for (i in 2:(length(changeTimes) + 1)) {
      ti <- changeTimes[i-1]

      ABCDE$A[i] <- sqrt((lambda[i] - mu[i] - phi[i])^2 + 4 * lambda[i] * phi[i])
      ABCDE$C[i] <- (1 - eventSamplingProbabilities[i]) * (((1 - burstBirthProbabilities[i])*(1 - massDeathProbabilities[i])*ABCDE$E_minus[i]) + ((1 - massDeathProbabilities[i])*ABCDE$E_minus[i]^2) + ((1 - burstBirthProbabilities[i])*massDeathProbabilities[i]))
      ABCDE$B[i] <- ((1 - 2 * ABCDE$C[i]) * lambda[i] + mu[i] + phi[i]) / ABCDE$A[i]

      ABCDE$E_minus[i] <- E.ebdstp(i-1, ti, lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE)
      ABCDE$D_minus[i] <- D.ebdstp(i-1, ti, lambda, mu, phi, r, burstBirthProbabilities, massDeathProbabilities, eventSamplingProbabilities, changeTimes, ABCDE)
    }

  }

  return(ABCDE)
}

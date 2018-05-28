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
# @param    times                                         vector        vector of branching times
# @param    times                                         vector        branching times
# @param    lambda                                        vector        speciation rates
# @param    mu                                            vector        extinction rates
# @param    rateChangeTimesLambda                         vector        speciation rates
# @param    rateChangeTimesMu                             vector        extinction rates
# @param    massExtinctionTimes                           vector        time at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: uniform|diversified|age
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @param    log                                           boolean       likelhood in log-scale?

# @return                                                 scalar        probability of the speciation times
#
################################################################################

tess.likelihood.efbd <- function( nodes,
                                  lambda,
                                  mu,
                                  phi,
                                  rateChangeTimesLambda = c(),
                                  rateChangeTimesMu = c(),
                                  rateChangeTimesPhi = c(),
                                  massExtinctionTimes = c(),
                                  massExtinctionSurvivalProbabilities = c(),
                                  samplingStrategy = "uniform",
                                  samplingProbability = 1.0,
                                  MRCA=TRUE,
                                  CONDITION="survival",
                                  log=TRUE) {

   if ( length(lambda) != (length(rateChangeTimesLambda)+1) || length(mu) != (length(rateChangeTimesMu)+1) || length(phi) != (length(rateChangeTimesPhi)+1) ) {
      stop("Number of rate-change times needs to be one less than the number of rates!")
   }

   if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
      stop("Number of mass-extinction times needs to equal the number of mass-extinction survival probabilities!")
   }

   if ( CONDITION != "time" && CONDITION != "survival" && CONDITION != "taxa" ) {
      stop("Wrong choice of argument for \"CONDITION\". Possible option are time|survival|taxa.")
   }

   if ( samplingStrategy != "uniform" && samplingStrategy != "diversified") {
      stop("Wrong choice of argument for \"samplingStrategy\". Possible option are uniform|diversified.")
   }

   # make sure the times and values are sorted
   if ( length(rateChangeTimesLambda) > 0 ) {
      sortedRateChangeTimesLambda <- sort( rateChangeTimesLambda )
      lambda <- c(lambda[1], lambda[ match(sortedRateChangeTimesLambda,rateChangeTimesLambda)+1 ] )
      rateChangeTimesLambda <- sortedRateChangeTimesLambda
   }
   if ( length(rateChangeTimesMu) > 0 ) {
      sortedRateChangeTimesMu <- sort( rateChangeTimesMu )
      mu <- c(mu[1], mu[ match(sortedRateChangeTimesMu,rateChangeTimesMu)+1 ] )
      rateChangeTimesMu <- sortedRateChangeTimesMu
   }
   if ( length(rateChangeTimesPhi) > 0 ) {
      sortedRateChangeTimesPhi <- sort( rateChangeTimesPhi )
      phi <- c(phi[1], phi[ match(sortedRateChangeTimesPhi,rateChangeTimesPhi)+1 ] )
      rateChangeTimesPhi <- sortedRateChangeTimesPhi
   }
   if ( length(massExtinctionTimes) > 0 ) {
      sortedMassExtinctionTimes <- sort( massExtinctionTimes )
      massExtinctionSurvivalProbabilities <- massExtinctionSurvivalProbabilities[ match(sortedMassExtinctionTimes,massExtinctionTimes) ]
      massExtinctionTimes <- sortedMassExtinctionTimes
   }

   # join the times of the rate changes and the mass-extinction events
   if ( length( rateChangeTimesLambda ) > 0 ||  length( rateChangeTimesMu ) > 0 ||  length( rateChangeTimesPhi ) > 0 || length( massExtinctionTimes ) > 0 ) {
      changeTimes <- sort( unique( c( rateChangeTimesLambda, rateChangeTimesMu, rateChangeTimesPhi, massExtinctionTimes ) ) )
   } else {
      changeTimes <- c()
   }
   speciation    <- rep(NaN,length(changeTimes)+1)
   extinction    <- rep(NaN,length(changeTimes)+1)
   fossilization <- rep(NaN,length(changeTimes)+1)
   mep <- rep(NaN,length(changeTimes))
   speciation[1] <- lambda[1]
   if ( length(lambda) > 1 ) {
      speciation[ match(rateChangeTimesLambda,changeTimes)+1 ] <- lambda[ 2:length(lambda) ]
   }
   extinction[1] <- mu[1]
   if ( length(mu) > 1 ) {
      extinction[ match(rateChangeTimesMu,changeTimes)+1 ] <- mu[ 2:length(mu) ]
   }
   fossilization[1] <- phi[1]
   if ( length(phi) > 1 ) {
      fossilization[ match(rateChangeTimesPhi,changeTimes)+1 ] <- phi[ 2:length(phi) ]
   }
   if ( length( massExtinctionSurvivalProbabilities ) > 0 ) {
      mep[ match(massExtinctionTimes,changeTimes) ] <- massExtinctionSurvivalProbabilities[ 1:length(massExtinctionSurvivalProbabilities) ]
   }
   for ( i in seq_len(length(changeTimes)) ) {
      if ( is.null(speciation[i+1]) || !is.finite(speciation[i+1]) ) {
         speciation[i+1] <- speciation[i]
      }
      if ( is.null(extinction[i+1]) || !is.finite(extinction[i+1]) ) {
         extinction[i+1] <- extinction[i]
      }
      if ( is.null(fossilization[i+1]) || !is.finite(fossilization[i+1]) ) {
         fossilization[i+1] <- fossilization[i]
      }
      if ( is.null(mep[i]) || !is.finite(mep[i]) ) {
         mep[i] <- 1.0
      }
   }

   # set the uniform taxon sampling probability
   if (samplingStrategy == "uniform") {
      rho <- samplingProbability
   } else {
      rho <- 1.0
   }
  
   rateChangeTimes <- changeTimes
   massExtinctionTimes <- changeTimes

   lambda <- rev(speciation)
   mu     <- rev(extinction)
   phi    <- rev(fossilization)
   massExtinctionSurvivalProbabilities <- c(rho,rev(mep))



   # initialize the log likelihood
   lnl <- 0

   if ( length(rateChangeTimes) > 0 ) {
      speciation.rate <- function(t) {
         idx <- findInterval(t,rateChangeTimes)+1
         idx[ idx > length(lambda) ] <- length(lambda)
         return ( lambda[idx] )
      }
      fossilization.rate <- function(t) {
         idx <- findInterval(t,rateChangeTimes)+1
         idx[ idx > length(phi) ] <- length(phi)
         return ( phi[idx] )
      }
   } else {
      speciation.rate     <- function(times) rep(lambda[1],length(times))
      fossilization.rate  <- function(times) rep(phi[1],length(times))
   }
   
   # add the serial tip age terms
   if ( sum( nodes$fossil_tip ) > 0 ) {
      t <- nodes$ages[ nodes$fossil_tip ]
      lnl <- lnl + sum( log( fossilization.rate(t) ) + log( E(t,lambda,mu,phi,massExtinctionSurvivalProbabilities,rateChangeTimes) ) )
   }


   # add the extant tip age term
   lnl <- lnl + sum( nodes$tip ) * log( rho )

    
   # add the sampled ancestor age terms
   if ( sum( nodes$sampled_ancestor ) > 0 ) {
      t <- nodes$ages[ nodes$sampled_ancestor ]
      lnl <- lnl + sum (log( fossilization.rate(t) ) )
   }
   
   # add the single lineage propagation terms
   if ( length(rateChangeTimes) >= 1 ) {
      f <- function(t) nodes$age < t & nodes$age_parent > t & is.finite(nodes$age_parent)
      survivors <- Vectorize( f )
      div <- colSums( survivors(rateChangeTimes) )
      survival_prob <- massExtinctionSurvivalProbabilities[-1]
      lnl <- lnl + sum( div * log(survival_prob) ) 
      lnl <- lnl + sum( div * log( D(rateChangeTimes,lambda,mu,phi,massExtinctionSurvivalProbabilities,rateChangeTimes) ) )
   }
   

   # add the bifurcation age terms
   t <- nodes$age[ nodes$fossil_tip == FALSE & nodes$tip == FALSE & is.finite(nodes$age_parent) ]
   lnl <- lnl + sum( log( speciation.rate(t) ) + log( D(t,lambda,mu,phi,massExtinctionSurvivalProbabilities,rateChangeTimes) ) )
   

   # add the initial age term
   num_initial_lineages <- ifelse(MRCA == TRUE,2,1)
   lnl <- lnl + num_initial_lineages * log( D( max(nodes$age),lambda,mu,phi,massExtinctionSurvivalProbabilities,rateChangeTimes ) )

    # condition on survival
    if ( CONDITION == "survival" )
    {
       lnl <- lnl - num_initial_lineages * log( 1.0 - E( max(nodes$age),lambda,mu,phi,massExtinctionSurvivalProbabilities,rateChangeTimes ) )
    }
#    // condition on nTaxa
#    else if ( condition == "nTaxa" )
#    {
#        lnProbTimes -= lnProbNumTaxa( value->getNumberOfTips(), 0, process_time, true );
#    }
#
#    if ( RbMath::isFinite(lnProbTimes) == false )
#    {
#        return RbConstants::Double::nan;
#    }








  # what do we condition on?
  # did we condition on survival?
#  if ( CONDITION == "survival" || CONDITION == "taxa" )    lnl <- - tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,rho,0,PRESENT,PRESENT,log=TRUE)


  # did we condition on observing n species today
#  if ( CONDITION == "taxa" )    lnl <- lnl + tess.equations.pN.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,rho,nTaxa,0,PRESENT,SURVIVAL=TRUE,MRCA,log=TRUE)


  if (is.nan(lnl)) lnl <- -Inf

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }

  return (lnl)
}


E <- function( t, lambda, mu, phi, rho, rateChangeTimes ) {

   idx <- findInterval(t,rateChangeTimes,left.open=TRUE)+1
   
   # get the parameters
   b   <- lambda[idx]
   d   <- mu[idx]
   f   <- phi[idx]
   r   <- rho[idx]
   ti  <- ifelse( idx <= 1, 0.0, rateChangeTimes[idx-1] )

   diff <- b - d - f
   dt   <- t - ti   
    
   A <- sqrt( diff*diff + 4.0*b*f)
   B <- ( (1.0 - 2.0*((1-r)+r*ifelse(ti==0.0,0,E(ti, lambda, mu, phi, rho, rateChangeTimes))) )*b + d + f ) / A

   e <- exp(-A*dt)
   tmp <- b + d + f - A * ((1.0+B)-e*(1.0-B))/((1.0+B)+e*(1.0-B))
   
   return ( tmp / (2.0*b) )
}



D <- function( t, lambda, mu, phi, rho, rateChangeTimes ) {
    
   idx <- findInterval(t,rateChangeTimes,left.open=TRUE)+1
   
   # get the parameters
   b   <- lambda[idx]
   d   <- mu[idx]
   f   <- phi[idx]
   r   <- rho[idx]
   tmp_ti <- c(0.0,rateChangeTimes)
   ti  <- tmp_ti[idx]
   
   diff <- b - d - f
   bp   <- b*f
   dt   <- ti - t   
   
   A <- sqrt( diff*diff + 4.0*b*f)
   B <- ( (1.0 - 2.0*((1-r)+r*ifelse(ti==0.0,0,E(ti, lambda, mu, phi, rho, rateChangeTimes))) )*b + d + f ) / A
  
   e   <- exp(A*dt)
   tmp <- (1.0+B) + e*(1.0-B)
     
   return ( 4.0*e / (tmp*tmp) )
}




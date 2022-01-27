test_that("not failing", {
  data(conifers)
  
  nodes <- tess.branching.times(conifers)

  lambda <- c(0.2, 0.1, 0.3)
  mu <- c(0.1, 0.05, 0.25)
  phi <- c(0.1, 0.2, 0.05)
  r <- c(0.0, 0.0, 0.0)
  
  changetimes <- c(100, 200)
  
  tess.likelihood.ebdstp(nodes,
                        lambda = lambda, 
                        mu = mu,
                        phi = phi,
                        r = 0.0,
                        rateChangeTimesLambda = changetimes,
                        rateChangeTimesMu = changetimes,
                        rateChangeTimesPhi = changetimes,
                        samplingProbability = 1.0
  )
  
  testthat::expect_false(FALSE)
})

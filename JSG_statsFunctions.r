# Load any custom functions called in this file####

#one Regression does a regression based on an azimuth (e.g. how well does northeasting (045) explain variation in my data?) This kernal geodetic regression is analogous to a best-fit line, and thus describes a steady variation, not a sudden change in orientation.
oneRegression <-
  function(follins,
           domain,
           pValuePerms,
           directionality) {
    regression <-
      oriGeodesicRegression(
        cos((directionality) * pi / 180) * follins$northing[domain] + sin((directionality) *
                                                                            pi / 180) * follins$easting[domain],
        follins$rotation[domain],
        oriLineInPlaneGroup,
        numSteps = 10000
      )
    
    if (pValuePerms > 0) {
      RSquareds <-
        oriGeodesicRegressionPermutations(
          cos((directionality) * pi / 180) * follins$northing[domain] + sin((directionality) *
                                                                              pi / 180) * follins$easting[domain],
          follins$rotation[domain],
          numPerms = pValuePerms,
          group = oriLineInPlaneGroup
        )
      sum(RSquareds > regression$rSquared)
      p <- sum(RSquareds > regression$rSquared) / length(RSquareds)
    } else {
      p <- 0
    }
    regressionStats <- zeros(1, 5)
    regressionStats[1, 1] = directionality
    regressionStats[1, 2] = regression$error
    regressionStats[1, 3] = regression$minEigenvalue
    regressionStats[1, 4] = regression$rSquared
    regressionStats[1, 5] = p
    names(regressionStats) = c("azimuth", "error", "minEigenValue", "R^2", "P")
    regression$stats <- regressionStats
    return(regression)
  }

#regressionSweep does a series of directional regressions from 0-180°, given an increment, and also calculates the pValue for each regression. BEWARE: doing this with p-values can take days-weeks-or-months of computing time. I suggest first running the function with degreeIncrement = 10 and pValuePerms=0. This may still take a couple hours, but will give you a first order picture of whether it is interesting to proceed.
regressionSweep <-
  function(follins,
           degreeIncrement,
           domain,
           pValuePerms) {
    v <- rbind(1)
    intervals <- 180 / degreeIncrement
    for (x in 2:intervals) {
      vTemp <- x
      v <- rbind(v, vTemp)
    }
    v <- as.matrix(as.vector(v))
    # Geodesic regression of pole vs. northing in domain 4. Check that error is 0 and minEigenvalue is positive.
    regressions <- list()
    regressionStats <- zeros(nrow(v), 5)
    for (i in 1:nrow(v)) {
      regressionTemp <-
        oriGeodesicRegression(
          cos((v[i, 1] - 1) * degreeIncrement * pi / 180) * follins$northing[domain] + sin((v[i, 1] -
                                                                                              1) * degreeIncrement * pi / 180) * follins$easting[domain],
          follins$rotation[domain],
          oriLineInPlaneGroup,
          numSteps = 10000
        )
      regressions[[i]] <- regressionTemp
      regressionStats[i, 1] = (i - 1) * degreeIncrement
      regressionStats[i, 2] = regressionTemp$error
      regressionStats[i, 3] = regressionTemp$minEigenvalue
      regressionStats[i, 4] = regressionTemp$rSquared
      if (pValuePerms > 0) {
        RSquareds <-
          oriGeodesicRegressionPermutations(
            cos((v[i, 1] - 1) * degreeIncrement * pi / 180) * follins$northing[domain] + sin((v[i, 1] -1) * degreeIncrement * pi / 180) * follins$easting[domain],
            follins$rotation[domain],
            numPerms = pValuePerms,
            group = oriLineInPlaneGroup)
        length(RSquareds)
        sum(RSquareds > regressionTemp$rSquared)
        p <-
          sum(RSquareds > regressionTemp$rSquared) / length(RSquareds)
      }
      else {
        p <- "nan"
      }
      
      regressionStats[i, 5] = p
    }
    regressions[[i + 1]] <- regressionStats
    return(regressions)
  }

#take two bootstrapped mean clouds and the number of data points in each, and calculate a distribution of angular differences. Probably not the most statistically robust.
distHist <- function(bootsOne, bootsTwo, numPerms, group) {
  distances <- list()
  boots1 <- bootsOne$bootstraps
  boots2 <- bootsTwo$bootstraps
  for (i in 1:numPerms) {
    a <- round(runif(1, 1, length(boots1)), 0)
    b <- round(runif(1, 1, length(boots2)), 0)
    distances[[i]] <- oriDistance(boots1[[a]], boots2[[b]], group)
  }
  
  return(distances)
}

#Using one two-sample bootstrapped mean difference cloud, calculate the distribution of angular differences. Probably not the most statistically robust.
diffDistHist <- function(bootsDiff, numPerms, group) {
  boots <- bootsDiff$bootstraps
  distances <-
    lapply(boots, function(s)
      oriDistance(s, diag(3), group) * 180 / pi)
  return(distances)
}

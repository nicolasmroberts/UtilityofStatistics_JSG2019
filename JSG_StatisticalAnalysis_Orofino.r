# # This script was written by nicolas m. roberts, primarily using the geologyGeometry library of functions, written by Joshua R. Davis. 
# # if you incounter bugs in the code, please inform nick at nmroberts@wisc.edu

### PRELIMINARY WORK ####

# Set the working directory.
setwd("~/Desktop/20170620geologyGeometry")

# [!!!] FOR WINDOWS: Set the working directory.
# setwd("C:/users/[Insert user name]/Desktop/20170620geologyGeometry")

# Load the necessary R libraries.
source("library/all.R")

# Load some custom functions
source("JSG_statsFunctions.r")

# #[!!!] Uncomment only if you have compiled C on your system and have the C libraries loaded at the top of this script. Uncomment and run the following line only if you have compiled C. # [!!!] Markov chain Monte Carlo and Kamb contouring in equal volume plots require C compiler. Skip MCMC and equal volume Kamb lines if you do not wish to install C. Load the necessary library
# source("libraryC/all.R")

#==========================LOAD THE DATA========================================

# Load the foliation-lineation data
Follins <- geoDataFromFile("data/Follins_Ahs.csv")

# Check how many measurements there are
nrow(Follins)

#==========================PLOT THE DATA========================================

# Plot foliation-lineation locations in map view
plot(
    Follins$easting,
    Follins$northing,
    xlab = "Easting (meters)",
    ylab = "Northing (meters)"
)

# Plot the data in Equal area plot. Pole to foliation (red), Lineation (cyan)
lineEqualAreaPlotTwo(lapply(Follins$rotation, function(s)
    s[1, ]),
    lapply(Follins$rotation, function(s)
        s[2, ]))

# Plot Follins in Equal Area plots with Kamb contours (numbers are 2 sigma values)
# Poles to foliation (circles)
lineKambPlot(lapply(Follins$rotation, function(s)
    s[1, ]), c(2, 6, 10, 14, 18))
# Lineations (squares)
lineKambPlot(lapply(Follins$rotation, function(s)
    s[2, ]),
    c(2, 6, 10, 14, 18),
    shapes = c("s"))

# Plot the Follins in an Equal Volume plot (after Davis and Titus, 2017)
oriEqualVolumePlot(
    Follins$rotation,
    group = oriLineInPlaneGroup,
    backgroundColor = "white",
    curveColor = "black",
    boundaryAlpha = 0.1,
    colors = "black" ,
    axesColors = c("black", "black", "black"),
    fogStyle = "none"
)

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow("Ahsahka_EqualVol1.png", "Ahsahka_EqualVol2.png")

# [!!!] Uncomment only if you have compiled C on your system and have the C libraries loaded at the top of this script.
# Skip if you have not compiled C. Plot 6-sigma Kamb contours for the data in an equal volume plot (after Davis and Titus, 2017)
# oricKambPlot(
#     Follins$rotation,
#     group = oriLineInPlaneGroup,
#     multiple = 6,
#     backgroundColor = "white",
#     curveColor = "black",
#     boundaryAlpha = 0.1,
#     colors = "black",
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none"
# )
# 
# # If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
# afterMaximizingWindow("Ahsahka_Kamb_equalVol1.png", "Ahsahka_Kamb_EqualVol2.png")


#==========================IDENTIFY GEOGRAPHIC TRENDS========================================

# 1) Visually inspect the possibility of geographic trends in the data

# Plot foliation-lineation locations in map view, colored in gray-scale by northing. This shading will be used in the following plots.
plot(
    Follins$easting,
    Follins$northing,
    col = shades(Follins$northing),
    pch = 19
)

# Plot foliation-lineation orientations in an equal volume plot (after Davis and Titus, 2017)
oriEqualVolumePlot(
    Follins$rotation,
    group = oriLineInPlaneGroup,
    simplePoints = FALSE,
    backgroundColor = "white",
    curveColor = "black",
    boundaryAlpha = 0.1,
    colors = shades(Follins$northing),
    axesColors = c("black", "black", "black"),
    fogStyle = "none"
)

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow("Ahsahka_shadesNorthing_equalVol1.png",
                      "Ahsahka_shadesNorthing_EqualVol2.png")

# Plot foliation-lineation orientations colored by northing. Poles to foliation (circles), Lineation (squares)
lineEqualAreaPlot(
    c(
        lapply(Follins$rotation, function(s)
            s[1, ]),
        lapply(Follins$rotation, function(s)
            s[2, ])
    ),
    col = shades(Follins$northing),
    shapes = c(replicate(length(Follins$rotation), "c"), replicate(length(Follins$rotation), "s"))
)

#————————————————————————————————————————————————————————————

# 2) Run geodesic regressions to quantify any extant trends.

# Create a temporary version of the foliation lineation data and establish a geographic criteria that encompasses all the data.
FollinsAll <- Follins
AhsAllCrit <-
    FollinsAll$easting > 550000 & Follins$northing > 5140000
FollinsAll$domain <- replicate(nrow(FollinsAll), 1)
Follins$domain[AhsAllCrit] <- 1

# 2A) Perform a geodesic regression with respect to azimuths every 10 degrees on all data, using the geographic criteria defined above. Each regression asks the question: "does orientation change linearly with respect to an azimuth towards x degrees"
regressions <- regressionSweep(FollinsAll, 10, (AhsAllCrit), 0)

# Create a data frame that contains the summary information for the regressions
regressionsSum <- data.frame(regressions[[19]])

# Name the columns of the data frame
names(regressionsSum) = c("Azimuth", "Error", "Min_Eigen", "R_squared", "Pvalue")

# View the summary table for the output.
regressionsSum

# Plot the R^2 value as a function of Azimuth
plot(
    as.vector(regressionsSum$Azimuth),
    as.vector(regressionsSum$R_squared),
    xlab = "Azimuth",
    ylab = "R-squared"
)

#————————————————————————————————————————————————————————————

# 3) Run kernal regression to quantify any extant trends.

# Recategorize the symmetric copies of the data.
mu <- oriFrechetMean(Follins$rotation, group=oriLineInPlaneGroup)
Follins$rotation <- oriNearestRepresentatives(Follins$rotation, mu, group=oriLineInPlaneGroup)

# Define the bandwith for a kernal regression. If desired, uncomment the next line (it may take ~20 minutes to run). Otherwise, use the pre-calculated value, "bandwidth".
# bandwidth <- rotBandwidthForKernelRegression(Follins$northing,Follins$rotation)
bandwidth <- 0.02314421

# Perform a kernel regression with respect to northing. The kernel function is analogous to finding a best fit curve in 2D.
kernelReg <- lapply(
    seq(
        from = min(Follins$northing),
        to = max(Follins$northing),
        by = 100
    ),
    rotKernelRegression,
    Follins$northing,
    Follins$rotation,
    bandwidth,
    numSteps = 1000
)

# Filter the regression results to only keep those that give an error of 0 and minimum eigen value > 0
sapply(kernelReg, function(regr)
    regr$error)
sapply(kernelReg, function(regr) {
    regr$minEigenvalue > 0
})


# Plot the regression curve in an equal volume plot.
kernelRegCurve <- lapply(kernelReg, function(regr)
    regr$r)
oriEqualVolumePlot(
    kernelRegCurve,
    group=oriLineInPlaneGroup,
    backgroundColor = "white",
    curveColor = "black",
    boundaryAlpha = 0.1,
    colors = "black",
    axesColors = c("black", "black", "black"),
    fogStyle = "none",
    simplePoints = TRUE
)

# Compute the R^2 for the kernel regression.
KernRsquared <-
    rotRsquaredForKernelRegression(Follins$northing, Follins$rotation, bandwidth, numSteps =
                                       1000)

# View the R^2 value
KernRsquared

# Do a permutation test for significance. This takes a significant amount of time (>>1 hr)
rSquareds <-
    rotKernelRegressionPermutations(Follins$northing, Follins$rotation, bandwidth, numPerms =
                                        1000)
length(rSquareds)
p <- sum(rSquareds > KernRsquared$rSquared)
p


#==========================DEFINE GEOGRAPHIC DOMAINS========================================

# Define geographic criteria, based on sampling area, with the proposed shear zone boundary into account.
domain1Crit <- Follins$easting < 560000 & Follins$northing < 5158000
domain2Crit <- (Follins$easting >= 560000 & Follins$northing < 5163000) | Follins$easting >= 565000
domain3Crit <- !(domain1Crit | domain2Crit)

# # Create a new column in the dataframe in which to store the domain information
Follins$domain <- replicate(nrow(Follins), 1)

# Classify the foliation-only dataset by domain
Follins$domain[domain1Crit] <- 1
Follins$domain[domain2Crit] <- 2
Follins$domain[domain3Crit] <- 3

# Plot the locations of the foliation-lineation data in map view, each domain a different color. Domain 1 (Red), Domain 2 (Green), Domain 3 (Blue)
plot(
    x = Follins$easting,
    y = Follins$northing,
    xlab = "Easting (meters)",
    ylab = "Northing (meters)",
    col = hues(Follins$domain),
    pch = 19
)

#==========================IDENTIFY GEOGRAPHIC TRENDS WITHIN DOMAINS========================================

# 1) Make sure each domain is roughly unimodal
#DOMAIN 1, n = 23.
oriEqualVolumePlot(
    Follins$rotation[domain1Crit],
    oriLineInPlaneGroup,
    backgroundColor = "white",
    curveColor = "black",
    boundaryAlpha = 0.1,
    colors = "black",
    axesColors = c("black", "black", "black"),
    fogStyle = "none"
)

#[!!!] Uncomment only if you have compiled C on your system and have the C libraries loaded at the top of this script.
# oricKambPlot(
#     Follins$rotation[domain1Crit],
#     oriLineInPlaneGroup,
#     backgroundColor = "white",
#     curveColor = "black",
#     boundaryAlpha = 0.1,
#     colors = "black",
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none"
# )

#DOMAIN 2, n = 14.
oriEqualVolumePlot(
    Follins$rotation[domain2Crit],
    oriLineInPlaneGroup,
    backgroundColor = "white",
    curveColor = "black",
    boundaryAlpha = 0.1,
    colors = "black",
    axesColors = c("black", "black", "black"),
    fogStyle = "none"
)

#[!!!] Uncomment only if you have compiled C on your system and have the C libraries loaded at the top of this script.
# oricKambPlot(
#     Follins$rotation[domain2Crit],
#     oriLineInPlaneGroup,
#     backgroundColor = "white",
#     curveColor = "black",
#     boundaryAlpha = 0.1,
#     colors = "black",
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none"
# )

#DOMAIN 3, n = 32.
oriEqualVolumePlot(
    Follins$rotation[domain2Crit],
    oriLineInPlaneGroup,
    backgroundColor = "white",
    curveColor = "black",
    boundaryAlpha = 0.1,
    colors = "black",
    axesColors = c("black", "black", "black"),
    fogStyle = "none"
)
oriEqualVolumePlot(Follins$rotation[domain2Crit],
                   oriLineInPlaneGroup,
                   col = hues(Follins$easting[domain2Crit]))

#[!!!] Uncomment only if you have compiled C on your system and have the C libraries loaded at the top of this script.
# oricKambPlot(
#     Follins$rotation[domain2Crit],
#     oriLineInPlaneGroup,
#     backgroundColor = "white",
#     curveColor = "black",
#     boundaryAlpha = 0.1,
#     colors = "black",
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none"
# )


# Plot all the foliation data in an equal volume plot (after Davis and Titus, 2017), colored by domain. Domain 1 (red), Domain 2 (green), Domain 3 (blue)
oriEqualVolumePlot(
    Follins$rotation,
    oriLineInPlaneGroup,
    col = hues(Follins$domain),
    backgroundColor = "white",
    curveColor = "black",
    boundaryAlpha = 0.2,
    axesColors = c("black", "black", "black"),
    fogStyle = "none"
)

#————————————————————————————————————————————————————————————

# 2) Look for geographic patterns in plots.
#DOMAIN 1, n = 23. maybe a very faint rainbow--more like clumping of the dark blues together.
oriEqualVolumePlot(Follins$rotation[domain1Crit],
                   oriLineInPlaneGroup,
                   col = shades(Follins$northing[domain1Crit]))

#DOMAIN 2, n = 14. Part of the problem with small sample sizes--you can't really know whether the datapoints are independent...
oriEqualVolumePlot(Follins$rotation[domain2Crit],
                   oriLineInPlaneGroup,
                   col = shades(Follins$northing[domain2Crit]))

#DOMAIN 3, n = 32. No rainbow.
oriEqualVolumePlot(Follins$rotation[domain2Crit],
                   oriLineInPlaneGroup,
                   col = shades(Follins$easting[domain2Crit]))

#————————————————————————————————————————————————————————————

# 3) From visual inspection of the above plots, Domain 1 and Domain 3 potentially have geographic dependency. Perform regressions on these two domains.

#DOMAIN 1, n = 23. Perform a series of Geodesic Regressions for every 10° azimuth, with p-values not calculated.
regressionsD1 <- regressionSweep(Follins, 10, domain1Crit, 0)

# Create a dataframe with the summary information for the regressions
regressionsD1Sum <- data.frame(regressionsD1[[19]])

# Name the columns of the data frame
names(regressionsD1Sum) = c("Azimuth", "Error", "Min_Eigen", "R_squared", "Pvalue")

# View the summary table for the output.
regressionsD1Sum

# Plot the R^2 value as a function of Azimuth
plot(
    as.vector(regressionsD1Sum$Azimuth),
    as.vector(regressionsD1Sum$R_squared),
    xlab = "Azimuth",
    ylab = "R-squared"
)

#————————————————————————————————————————————————————————————

# Define the Azimuth with the highest R-squared value. (In this case, 130)
azimuthD1 <-
    Follins$easting[domain1Crit] * sin(130 * degree) + Follins$northing[domain1Crit] *
    cos(130 * degree)

# DOMAIN 1, n = 23. Define the bandwith for a kernal regression. If desired, uncomment the next line (it may take ~20 minutes to run). Otherwise, use the pre-calculated value, "bandwidthD1".
# bandwidthD1 <- rotBandwidthForKernelRegression(azimuthD1,Follins$rotation[domain1Crit])
bandwidthD1 <- 0.02202349

#Run the kernel regression
kernelRegD1 <- lapply(
    seq(
        from = min(azimuthD1),
        to = max(azimuthD1),
        by = 100
    ),
    rotKernelRegression,
    azimuthD1,
    Follins$rotation[domain1Crit],
    bandwidthD1,
    numSteps = 1000
)
sapply(kernelRegD1, function(regr)
    regr$error)
sapply(kernelRegD1, function(regr) {
    regr$minEigenvalue > 0
})

# Plot the regression curve
kernelRegCurveD1 <- lapply(kernelRegD1, function(regr)
    regr$r)
oriEqualVolumePlot(Follins$rotation[domain1Crit],
                   curves = list(kernelRegCurveD1),
                   simplePoints = TRUE, oriLineInPlaneGroup)

# Compute the R^2.
KernRsquaredD1 <-
    rotRsquaredForKernelRegression(azimuthD1, Follins$rotation[domain1Crit], bandwidthD1, numSteps =
                                       1000)
KernRsquaredD1

# Do a permutation test for significance.
rSquaredsD1 <-
    rotKernelRegressionPermutations(azimuthD1, Follins$rotation[domain1Crit], bandwidthD1, numPerms =
                                        1000)
length(rSquareds)
pD1 <- sum(rSquareds > KernRsquared$rSquared)

#the p-value for this regression.
pD1

#————————————————————————————————————————————————————————————

#Now do the same thing for Domain 3

#DOMAIN 3, n = 32. A series of Geodesic Regressions for every 10° azimuth, with p-values not calculated.
regressionsD3 <- regressionSweep(Follins, 10, domain3Crit, 0)

# Create a dataframe with the summary information for the regressions
regressionsD3Sum <- data.frame(regressionsD3[[19]])

# Name the columns of the data frame
names(regressionsD3Sum) = c("Azimuth", "Error", "Min_Eigen", "R_squared", "Pvalue")

# View the summary table for the output.
regressionsD3Sum

# Plot the R^2 value as a function of Azimuth
plot(
    as.vector(regressionsD3Sum$Azimuth),
    as.vector(regressionsD3Sum$R_squared),
    xlab = "Azimuth",
    ylab = "R-squared"
)

#————————————————————————————————————————————————————————————

# Define the Azimuth with the highest R-squared value. (In this case, 60)
azimuthD3 <-
    Follins$easting[domain3Crit] * sin(60 * degree) + Follins$northing[domain3Crit] *
    cos(60 * degree)

#DOMAIN 3, n = 32. Define the bandwith for a kernal regression. If desired, uncomment the next line (it may take ~20 minutes to run). Otherwise, use the pre-calculated value, "bandwidthD3".
# bandwidthD3 <-
#     rotBandwidthForKernelRegression(Follins$northing[domain3Crit], Follins$rotation[domain3Crit])
bandwidthD3 <- 0.1331092

#Run the kernel regression
kernelRegD3 <- lapply(
    seq(
        from = min(azimuthD3),
        to = max(azimuthD3),
        by = 100
    ),
    rotKernelRegression,
    azimuthD3,
    Follins$rotation[domain3Crit],
    bandwidthD3,
    numSteps = 1000
)
sapply(kernelRegD3, function(regr)
    regr$error)
sapply(kernelRegD3, function(regr) {
    regr$minEigenvalue > 0
})

# Plot the regression curve
kernelRegCurveD3 <- lapply(kernelRegD3, function(regr)
    regr$r)
oriEqualVolumePlot(Follins$rotation[domain3Crit],
                   curves = list(kernelRegCurveD3),
                   simplePoints = TRUE,
                   oriLineInPlaneGroup)

# Compute the R^2.
KernRsquaredD3 <-
    rotRsquaredForKernelRegression(azimuthD3, Follins$rotation[domain3Crit], bandwidthD3, numSteps = 1000)
KernRsquaredD3

# Do a permutation test for significance.
rSquaredsD3 <-
    rotKernelRegressionPermutations(Follins$northing[domain3Crit], Follins$rotation[domain3Crit], bandwidthD3, numPerms = 1000)
length(rSquareds)
pD3 <- sum(rSquareds > KernRsquaredD3$rSquared)

#the p-value for this regression.
pD3

##==========================STATISTICAL DESCRIPTORS========================================

# 1) Domain 1, n = 23

# Frechet mean (minimizes the Frechet variance).
FrechetMeanDom1 <-
    oriFrechetMean(Follins$rotation[domain1Crit], oriLineInPlaneGroup)
geoStrikeDipRakeDegFromRotation(FrechetMeanDom1)
FrechetVarDom1 <-
    oriVariance(Follins$rotation[domain1Crit], FrechetMeanDom1, oriLineInPlaneGroup)

# Plot the FrechetMean in an Equal Angle plot
FrechetCurvesDom1 <-
    lapply(Follins$rotation[domain1Crit], function(r)
        rotGeodesicPoints(FrechetMeanDom1, r, 10))
oriEqualAnglePlot(points = Follins$rotation[domain1Crit], curves = FrechetCurvesDom1, group = oriLineInPlaneGroup)

# Print the Strike, Dip, Rake of the Frechet mean.
geoStrikeDipRakeDegFromRotation(FrechetMeanDom1)

# Print the Frechet variance.
FrechetVarDom1

#————————————————————————————————————————————————————————————

# 2) Domain 2, n = 14
#Frechet mean minimizes the Frechet variance.
FrechetMeanDom2 <-
    oriFrechetMean(Follins$rotation[domain2Crit], oriLineInPlaneGroup)
geoStrikeDipRakeDegFromRotation(FrechetMeanDom2)
FrechetVarDom2 <-
    oriVariance(Follins$rotation[domain2Crit], FrechetMeanDom2, oriLineInPlaneGroup)

# Plot the FrechetMean in an Equal Angle plot
FrechetCurvesDom2 <-
    lapply(Follins$rotation[domain2Crit], function(r)
        rotGeodesicPoints(FrechetMeanDom2, r, 10))
oriEqualAnglePlot(points = Follins$rotation[domain2Crit], curves = FrechetCurvesDom2, oriLineInPlaneGroup)

# Print the Strike, Dip, Rake of the Frechet mean.
geoStrikeDipRakeDegFromRotation(FrechetMeanDom2)

# Print the Frechet variance
FrechetVarDom2

#————————————————————————————————————————————————————————————

# 3) Domain 3, n = 32
#Frechet mean minimizes the Frechet variance.
FrechetMeanDom3 <-
    oriFrechetMean(Follins$rotation[domain3Crit], oriLineInPlaneGroup)
geoStrikeDipRakeDegFromRotation(FrechetMeanDom3)
FrechetVarDom3 <-
    oriVariance(Follins$rotation[domain3Crit], FrechetMeanDom3, oriLineInPlaneGroup)

# Plot the FrechetMean in an Equal Angle plot
FrechetCurvesDom3 <-
    lapply(Follins$rotation[domain3Crit], function(r)
        rotGeodesicPoints(FrechetMeanDom3, r, 10))
oriEqualAnglePlot(points = Follins$rotation[domain3Crit], curves = FrechetCurvesDom3, oriLineInPlaneGroup)

# Print the Strike, Dip, Rake of the Frechet mean.
geoStrikeDipRakeDegFromRotation(FrechetMeanDom3)

# Print the Frechet variance.
FrechetVarDom3

#————————————————————————————————————————————————————————————

# 4) Dispersion of the data

#The standard way to get at dispersion is to compute the maximum matrix fisher likelihood. The matrix fisher distribution comprises a kind of "mean" and a positive definite symmetric matrix, which characterises the anisotropic dispersion.

# Redefine the rotations to be within one of the symmetric copies
mu <- oriFrechetMean(Follins$rotation, group = oriLineInPlaneGroup)
Follins$rotation <-
    oriNearestRepresentatives(Follins$rotation, mu, oriLineInPlaneGroup)

# Northern Domain, n = 16. Fisher maximum likelihood
mleDom1 <- rotFisherMLE(Follins$rotation[domain1Crit])
mleDom1
eigen(mleDom1$kHat, symmetric = TRUE, only.value = TRUE)$values

#Central Domain, n = 34. Fisher maximum likelihood
mleDom2 <- rotFisherMLE(Follins$rotation[domain2Crit])
mleDom2
eigen(mleDom2$kHat, symmetric = TRUE, only.value = TRUE)$values

#Southern Domain, n = 79. Fisher maximum likelihood
mleDom3 <- rotFisherMLE(Follins$rotation[domain3Crit])
mleDom3
eigen(mleDom3$kHat, symmetric = TRUE, only.value = TRUE)$values

#==========================INFERENCE========================================

# Here, we can perform some statistical techniques to use information about the samples (each domain) to compute a probability cloud of the mean of the population(s) from which those samples were taken.
# From numerical work in Davis and Titus (2017), MCMC will work the best for small sample sizes that have the matrix fisher anisotropy (eigenvalues) of (large, large, small). All four domains have this anisotropy, and have too few data points to use the (computationally quicker) bootstrapping method.

#We'll do both and compare.

# Compute the bootstrap cloud of means
Follins_domain1Boots <-
    oriBootstrapInference(Follins$rotation[domain1Crit], 10000, oriLineInPlaneGroup)
Follins_domain2Boots <-
    oriBootstrapInference(Follins$rotation[domain2Crit], 10000, oriLineInPlaneGroup)
Follins_domain3Boots <-
    oriBootstrapInference(Follins$rotation[domain3Crit], 10000, oriLineInPlaneGroup)

# [!!!] Requires C libraries to run. 
# [!!!] Compute the Markov chain Monte Carlo cloud of means
# [!!!]Uncomment lines if you have compiled C on your computer and have loaded the C libraries at the top of this script.

# Follins_domain1MCMC <-
#     oricWrappedTrivariateNormalMCMCInference(Follins$rotation[domain1Crit],
#                                              group = oriLineInPlaneGroup,
#                                              numCollection = 100)
# Follins_domain2MCMC <-
#     oricWrappedTrivariateNormalMCMCInference(Follins$rotation[domain2Crit],
#                                              group = oriLineInPlaneGroup,
#                                              numCollection = 100)
# Follins_domain3MCMC <-
#     oricWrappedTrivariateNormalMCMCInference(Follins$rotation[domain3Crit],
#                                              group = oriLineInPlaneGroup,
#                                              numCollection = 100)

#————————————————————————————————————————————————————————————

# 1) Northern v. Central domains

# A) Using MCMC

# [!!!] Requires C libraries to run. 
# Construct the 95% confidence ellipsoids from small triangles
# tris1_MCMC <-
#     rotEllipsoidTriangles(
#         Follins_domain1MCMC$mBar,
#         Follins_domain1MCMC$leftCovarInv,
#         Follins_domain1MCMC$q095,
#         numNonAdapt = 4
#     )
# tris2_MCMC <-
#     rotEllipsoidTriangles(
#         Follins_domain2MCMC$mBar,
#         Follins_domain2MCMC$leftCovarInv,
#         Follins_domain2MCMC$q095,
#         numNonAdapt = 4
#     )
# 
# # Plot the MCMC comparison in an equal Volume plot, with 95% confidence ellipsoids
# oriEqualAnglePlot(
#     points = c(Follins_domain1MCMC$ms, Follins_domain2MCMC$ms),
#     boundaryAlpha = .1,
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none",
#     background = "white",
#     triangles = c(tris1_MCMC, tris2_MCMC),
#     simplePoints = TRUE,
#     colors = c(replicate(length(
#         Follins_domain1MCMC$ms
#     ), "black"), replicate(
#         length(Follins_domain2MCMC$ms), "orange"
#     )),
#     group = oriTrivialGroup
# )
# 
# # If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
# afterMaximizingWindow("MCMC_Ahs_d1d2_1.png", "MCMC_WestMt_d1d2_2.png")
# 
# # Plot the MCMC comparison in an Equal Area plot
# lineEqualAreaPlotTwo(c(
#     lapply(Follins_domain1MCMC$ms, function(s)
#         s[1, ]),
#     lapply(Follins_domain1MCMC$ms, function(s)
#         s[2, ])
# ),
# c(
#     lapply(Follins_domain2MCMC$ms, function(s)
#         s[1, ]),
#     lapply(Follins_domain2MCMC$ms, function(s)
#         s[2, ])
# ))

#....................................................................

# B) Using bootstrapping

# Construct the 95% confidence ellipsoids from small triangles
tris1_Boot <-
    rotEllipsoidTriangles(
        Follins_domain1Boots$center,
        Follins_domain1Boots$leftCovarInv,
        Follins_domain1Boots$q095,
        numNonAdapt = 4
    )
tris2_Boot <-
    rotEllipsoidTriangles(
        Follins_domain2Boots$center,
        Follins_domain2Boots$leftCovarInv,
        Follins_domain2Boots$q095,
        numNonAdapt = 4
    )

# Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids
rotEqualAnglePlot(
    points = c(
        Follins_domain1Boots$bootstraps,
        Follins_domain2Boots$bootstraps
    ),
    triangles = c(tris1_Boot, tris2_Boot),
    boundaryAlpha = .1,
    axesColors = c("black", "black", "black"),
    fogStyle = "none",
    background = "white",
    simplePoints = TRUE,
    colors = c(replicate(
        length(Follins_domain1Boots$bootstraps), "black"
    ), replicate(
        length(Follins_domain2Boots$bootstraps), "orange"
    ))
)

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow("Boots_Ahs_d1d2_1.png", "Boots_Ahs_d1d2_2.png")

# Plot the bootstrap comparison in an Equal Area plot
lineEqualAreaPlotTwo(c(
    lapply(Follins_domain1Boots$bootstraps, function(s)
        s[1, ]),
    lapply(Follins_domain1Boots$bootstraps, function(s)
        s[2, ])
),
c(
    lapply(Follins_domain2Boots$bootstraps, function(s)
        s[1, ]),
    lapply(Follins_domain2Boots$bootstraps, function(s)
        s[2, ])
))

#————————————————————————————————————————————————————————————

# 2) Northern vs. Southern domains

# A) Using MCMC

# [!!!] Requires C libraries to run. 
# Construct the 95% confidence ellipsoids from small triangles
# tris1_MCMC <-
#     rotEllipsoidTriangles(
#         Follins_domain1MCMC$mBar,
#         Follins_domain1MCMC$leftCovarInv,
#         Follins_domain1MCMC$q095,
#         numNonAdapt = 5
#     )
# tris3_MCMC <-
#     rotEllipsoidTriangles(
#         Follins_domain3MCMC$mBar,
#         Follins_domain3MCMC$leftCovarInv,
#         Follins_domain3MCMC$q095,
#         numNonAdapt = 5
#     )
# 
# # Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids
# oriEqualAnglePlot(
#     points = c(Follins_domain1MCMC$ms, Follins_domain3MCMC$ms),
#     triangles = c(tris1_MCMC, tris3_MCMC),
#     boundaryAlpha = 0.1,
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none",
#     background = "white",
#     simplePoints = TRUE,
#     colors = c(replicate(length(
#         Follins_domain1MCMC$ms
#     ), "black"), replicate(length(
#         Follins_domain3MCMC$ms
#     ), "blue")),
#     group = oriTrivialGroup
# )
# 
# # If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
# afterMaximizingWindow("MCMC_Ahs_d1d3_1.png", "MCMC_Ahsd1d3_2.png")
# 
# # Plot the MCMC comparison in an Equal Area plot
# lineEqualAreaPlotTwo(c(
#     lapply(Follins_domain1MCMC$ms, function(s)
#         s[1, ]),
#     lapply(Follins_domain1MCMC$ms, function(s)
#         s[2, ])
# ),
# c(
#     lapply(Follins_domain3MCMC$ms, function(s)
#         s[1, ]),
#     lapply(Follins_domain3MCMC$ms, function(s)
#         s[2, ])
# ))

#....................................................................

# B) Using bootstrapping
# Construct the 95% confidence ellipsoids from small triangles
tris1_Boot <-
    rotEllipsoidTriangles(
        Follins_domain1Boots$center,
        Follins_domain1Boots$leftCovarInv,
        Follins_domain1Boots$q095,
        numNonAdapt = 4
    )
tris3_Boot <-
    rotEllipsoidTriangles(
        Follins_domain3Boots$center,
        Follins_domain3Boots$leftCovarInv,
        Follins_domain3Boots$q095,
        numNonAdapt = 4
    )

# Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids.
rotEqualAnglePlot(
    points = c(
        Follins_domain1Boots$bootstraps,
        Follins_domain3Boots$bootstraps
    ),
    triangles = c(tris1_Boot, tris3_Boot),
    boundaryAlpha = .1,
    axesColors = c("black", "black", "black"),
    fogStyle = "none",
    background = "white",
    simplePoints = TRUE,
    colors = c(replicate(
        length(Follins_domain1Boots$bootstraps), "black"
    ),
    replicate(
        length(Follins_domain3Boots$bootstraps), "blue"
    ))
)

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow("Boots_Ahs_d1d3_1.png", "Boots_Ahs_d1d3_2.png")

# Plot the bootstrap comparison in an Equal Area plot
lineEqualAreaPlotTwo(c(
    lapply(Follins_domain1Boots$bootstraps, function(s)
        s[1, ]),
    lapply(Follins_domain1Boots$bootstraps, function(s)
        s[2, ])
),
c(
    lapply(Follins_domain3Boots$bootstraps, function(s)
        s[1, ]),
    lapply(Follins_domain3Boots$bootstraps, function(s)
        s[2, ])
))



#————————————————————————————————————————————————————————————

# 3) Central vs. Southern domains

# A) Using MCMC

# [!!!] Requires C libraries to run. 
# Construct the 95% confidence ellipsoids from small triangles
# tris2_MCMC <-
#     rotEllipsoidTriangles(
#         Follins_domain2MCMC$mBar,
#         Follins_domain2MCMC$leftCovarInv,
#         Follins_domain2MCMC$q095,
#         numNonAdapt = 4
#     )
# tris3_MCMC <-
#     rotEllipsoidTriangles(
#         Follins_domain3MCMC$mBar,
#         Follins_domain3MCMC$leftCovarInv,
#         Follins_domain3MCMC$q095,
#         numNonAdapt = 4
#     )
# 
# # Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids
# oriEqualAnglePlot(
#     points = c(Follins_domain2MCMC$ms, Follins_domain3MCMC$ms),
#     triangles = c(tris2_MCMC, tris3_MCMC),
#     boundaryAlpha = 0.1,
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none",
#     background = "white",
#     simplePoints = TRUE,
#     colors = c(replicate(
#         length(Follins_domain2MCMC$ms), "orange"
#     ), replicate(length(
#         Follins_domain3MCMC$ms
#     ), "blue")),
#     group = oriTrivialGroup
# )
# 
# # If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
# afterMaximizingWindow("MCMC_Ahs_d2d3_1.png", "MCMC_Ahs_d2d3_2.png")
# 
# # Plot the MCMC comparison in an Equal Area plot
# lineEqualAreaPlotTwo(c(
#     lapply(Follins_domain2MCMC$ms, function(s)
#         s[1, ]),
#     lapply(Follins_domain2MCMC$ms, function(s)
#         s[2, ])
# ),
# c(
#     lapply(Follins_domain3MCMC$ms, function(s)
#         s[1, ]),
#     lapply(Follins_domain3MCMC$ms, function(s)
#         s[2, ])
# ))

#....................................................................


# B) Using bootstrapping

# Construct the 95% confidence ellipsoids from small triangles
tris2_Boot <-
    rotEllipsoidTriangles(
        Follins_domain2Boots$center,
        Follins_domain2Boots$leftCovarInv,
        Follins_domain2Boots$q095,
        numNonAdapt = 4
    )
tris3_Boot <-
    rotEllipsoidTriangles(
        Follins_domain3Boots$center,
        Follins_domain3Boots$leftCovarInv,
        Follins_domain3Boots$q095,
        numNonAdapt = 4
    )

# Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids
oriEqualAnglePlot(
    points = c(
        Follins_domain2Boots$bootstraps,
        Follins_domain3Boots$bootstraps
    ),
    triangles = c(tris2_Boot, tris3_Boot),
    boundaryAlpha = 0.1,
    axesColors = c("black", "black", "black"),
    fogStyle = "none",
    background = "white",
    simplePoints = TRUE,
    colors = c(replicate(
        length(Follins_domain1Boots$bootstraps), "orange"
    ), replicate(
        length(Follins_domain3Boots$bootstraps), "blue"
    )),
    group = oriTrivialGroup
)

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow("Boots_Ahs_d2d3_1.png", "Boots_Ahs_d2d3_2.png")

# Plot the bootstrap comparison in an Equal Area plot
lineEqualAreaPlotTwo(c(
    lapply(Follins_domain2Boots$bootstraps, function(s)
        s[1, ]),
    lapply(Follins_domain2Boots$bootstraps, function(s)
        s[2, ])
),
c(
    lapply(Follins_domain3Boots$bootstraps, function(s)
        s[1, ]),
    lapply(Follins_domain3Boots$bootstraps, function(s)
        s[2, ])
))

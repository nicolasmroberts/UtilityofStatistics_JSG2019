# # This script was written by nicolas m. roberts, primarily using the geologyGeometry library of functions, written by Joshua R. Davis. 
# # if you incounter bugs in the code, please inform nick at nmroberts@wisc.edu

### PRELIMINARY WORK ####

# FOR MAC: Set the working directory.
setwd("~/Desktop/20170620geologyGeometry")

# [!!!] FOR WINDOWS: Set the working directory.
# setwd("C:/users/[Insert user name]/Desktop/20170620geologyGeometry")

# Load the necessary R libraries.
source("library/all.R")

# Load some custom functions.
source("JSG_statsFunctions.r")

# [!!!] Markov chain Monte Carlo and Kamb contouring in equal volume plots require C compiler. Skip MCMC and equal volume Kamb lines if you do not wish to install C. Load the necessary library
# source("libraryC/all.R")

#==========================LOAD THE DATA========================================

# 1) Foliation only
Fols <- geoDataFromFile("data/Fols_WestMt.csv")

# Plot foliation locations in map view
plot(Fols$easting, Fols$northing, xlab = "Easting (meters)", ylab = "Northing (meters)")

# Check how many measurments there are
nrow(Fols)

#————————————————————————————————————————————————————————————

# 2) Foliation-lineation pairs
Follins <- geoDataFromFile("data/Follins_WestMt.csv")

# Plot foliation-lineation locations in map view
plot(
    Follins$easting,
    Follins$northing,
    xlab = "Easting (meters)",
    ylab = "Northing (meters)"
)

# Check how many measurements there are
nrow(Follins)



#==========================PART I: DIRECTIONAL STATISTICS ON FOLIATIONS========================================
#==========================DEFINE GEOGRAPHIC DOMAINS FOR FOLIATION-ONLY DATA=================================

#(After Braudy et al., 2017)

# Define the geographic criteria that defines the different domains. The numbers are UTM coordinates
Fols_northCrit <- Fols$northing > 4922736
Fols_centerCrit <- Fols$northing < 4922736 & Fols$northing > 4919205
Fols_southCrit <- Fols$northing < 4919205

# Create a new column in the dataframe in which to store the domain information
Fols$domain <- replicate(nrow(Fols), 1)

# Classify the foliation-only dataset by domain
Fols$domain[Fols_northCrit] <- 1
Fols$domain[Fols_centerCrit] <- 2
Fols$domain[Fols_southCrit] <- 3

# Plot the locations of the foliation-only data in map view, each domain a different color.
plot(
    x = Fols$easting,
    y = Fols$northing,
    xlab = "Easting (meters)",
    ylab = "Northing (meters)",
    col = hues(Fols$domain),
    pch = 19
)

# Return the number of datapoints in each domain.
length(Fols$domain[Fols_northCrit])
length(Fols$domain[Fols_centerCrit])
length(Fols$domain[Fols_southCrit])

#==========================DEFINE GEOGRAPHIC DOMAINS FOR FOLIATION-LINEATION DATA========================================

# (After Braudy et al., 2017)

# Define the geographic criteria that defines the different domains
Follins_northCrit <-  Follins$northing > 4922736
Follins_centerCrit <- Follins$northing < 4922736 & Follins$northing > 4919205
Follins_southCrit <-  Follins$northing < 4919205

# Create a new column in the dataframe in which to store the domain information
Follins$domain <- replicate(nrow(Follins), 1)

# Classify the foliation-lineation dataset by domain
Follins$domain[Follins_northCrit] <- 1
Follins$domain[Follins_centerCrit] <- 2
Follins$domain[Follins_southCrit] <- 3

# Plot the locations of the foliation-lineation data in map view, each domain a different color.
plot(
    x = Follins$easting,
    y = Follins$northing,
    xlab = "Easting (meters)",
    ylab = "Northing (meters)",
    col = hues(Follins$domain),
    pch = 19
)

# Return the number of datapoints in each domain.
length(Follins$domain[Follins_northCrit])
length(Follins$domain[Follins_centerCrit])
length(Follins$domain[Follins_southCrit])

#==========================FOLIATION-ONLY STATISTCAL TREATMENT========================================

# 1) Some parametric two-sample tests. The null hypothesis for all tests is that the two domains being tested come from the same population.

# Three Wellner tests (Wellner, 1979), one for each pair of domains. Each test is based on 10,000 permutations. This tests not only if the samples come from the populations with the same mean, but also with the same dispersion.
lineWellnerInference(Fols$pole[Fols_northCrit], Fols$pole[Fols_southCrit], 10000)
lineWellnerInference(Fols$pole[Fols_northCrit], Fols$pole[Fols_centerCrit], 10000)
lineWellnerInference(Fols$pole[Fols_centerCrit], Fols$pole[Fols_southCrit], 10000)

# Three Watson tests that assume large sample size (Mardia and Jupp, 2000), one for each pair of domains.
lineLargeMultiSampleWatsonInference(list(Fols$pole[Fols_northCrit],  Fols$pole[Fols_southCrit]))
lineLargeMultiSampleWatsonInference(list(Fols$pole[Fols_northCrit],  Fols$pole[Fols_centerCrit]))
lineLargeMultiSampleWatsonInference(list(Fols$pole[Fols_centerCrit], Fols$pole[Fols_southCrit]))

# Three Watson tests that assume tightly concentrated datasets (Mardia and Jupp, 2000), one for each pair of domains.
lineConcentratedMultiSampleWatsonInference(list(Fols$pole[Fols_northCrit],  Fols$pole[Fols_southCrit]))
lineConcentratedMultiSampleWatsonInference(list(Fols$pole[Fols_northCrit],  Fols$pole[Fols_centerCrit]))
lineConcentratedMultiSampleWatsonInference(list(Fols$pole[Fols_centerCrit], Fols$pole[Fols_southCrit]))

#————————————————————————————————————————————————————————————

# 2) Non-parametric bootstrapping

# Perform the bootstrapping routine for each domain. Each bootstrapped dataset is based on 10,000 iterations.
Fols_northBoots <-  lineBootstrapInference(Fols$pole[Fols_northCrit], 10000, numPoints = 50)
Fols_centerBoots <- lineBootstrapInference(Fols$pole[Fols_centerCrit], 10000, numPoints = 50)
Fols_southBoots <-  lineBootstrapInference(Fols$pole[Fols_southCrit], 10000, numPoints = 50)

# Plot data for each domain. Northern (red), central (green), southern (blue)
lineEqualAreaPlotThree(
                       Fols$pole[Fols_northCrit],
                       Fols$pole[Fols_centerCrit],
                       Fols$pole[Fols_southCrit]
)

# Plot the bootstrapped mean clouds for each domain.  Northern (red), central (green), southern (blue)
lineEqualAreaPlotThree(
    Fols_northBoots$bootstraps,
    Fols_centerBoots$bootstraps,
    Fols_southBoots$bootstraps
)

# Plot the bootstrapped mean with 95% confidence ellipses superimposed.  Northern (red), central (green), southern (blue)
lineEqualAreaPlotThree(
    list(Fols_northBoots$center),
    list(Fols_centerBoots$center),
    list(Fols_southBoots$center),
    curves = list(
        Fols_northBoots$points,
        Fols_centerBoots$points,
        Fols_southBoots$points
    )
)

#————————————————————————————————————————————————————————————

# 3) Compute the rotation between the north and south domains.

# Create empty lists to store data
northSouthDiff <- list()

# Assign variables for the orientations of all bootstrapped means.
northB <- Fols_northBoots$bootstraps
southB <- Fols_southBoots$bootstraps

# Calculate the smallest possible rotations that bring means from the northern bootstrap cloud to the southern bootstrap clouds.
count = 1
for (i in 1:300) {
    for (n in i:300) {
        northSouthDiff$rotation[[count]] = rotSmallestRotationFromTwoLines(northB[[i]], southB[[n]])
        if (n < 90000) {
            count = count + 1
        }
    }
}

# Display the raw results of the preceting for loop
northSouthDiff

# Plot these rotations in an equal volume plot
oriEqualVolumePlot(
                   northSouthDiff$rotation,
                   group = oriLineInPlaneGroup,
                   simplePoints = TRUE
)

# Compute the axis and rotation amount from all the rotations
Fols_angleAxis <- lapply(northSouthDiff$rotation, function(s) rotAxisAngleFromMatrix(s))


#Identify the rotation axes
rotationAxes <- lapply(Fols_angleAxis, function(s) c(s[1], s[2], s[3]))

rotationAxisMean <- lineProjectedMean(rotationAxes)
inf <- rayMahalanobisInference(rotationAxes, rotationAxisMean, numPoints=100)

# Plot the rotation axes in an equal area plot
lineEqualAreaPlot(
                  rotationAxes, 
                  shapes = '.',
)

# Plot the rotation axis mean orientation with the 95% confidence ellipse
lineEqualAreaPlot(
                  list(inf$center),
                  curves = list(inf$points)
)


# Plot a histogram of rotation amount between the northern domain and the southern domain, in degrees.
hist(
     as.numeric(lapply(Fols_angleAxis, function(s) s[4] * 180 / pi)),
     30,
     xlab = "Angular Distance (degrees)",
     ylab = "Frequency",
     main = "Angular distances, Northern to Southern"
)

# Compute the mean and 1-sigma standard deviation of the rotation amount, in degrees
mean(sapply(Fols_angleAxis, function(s) s[4] * 180 / pi))
sd(sapply(Fols_angleAxis, function(s) s[4] * 180 / pi))

# Compute the mean trend and plunge in degrees
Fols_axes <- lapply(Fols_angleAxis, function(s) c(s[1], s[2], s[3]))
Fols_meanAxisDeg <- geoTrendPlungeDegFromCartesian(lower(lineProjectedMean(Fols_axes)))

#Print the mean axis of rotation, in trend and plunge
Fols_meanAxisDeg

#==========================FOLIATIONS FROM FOLIATION-LINEATION PAIRS STATISTICAL TREATMENT========================================

#### FOLIATIONS FROM FOLIATION-LINEATION PAIRS STATISTICAL TREATMENT

# 1) Some parametric two-sample tests.

# Three Wellner tests (Wellner, 1979), one for each pair of domains. Each test is based on 10,000 permutations.
lineWellnerInference(Follins$pole[Follins_northCrit],  Follins$pole[Follins_southCrit], 10000)
lineWellnerInference(Follins$pole[Follins_northCrit],  Follins$pole[Follins_centerCrit], 10000)
lineWellnerInference(Follins$pole[Follins_centerCrit], Follins$pole[Follins_southCrit], 10000)

# Three Watson tests that assume large sample size (Mardia and Jupp, 2000), one for each pair of domains.
lineLargeMultiSampleWatsonInference(list(Follins$pole[Follins_northCrit], Follins$pole[Follins_southCrit]))
lineLargeMultiSampleWatsonInference(list(Follins$pole[Follins_northCrit], Follins$pole[Follins_centerCrit]))
lineLargeMultiSampleWatsonInference(list(Follins$pole[Follins_centerCrit], Follins$pole[Follins_southCrit]))

# Three Watson tests that assume tightly concentrated datasets (Mardia and Jupp, 2000), one for each pair of domains.
lineConcentratedMultiSampleWatsonInference(list(Follins$pole[Follins_northCrit], Follins$pole[Follins_southCrit]))
lineConcentratedMultiSampleWatsonInference(list(Follins$pole[Follins_northCrit], Follins$pole[Follins_centerCrit]))
lineConcentratedMultiSampleWatsonInference(list(Follins$pole[Follins_centerCrit], Follins$pole[Follins_southCrit]))

#————————————————————————————————————————————————————————————

# 2) Non-parametric bootstrapping

# Perform the bootstrapping routine for each domain. Each bootstrapped dataset is based on 10,000 iterations.
FollinFols_northBoots <-  lineBootstrapInference(Follins$pole[Follins_northCrit], 10000, numPoints = 50)
FollinFols_centerBoots <- lineBootstrapInference(Follins$pole[Follins_centerCrit], 10000, numPoints = 50)
FollinFols_southBoots <-  lineBootstrapInference(Follins$pole[Follins_southCrit], 10000, numPoints = 50)


# Plot data for each domain. Northern (red), central (green), southern (blue)
lineEqualAreaPlotThree(
                       Follins$pole[Follins_northCrit], 
                       Follins$pole[Follins_centerCrit], 
                       Follins$pole[Follins_southCrit]
)

# Plot the bootstrapped mean clouds for each domain. Northern (red), central (green), southern (blue)
lineEqualAreaPlotThree(
                       FollinFols_northBoots$bootstraps,
                       FollinFols_centerBoots$bootstraps,
                       FollinFols_southBoots$bootstraps
)

# Plot the bootstrapped mean clouds with 95% confidence ellipses superimposed. Northern (red), central (green), southern (blue)
lineEqualAreaPlotThree(
                       list(FollinFols_northBoots$center),
                       list(FollinFols_centerBoots$center),
                       list(FollinFols_southBoots$center),
                       curves = list(
                                     FollinFols_northBoots$points,
                                     FollinFols_centerBoots$points,
                                     FollinFols_southBoots$points
                                     )
)

#————————————————————————————————————————————————————————————

# 3) Compute the rotation between the north and south domains.

# Create empty lists to store data
Follins_northSouthDif <- list()

# Assign variables for the orientations of all bootstrapped means.
Follins_northB <- FollinFols_northBoots$bootstraps
Follins_southB <- FollinFols_southBoots$bootstraps

# Calculate the smallest possible rotations that bring means from the northern bootstrap cloud to the southern bootstrap clouds.
count = 1
for (i in 1:100) {
    for (n in i:100) {
        Follins_northSouthDif$rotation[[count]] = rotSmallestRotationFromTwoLines(Follins_northB[[i]], Follins_southB[[n]])
        if (n < 1000) {
            count = count + 1
        }
    }
}

# Display the raw results of the preceting for loop
Follins_northSouthDif

# Plot these rotations in an equal volume plot
oriEqualVolumePlot(
                   Follins_northSouthDif$rotation,
                   group = oriLineInPlaneGroup,
                   simplePoints = TRUE
)

# Compute the axis and rotation amount from all the rotations
Follins_angleAxis <- lapply(Follins_northSouthDif$rotation, function(s) rotAxisAngleFromMatrix(s))

# Plot the rotation axes in an equal area plot
lineEqualAreaPlot(lapply(Follins_angleAxis, function(s) c(s[1], s[2], s[3])), shapes = '.')

# Plot a histogram of rotation amount between the northern domain and the southern domain, in degrees.
hist(
     as.numeric(lapply(Follins_angleAxis, function(s)
        s[4] * 180 / pi)),
     30,
     xlab = "Angular Distance (degrees)",
     ylab = "Frequency",
     main = "Angular distances, Northern to Southern"
)

# Compute the mean and 1-sigma standard deviation of the rotation amount, in Follins_angleAxisrees
mean(sapply(Follins_angleAxis, function(s) s[4] * 180 / pi))
sd(sapply(Follins_angleAxis, function(s) s[4] * 180 / pi))

# Compute the mean trend and plunge in Follins_angleAxisrees
Follins_axes <- lapply(Follins_angleAxis, function(s) c(s[1], s[2], s[3]))
Follins_meanAxisDeg <- geoTrendPlungeDegFromCartesian(lower(lineProjectedMean(Follins_axes)))

#Print the mean axis of rotation, in trend and plunge
Follins_meanAxisDeg

#==========================COMPARE FOLIATIONS FROM THE TWO DATASETS========================================

# (Foliation-only vs. Foliation-lineation datasets)

# 1) Three plots to compare the foliations of the Fols (Foliation only data set) and Follins (Foliation-lineation dataset) in the northern domain

# Equal area plot of the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(Follins$pole[Follins_northCrit], Fols$pole[Fols_northCrit])

# Equal area plot of the bootstrapped means for the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(
                     FollinFols_northBoots$bootstraps,
                     Fols_northBoots$bootstraps
)

# Equal area plot of 95% confidence ellipses from bootstrapping for the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(
                     list(FollinFols_northBoots$center),
                     list(Fols_northBoots$center),
                     curves = list(FollinFols_northBoots$points,
                                   Fols_northBoots$points
                     )
)

#————————————————————————————————————————————————————————————

# 2) Three plots to compare the foliations of the Fols (Foliation only data set) and Follins (Foliation-lineation dataset) in the central domain

# Equal area plot of the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(
                     Follins$pole[Follins_centerCrit], 
                     Fols$pole[Fols_centerCrit]
)

# Equal area plot of the bootstrapped means for the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(
                     FollinFols_centerBoots$bootstraps,
                     Fols_centerBoots$bootstraps
)

# Equal area plot of 95% confidence ellipses from bootstrapping for the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(
                     list(FollinFols_centerBoots$center),
                     list(Fols_centerBoots$center),
                     curves = list(FollinFols_centerBoots$points,
                                   Fols_centerBoots$points
                     )
)

#————————————————————————————————————————————————————————————

# 3) Three plots to compare the foliations of the Fols (Foliation only data set) and Follins (Foliation-lineation dataset) in the southern domain

# Equal area plot of the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(
                     Follins$pole[Follins_southCrit], 
                     Fols$pole[Fols_southCrit]
)

# Equal area plot of the bootstrapped means for the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(
                     FollinFols_southBoots$bootstraps,
                     Fols_southBoots$bootstraps
)

# Equal area plot of 95% confidence ellipses from bootstrapping for the Fols (cyan) and Follins (red)
lineEqualAreaPlotTwo(
                     list(FollinFols_southBoots$center),
                     list(Fols_southBoots$center),
                     curves = list(FollinFols_southBoots$points,
                                   Fols_southBoots$points
                     )
)


#==================================================================
#==================================================================



#==========================PART II: ORIENTATION STATISTICS ON FOLIATION-LINEATION PAIRS========================================

#These are all foliation-lineation pairs, so they are orientation data. We will proced with methods outlined in Davis and Titus, 2017.

#Since lines on foliation are bidirectional, these are line-in-plane data, and have a four fold symmetry (see Davis and Titus, 2017).

#==========================PLOT THE DATA========================================

# Plot the data in Equal area plot. Lineation (red), Pole to foliation (cyan)
lineEqualAreaPlotTwo(
                     Follins$direction, 
                     Follins$pole
)

# Plot the data in an Equal Volume Plot (after Davis and Titus, 2017). Each point represents a foliation-lineation pair. (There are four copies of the data due to mathematical symmetry)
oriEqualVolumePlot(
                   Follins$rotation, 
                   oriLineInPlaneGroup
)

# We can also do some basic directional kamb contouring for poles and lines. The numbers are the kamb intervals
lineKambPlot(
             Follins$pole, 
             c(2, 4, 6, 8, 10, 12, 14, 18)
)
lineKambPlot(
             Follins$direction, 
             c(2, 4, 6, 8, 10, 12, 14, 18)
)

# [!!!]NOTE: This line requires c library to run--skip if you have not compiled c.
# [!!!]Uncomment lines below to Plot 6-sigma kamb contour in an Equal volume plot of Follins.

# oricKambPlot(Follins$rotation,
#              group=oriLineInPlaneGroup,
#              multiple = 6,
#              simplePoints = TRUE,
#              backgroundColor="white", curveColor="black",
#              boundaryAlpha=0,
#              colors="black", axesColors=c("black", "black", "black"),
#              fogStyle="exp")

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.

# afterMaximizingWindow("WestMt_kamb_EqualVol1.png", "WestMt_kamb_EqualVol2.png")


#==========================IDENTIFY GEOGRAPHIC TRENDS========================================

# 1) Visually inspect the possibility of geographic trends in the data

# Plot an Equal Volume plot (after Davis and Titus, 2017), colored by the three domains. Northern (red); Central (green); Southern (blue).
oriEqualVolumePlot(
                   Follins$rotation,
                   group = oriLineInPlaneGroup,
                   backgroundColor = "white",
                   curveColor = "black",
                   boundaryAlpha = 0.2,
                   colors = hues(Follins$domain),
                   axesColors = c("black", "black", "black")
)

# Plot an Equal volume plot of the data, colored by northing
oriEqualVolumePlot(
                   Follins$rotation,
                   group = oriLineInPlaneGroup,
                   backgroundColor = "white",
                   curveColor = "black",
                   boundaryAlpha = 0.2,
                   colors = hues(Follins$northing),
                   axesColors = c("black", "black", "black")
)

# Plot an Equal volume plot of the data, colored by easting
oriEqualVolumePlot(
                   Follins$rotation,
                   group = oriLineInPlaneGroup,
                   backgroundColor = "white",
                   curveColor = "black",
                   boundaryAlpha = 0.2,
                   colors = hues(Follins$easting),
                   axesColors = c("black", "black", "black")
)

# Save EqualVolumePlot figure. If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow(
                      "WestMt_domains_EqualVol.png",
                      "WestMt_domains_EqualVol.png"
)

# Plot Follins in an Equal area plot, colored by domain. Lineation (squares); Foliation (circles)
lineEqualAreaPlot(
                  c(Follins$pole, Follins$direction),
                  col = hues(Follins$domain),
                  shapes = c(replicate(length(Follins$pole), "c"),
                             replicate(length(Follins$direction), "s"))
)

#————————————————————————————————————————————————————————————

# 2) Run regressions to quantify any extant trends.

# 2A) Geodesic regression on all data. Perform a geodesic regression with respect to azimuths every 10 degrees. Each regression asks the question: "does orientation change linearly with respect to an azimuth towards x degrees"
regressionSWestAll <- regressionSweep(Follins, 10, Follins_northCrit | Follins_centerCrit | Follins_southCrit , 0)
westAllRegSum <- data.frame(regressionSWestAll[[19]])
names(westAllRegSum) <- c("Azimuth", "Error", "Min_Eigen", "R_squared", "Pvalue")


# Plot the R^2 value (y-axis) as a function of the azimuth in degrees (x-axis)
plot(
     as.vector(westAllRegSum[, 1]),
     as.vector(westAllRegSum[, 4]),
     main = "Geodesic regressions (all Domains)",
     xlab = "Azimuth in degrees",
     ylab = "R-squared"
)   

#————————————————————————————————————————————————————————————

#2B) Perform geodesic regressions on each domain

# Perform a geodesic regression for the northern domain with respect to azimuths every 10 degrees.
regressionSNorth <- regressionSweep(Follins, 10, Follins_northCrit, 0)
NorthRegSum <- data.frame(regressionSNorth[[19]])
names(NorthRegSum) = c("Azimuth", "Error", "Min_Eigen", "R_squared", "Pvalue")

# Check that "Error" = 0, "Min_Eigen" always > 0.
NorthRegSum

# Plot the R^2 value (y-axis) as a function of the azimuth in degrees (x-axis)
plot(
     as.vector(NorthRegSum$Azimuth),
     as.vector(NorthRegSum$R_squared),
     main = "Geodesic regressions (Northern Domain)",
     xlab = "Azimuth in degrees",
     ylab = "R-squared"
)

#————————————————————————————————————————————————————————————

# Perform a geodesic regression for the southern domain with respect to azimuths every 10 degrees.
regressionSCenter <- regressionSweep(Follins, 10, Follins_centerCrit, 0)
CenterRegSum <- data.frame(regressionSCenter[[19]])
names(CenterRegSum) = c("Azimuth", "Error", "Min_Eigen", "R_squared", "Pvalue")

# Check that "Error" = 0, "Min_Eigen" always > 0.
CenterRegSum

# Plot the R^2 value (y-axis) as a function of the azimuth in degrees (x-axis)
plot(
     as.vector(CenterRegSum$Azimuth),
     as.vector(CenterRegSum$R_squared),
     main = "Geodesic regressions (Central Domain)",
     xlab = "Azimuth in degrees",
     ylab = "R-squared"
) 

#————————————————————————————————————————————————————————————

# Perform a geodesic regression for the southern domain with respect to azimuths every 10 degrees.
regressionSSouth <- regressionSweep(Follins, 10, Follins_southCrit, 0)
SouthRegSum <- data.frame(regressionSSouth[[19]])
names(SouthRegSum) = c("Azimuth", "Error", "Min_Eigen", "R_squared", "Pvalue")

# Check that "Error" = 0, "Min_Eigen" always > 0.
SouthRegSum

# Plot the R^2 value (y-axis) as a function of the azimuth in degrees (x-axis)
plot(
     as.vector(SouthRegSum$Azimuth),
     as.vector(SouthRegSum$R_squared),
     main = "Geodesic regressions (Southern Domain)",
     xlab = "Azimuth in degrees",
     ylab = "R-squared"
)

#————————————————————————————————————————————————————————————

# 3) Perform a Kernal regression for all Follins.

# Recategorize the symmetric copies of the data.
mu <- oriFrechetMean(Follins$rotation, group=oriLineInPlaneGroup)
Follins$rotation <- oriNearestRepresentatives(Follins$rotation, mu, group=oriLineInPlaneGroup)

# A precomputed value for bandwidth--use this value to save time.

bandwidth <- 0.5879134

# Uncomment line below to compute the appropriate bandwidth for Kernal regression with respect to easting.
bandwidth <- rotBandwidthForKernelRegression(Follins$easting,Follins$rotation, dnorm)

#Run the kernel regression to generate a bunch of points.There were no errors and all minEigenvalues were >0
kernelReg <- lapply(
                    seq(from = min(Follins$easting), to = max(Follins$easting), by = 100),
                    rotKernelRegression,
                    Follins$easting,
                    Follins$rotation,
                    bandwidth,
                    numSteps = 1000
)
sapply(kernelReg, function(regr) regr$error)
sapply(kernelReg, function(regr) {regr$minEigenvalue > 0})

# Plot the regression curve.
kernelRegCurve <- lapply(kernelReg, function(regr) regr$r)
oriEqualVolumePlot(
                   kernelRegCurve, 
                   simplePoints = TRUE,
                   oriLineInPlaneGroup
)

# Compute the R^2.
KernRsquared <- rotRsquaredForKernelRegression(Follins$easting, Follins$rotation, bandwidth, numSteps =
                                       1000)
KernRsquared

# Do a permutation test for significance. For the sake of time, we do only 10 permutations, although that is far too few to tell us anything.
rSquareds <- rotKernelRegressionPermutations(Follins$easting, Follins$rotation, bandwidth, numPerms =
                                        1000)
length(rSquareds)
p <- sum(rSquareds > KernRsquared$rSquared)
p


#==========================STATISTICAL DESCRIPTORS========================================

# 1) Northern Domain, n = 16

# Frechet mean (minimizes the Frechet variance).
FrechetMeanNorth <- oriFrechetMean(Follins$rotation[Follins_northCrit], oriLineInPlaneGroup)
FrechetVarNorth <-  oriVariance(Follins$rotation[Follins_northCrit], FrechetMeanNorth, oriLineInPlaneGroup)

# Plot the FrechetMean in an Equal Volume plot
FrechetCurvesNorth <- lapply(Follins$rotation[Follins_northCrit], function(r) rotGeodesicPoints(FrechetMeanNorth, r, 10))
rotEqualAnglePlot(points = Follins$rotation[Follins_northCrit], curves = FrechetCurvesNorth)

# Print the Strike, Dip, Rake of the Frechet mean.
geoStrikeDipRakeDegFromRotation(FrechetMeanNorth)

# Print the Frechet variance.
FrechetVarNorth

#————————————————————————————————————————————————————————————

# 2) Central Domain, n = 34
#Frechet mean minimizes the Frechet variance.
FrechetMeanCenter <- oriFrechetMean(Follins$rotation[Follins_centerCrit], oriLineInPlaneGroup)
FrechetVarCenter <- oriVariance(Follins$rotation[Follins_centerCrit], FrechetMeanCenter, oriLineInPlaneGroup)

#plot the FrechetMean
FrechetCurvesCenter <- lapply(Follins$rotation[Follins_centerCrit], function(r) rotGeodesicPoints(FrechetMeanCenter, r, 10))
rotEqualAnglePlot(
                  points = Follins$rotation[Follins_centerCrit], 
                  curves = FrechetCurvesCenter
)

# Print the Strike, Dip, Rake of the Frechet mean.
geoStrikeDipRakeDegFromRotation(FrechetMeanCenter)

# Print the Frechet variance.
FrechetVarCenter

#————————————————————————————————————————————————————————————

# 3) Southern Domain, n = 79
#Frechet mean minimizes the Frechet variance.
FrechetMeanSouth <- oriFrechetMean(Follins$rotation[Follins_southCrit], oriLineInPlaneGroup)
FrechetVarSouth <-  oriVariance(Follins$rotation[Follins_southCrit], FrechetMeanSouth, oriLineInPlaneGroup)

#plot the FrechetMean
FrechetCurvesSouth <- lapply(Follins$rotation[Follins_southCrit], function(r) rotGeodesicPoints(FrechetMeanSouth, r, 10))
rotEqualAnglePlot(points = Follins$rotation[Follins_southCrit], curves = FrechetCurvesSouth)

# Print the Strike, Dip, Rake of the Frechet mean.
geoStrikeDipRakeDegFromRotation(FrechetMeanSouth)

# Print the Frechet variance.
FrechetVarSouth

#————————————————————————————————————————————————————————————

# 4) Dispersion of the data

#The standard way to get at dispersion is to compute the maximum matrix fisher likelihood. The matrix fisher distribution comprises a kind of "mean" and a positive definite symmetric matrix, which characterises the anisotropic dispersion.

# Northern Domain, n = 16. Fisher maximum likelihood
mleNorth <- rotFisherMLE(Follins$rotation[Follins_northCrit])
mleNorth
eigen(mleNorth$kHat, symmetric = TRUE, only.value = TRUE)$values

#Central Domain, n = 34. Fisher maximum likelihood
mleCenter <- rotFisherMLE(Follins$rotation[Follins_centerCrit])
mleCenter
eigen(mleCenter$kHat, symmetric = TRUE, only.value = TRUE)$values

#Southern Domain, n = 79. Fisher maximum likelihood
mleSouth <- rotFisherMLE(Follins$rotation[Follins_southCrit])
mleSouth
eigen(mleSouth$kHat, symmetric = TRUE, only.value = TRUE)$values



#==========================INFERENCE + HYPOTHESIS TESTING========================================

# Here, we can perform some statistical techniques to use information about the samples (each domain) to compute a probability cloud of the mean of the population(s) from which those samples were taken.
# From numerical work in Davis and Titus (2017), MCMC will work the best for small sample sizes that have the matrix fisher anisotropy (eigenvalues) of (large, large, small). All four domains have this anisotropy, and have too few data points to use the (computationally quicker) bootstrapping method.

#We'll do both and compare.

# Compute the bootstrap cloud of means
Follins_northBoots <-  oriBootstrapInference(Follins$rotation[Follins_northCrit], 10000, oriLineInPlaneGroup)
Follins_centerBoots <- oriBootstrapInference(Follins$rotation[Follins_centerCrit], 10000, oriLineInPlaneGroup)
Follins_southBoots <-  oriBootstrapInference(Follins$rotation[Follins_southCrit], 10000, oriLineInPlaneGroup)

# [!!!] Requires C libraries to run. 
# [!!!] Compute the Markov chain Monte Carlo cloud of means
# [!!!]Uncomment lines if you have compiled C on your computer and have loaded the C libraries at the top of this script.
# 
# northMCMC <-  oricWrappedTrivariateNormalMCMCInference(Follins$rotation[Follins_northCrit],
#                                                       group = oriLineInPlaneGroup,
#                                                       numCollection = 100
#               )
# centerMCMC <- oricWrappedTrivariateNormalMCMCInference(Follins$rotation[Follins_centerCrit],
#                                                        group = oriLineInPlaneGroup,
#                                                        numCollection = 100
#               )
# southMCMC <-  oricWrappedTrivariateNormalMCMCInference(Follins$rotation[Follins_southCrit],
#                                                        group = oriLineInPlaneGroup,
#                                                        numCollection = 100
#               )

#————————————————————————————————————————————————————————————

# 1) Northern v. Central domains

# A) Using MCMC

# [!!!] Uncomment if you ran the MCMC lines above, which required C to be compiled on your system.
# Construct the 95% confidence ellipsoids from small triangles
# trisNorthMCMC <-  rotEllipsoidTriangles(northMCMC$mBar,
#                                         northMCMC$leftCovarInv,
#                                         northMCMC$q095,
#                                         numNonAdapt = 4)
# trisCenterMCMC <- rotEllipsoidTriangles(centerMCMC$mBar,
#                                         centerMCMC$leftCovarInv,
#                                         centerMCMC$q095,
#                                         numNonAdapt = 4)
# Plot the MCMC comparison in an equal Volume plot, with 95% confidence ellipsoids
# oriEqualAnglePlot(
#                   points = c(northMCMC$ms, centerMCMC$ms),
#                   boundaryAlpha = .1,
#                   axesColors = c("black", "black", "black"),
#                   fogStyle = "none",
#                   background = "white",
#                   triangles = c(trisNorthMCMC, trisCenterMCMC),
#                   simplePoints = TRUE,
#                   colors = c(replicate(length(northMCMC$ms), "black"), replicate(length(centerMCMC$ms), "orange")),
#                   group = oriTrivialGroup
# )
# # If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
# afterMaximizingWindow("MCMC_WestMt_NC_1.png", "MCMC_WestMt_NC_2.png")
# 
# # Plot the MCMC comparison in an Equal Area plot
# lineEqualAreaPlotTwo(c(
#                        lapply(northMCMC$ms, function(s) s[1,]),
#                        lapply(northMCMC$ms, function(s) s[2,])
#                        ),
#                      c(
#                        lapply(centerMCMC$ms, function(s) s[1,]),
#                        lapply(centerMCMC$ms, function(s) s[2,])
# ))

#....................................................................

# B) Using bootstrapping

# Construct the 95% confidence ellipsoids from small triangles
trisNorthBoot <-
    rotEllipsoidTriangles(
        Follins_northBoots$center,
        Follins_northBoots$leftCovarInv,
        Follins_northBoots$q095,
        numNonAdapt = 4
    )
trisCenterBoot <-
    rotEllipsoidTriangles(
        Follins_centerBoots$center,
        Follins_centerBoots$leftCovarInv,
        Follins_centerBoots$q095,
        numNonAdapt = 4
    )

# Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids
rotEqualAnglePlot(
    points = c(
        Follins_northBoots$bootstraps,
        Follins_centerBoots$bootstraps
    ),
    triangles = c(trisNorthBoot, trisCenterBoot),
    boundaryAlpha = .1,
    axesColors = c("black", "black", "black"),
    fogStyle = "none",
    background = "white",
    simplePoints = TRUE,
    colors = c(replicate(
        length(Follins_northBoots$bootstraps), "black"
    ), replicate(
        length(Follins_centerBoots$bootstraps), "orange"
    ))
)

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow("Boots_WestMt_NC_1.png", "Boots_WestMt_NC_2.png")

# Plot the bootstrap comparison in an Equal Area plot
lineEqualAreaPlotTwo(c(
                       lapply(Follins_northBoots$bootstraps, function(s) s[1,]),
                       lapply(Follins_northBoots$bootstraps, function(s) s[2,])
                       ),
                     c(
                       lapply(Follins_centerBoots$bootstraps, function(s) s[1,]),
                       lapply(Follins_centerBoots$bootstraps, function(s) s[2,])
))

#————————————————————————————————————————————————————————————

# 2) Northern vs. Southern domains

# A) Using MCMC

# [!!!] Uncomment if you ran the MCMC lines above, which required C to be compiled on your system.
# Construct the 95% confidence ellipsoids from small triangles
# trisNorthMCMC <-
#     rotEllipsoidTriangles(northMCMC$mBar,
#                           northMCMC$leftCovarInv,
#                           northMCMC$q095,
#                           numNonAdapt = 5)
# trisSouthMCMC <-
#     rotEllipsoidTriangles(southMCMC$mBar,
#                           southMCMC$leftCovarInv,
#                           southMCMC$q095,
#                           numNonAdapt = 5)
# 
# # Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids
# oriEqualAnglePlot(
#     points = c(northMCMC$ms, southMCMC$ms),
#     triangles = c(trisNorthMCMC, trisSouthMCMC),
#     boundaryAlpha = 0.1,
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none",
#     background = "white",
#     simplePoints = TRUE,
#     colors = c(replicate(length(northMCMC$ms), "black"), replicate(length(southMCMC$ms), "blue")),
#     group = oriTrivialGroup
# )
# 
# # If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
# afterMaximizingWindow("MCMC_WestMt_NS_1.png", "MCMC_WestMt_NS_2.png")
# 
# # Plot the MCMC comparison in an Equal Area plot
# lineEqualAreaPlotTwo(c(
#     lapply(northMCMC$ms, function(s)
#         s[1,]),
#     lapply(northMCMC$ms, function(s)
#         s[2,])
# ),
# c(
#     lapply(southMCMC$ms, function(s)
#         s[1,]),
#     lapply(southMCMC$ms, function(s)
#         s[2,])
# ))

#....................................................................

# B) Using bootstrapping
# Construct the 95% confidence ellipsoids from small triangles
trisNorthBoot <-
    rotEllipsoidTriangles(
        Follins_northBoots$center,
        Follins_northBoots$leftCovarInv,
        Follins_northBoots$q095,
        numNonAdapt = 4
    )
trisSouthBoot <-
    rotEllipsoidTriangles(
        Follins_southBoots$center,
        Follins_southBoots$leftCovarInv,
        Follins_southBoots$q095,
        numNonAdapt = 4
    )

# Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids.
rotEqualAnglePlot(
    points = c(
        Follins_northBoots$bootstraps,
        Follins_southBoots$bootstraps
    ),
    triangles = c(trisNorthBoot, trisSouthBoot),
    boundaryAlpha = .1,
    axesColors = c("black", "black", "black"),
    fogStyle = "none",
    background = "white",
    simplePoints = TRUE,
    colors = c(replicate(
        length(Follins_northBoots$bootstraps), "black"
    ),
    replicate(
        length(Follins_southBoots$bootstraps), "blue"
    ))
)

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow("Boots_WestMt_NS_1.png", "Boots_WestMt_NS_2.png")

# Plot the bootstrap comparison in an Equal Area plot
lineEqualAreaPlotTwo(c(
                       lapply(Follins_northBoots$bootstraps, function(s) s[1,]),
                       lapply(Follins_northBoots$bootstraps, function(s) s[2,])
                       ),
                     c(
                       lapply(Follins_southBoots$bootstraps, function(s) s[1,]),
                       lapply(Follins_southBoots$bootstraps, function(s) s[2,])
))



#————————————————————————————————————————————————————————————

# 3) Central vs. Southern domains

# A) Using MCMC
# [!!!] Uncomment if you ran the MCMC lines above, which required C to be compiled on your system.
# # Construct the 95% confidence ellipsoids from small triangles
# trisCenterMCMC <-
#     rotEllipsoidTriangles(centerMCMC$mBar,
#                           centerMCMC$leftCovarInv,
#                           centerMCMC$q095,
#                           numNonAdapt = 4)
# trisSouthMCMC <-
#     rotEllipsoidTriangles(southMCMC$mBar,
#                           southMCMC$leftCovarInv,
#                           southMCMC$q095,
#                           numNonAdapt = 4)
# 
# # Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids
# oriEqualAnglePlot(
#     points = c(centerMCMC$ms, southMCMC$ms),
#     triangles = c(trisCenterMCMC, trisSouthMCMC),
#     boundaryAlpha = 0.1,
#     axesColors = c("black", "black", "black"),
#     fogStyle = "none",
#     background = "white",
#     simplePoints = TRUE,
#     colors = c(replicate(length(centerMCMC$ms), "orange"), replicate(length(southMCMC$ms), "blue")),
#     group = oriTrivialGroup
# )
# 
# # If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
# afterMaximizingWindow("MCMC_WestMt_CS_1.png", "MCMC_WestMt_CS_2.png")
# 
# # Plot the MCMC comparison in an Equal Area plot
# lineEqualAreaPlotTwo(c(
#     lapply(centerMCMC$ms, function(s)
#         s[1,]),
#     lapply(centerMCMC$ms, function(s)
#         s[2,])
# ),
# c(
#     lapply(southMCMC$ms, function(s)
#         s[1,]),
#     lapply(southMCMC$ms, function(s)
#         s[2,])
# ))

#....................................................................


# B) Using bootstrapping

# Construct the 95% confidence ellipsoids from small triangles
trisCenterBoot <-
    rotEllipsoidTriangles(
        Follins_centerBoots$center,
        Follins_centerBoots$leftCovarInv,
        Follins_centerBoots$q095,
        numNonAdapt = 4
    )
trisSouthBoot <-
    rotEllipsoidTriangles(
        Follins_southBoots$center,
        Follins_southBoots$leftCovarInv,
        Follins_southBoots$q095,
        numNonAdapt = 4
    )

# Plot the bootstrap comparison in an equal Volume plot, with 95% confidence ellipsoids
oriEqualAnglePlot(
    points = c(
        Follins_centerBoots$bootstraps,
        Follins_southBoots$bootstraps
    ),
    triangles = c(trisCenterBoot, trisSouthBoot),
    boundaryAlpha = 0.1,
    axesColors = c("black", "black", "black"),
    fogStyle = "none",
    background = "white",
    simplePoints = TRUE,
    colors = c(replicate(
        length(Follins_northBoots$bootstraps), "orange"
    ), replicate(
        length(Follins_southBoots$bootstraps), "blue"
    )),
    group = oriTrivialGroup
)

# If you wish to save this figure, maximize the plot window on your screen before running this line. It will save to your working directory folder.
afterMaximizingWindow("Boots_WestMt_CS_1.png", "Boots_WestMt_CS_2.png")

# Plot the bootstrap comparison in an Equal Area plot
lineEqualAreaPlotTwo(c(
                       lapply(Follins_centerBoots$bootstraps, function(s) s[1,]),
                       lapply(Follins_centerBoots$bootstraps, function(s) s[2,])
                       ),
                     c(
                       lapply(Follins_southBoots$bootstraps, function(s) s[1,]),
                       lapply(Follins_southBoots$bootstraps, function(s) s[2,])
))

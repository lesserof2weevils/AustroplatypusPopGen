# Title: Population genetics and phylogeny of Austroplatypus incompertus
# Script purpose: Descriptive statistics
# Authors: James Bickerstaff & Markus Riegler
# Date: 22/2/2021

# Packages to call ####
libs <- c("dartR", "adegenet", "LEA", "conStruct", "parallel", "fields",
          "maps", "ecodist", "hierfstat", "StAMPP", "ade4", "poppr",
          "netview", "visreg", "devtools", "ggplot2", "tidyverse",
          "doBy", "plotly", "reshape", "ggspatial", "ozmaps", "sf",
          "ggrepel", "patchwork", 'viridis', 'raster', 'ggnewscale',
          'patchwork', 'ggrepel', 'scatterpie', 'PopGenome')  # List of required packages

# Install any packages that are not already installed
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == FALSE)) {
  install.packages(libs[!installed_libs])
}

# Load all required libraries
invisible(lapply(libs, library, character.only = TRUE))

# Set the working directory to the project folder
setwd("~/Documents/GitHub/AustroplatypusPopGen/")

# 1 - Split the dataset according to conStruct and filter ####

## 1.1 Splitting the dataset into Northern (Nth) vs Southern (Sth) populations #####
# Read in DArT SNP data and individual metadata
gl <- gl.read.dart(
  filename = "data/Report_DAmb20-5632_1_moreOrders_SNP_2_rename_ordered.csv",
  ind.metafile = "data/Ai_ind_metrics.csv"
)

# Keep only the northern populations: Dorrigo, Kerewong, NthBarrington
glNth <- gl.keep.pop(gl, pop.list = c('Dorrigo', 'Kerewong', 'NthBarrington'), recalc = TRUE)
nInd(glNth)               # Number of individuals in the northern dataset
nLoc(glNth)               # Number of loci in the northern dataset
levels(pop(glNth))        # Population levels in the northern dataset
glNth$ind.names           # Names of individuals in the northern dataset

# Keep only the southern populations by dropping northern ones and WilliWilli
glSth <- gl.drop.pop(gl, pop.list = c('Dorrigo', 'Kerewong', 'NthBarrington', 'WilliWilli'), recalc = TRUE)
nInd(glSth)               # Number of individuals in the southern dataset
nLoc(glSth)               # Number of loci in the southern dataset
levels(pop(glSth))        # Population levels in the southern dataset
glSth$ind.names           # Names of individuals in the southern dataset

## 1.2 Filtering #####
# Apply filters to the full dataset
glFilt <- gl %>%
  gl.filter.reproducibility(threshold = 0.99) %>%  # Remove loci with reproducibility below 99%
  gl.filter.monomorphs() %>%                       # Remove monomorphic loci
  gl.filter.callrate(threshold = 0.95, method = 'loc') %>%  # Remove loci with call rate below 95%
  gl.filter.callrate(method = 'ind', threshold = 0.9, recalc = TRUE)  # Remove individuals with call rate below 90%

nInd(glFilt)             # Number of individuals after filtering
nLoc(glFilt)             # Number of loci after filtering
levels(pop(glFilt))      # Population levels after filtering
levels(gallery(glFilt))

glFilt@other$ind.metrics$gallery

# Apply the same filters to the northern dataset
glNthFilt <- glNth %>%
  gl.filter.reproducibility(threshold = 0.99) %>%
  gl.filter.monomorphs() %>%
  gl.filter.callrate(threshold = 0.95, method = 'loc') %>%
  gl.filter.callrate(method = 'ind', threshold = 0.9, recalc = TRUE)

nInd(glNthFilt)          # Number of individuals in the filtered northern dataset
nLoc(glNthFilt)          # Number of loci in the filtered northern dataset
levels(pop(glNthFilt))   # Population levels in the filtered northern dataset
namesNth <- glNthFilt$ind.names  # Names of individuals
write.csv(namesNth, file = 'data/derived/namesNth.csv')  # Save the names to a CSV file
ind.callrate <- NA.posi(glNth)  # Identify positions with missing data

# Read in hierarchical information for the northern populations
Nthhierarchy <- read.csv('data/derived/heirarchy_Nth.csv')
glNth$strata <- Nthhierarchy  # Assign strata (e.g., population or gallery) to the genlight object
head(strata(glNth, ~ pop/gallery, combine = FALSE))  # Preview strata assignments

# Convert the genlight object to a genind object for analysis
glNthi <- gl2gi(glNthFilt)
glNthi$strata <- Nthhierarchy  # Assign strata to the genind object
names(other(glNthi))           # List additional data stored in the genind object
head(strata(glNthi, ~pop/gallery, combine = FALSE))  # Preview strata
glNthclone <- as.genclone(glNthi)  # Convert to a genclone object for clone-aware analyses

# Apply filters to the southern dataset
glSthFilt <- glSth %>%
  gl.filter.reproducibility(threshold = 0.99) %>%
  gl.filter.monomorphs() %>%
  gl.filter.callrate(threshold = 0.95, method = 'loc') %>%
  gl.filter.callrate(method = 'ind', threshold = 0.9, recalc = TRUE)

nInd(glSthFilt)          # Number of individuals in the filtered southern dataset
nLoc(glSthFilt)          # Number of loci in the filtered southern dataset
levels(pop(glSthFilt))   # Population levels in the filtered southern dataset
namesSth <- glSthFilt$ind.names  # Names of individuals
write.csv(namesSth, file = 'data/derived/namesSth.csv')  # Save the names to a CSV file
ind.callrate <- NA.posi(glSth)  # Identify positions with missing data

# Read in hierarchical information for the all populations
glFilt@strata <- data.frame(gallery = as.character(glFilt@other$ind.metrics$gallery)) # Convert the factor to a data frame
glFilt@strata$pop <- glFilt@pop # Add in population data into strata assignment
head(strata(glFilt, ~ pop/gallery, combine = FALSE))  # Preview strata assignments

# Read in hierarchical information for the southern populations
Sthhierarchy <- read.csv('data/derived/hierarchy_Sth.csv')
glSthFilt$strata <- Sthhierarchy  # Assign strata to the genlight object
head(strata(glSthFilt, ~ pop/gallery, combine = FALSE))  # Preview strata assignments

# Convert the genlight object to a genind object for analysis
glSthi <- gl2gi(glSthFilt)
glSthi$strata <- Sthhierarchy  # Assign strata to the genind object
names(other(glSthi))           # List additional data stored in the genind object
head(strata(glSthi, ~pop/gallery, combine = FALSE))  # Preview strata
glSthclone <- as.genclone(glSthi)  # Convert to a genclone object

# 2 - Generate report of expected and observed heterozygosity ####
# Calculate heterozygosity for all datasets at the population level
HeHoAll <- gl.report.heterozygosity(glFilt, method = 'pop')
HeHoNth <- gl.report.heterozygosity(glNthFilt, method = 'pop', verbose = 2)
HeHoSth <- gl.report.heterozygosity(glSthFilt, method = 'pop', verbose = 2)
HeHo <- gl.report.heterozygosity(gl, method = 'pop', verbose = 2)
par(mfrow = c(1,1))  # Reset plotting parameters

# Calculate Shannon's diversity index for all datasets
AllDiv <- gl.report.diversity(glFilt)
NthDiv <- gl.report.diversity(glNthFilt)
SthDiv <- gl.report.diversity(glSthFilt)

# Calculate Tajima's D for all datasets
AllTajima <- get_tajima_D(glFilt)
NthTajima <- get_tajima_D(glNthFilt)
SthTajima <- get_tajima_D(glSthFilt)

## 2.1 Calculation of Heterozygosity measures and Fis #####
# Compute basic statistics for the northern dataset (e.g., Ho, He, Fis)
DSNth <- glNthi %>%
  basic.stats(diploid = TRUE)

# Compute basic statistics for the southern dataset
DSSth <- glSthi %>%
  basic.stats(diploid = TRUE)

# Custom function to calculate the standard error of the mean for Ho, He, and Fis
SEmean <- function(x){
  x <- x[!is.na(x)]               # Remove NA values
  se <- sd(x) / sqrt(length(x))   # Calculate standard error
  return(se)
}

## 2.2 Generate a summary statistics table for loci, heterozygosity, and inbreeding #####

# Prepare summary statistics for all populations
AllPopStats <- HeHoAll[, c("pop", "nLoc", "polyLoc", "monoLoc", "Ho", "HoSD", "He", "HeSD")]
AllPopStats$ShannonDiversity <- AllDiv$one_H_alpha         # Add Shannon's diversity index
AllPopStats$ShannonDiversitySD <- AllDiv$one_H_alpha_sd    # Add standard deviation
AllPopStats <- merge(AllPopStats, AllTajima, by.x = "pop", by.y = "population", all.x = TRUE)  # Merge with Tajima's D
desired_order <- c('Dorrigo', 'Kerewong', 'NthBarrington', 'SthBarrington', 
                   'Olney', 'Cumberland', 'MtWilson', 'MaresForest', 
                   'Termeil', 'Dampier', 'Broadwater')  # Desired population order
AllPopStats <- AllPopStats %>%
  mutate(pop = factor(pop, levels = desired_order)) %>%  # Set factor levels
  arrange(pop)  # Arrange by population

# Prepare summary statistics for northern populations
NthPopStats <- HeHoNth[, c("pop", "nLoc", "polyLoc", "monoLoc", "Ho", "HoSD", "He", "HeSD")]
NthPopStats$ShannonDiversity <- NthDiv$one_H_alpha
NthPopStats$ShannonDiversitySD <- NthDiv$one_H_alpha_sd
NthPopStats <- merge(NthPopStats, NthTajima, by.x = "pop", by.y = "population", all.x = TRUE)
desired_order <- c('Dorrigo', 'Kerewong', 'NthBarrington')
NthPopStats <- NthPopStats %>%
  mutate(pop = factor(pop, levels = desired_order)) %>%
  arrange(pop)

# Prepare summary statistics for southern populations
SthPopStats <- HeHoSth[, c("pop", "nLoc", "polyLoc", "monoLoc", "Ho", "HoSD", "He", "HeSD")]
SthPopStats$ShannonDiversity <- SthDiv$one_H_alpha
SthPopStats$ShannonDiversitySD <- SthDiv$one_H_alpha_sd
SthPopStats <- merge(SthPopStats, SthTajima, by.x = "pop", by.y = "population", all.x = TRUE)
desired_order <- c('SthBarrington', 'Olney', 'Cumberland', 'MtWilson', 'MaresForest', 
                   'Termeil', 'Dampier', 'Broadwater')
SthPopStats <- SthPopStats %>%
  mutate(pop = factor(pop, levels = desired_order)) %>%
  arrange(pop)

# Save the summary statistics tables to CSV files
write.csv(AllPopStats, file = 'outputs/AllPopStats.csv')
write.csv(NthPopStats, file = 'outputs/NthPopStats.csv')
write.csv(SthPopStats, file = 'outputs/SthPopStats.csv')

# 3 - AMOVA tests across populations ####
# Perform AMOVA for northern populations considering hierarchy of pop/gallery
NthAMOVA <- poppr.amova(glNthclone, ~ pop/gallery, cutoff = 0.6, nperm = 999)
NthAMOVA                    # Display AMOVA results
Nthamova.test <- randtest(NthAMOVA)  # Perform randomization test
Nthamova.test               # Display test results
plot(Nthamova.test)         # Plot the distribution of simulated statistics

# Perform AMOVA for southern populations
SthAMOVA <- poppr.amova(glSthclone, ~ pop/gallery, cutoff = 0.6, nperm = 999)
SthAMOVA                    # Display AMOVA results
Sthamova.test <- randtest(SthAMOVA)  # Perform randomization test
Sthamova.test               # Display test results
plot(Sthamova.test)         # Plot the distribution of simulated statistics

# Perform AMOVA for all populations
AMOVA <- poppr.amova(glFilt, ~ pop/gallery, cutoff = 0.6, nperm = 999)
AMOVA                    # Display AMOVA results
amova.test <- randtest(AMOVA)  # Perform randomization test
amova.test               # Display test results
plot(amova.test)         # Plot the distribution of simulated statistics



#Title: Population genetics and phylogeny of Austroplatypus incompertus
#Script purpose: Phylogenetic analyses
#Author: James Bickerstaff & Markus Riegler
#Date: 22/2/2021

#Packages to call####
libs <- c("dartR", "adegenet", "LEA", "conStruct", "parallel", "fields",
          "maps", "ecodist", "hierfstat", "StAMPP", "ade4", "poppr",
          "netview", "visreg", "devtools", "ggplot2", "tidyverse",
          "doBy", "plotly", "reshape", "ggspatial", "ozmaps", "sf",
          "ggrepel", "patchwork", 'viridis', 'raster', 'ggnewscale',
          'patchwork', 'ggrepel', 'scatterpie') # Define a list of required packages

# libraries we need
installed_libs <- libs %in% rownames(installed.packages()) # Check which packages are already installed
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs]) # Install packages that are not yet installed
}

# load libraries
invisible(lapply(libs, library, character.only = T)) # Load all required libraries silently

setwd("~/Documents/GitHub/AustroplatypusPopGen/") # Set the working directory

#1 - SNAPP ####
##1.1 Filtering####
gl <- gl.read.dart(filename = "data/Report_DAmb20-5632_1_moreOrders_SNP_2_rename_ordered.csv", ind.metafile = "data/Ai_ind_metrics.csv") # Read SNP data and individual metadata into a genlight object

gl.report.rdepth(gl) # Generate a report on read depth

indlist <- gl$ind.names  # Extract individual names from the genlight object
write.csv(indlist, file='glSNAPP.csv') # Save the list of individual names to a CSV file
drop.list <- as.matrix(read_csv("glSNAPP_keep.csv")) # Read a list of individuals to keep

glDrop <- gl.keep.ind(gl, drop.list, recalc = T, mono.rm = T) # Keep specified individuals and remove monomorphic loci
nLoc(glDrop) # Display the number of loci after filtering
nInd(glDrop) # Display the number of individuals after filtering
levels(pop(glDrop)) # Display the population retained after filtering

glSNAPP <- glDrop %>% 
  gl.filter.reproducibility(., threshold = .99) %>%  # Filter loci with reproducibility below 99%
  gl.filter.monomorphs(.) %>%  # Remove monomorphic loci
  gl.filter.rdepth(., 7, 50) %>% # Filter loci with read depth outside 7-50
  gl.filter.callrate(., threshold = 0.9, method = 'loc') # Filter loci with call rate below 90%

nLoc(glSNAPP) # Display the number of loci after further filtering
nInd(glSNAPP) # Display the number of individuals after further filtering
levels(pop(glSNAPP))  #Display the population retained after filtering

##1.2 - Check to make sure you haven't removed EVERYTHING####
pc <- gl.pcoa(glSNAPP, nfactors = 5) # Perform Principal Coordinate Analysis (PCoA)
gl.pcoa.plot(pc, glSNAPP, pop = 'pop', xaxis = 1, yaxis = 2) # Plot the first two principal coordinates

##1.3 Export the nexus file for beauti and beast####
gl2snapp(glSNAPP, 'Ai.nex', outpath = 'data/derived', verbose = 2) # Export data to a Nexus file for SNAPP analysis

#2 - IQtree####
##2.1 Filtering####
glML <- gl %>% 
  gl.filter.reproducibility(., threshold = .99) %>% # Filter loci with reproducibility below 99%
  gl.filter.monomorphs(.) %>% # Remove monomorphic loci
  gl.filter.rdepth(., 7, 50) %>%  # Filter loci with read depth outside 7-50
  gl.filter.callrate(., threshold = 0.9, method = 'loc') # Filter loci with call rate below 90%

nLoc(glML) # Display the number of loci after filtering
nInd(glML) # Display the number of individuals after filtering
levels(pop(glML)) #Display the population retained after filtering

gl2fasta(glSNAPP, method = 3, outfile = 'gl_drop_indivs.fasta', outpath = './', probar = T, verbose = T) # Export filtered data to a FASTA file
gl2fasta(glML, method = 3, outfile = 'gl_all_indivs.fasta', outpath = './', probar = T, verbose = T) # Export all data to a FASTA file

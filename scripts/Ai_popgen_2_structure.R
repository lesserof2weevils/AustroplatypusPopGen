#Title: Population genetics and phylogeny of Austroplatypus incompertus
#Script purpose: Structure analyses
#Author: James Bickerstaff & Markus Riegler
#Date: 22/2/2021

#Packages to call####
libs <- c("dartR", "adegenet", "LEA", "conStruct", "parallel", "fields",
          "maps", "ecodist", "hierfstat", "StAMPP", "ade4", "poppr",
          "netview", "visreg", "devtools", "ggplot2", "tidyverse",
          "doBy", "plotly", "reshape", "ggspatial", "ozmaps", "sf",
          "ggrepel", "patchwork", 'viridis', 'raster', 'ggnewscale',
          'patchwork', 'ggrepel', 'scatterpie', 'tidyverse') # Define a list of required packages

# libraries we need
installed_libs <- libs %in% rownames(installed.packages()) # Check which packages are installed
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs]) # Install missing packages
}

# load libraries
invisible(lapply(libs, library, character.only = T)) # Load all required libraries silently

setwd("~/Documents/GitHub/AustroplatypusPopGen/") # Set the working directory

#1 - Filtering SNPs with DartR####

##1.1. Read in dartseq SNPs and metafile data####
gl <- gl.read.dart(filename = "data/Report_DAmb20-5632_1_moreOrders_SNP_2_rename_ordered.csv", 
                   ind.metafile = "data/Ai_ind_metrics.csv") # Read SNP data and individual metadata into a genlight object
gl # Display the genlight object

##1.2. Report loci, individuals, and populations of the dartseq dataset####
nLoc(gl) # Number of loci
nInd(gl) # Number of individuals
nPop(gl) # Number of populations
levels(pop(gl)) # List populations

##1.2. Filtering SNPs that may skew later analyses####
glx <- gl %>% 
  gl.filter.reproducibility(threshold = 0.95) %>%  # Filter loci with reproducibility below 95%
  gl.filter.monomorphs() %>%  # Remove monomorphic loci
  # gl.filter.secondaries() %>%  # (Commented out) Remove secondary SNPs within the same tag
  gl.filter.maf(threshold = 0.02, verbose = 2) %>%  # Filter loci with minor allele frequency below 2%
  # gl.filter.hwe(alpha = .05, basis = 'any', bon = TRUE, verbose = 2) %>%  # (Commented out) Filter loci deviating from Hardy-Weinberg equilibrium
  gl.filter.hamming(threshold = 0.2, rs = 5, pb = TRUE, verbose = 2) %>%  # Filter based on Hamming distance
  gl.filter.callrate(threshold = 0.7, method = 'loc') %>%  # Filter loci with call rate below 70%
  gl.filter.callrate(method = 'ind', threshold = 0.4, recalc = TRUE)  # Filter individuals with call rate below 40%

gl.report.sexlinked(glx)  # Generate a report on sex-linked loci

ind.callrate <- NA.posi(glx)  # Identify positions of missing data

NA.posi(glx)  # Display positions of missing data

nInd(glx)  # Number of individuals after filtering
gl_allfilt <- nLoc(glx)  # Number of loci after filtering
nPop(glx)  # Number of populations after filtering
levels(pop(glx))  # List populations after filtering
glx$ind.names  # Names of individuals

saveRDS(glx, file = 'glx.rds')  # Save the filtered genlight object

##1.3 Import hierarchy into genlight file, and convert into genind and genclone for downstream analysis####
Aihierarchy <- read.csv('data/derived/hierarchy_inds.csv')  # Read hierarchy data
glx$strata <- Aihierarchy  # Assign strata to genlight object
head(strata(glx, ~ region/pop, combine = FALSE))  # Preview strata

glxi <- gl2gi(glx)  # Convert genlight object to genind object
glxi$strata <- Aihierarchy  # Assign strata to genind object
names(other(glxi))  # List additional data in genind object
head(strata(glxi, ~region/pop, combine = FALSE))  # Preview strata

saveRDS(glxi, file = 'glxi.rds')  # Save the genind object

glxclone <- as.genclone(glxi)  # Convert genind object to genclone object

##1.4. Various visualisations of data without model constraints####
pc <- gl.pcoa(glx, nfactors = 5)  # Perform Principal Coordinate Analysis (PCoA)
p1 <- gl.pcoa.plot(pc, glx, pop.labels = 'pop', xaxis = 1, yaxis = 2)  # Plot PCoA
gl.pcoa.plot(pc, glx, interactive = TRUE, xaxis = 1, yaxis = 2)  # Interactive PCoA plot

gl.tree.nj(glx, type = 'unrooted')  # Generate an unrooted Neighbor-Joining tree

#2 - PCAs, Discriminant Analysis of Variance & Fst####
##2.1 PCA####
x <- tab(glx, NA.method = 'mean')  # Create a genotype matrix, imputing missing values with mean allele frequencies
pca1 <- dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)  # Perform PCA
percent = pca1$eig / sum(pca1$eig) * 100  # Calculate variance explained by each principal component
barplot(percent, ylab = 'Genetic variance explained by eigenvectors (%)', ylim = c(0, 12),
        names.arg = round(percent, 1))  # Plot the variance explained

ind_coords = as.data.frame(pca1$li)  # Extract individual coordinates
colnames(ind_coords) <- c('Axis1', 'Axis2', 'Axis3')  # Rename columns
ind_coords$Ind = indNames(glx)  # Add individual names
ind_coords$population = glx$pop  # Add population information

# Define the order of populations
ord <- c('Dorrigo', 'Kerewong', 'NthBarrington', 'SthBarrington', 'Olney',
         'Cumberland', 'MtWilson', 'MaresForest', 'Termeil', 'Dampier', 'Broadwater')
ind_coords$population <- factor(ind_coords$population, ord)  # Set factor levels

# Assign regions based on populations
ind_coords <- mutate(ind_coords, region = case_when(
  population %in% c('Kerewong', 'NthBarrington', 'Dorrigo') ~ 'North',
  population %in% c('SthBarrington', 'Olney', 'Cumberland', 'MtWilson', 'MaresForest') ~ 'Central',
  population %in% c('Termeil', 'Dampier', 'Broadwater') ~ 'South'
))

# Calculate centroids for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ population, data = ind_coords, FUN = mean)
ind_coords = left_join(ind_coords, centroid, by = "population", suffix = c("", ".cen"))  # Merge centroids

# Create labels for axes
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall = 1), " %)", sep = "")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall = 1), " %)", sep = "")

ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, linewidth =1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15))

# Define colors for populations (optional)
population_colors <- c(
  # North – cool blues & teals
  'Dorrigo'        = '#1b9e77',  # teal green
  'Kerewong'       = '#377eb8',  # mid blue
  'NthBarrington'  = '#a6cee3',  # pale blue
  
  # Central – warm reds, purples, and browns
  'SthBarrington'  = '#e41a1c',  # bright red
  'Olney'          = '#984ea3',  # purple
  'Cumberland'     = '#ff7f00',  # orange
  'MtWilson'       = '#a65628',  # dark brown
  'MaresForest'    = '#f781bf',  # pink
  
  # South – green/yellow tones
  'Termeil'        = '#4daf4a',  # green
  'Dampier'        = '#cce61d',  # lime
  'Broadwater'     = '#006d2c'   # dark green
)

# Scatter plot for PCA
p1 <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = population), shape = 21, size = 3, show.legend = TRUE) +  # Plot individuals
  stat_ellipse(aes(group = region), level = 0.99, color = "black", linetype = "solid", type = 'norm') +  # Add ellipses
  scale_fill_viridis_d(direction = -1) +  # Color scale
  labs(x = xlab, y = ylab) +
  theme_classic(base_size = 17) + 
  guides(fill = 'none') + 
  labs(title = 'A)')

##2.2. Run the dapc####
glx <- readRDS(file = 'glx.rds')  # Load genlight object
glxi <- readRDS(file = 'glxi.rds')  # Load genind object

set.seed(999)  # Set seed for reproducibility
Aix <- xvalDapc(tab(glxi, NA.method = 'mean'), pop(glxi))  # Cross-validation for DAPC
system.time(Aix <- xvalDapc(tab(glxi, NA.method = 'mean'), pop(glxi),
                            n.pca = 10:20, n.rep = 1000,
                            parallel = 'multicore', ncpus = 4L))  # Run DAPC with multiple PCs and replicates
Aix[-1]  # Display results excluding the first element

##2.3 Plot the DAPC####
# Reorder factors from north to south
ord <- c('Dorrigo', 'Kerewong', 'NthBarrington', 'SthBarrington', 'Olney',
         'Cumberland', 'MtWilson', 'MaresForest', 'Termeil', 'Dampier', 'Broadwater')
Aix$DAPC$assign <- factor(Aix$DAPC$assign, levels = ord)
Aix$DAPC$grp <- factor(Aix$DAPC$grp, levels = ord)

# Create a dataframe for plotting
dapc_data_df <- 
  as_tibble(Aix$DAPC$ind.coord, rownames = "individual") %>% 
  mutate(population = factor(glx@other$ind.metrics$pop, levels = ord),
         region = case_when(
           population %in% c('Kerewong', 'NthBarrington', 'Dorrigo') ~ 'North',
           population %in% c('SthBarrington', 'Olney', 'Cumberland', 'MtWilson', 'MaresForest') ~ 'Central',
           population %in% c('Termeil', 'Dampier', 'Broadwater') ~ 'South'
         ))

# Calculate centroids for regions
centroids <- dapc_data_df %>%
  group_by(region) %>%
  summarize(LD1 = mean(LD1), LD2 = mean(LD2))

# Create the DAPC plot
dapc_plot <- ggplot(dapc_data_df, aes(x = LD1, y = LD2, fill = population)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_viridis_d(direction = -1) +  # Color scale
  stat_ellipse(aes(group = region), level = 0.99, color = "black", linetype = "solid", type = 'norm') +
  theme_classic(base_size = 16) + labs(title = "B)")

## 2.4 Print the PCA & DAPC plot
print(dapc_plot)  # Display the DAPC plot
p2 <- scatter(Aix$DAPC, cex = 1, legend = TRUE, clabel = FALSE, posi.leg = 'bottomright',
              scree.pca = FALSE, cleg = 0.75, xax = 1, yax = 2, inset.solid = 1, posi.da = 'topright')  # Alternative scatter plot

svg("outputs/Figure_1_PCA_DAPC.svg")  # Open SVG device
p1 + dapc_plot  # Combine PCA and DAPC plots
dev.off()  # Close the device

plot <- p1 + dapc_plot  # Combine plots
plot  # Display the combined plot

ggsave("Figure_1_PCA_DAPC.svg", plot = plot, device = 'svg', path = './outputs', dpi = 300)  # Save the plot
dev.off()  # Close any open graphic devices

## 2.5 Pairwise Fst across populations
Fst <- gl.fst.pop(glx, nboots = 1000, percent = 95, nclusters = 1, verbose = NULL)  # Calculate pairwise Fst values
fst_pair <- as.matrix(as.dist(Fst$Fsts))  # Convert Fst values to a matrix

# Desired order for populations
ord <- c('Dorrigo', 'Kerewong', 'NthBarrington', 'SthBarrington', 'Olney',
         'Cumberland', 'MtWilson', 'MaresForest', 'Termeil', 'Dampier', 'Broadwater')
ord_inverted <- rev(ord) # Invert the order for bottom triangle

# Reorder the matrix
Fst_ordered <- fst_pair[ord, ord]

# Convert the matrix to a long format
Fst_long <- as.data.frame(fst_pair) %>%
  mutate(X1 = rownames(.)) %>%
  pivot_longer(
    cols = -X1,
    names_to = "X2",
    values_to = "value"
  )
# Convert Var1 and Var2 to factors with the specified order
Fst_long$X1 <- factor(Fst_long$X1, levels = ord_inverted)
Fst_long$X2 <- factor(Fst_long$X2, levels = ord_inverted)

# Rename populations for clarity
Fst_long <- Fst_long %>%
  mutate(X1 = recode(X1,
                     "Dorrigo" = "Dorrigo (n)",
                     "Kerewong" = "Kerewong (n)",
                     "NthBarrington" = "NthBarrington (n)",
                     "SthBarrington" = "SthBarrington (c)",
                     "Olney" = "Olney (c)",
                     "Cumberland" = "Cumberland (c)",
                     "MtWilson" = "MtWilson (c)",
                     "MaresForest" = "MaresForest (c)",
                     "Termeil" = "Termeil (s)",
                     "Dampier" = "Dampier (s)",
                     "Broadwater" = "Broadwater (s)"),
         X2 = recode(X2,
                     "Dorrigo" = "Dorrigo (n)",
                     "Kerewong" = "Kerewong (n)",
                     "NthBarrington" = "NthBarrington (n)",
                     "SthBarrington" = "SthBarrington (c)",
                     "Olney" = "Olney (c)",
                     "Cumberland" = "Cumberland (c)",
                     "MtWilson" = "MtWilson (c)",
                     "MaresForest" = "MaresForest (c)",
                     "Termeil" = "Termeil (s)",
                     "Dampier" = "Dampier (s)",
                     "Broadwater" = "Broadwater (s)"))

# Create the heatmap of Fst values
fst_heatmap <- ggplot(Fst_long, aes(X1, X2, fill = value)) + 
  geom_tile() +
  geom_text(aes(label = round(value, 3)), color = "white") +  # Add Fst values as text
  scale_fill_viridis(name = "Fst", direction = -1) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Population", y = "Population")

print(fst_heatmap)  # Display the heatmap

ggsave("Figure_2_Fst_heatmap.svg", plot = fst_heatmap, device = 'svg', path = './outputs', dpi = 300)  # Save the heatmap
dev.off()  # Close any open graphic devices

#3 - Testing for Isolation by Distance####
# 3 - Testing for Isolation by Distance
## 3.1. Make geographic and genetic matrices
popDist <- as.data.frame(glx$other$latlon, stringsAsFactors = FALSE)  # Extract latitude and longitude data from the genlight object and convert to a data frame
popDist <- rdist.earth(popDist, popDist, miles = FALSE, R = NULL)     # Compute pairwise geographic distances (in kilometers) between all individuals
popDist[1:5, 1:5]                                                     # Display the first 5 rows and columns of the geographic distance matrix

glxi <- gl2gi(glx)                                                    # Convert the genlight object to a genind object
genDist <- diss.dist(glxi)                                            # Compute pairwise genetic distances (Nei's distance) between individuals
genDist                                                               # Display the genetic distance matrix
genDist[1:5, 1:5]                                                     # Display the first 5 rows and columns of the genetic distance matrix
Gdis.vec <- lower(genDist)                                            # Extract the lower triangle of the genetic distance matrix as a vector
Edis.vec <- lower(popDist)                                            # Extract the lower triangle of the geographic distance matrix as a vector
plot(Gdis.vec ~ Edis.vec,
     xlab = 'Geographic distance (km)',
     ylab = 'Genetic distance')                                       # Plot genetic distance against geographic distance

## 3.2. Testing for Spatial Autocorrelation using Mantel's R
SA <- lm(Gdis.vec ~ Edis.vec)                                         # Perform a linear regression of genetic distance on geographic distance
summary(SA)                                                           # Display a summary of the linear model
par(mfrow = c(1,1))                                                   # Set the plotting area to 1 row and 1 column
visreg(SA, xvar = 'Edis.vec', overlay = TRUE,
       xlab = "Geographic distance (km)",
       ylab = "Genetic distance (Nei's)")                             # Visualize the regression with confidence intervals

Mantelr <- ecodist::mantel(Gdis.vec ~ Edis.vec, nperm = 1e6, nboot = 5e5, 
                           pboot = 0.9, cboot = 0.95)                 # Perform a Mantel test with permutations and bootstrapping
Mantelr                                                              # Display the results of the Mantel test

dens <- kde2d(Edis.vec, Gdis.vec, n = 300)                            # Compute a 2D kernel density estimate
myPal <- colorRampPalette(c('white', 'blue', 'gold', 'orange', 'red'))# Create a custom color palette
plot(Edis.vec, Gdis.vec, pch = 20, cex = 0.5)                         # Plot the scatter plot of distances
image(dens, col = transp(myPal(300), 0.5), add = TRUE)                # Overlay the density estimate on the scatter plot
abline(SA)                                                            # Add the linear regression line to the plot

IBD_plot <- ggplot(SA, aes(x = Edis.vec, y = Gdis.vec)) +
  geom_jitter(alpha = 0.5, width = 0.2) +                             # Add jittered points to reduce overplotting
  geom_smooth(method = 'lm') +                                        # Add a linear regression line with confidence intervals
  theme_classic(base_size = 16) +
  ylab("Genetic distance (Nei's)") +
  xlab('Geographic distance (km)')                                    # Customize the plot appearance

ggsave("IBD.svg", plot = IBD_plot, device = 'svg', path = './', dpi = 400)  # Save the plot as an SVG file

# 4 - Inferring Continuous and Discrete Population Genetic Structure Across Space with conStruct ####

## 4.1 Set up population and name lines; STRUCTURE takes numeric population IDs #####
drop_sibs <- as.matrix(read.csv('data/drop_sibs.csv'))  # Read a CSV file containing individuals to drop (e.g., siblings) and convert it to a matrix

glxunique <- gl.drop.ind(glxunique, ind.list = drop_sibs, recalc = TRUE)  # Remove specified individuals from the genlight object and recalculate statistics

glxunique <- gl.filter.callrate(glxunique, method = 'ind', threshold = 0.8, recalc = TRUE)  # Filter out individuals with a call rate below 80%
glxunique <- gl.filter.callrate(glxunique, threshold = 0.75, method = 'loc')  # Filter out loci with a call rate below 75%

save(glxunique, file = 'data/derived/glxunique.RData')  # Save the filtered genlight object for future use

load('data/derived/glxunique.RData')  # Load the genlight object from the saved file

levels(pop(glxunique))  # Display the levels of the population factor in the genlight object
glxunique$ind.names     # Display the names of individuals in the genlight object

pop_lvls <- levels(glxunique@pop)                         # Get the population levels
Pp <- glxunique@pop                                       # Get the population assignments for each individual
POP <- as.numeric(factor(Pp, levels = pop_lvls))          # Convert population factors to numeric IDs
ind.names <- glxunique$ind.names                          # Retrieve individual names

xy <- as.matrix(glxunique$other$latlong)                  # Extract latitude and longitude coordinates as a matrix
popDist <- rdist.earth(xy, xy, miles = FALSE, R = NULL)   # Calculate pairwise geographic distances between individuals

## 4.8 Cross Validation ####

# Cross validation for the non-spatial models
library(LEA)
# Make a .geno file from genlight object
gl2geno(glxunique, outfile = 'glb_2_geno', outpath = './', verbose = NULL)

# Run non-spatial cross-validation
glb_2_nonspatial = snmf('glb_2_geno.lfmm',
                        K = 1:11,
                        entropy = TRUE,
                        repetitions = 5,
                        project = 'new',
                        alpha = 10)

glb_2_lfmm <- read.lfmm('./glb_2_geno.lfmm')
replace_9_with_NA <- function(x) {
  x[x == 9] <- NA
  return(x)
}
glb_2_lfmm_NA <- apply(glb_2_lfmm, 2, replace_9_with_NA)

all_glx.snmf = snmf("glb_2_geno.lfmm", 
                    K = 1:11, 
                    entropy = T, 
                    ploidy = 2, 
                    project="new", 
                    repetitions = 100)  

# Save and reload to be safe
saveRDS(all_glx.snmf, file = 'glb_2_nonspatial.rds')
all_glxunique.nonspatial <- readRDS('glb_2_nonspatial.rds')

# Cross-validation for the spatial models
library(tess3r)

# Run tess3r (spatial)
project <- tess3(X = glxunique, coord = xy, K = 1:11, method = 'projected.ls', ploidy = 2, max.iteration = 100, rep = 100)

# Save spatial project
saveRDS(project, file = "glb_2_spatial.rds")
project <- readRDS("glb_2_spatial.rds")

# Combine cross-validation results and plot with ggplot2
library(ggplot2)
library(dplyr)

## Define range of K and number of repetitions
K_vals <- 1:11
n_reps <- 5  # Change to however many you ran

## Extract from SNMF (non-spatial)
cv_nonspatial <- expand.grid(K = K_vals, rep = 1:n_reps) |>
  rowwise() |>
  mutate(value = cross.entropy(all_glxunique.nonspatial, K = K, run = rep),
         method = "sNMF") |>
  ungroup()

## Extract from tess3r (spatial)
K_spatial <- 1:11  # Or whatever range you ran
n_rep_spatial <- 5  # Again, adjust if needed

cv_spatial <- expand.grid(K = K_spatial, rep = 1:n_rep_spatial) |>
  rowwise() |>
  mutate(value = list(cross.entropy(project, K = K, run = rep)),
         method = "spatial") |>
  ungroup() |>
  mutate(value = unlist(value))

# Make sure both 'value' columns are numeric (not list-columns)
cv_nonspatial <- cv_nonspatial |> 
  mutate(value = as.numeric(value))

# Extract crossentropy values for each K
cv_spatial <- purrr::map_dfr(
  .x = seq_along(project),
  .f = function(k) {
    data.frame(
      K = k,
      rep = 1:length(project[[k]]$crossentropy),
      value = project[[k]]$crossentropy,
      method = "tess3r"
    )
  }
)

# Now combine
cv_all <- bind_rows(cv_nonspatial, cv_spatial)

# Plot with ggplot
ggplot(cv_all, aes(x = K, y = value, color = method)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  stat_summary(fun = mean, geom = "line", aes(group = method), size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 2.5) +
  theme_minimal() +
  labs(title = "Cross-validation scores",
       subtitle = "Spatial vs Non-Spatial models",
       x = "Number of ancestral populations (K)",
       y = "Cross-entropy",
       color = "Model type") +
  scale_x_continuous(breaks = unique(cv_all$K))

ggplot(cv_all, aes(x = K, y = value, color = method)) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = FALSE, method = "loess", span = 0.75) +
  labs(title = "Cross-validation comparison",
       x = "K-values",
       y = "Cross-entropy",
       color = "Analysis") +
  theme_minimal()

# Save as SVG
ggsave("cross_validation_comparison.svg", width = 6, height = 5)

# retrieve tess3 Q matrix for K = 4 clusters 
q.matrix <- qmatrix(project, K = 4)

k4_tess <- as.data.frame(q.matrix)
k4_tess$ind <- glxunique$ind.names
k4_tess$locality <- glxunique$pop    
k4_tess$lat <- glxunique$other$latlon$lat
k4_tess$lat<- as.numeric(gsub("−", "-", as.character(glxunique$other$latlon$lat)))
k4_tess$lon <- glxunique$other$latlon$lon
k4_tess <- k4_tess[,c(5:6, 1:4, 7:8)]
colnames(k4_tess) <- c("ind", "pop", "P1", "P2", "P3", "P4", "lat", "lon")
write.csv(k4_tess, "all_k4_tess_r68_latlong.csv", row.names = F)

k4_tess_df <- k4_tess[3:7] %>% 
  as_tibble() %>% 
  mutate(individual = k4_tess$ind,
         locality = k4_tess$pop)

k4_tess_df_long <- k4_tess_df %>% 
  pivot_longer(cols = starts_with("P"),
               names_to = 'pop',
               values_to = 'q')

k4_tess_df_ordered <- k4_tess_df_long %>% 
  group_by(individual) %>% 
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>% 
  mutate(lat = -lat) %>% 
  arrange(lat, likely_assignment, assignment_prob) %>% 
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))

k4_tess_df_ordered

k4_q_palette <- c("blue","red","yellow","darkgreen")

tess3rk4 <- k4_tess_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = k4_q_palette) +
  labs(fill = 'ancestry assignment') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = 'none') +
  ylab('Ancestry coefficient') + xlab('Individuals') + ggtitle('tess3r K=4')

# retrieve tess3 Q matrix for K = 5 clusters 
q.matrix5 <- qmatrix(project, K = 5)

k5_tess <- as.data.frame(q.matrix5)
k5_tess$ind <- glxunique$ind.names
k5_tess$pop <- glxunique$pop    
k5_tess$lat <- glxunique$other$latlon$lat
k5_tess$lat<- as.numeric(gsub("−", "-", as.character(glxunique$other$latlon$lat)))
k5_tess$lon <- glxunique$other$latlon$lon
k5_tess <- k5_tess[,c(6:7, 1:5, 8:9)]
colnames(k5_tess) <- c("ind", "pop", "P1", "P2", "P3", "P4", "P5", "lat", "lon")
write.csv(k5_tess, "all_k4_tess_r68_latlong.csv", row.names = F)

k5_tess_df <- k5_tess[3:8] %>% 
  as_tibble() %>% 
  mutate(individual = k5_tess$ind,
         locality = k5_tess$pop)

k5_tess_df_long <- k5_tess_df %>% 
  pivot_longer(cols = starts_with("P"),
               names_to = 'pop',
               values_to = 'q')

k5_tess_df_ordered <- k5_tess_df_long %>% 
  group_by(individual) %>% 
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>% 
  mutate(lat = -lat) %>% 
  arrange(lat, likely_assignment, assignment_prob) %>% 
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))

k5_tess_df_ordered

k5_q_palette <- c("blue","red","darkgreen","orange" ,"yellow")

tess3rk5 <- k5_tess_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = k5_q_palette) +
  labs(fill = 'ancestry assignment') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = 'none') +
  ylab('Ancestry coefficient') + xlab('Individuals') + ggtitle('tess3r K=5')

ggsave("k5_tess3r_plot.svg", plot = last_plot(), device = 'svg', path = './', dpi = 400)

#Draw the non-spatial sNMF plots
ce_4 = cross.entropy(all_glxunique.nonspatial, K = 4)
best_4 <- which.min(ce_4)
ce_5 = cross.entropy(all_glxunique.nonspatial, K = 5)
best_5 <- which.min(ce_5)

#4.3 Draw the plots for K= 4 and export the files
k4 <- as.data.frame(Q(all_glxunique.nonspatial, K = 4, run = best_4))
k4$ind <- glxunique$ind.names
k4$pop <- glxunique$pop    
k4$lat <- glxunique$other$latlon$lat
k4$lat<- as.numeric(gsub("−", "-", as.character(glxunique$other$latlon$lat)))
k4$lon <- glxunique$other$latlon$lon
k4 <- k4[,c(5:6, 1:4, 7:8)]
colnames(k4) <- c("ind", "pop", "P1", "P2", "P3", "P4", "lat", "lon")
write.csv(k4, "all_k4_r68_latlong.csv", row.names = F)

k4_df <- k4[3:7] %>% 
  as_tibble() %>% 
  mutate(individual = k4$ind,
         locality = k4$pop)

k4_df_long <- k4_df %>% 
  pivot_longer(cols = starts_with("P"),
               names_to = 'pop',
               values_to = 'q')
k4_df_ordered <- k4_df_long %>% 
  group_by(individual) %>% 
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>% 
  mutate(lat = -lat) %>% 
  arrange(lat, likely_assignment, assignment_prob) %>% 
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))

k4_df_ordered

k4_q_palette <- c("darkgreen","red","blue","yellow")

snmfk4 <- k4_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = k4_q_palette) +
  labs(fill = 'ancestry assignment') +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = 'none') +
  ylab('Ancestry coefficient') + xlab('Individuals') + ggtitle('snmf K=4')

ggsave("k4_snmf_plot.svg", plot = last_plot(), device = 'svg', path = './', dpi = 400)

#Draw the plot for K = 5
k5 <- as.data.frame(Q(all_glxunique.nonspatial, K = 5, run = best_5))
k5$ind <- glxunique$ind.names
k5$pop <- glxunique$pop    
k5$lat <- glxunique$other$latlon$lat
k5$lat<- as.numeric(gsub("−", "-", as.character(glxunique$other$latlon$lat)))
k5$lon <- glxunique$other$latlon$lon
k5 <- k5[,c(6:7, 1:5, 8:9)]
colnames(k5) <- c("ind", "pop", "P1", "P2", "P3", "P4", "P5", "lat", "lon")
write.csv(k5, "all_k5_latlong.csv", row.names = F)

k5_df <- k5[3:8] %>% 
  as_tibble() %>% 
  mutate(individual = k5$ind,
         locality = k5$pop)

k5_df_long <- k5_df %>% 
  pivot_longer(cols = starts_with("P"),
               names_to = 'pop',
               values_to = 'q')

k5_df_ordered <- k5_df_long %>% 
  group_by(individual) %>% 
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>% 
  mutate(lat = -lat) %>% 
  arrange(lat, likely_assignment, assignment_prob) %>% 
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))

k5_df_ordered

k5_q_palette <- c("orange","darkgreen","yellow","blue","red")

snmfk5 <- k5_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = k5_q_palette) +
  labs(fill = 'ancestry assignment') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = 'none') +
  ylab('Ancestry coefficient') + xlab('Individuals') + ggtitle('snmf K=5')

ggsave("k5_snmf_plot.svg", plot = last_plot(), device = 'svg', path = './', dpi = 400)

library(patchwork)

# Remove legends and x-axis text/labels from all four plots
tess3rk4_clean <- tess3rk4

tess3rk5_clean <- tess3rk5 +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank())

snmfk4_clean <- snmfk4 +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank())

snmfk5_clean <- snmfk5 +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank())

# Combine vertically
combined_plot <-snmfk5_clean / snmfk4_clean/ tess3rk5_clean / tess3rk4_clean  +
  plot_layout(ncol = 1)

# Save the plot
ggsave("combined_ancestry_barplots_clean.svg",
       plot = combined_plot,
       device = "svg",
       width = 10, height = 12, dpi = 300)

## 4.2 Convert genlight file into a STRUCTURE file to be used by conStruct ####
gl2structure(glxunique, indNames = ind.names, addcolumns = c(POP), ploidy = 2, exportMarkerNames = FALSE, 
             outfile = 'glxunique.str', outpath = 'data/derived')  # Export the genlight object to a STRUCTURE format file for conStruct analysis

## 4.3 Construct the conStruct file from the STRUCTURE file ####
Ai.construct <- structure2conStruct(
  infile = 'data/derived/glxunique.str',    # Input STRUCTURE format file
  onerowperind = FALSE,                      # Data is formatted with two rows per individual
  start.loci = 3,                            # Column index where loci data starts
  missing.datum = -9,                        # Indicator for missing data in the file
  outfile = 'Ai_construct'      # Output file name for conStruct
)  # Convert STRUCTURE file to conStruct format

Ai.construct <- readRDS('data/derived/Ai_construct.RData')

Ai.construct[1:5, 2:10]  # Display a subset of the converted data for verification

# Set options for parallel processing and Stan configuration
options(mc.cores = parallel::detectCores())  # Use all available cores for parallel processing
rstan_options(auto_write = TRUE)             # Allow Stan to cache compiled models on disk
options(mc.cores = 5)                        # Limit the number of cores to 5
options()                                    # Display current global options

## 4.4 Run spatial analyses for K = 1 to 5 ####
# Run conStruct spatial models for different numbers of layers (K)
Aispk1 <- conStruct(
  spatial = TRUE,                            # Enable spatial modeling
  K = 1,                                     # Number of ancestral layers
  freqs = Ai.construct,                      # Allele frequency data
  geoDist = popDist,                         # Geographic distance matrix
  coords = xy,                               # Coordinates of individuals
  n.iter = 100000,                           # Number of MCMC iterations
  n.chains = 10,                             # Number of MCMC chains
  make.figs = TRUE,                          # Generate figures during the run
  save.files = TRUE,                         # Save output files
  prefix = 'spk1',                           # Prefix for output files
  control = setNames(list(0.95), "adapt_delta")  # Adjust Stan's adapt_delta parameter for convergence
)

Aispk2 <- conStruct(
  spatial = TRUE,
  K = 2,
  freqs = Ai.construct,
  geoDist = popDist,
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  make.figs = TRUE,
  save.files = TRUE,
  prefix = 'spk2',
  control = setNames(list(0.95), "adapt_delta")
)

Aispk3 <- conStruct(
  spatial = TRUE,
  K = 3,
  freqs = Ai.construct,
  geoDist = popDist,
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  make.figs = TRUE,
  save.files = TRUE,
  prefix = 'spk3'
)

Aispk4 <- conStruct(
  spatial = TRUE,
  K = 4,
  freqs = Ai.construct,
  geoDist = popDist,
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  make.figs = TRUE,
  save.files = TRUE,
  prefix = 'spk4'
)

Aispk5 <- conStruct(
  spatial = TRUE,
  K = 5,
  freqs = Ai.construct,
  geoDist = popDist,
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  make.figs = TRUE,
  save.files = TRUE,
  prefix = 'spk5'
)

## 4.5 Run non-spatial analyses for K = 1 to 5 ####
# Run conStruct non-spatial models (without geographic data)
nspk1 <- conStruct(
  spatial = FALSE,                           # Disable spatial modeling
  K = 1,
  freqs = Ai.construct,
  geoDist = NULL,                            # No geographic distances used
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  prefix = 'nspk1'
)

nspk2 <- conStruct(
  spatial = FALSE,
  K = 2,
  freqs = Ai.construct,
  geoDist = NULL,
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  prefix = 'nspk2',
  make.figs = TRUE,
  save.files = TRUE
)

nspk3 <- conStruct(
  spatial = FALSE,
  K = 3,
  freqs = Ai.construct,
  geoDist = NULL,
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  prefix = 'nspk3',
  make.figs = TRUE,
  save.files = TRUE
)

nspk4 <- conStruct(
  spatial = FALSE,
  K = 4,
  freqs = Ai.construct,
  geoDist = NULL,
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  prefix = 'nspk4',
  make.figs = TRUE,
  save.files = TRUE
)

nspk5 <- conStruct(
  spatial = FALSE,
  K = 5,
  freqs = Ai.construct,
  geoDist = NULL,
  coords = xy,
  n.iter = 100000,
  n.chains = 10,
  prefix = 'nspk5',
  make.figs = TRUE,
  save.files = TRUE
)

## 4.6 Model Contributions ####
### 4.6.1 Spatial model #####
Ai_construct <- load("data/derived/Ai_construct.RData")  # Load the Ai_construct data from an RData file

layer.contributionsS <- matrix(NA, nrow = 5, ncol = 5)   # Initialize a 5x5 matrix to store layer contributions

load("outputs/conStruct/spk1_conStruct.results.Robj")    # Load conStruct results for spatial model K=1
load("outputs/conStruct/spk1_data.block.Robj")           # Load the corresponding data block

# Calculate layer contributions for K=1 and store in the first column
layer.contributionsS[,1] <- c(
  calculate.layer.contribution(conStruct.results[[1]], data.block),
  rep(0, 4)  # Fill remaining rows with zeros to match matrix dimensions
)

tmpS <- conStruct.results[[1]]$MAP$admix.proportions     # Extract admixture proportions for matching layers

# Loop over K values from 2 to 5 to calculate layer contributions
for (i in 2:5) {
  load(sprintf("outputs/conStruct/spk%s_conStruct.results.Robj", i))  # Load conStruct results for spatial model K=i
  load(sprintf("outputs/conStruct/spk%s_data.block.Robj", i))         # Load the corresponding data block
  
  # Match layers between runs to ensure consistency
  tmpS.order <- match.layers.x.runs(tmpS, conStruct.results[[1]]$MAP$admix.proportions)
  
  # Calculate layer contributions for current K and store in matrix
  layer.contributionsS[, i] <- c(
    calculate.layer.contribution(
      conStruct.results = conStruct.results[[1]],
      data.block = data.block,
      layer.order = tmpS.order
    ),
    rep(0, 5 - i)  # Pad with zeros to match matrix dimensions
  )
  
  # Update tmpS with ordered admixture proportions for the next iteration
  tmpS <- conStruct.results[[1]]$MAP$admix.proportions[, tmpS.order]
}

# Create a barplot of layer contributions for spatial models
p1 <- barplot(
  layer.contributionsS,
  col = c("blue", "red", "gold", "forestgreen", "darkorchid1"),  # Define colors for layers
  xlab = "Spatial Model K-values",
  ylab = "Layer Contributions",
  names.arg = paste0("K=", 1:5)  # Label x-axis with K values
)

### 4.6.2 Non-spatial model #####
layer.contributionsS <- matrix(NA, nrow = 5, ncol = 5)  # Initialize a 5x5 matrix for non-spatial layer contributions

load("outputs/conStruct/nspk1_conStruct.results.Robj")   # Load conStruct results for non-spatial model K=1
load("outputs/conStruct/nspk1_data.block.Robj")          # Load the corresponding data block

# Calculate layer contributions for K=1 and store in the first column
layer.contributionsS[,1] <- c(
  calculate.layer.contribution(conStruct.results[[1]], data.block),
  rep(0, 4)
)

tmpS <- conStruct.results[[1]]$MAP$admix.proportions     # Extract admixture proportions

# Loop over K values from 2 to 5 for non-spatial models
for (i in 2:5) {
  load(sprintf("outputs/conStruct/nspk%s_conStruct.results.Robj", i))  # Load conStruct results for non-spatial model K=i
  load(sprintf("outputs/conStruct/nspk%s_data.block.Robj", i))         # Load the corresponding data block
  
  # Match layers between runs
  tmpS.order <- match.layers.x.runs(tmpS, conStruct.results[[1]]$MAP$admix.proportions)
  
  # Calculate layer contributions and store in matrix
  layer.contributionsS[, i] <- c(
    calculate.layer.contribution(
      conStruct.results = conStruct.results[[1]],
      data.block = data.block,
      layer.order = tmpS.order
    ),
    rep(0, 5 - i)
  )
  
  # Update tmpS for next iteration
  tmpS <- conStruct.results[[1]]$MAP$admix.proportions[, tmpS.order]
}

# Create a barplot of layer contributions for non-spatial models
p2 <- barplot(
  layer.contributionsS,
  col = c("blue", "red", "gold", "forestgreen", "darkorchid1"),
  xlab = "Non-spatial Model K-values",
  ylab = "Layer Contributions",
  names.arg = paste0("K=", 1:5)
)

## 4.7 Making the admixture plots ####
### 4.7.1 Spatial Models #####

# Load the conStruct results and data block for spatial model K=2
load('outputs/conStruct/spk2_conStruct.results.Robj')
load('outputs/conStruct/spk2_data.block.Robj')

# Extract the admixture proportions from the MAP estimates of the first MCMC chain
admix.props <- conStruct.results$chain_1$MAP$admix.proportions

# Match layers across runs (here, matching within the same run)
match.layers.x.runs(admix.props, admix.props)

# Generate a structure-like plot of the admixture proportions
make.structure.plot(
  admix.proportions = admix.props,
  sample.order = order(data.block$coords[,1]),       # Order samples based on their coordinates (latitude)
  sample.names = row.names(data.block$coords),       # Use sample names from the data block
  mar = c(4.5, 4, 2, 2)                              # Set plot margins
)

# Convert admixture proportions to a data frame
spk4 <- as.data.frame(admix.props)

# Add individual names to the data frame
spk4$ind <- rownames(data.block$coords)

# Create a data frame with individual names and their coordinates
coords <- data.frame(ind = rownames(data.block$coords), data.block$coords)

# Create a data frame with population assignments
pops <- data.frame(pop = glx$pop)
pops$ind <- glx$ind.names

# Merge admixture proportions with coordinates based on individual names
spk4 <- merge(spk4, coords, by = "ind")

# Merge the result with population data
spk4 <- merge(spk4, pops, by = "ind")

# Reorder columns for clarity
spk4 <- spk4[, c(1, 8, 2:5, 6:7)]

# Rename columns for easier interpretation
colnames(spk4) <- c("ind", "pop", "P1", "P2","P3","P4", "lat", "lon")

# Save the combined data frame to a CSV file
write.csv(spk4, "spk2_latlong.csv", row.names = FALSE)

# Select relevant columns and convert to a tibble
spk4_df <- spk4[3:7] %>% 
  as_tibble() %>% 
  mutate(
    individual = spk4$ind,    # Add individual names
    locality = spk4$pop       # Add population names
  )

# Reshape data from wide to long format for plotting
spk4_df_long <- spk4_df %>% 
  pivot_longer(
    cols = starts_with("P"),  # Columns starting with "P" (P1, P2)
    names_to = 'pop',         # New column for population components
    values_to = 'q'           # Admixture proportions
  )

# Order the data based on latitude and assignment probabilities
spk4_df_ordered <- spk4_df_long %>% 
  group_by(individual) %>% 
  mutate(
    likely_assignment = pop[which.max(q)],  # Determine the population with the highest proportion
    assignment_prob = max(q)                # Get the highest assignment probability
  ) %>% 
  mutate(lat = -lat) %>%                    # Negate latitude to arrange from north to south
  arrange(lat, likely_assignment, assignment_prob) %>%  # Arrange data for plotting
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))  # Preserve the order of individuals

# Display the ordered data frame (optional)
spk4_df_ordered

# Define a custom color palette for the plot
spk4_q_palette <- c("red", "blue","gold" , "darkgreen")

# Plot the admixture proportions as a bar plot
spk4_bar <- spk4_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +    # Create stacked bars for each individual
  scale_fill_manual(values = spk4_q_palette) +          # Apply the custom color palette
  labs(fill = 'Ancestry Assignment') +                  # Label for the legend
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), legend.position = 'none') +
  ylab('Admixture') + xlab('Individuals') +             # Axis labels
  ggtitle('Spatial Model K=2')                          # Plot title

# Save the plot as an SVG file with high resolution
ggsave("spk2_barplot.svg", plot = spk2_bar, device = 'svg', path = './', dpi = 400)

# Repeat the process for spatial model K=3

# Load the conStruct results and data block for spatial model K=3
load('outputs/conStruct/spk5_conStruct.results.Robj')
load('outputs/conStruct/spk5_data.block.Robj')

# Extract the admixture proportions
admix.props <- conStruct.results$chain_1$MAP$admix.proportions

# Match layers across runs
match.layers.x.runs(admix.props, admix.props)

# Generate a structure-like plot
make.structure.plot(
  admix.proportions = admix.props,
  sample.order = order(data.block$coords[,1]),
  sample.names = row.names(data.block$coords),
  mar = c(4.5, 4, 2, 2)
)

# Convert admixture proportions to a data frame
spk5 <- as.data.frame(admix.props)

# Add individual names
spk5$ind <- rownames(data.block$coords)

# Create data frames for coordinates and populations
coords <- data.frame(ind = rownames(data.block$coords), data.block$coords)
pops <- data.frame(pop = glxunique$pop)
pops$ind <- glxunique$ind.names

# Merge data frames to combine all information
spk5 <- merge(spk5, coords, by = "ind")
spk5 <- merge(spk5, pops, by = "ind")

# Reorder columns
spk5 <- spk5[, c(1, 9, 2:6, 7:7)]

# Rename columns
colnames(spk5) <- c("ind", "pop", "P1", "P2", "P3","P4","P5", "lat", "lon")

# Save to CSV
write.csv(spk5, "spk5_latlong.csv", row.names = FALSE)

# Select relevant columns and convert to a tibble
spk5_df <- spk5[3:6] %>% 
  as_tibble() %>% 
  mutate(
    individual = spk3$ind,
    locality = spk3$pop
  )

# Reshape data for plotting
spk3_df_long <- spk3_df %>% 
  pivot_longer(
    cols = starts_with("P"),   # Columns P1, P2, P3
    names_to = 'pop',
    values_to = 'q'
  )

# Order the data based on latitude and assignment probabilities
spk3_df_ordered <- spk3_df_long %>% 
  group_by(individual) %>% 
  mutate(
    likely_assignment = pop[which.max(q)],
    assignment_prob = max(q)
  ) %>% 
  mutate(lat = -lat) %>% 
  arrange(lat, likely_assignment, assignment_prob) %>% 
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))

# Display the ordered data frame (optional)
spk3_df_ordered

# Define a color palette for K=3
spk3_q_palette <- c("blue", "red", "gold")

# Plot the admixture proportions as a bar plot
spk3_bar <- spk3_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = spk3_q_palette) +
  labs(fill = 'Ancestry Assignment') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels
    legend.position = 'none'
  ) +
  ylab('Admixture') + xlab('Individuals') + ggtitle('Spatial Model K=3')

# Save the plot as an SVG file
ggsave("spk3_barplot.svg", plot = spk3_bar, device = 'svg', path = './', dpi = 400)

# Now plot as a map

# Load the boundary of Australia using the geoboundaries package
aus_bound <- rgeoboundaries::geoboundaries('Australia')

# Load Australian state boundaries and transform to WGS84 coordinate system (EPSG:4326)
aus <- st_transform(ozmaps::ozmap_states, 4326)

# Get elevation data for Australia at zoom level 7 and clip to the country's boundaries
elev_data <- elevatr::get_elev_raster(locations = aus_bound, z = 7, clip = 'locations')

# Convert the elevation raster data to a data frame with x and y coordinates, and remove missing values
elev_data <- elev_data %>%
  as.data.frame(xy = TRUE) %>%
  na.omit()

# Rename the third column to 'elevation' for clarity
colnames(elev_data)[3] <- 'elevation'

# Filter the elevation data to include only non-negative elevations within specified longitude and latitude ranges
elev_data_filt <- elev_data %>% 
  filter(elevation >= 0) %>% 
  filter(between(x, 120, 160) & between(y, -44, -10))

# Create a data frame of sampling sites with their names and geographic coordinates
points <- tibble::tribble(
  ~site, ~lat, ~lon,
  'Broadwater SF', -37.0116, 149.9026333,
  'Dampier SF', -36.092783, 149.8720167,
  'Termeil SF', -35.439667, 150.3634167, 
  'Mares Forest SF', -34.2882, 149.9493167,
  'Mt Wilson', -33.495483, 150.4030667,
  'Cumberland SF', -33.743641, 151.040079, 
  'Olney SF', -33.0964, 151.3515,
  'South Barrington Tops NP', -32.156766, 151.324950,
  'North Barrington Tops NP', -31.911183, 151.6354167,
  'Willi Willi NP', -31.211832, 152.415809,
  'Kerewong SF', -31.60764, 152.560996,
  'Dorrigo NP', -30.363919, 152.79882
)
points <- as.data.frame(points)  # Convert the tibble to a data frame

# Create the map using ggplot2
spk3_map <- ggplot() +
  geom_raster(data = elev_data_filt, aes(x = x, y = y, fill = elevation)) +  # Add elevation data as a raster layer
  geom_sf(data = aus, colour = 'black', fill = NA) +  # Overlay Australian state boundaries
  coord_sf(xlim = c(141, 155), ylim = c(-28, -38.7)) + # Set the coordinate limits for the map (longitude and latitude)
  scale_fill_gradient(low = 'white', high = 'black') + # Define the color gradient for elevation
  new_scale_fill() + # Introduce a new fill scale for subsequent layers
  geom_point(data = points, mapping = aes(x = lon, y = lat), color = 'black') +  # Add points for the sampling sites
  geom_text_repel(data = points, aes(x = lon, y = lat, label = site)) +  # Add labels for each sampling site, adjusting to avoid overlap
  # geom_scatterpie(aes(x = lon, y = lat, r = 0.25), data = spk3, cols = c("P1", "P2", "P3")) + # (Optional) Add scatter pie charts representing admixture proportions at each site
  # scale_fill_manual(values = spk3_q_palette) +
  theme(legend.position = 'none') + # Remove the legend from the plot
  labs(x = 'Longitude', y = 'Latitude') +  # Add axis labels
  annotation_scale(location = 'br', width_hint = 0.3, height = unit(0.15, 'cm')) + # Add a scale bar to the bottom right of the map
  annotation_north_arrow(location = 'br', which_north = 'true',
                         pad_x = unit(0.75, 'cm'), pad_y = unit(0.5, 'cm'),
                         style = north_arrow_fancy_orienteering) +  # Add a north arrow to the bottom right of the map
  theme_minimal() # Apply a minimal theme to the plot

# Display the map
spk3_map

# Save the map as an SVG file with high resolution (400 dpi)
ggsave("spk3_map.svg", plot = spk3_map, device = 'svg', path = './', dpi = 400)

### 4.7.2 Non-Spatial Models ####

# Load the conStruct results and data block for the non-spatial model with K = 4
load('outputs/conStruct/nspk4_conStruct.results.Robj')
load('outputs/conStruct/nspk4_data.block.Robj')

# Extract the admixture proportions from the MAP estimates of chain 4
admix.props <- conStruct.results$chain_4$MAP$admix.proportions

# Match layers across runs (here, matching within the same run for consistency)
match.layers.x.runs(admix.props, admix.props)

# Generate a structure-like plot of the admixture proportions
make.structure.plot(
  admix.proportions = admix.props,
  sample.order = order(data.block$coords[, 1]),    # Order samples by latitude
  sample.names = row.names(data.block$coords),     # Use sample names from the data block
  mar = c(4.5, 4, 2, 2)                            # Set plot margins
)

# Convert admixture proportions to a data frame
nspk4 <- as.data.frame(admix.props)

# Add individual identifiers
nspk4$ind <- rownames(data.block$coords)

# Create a data frame with individual coordinates
coords <- data.frame(
  ind = rownames(data.block$coords),
  data.block$coords
)

# Create a data frame with population assignments from the original genlight object
pops <- data.frame(
  pop = glx$pop,
  ind = glx$ind.names
)

# Merge admixture proportions with coordinates and population data
nspk4 <- merge(nspk4, coords, by = "ind")
nspk4 <- merge(nspk4, pops, by = "ind")

# Reorder columns for clarity
nspk4 <- nspk4[, c(1, 8, 2:5, 6:7)]

# Rename columns for easier interpretation
colnames(nspk4) <- c("ind", "pop", "P1", "P2", "P3", "P4", "lat", "lon")

# Save the combined data frame to a CSV file
write.csv(nspk4, "nspk4_latlong.csv", row.names = FALSE)

# Prepare data for plotting by selecting admixture proportions and adding identifiers
nspk4_df <- nspk4[3:7] %>%
  as_tibble() %>%
  mutate(
    individual = nspk4$ind,
    locality = nspk4$pop
  )

glx <- readRDS('glx.rds')

#Cross validation for the non-spatial models
library(LEA)
#4.1 First we need to make a .geno file from out genlight object
gl2geno(glx, outfile = 'glb_2_geno', outpath = './', verbose = NULL)
glb_2_snmf = snmf('glb_2_geno.lfmm',
                  K=1:15,
                  entropy = T,
                  repetitions = 5,
                  project = 'new',
                  alpha = 1000)

glb_2_lfmm <- read.lfmm('./glb_2_geno.lfmm')

replace_9_with_NA <- function(x) {
  x[x == 9] <- NA
  return(x)
}
glb_2_lfmm_NA <- apply(glb_2_lfmm, 2, replace_9_with_NA)

#4.2 First we test with snmf as our null hypothesis (and to match Renee's
#mapping pipeline)
all_glx.snmf = snmf("glb_2_geno.lfmm", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 100)  

saveRDS(all_glx.snmf, file = 'all_glx.snmf.rds')
all_glx.snmf <- readRDS('all_glx.snmf.rds')

CE_plot <- plot(all_glx.snmf, pch = 19, col = 'black',
                main = 'smnf cross entropy')

#Cross validation for the spatial models
library(tess3r)

# Run tess3 for a range of K values
project <- tess3(X = glxunique,
                 coord = xy,
                 K = 1:6,
                 ploidy = 2,
                 rep = 5,             # repeat runs to assess stability
                 max.iteration = 500,
                 openMP.core.num = 1) # use more cores if desired

# Plot cross-entropy to choose K
# Save as SVG
svg("cross_entropy_per_K.svg", width = 6, height = 5)

# Your plot
plot(project, pch = 19, col = "blue", main = "Cross-entropy per K",
     xlab = "Values of K",
     ylab = "Cross-validation score")

# Close the device
dev.off()


# Extract ancestry coefficients for your chosen K
bestK <- which.min(cross.entropy(project))
Qmatrix <- qmatrix(project, K = bestK)

# to run a cross-validation analysis
# you have to specify:
#       the numbers of layers you want to compare (K)
#       the allele frequency data (freqs)
#       the geographic distance matrix (geoDist)
#       the sampling coordinates (coords)

# Transpose to get individuals as rows
Ai.construct <- t(Ai.construct)

# Now assign rownames (individual IDs)
rownames(Ai.construct) <- rownames(glx@other$latlon)[match(rownames(Ai.construct), rownames(glx@other$latlon))]

# Now transpose back (if needed)
Ai.construct <- t(Ai.construct)

# Set column names to the row names (individual IDs)
colnames(Ai.construct) <- rownames(glx@other$latlon)[match(colnames(Ai.construct), rownames(glx@other$latlon))]

# Subset coordinates to match these individuals
xy <- glx@other$latlon[match(colnames(Ai.construct), rownames(glx@other$latlon)), ]

# Check that everything matches now
ncol(Ai.construct) == nrow(xy)

# Ensure coords are in matrix format
xy <- as.matrix(xy)

# Step 1: Manually assign data partitions
parts <- conStruct:::make.data.partitions(
  freqs = Ai.construct,
  train.prop = 0.7,
  n.reps = 3
)
# Step 2: Run x.validation with pre-computed partitions
my.xvals <- x.validation(
  freqs = Ai.construct,
  data.partitions = parts,
  geoDist = popDist,
  coords = as.matrix(xy),
  train.prop = 0.7,
  n.reps = 3,
  K = 1:5,
  prefix = "xval_run",
  n.iter = 1000,
  make.figs = TRUE,
  save.files = TRUE,
  parallel = FALSE
)

for (prop in seq(0.7, 0.5, by = -0.05)) {
  cat("Trying train.prop =", prop, "\n")
  try({
    parts <- conStruct:::make.data.partitions(freqs = Ai.construct, train.prop = prop, n.reps = 3)
    my.xvals <- x.validation(
      freqs = Ai.construct,
      data.partitions = parts,
      geoDist = popDist,
      coords = as.matrix(xy),
      train.prop = prop,
      n.reps = 3,
      K = 1:3,  # use smaller K for testing
      prefix = paste0("xval_test_", prop),
      n.iter = 500,
      make.figs = FALSE,
      save.files = FALSE,
      parallel = FALSE
    )
    cat("✅ Success with train.prop =", prop, "\n")
    break
  }, silent = TRUE)
}

##4.8 Layer Contributions####

###4.8.1 Spatial Models####

# Initialise the matrix
layer.contributions_spatial <- matrix(NA, nrow = 5, ncol = 5)

# Load K = 1 data
load("outputs/conStruct/spk1_conStruct.results.Robj")
load("outputs/conStruct/spk1_data.block.Robj")

# Calculate contributions for K=1
layer.contributions_spatial[, 1] <- c(
  calculate.layer.contribution(conStruct.results[[1]], data.block),
  rep(0, 4)
)
tmp <- conStruct.results[[1]]$MAP$admix.proportions

# Loop from K = 2 to 5
for (i in 2:5) {
  load(sprintf("outputs/conStruct/spk%d_conStruct.results.Robj", i))
  load(sprintf("outputs/conStruct/spk%d_data.block.Robj", i))
  
  tmp.order <- match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
  
  layer.contributions_spatial[, i] <- c(
    calculate.layer.contribution(
      conStruct.results = conStruct.results[[1]],
      data.block = data.block,
      layer.order = tmp.order
    ),
    rep(0, 5 - i)
  )
  
  tmp <- conStruct.results[[1]]$MAP$admix.proportions[, tmp.order]
}

barplot(layer.contributions_spatial,
        col=c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1"),
        xlab="",
        ylab="layer contributions",
        names.arg=paste0("K=",1:5))

###4.8.2 Non-Spatial Models####

# Initialise the matrix
layer.contributions_nonspatial <- matrix(NA, nrow = 5, ncol = 5)

# Load K = 1 data
load("outputs/conStruct/nspk1_conStruct.results.Robj")
load("outputs/conStruct/nspk1_data.block.Robj")

# Calculate contributions for K=1
layer.contributions_nonspatial[, 1] <- c(
  calculate.layer.contribution(conStruct.results[[1]], data.block),
  rep(0, 4)
)
tmp <- conStruct.results[[1]]$MAP$admix.proportions

# Loop from K = 2 to 5
for (i in 2:5) {
  load(sprintf("outputs/conStruct/nspk%d_conStruct.results.Robj", i))
  load(sprintf("outputs/conStruct/nspk%d_data.block.Robj", i))
  
  tmp.order <- match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
  
  layer.contributions_nonspatial[, i] <- c(
    calculate.layer.contribution(
      conStruct.results = conStruct.results[[1]],
      data.block = data.block,
      layer.order = tmp.order
    ),
    rep(0, 5 - i)
  )
  
  tmp <- conStruct.results[[1]]$MAP$admix.proportions[, tmp.order]
}

barplot(layer.contributions_nonspatial,
        col=c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1"),
        xlab="",
        ylab="layer contributions",
        names.arg=paste0("K=",1:5))

# Reshape data from wide to long format for plotting
nspk4_df_long <- nspk4_df %>%
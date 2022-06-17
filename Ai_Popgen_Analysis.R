#Title: Population genetics and phylogeny of Austroplatypus incompertus
#Author: James Bickerstaff & Markus Riegler
#Date: 22/2/2021

#Packages to call####
library(dartR)
library(adegenet)
library(LEA)
library(conStruct)
library(rstan)
library(parallel)
library(fields)
library(maps)
library(ecodist)
library(hierfstat)
library(StAMPP)
library(ade4)
library(poppr)
library(netview)
library(visreg)
library(devtools)
library(ggplot2)
library(tidyverse)
library(doBy)
library(plotly)
library(reshape)
library(ggspatial)
library(ozmaps)
library(sf)
library(ggrepel)
library(patchwork)

##1 - Filtering SNPs with DartR####

#1.1. Read in dartseq SNPs and metafile data
gl <- gl.read.dart(filename = "data/Report_DAmb20-5632_1_moreOrders_SNP_2_rename_ordered.csv", ind.metafile = "data/Ai_ind_metrics.csv")
gl

#1.2. Report loci, individuals, and populations of the dartseq dataset
nLoc(gl)
nInd(gl)
nPop(gl)
levels(pop(gl))

drop_sibs <- as.matrix(read.csv('data/drop_sibs.csv'))

glunique <- gl.drop.ind(gl, ind.list = drop_sibs, recalc = T)

#1.2. Filtering SNPs that may skew later analyses
glx <- gl %>% 
  gl.filter.reproducibility(., threshold = .95) %>% 
  gl.filter.monomorphs(.) %>% 
  gl.filter.secondaries(.) %>% 
  gl.filter.maf(., threshold = 0.02, verbose = 2) %>% 
  gl.filter.hwe(., alpha = .05, basis = 'any', bon = T, verbose = 2) %>% 
  gl.filter.hamming(., threshold = 0.2, rs = 5, pb = TRUE, verbose = 2) %>% 
  gl.filter.callrate(., threshold = .7, method = 'loc') %>% 
  gl.filter.callrate(., method = 'ind', threshold = .4, recalc = T) 

gl.report.sexlinked(glx)

ind.callrate <- NA.posi(glx)

NA.posi(glx)

nInd(glx)
nLoc(glx)
nPop(glx)
levels(pop(glx))
glx$ind.names  

#1.x Import hierarchy into genlight file, and convert into genind and genclone for downstream analysis
Aihierarchy <- read.csv('data/derived/hierarchy_inds.csv')
glx$strata <- Aihierarchy
head(strata(glx, ~ region/pop, combine = F))

glxi <- gl2gi(glx)
glxi$strata <- Aihierarchy
names(other(glxi))
head(strata(glxi, ~region/pop, combine = F))

glxclone <- as.genclone(glxi)

#1.x Export the genlight, genind, and genclone files


#1.3. Various visualisations of data without model constraints
pc <- gl.pcoa(glx, nfactors = 5)
p1 <- gl.pcoa.plot(pc, glx, labels = 'pop', xaxis = 1, yaxis = 2)
gl.pcoa.plot(pc, glx, labels = 'interactive', xaxis = 1, yaxis=2)

gl.tree.nj(glx, type = 'unrooted')

#1.4 Let's make a map
oz_states <- ozmaps::ozmap_states
oz_states

points <- tibble::tribble(
  ~ site, ~Latitude, ~Longitude,
  'Broadwater SF',-37.0116, 149.9026333,
  'Dampier SF', -36.092783, 149.8720167,
  'Termeil SF', -35.439667, 150.3634167, 
  'Mares Forest SF', -34.2882, 149.9493167,
  'Mt Wilson', -33.495483, 150.4030667 ,
  'Cumberland SF', -33.743641, 151.040079, 
  'Olney SF', -33.0964, 151.3515,
  'South Barrington Tops NP', -32.156766, 151.324950,
  'North Barrington Tops NP', -31.911183, 151.6354167,
  'Willi Willi NP', -31.211832, 152.415809,
  'Kerewong SF', -31.60764, 152.560996,
  'Dorrigo NP', -30.363919, 152.79882
)
points <- as.data.frame(points)
sites <- st_as_sf(points, coords = c('Latitude', 'Longitude'), remove = F, crs = 4283, agr = 'constant')

ggplot() +
  geom_sf(data = oz_states, colour = 'grey', fill = 'white') +
  geom_point(data = sites, mapping = aes(x = Longitude, y = Latitude), color = 'red') +
  coord_sf(xlim = c(135, 160), ylim = c(-43, -27)) + 
  #geom_text_repel(data = sites, aes(x = Longitude, y = Latitude, label = site), 
                   #nudge_x = c(2, 1, 1.75, 1.75, -.5, 2.5, 3, 2.7, -1, -.4, 1.5, 1.5),
                   #nudge_y = c(0, 0, 0, 0, .25, -.2, .5, .1, 0, .5, 0, 0)) +
  annotation_scale(location = 'br', width_hint = 0.3, height = unit(.15, 'cm')) +
  annotation_north_arrow(location = 'br', which_north = 'true', pad_x = unit(.75, 'cm'), pad_y = unit(.5,'cm'),
                         style = north_arrow_fancy_orienteering) +
  theme_light() +
  theme(panel.grid = element_line(colour = 'white'))

#4 - Discriminant PCoAs#### ###Try making comparative PCAs in adegent
#4.1. Run the dapc (careful this makes the mac very hot)
set.seed(999)
Aix <- xvalDapc(tab(glxi, NA.method = 'mean'), pop(glxi))
system.time(Aix <- xvalDapc(tab(glxi, NA.method='mean'), pop(glxi),
                            n.pca = 10:20, n.rep = 1000,
                            parallel = 'multicore', ncpus = 4L))
Aix[-1]

#4.2 Plot the DAPC, first reorder factors in nth -> sth direction 
ord <- c('Dorrigo', 'Kerewong', 'NthBarrington', 'SthBarrington', 'Olney', 'Cumberland', 'MtWilson', 'MaresForest', 'Termeil', 'Dampier', 'Broadwater')
Aix$DAPC$assign <- factor(Aix$DAPC$assign, ord)
Aix$DAPC$grp <- factor(Aix$DAPC$grp, ord)
p1 <- scatter(Aix$DAPC, cex = 1, legend = TRUE, clabel = FALSE, posi.leg = 'bottomright',
              scree.pca = F, cleg = .75, xax = 1, yax=2, inset.solid = 1, posi.da = 'topright')

#3 - Testing for Isolation by Distance####
#3.1. Make geographic and genetic matricies
popDist <- as.data.frame(glx$other$latlong, stringsAsFactors = F)
popDist <- rdist.earth(popDist, popDist, miles = F, R = NULL)
popDist[1:5, 1:5]

genDist <- provesti.dist(glxi)
genDist
genDist[1:5, 1:5]
Gdis.vec <- lower(genDist)
Edis.vec <- lower(popDist)
plot(Gdis.vec~Edis.vec,
     xlab = 'geographic distance (km)',
     ylab = 'genetic distance')

#3.2. Testing for Spatial Autocorrelation using Mantel's R
SA <- lm(Gdis.vec ~ Edis.vec)
summary(SA)
par(mfrow =c(1,1))
visreg(SA, xvar = 'Edis.vec', overlay=T,
       xlab = "geographic distance (km)",
       ylab = "genetic distance (Nei's)")

Mantelr <- ecodist::mantel(Gdis.vec~Edis.vec, nperm = 1e6, nboot = 5e5, 
                           pboot = 0.9, cboot = 0.95)
Mantelr

dens <- kde2d(Edis.vec, Gdis.vec, n=300)
myPal <- colorRampPalette(c('white', 'blue', 'gold', 'orange', 'red'))
plot(Edis.vec, Gdis.vec, pch=20, cex=.5)
image(dens, col=transp(myPal(300), .5), add=T)
abline(SA)

ggplot(SA, aes(x = Edis.vec, y= Gdis.vec)) +
  geom_smooth(method = 'lm') +
  geom_jitter(alpha = 0.5, width = 0.2) +
  theme_classic() +
  ylab("Genetic distance (Nei's)") +
  xlab('Geographic distance (km)')

#5 - Inferring Continuous and Discrete Population Genetic Structure Across Space with conStruc####
#5.1 set up pop and name lines, STRUCTURE takes numeric pop-IDs
drop_sibs <- as.matrix(read.csv('data/drop_sibs.csv'))

glxunique <- gl.drop.ind(glxunique, ind.list = drop_sibs, recalc = T)

glxunique <- gl.filter.callrate(glxunique, method = 'ind', threshold = .8, recalc = T)
glxunique <- gl.filter.callrate(glxunique, threshold = .75, method = 'loc')

save(glxunique, file = 'data/derived/glxunique.RData')

load('data/derived/glxunique.RData')

levels(pop(glxunique))
glxunique$ind.names

pop_lvls <- levels(glxunique@pop)
Pp <- glxunique@pop
POP = as.numeric(factor(Pp, levels = pop_lvls))
ind.names <- glxunique$ind.names

xy <- as.matrix((glxunique$other$latlong))
popDist <- rdist.earth(xy, xy, miles = F, R = NULL)

#5.2 Convert Converting genlight file into a STRCUTURE file to be used by conStruct
gl2structure(glxunique, indNames = ind.names, addcolumns = c(POP), ploidy = 2, exportMarkerNames = F, 
             outfile = 'glxunique.str', outpath = 'data/derived')

#5.3 conStruct the conStruct file from the structure file
Ai.construct <- structure2conStruct(infile = 'data/derived/glxunique.str', 
                                    onerowperind = F,
                                    start.loci = 3,
                                    missing.datum = -9,
                                    outfile = 'data/derived/Ai_construct')
Ai.construct[1:5, 2:10]

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
options(mc.cores = 5)
options()

#5.4 Run Spatial analyses for K = 1:5
Aispk1 <- conStruct(spatial = T,
                    K = 1,
                    freqs = Ai.construct,
                    geoDist = popDist, 
                    coords = xy,
                    n.iter = 100000,
                    n.chains = 10,
                    make.figs = T,
                    save.files = T,
                    prefix = 'spk1',
                    control = setNames(list(0.95),"adapt_delta"))

Aispk2 <- conStruct(spatial = T,
                    K = 2,
                    freqs = Ai.construct,
                    geoDist = popDist, 
                    coords = xy,
                    n.iter = 100000,
                    n.chains = 10,
                    make.figs = T,
                    save.files = T,
                    prefix = 'spk2',
                    control = setNames(list(0.95),"adapt_delta"))

Aispk3 <- conStruct(spatial = T,
                    K = 3,
                    freqs = Ai.construct,
                    geoDist = popDist, 
                    coords = xy,
                    n.iter = 100000,
                    n.chains = 10,
                    make.figs = T,
                    save.files = T,
                    prefix = 'spk3')

Aispk4 <- conStruct(spatial = T,
                    K = 4,
                    freqs = Ai.construct,
                    geoDist = popDist, 
                    coords = xy,
                    n.iter = 100000,
                    n.chains = 10,
                    make.figs = T,
                    save.files = T,
                    prefix = 'spk4')

Aispk5 <- conStruct(spatial = T,
                    K = 5,
                    freqs = Ai.construct,
                    geoDist = popDist, 
                    coords = xy,
                    n.iter = 100000,
                    n.chains = 10,
                    make.figs = T,
                    save.files = T,
                    prefix = 'spk5')

#5.5 Run non-spatial analyses for K:1-5
nspk1 <- conStruct(spatial = FALSE,
                   K = 1,
                   freqs = Ai.construct,
                   geoDist = NULL,
                   coords = xy,
                   n.iter = 100000,
                   n.chains = 10,
                   prefix = 'nspk1')

nspk2 <- conStruct(spatial = FALSE,
                   K = 2,
                   freqs = Ai.construct,
                   geoDist = NULL,
                   coords = xy,
                   n.iter = 100000,
                   n.chains = 10,
                   prefix = 'nspk2', 
                   make.figs = T, 
                   save.files = T)

nspk3 <- conStruct(spatial = FALSE,
                   K = 3,
                   freqs = Ai.construct,
                   geoDist = NULL,
                   coords = xy,
                   n.iter = 100000,
                   n.chains = 10,
                   prefix = 'nspk3', 
                   make.figs = T, 
                   save.files = T)

nspk4 <- conStruct(spatial = FALSE,
                   K = 4,
                   freqs = Ai.construct,
                   geoDist = NULL,
                   coords = xy,
                   n.iter = 100000,
                   n.chains = 10,
                   prefix = 'nspk4', 
                   make.figs = T, 
                   save.files = T)

nspk5 <- conStruct(spatial = FALSE,
                   K = 5,
                   freqs = Ai.construct,
                   geoDist = NULL,
                   coords = xy,
                   n.iter = 100000,
                   n.chains = 10,
                   prefix = 'nspk5', 
                   make.figs = T, 
                   save.files = T)

##Model Contributions
#Spatial model
layer.contributionsS <- matrix(NA, nrow=5, ncol =5)
load("outputs/conStruct/spk1_conStruct.results.Robj")
load("outputs/conStruct/spk1_data.block.Robj")
layer.contributionsS[,1] <- c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(0,4))
tmpS <- conStruct.results[[1]]$MAP$admix.proportions

for(i in 2:5){
  load(sprintf("outputs/conStruct/spk%s_conStruct.results.Robj", i))
  load(sprintf("outputs/conStruct/spk%s_data.block.Robj", i))
  tmpS.order <- match.layers.x.runs(tmpS, conStruct.results[[1]]$MAP$admix.proportions)
  layer.contributionsS[, i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                              data.block = data.block,
                                                              layer.order = tmpS.order),
                                 rep(0,5-i))
  tmpS <- conStruct.results[[1]]$MAP$admix.proportions[,tmpS.order]
}

p1 <- barplot(layer.contributionsS,
        col=c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1"),
        xlab ="Spatial Model K-values",
        ylab = "layer contributions", 
        names.arg=paste0("K=",1:5))

#Non-spatial model
layer.contributionsS <- matrix(NA, nrow=5, ncol =5)
load("outputs/conStruct/nspk1_conStruct.results.Robj")
load("outputs/conStruct/nspk1_data.block.Robj")
layer.contributionsS[,1] <- c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(0,4))
tmpS <- conStruct.results[[1]]$MAP$admix.proportions

for(i in 2:5){
  load(sprintf("outputs/conStruct/nspk%s_conStruct.results.Robj", i))
  load(sprintf("outputs/conStruct/nspk%s_data.block.Robj", i))
  tmpS.order <- match.layers.x.runs(tmpS, conStruct.results[[1]]$MAP$admix.proportions)
  layer.contributionsS[, i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                              data.block = data.block,
                                                              layer.order = tmpS.order),
                                 rep(0,5-i))
  tmpS <- conStruct.results[[1]]$MAP$admix.proportions[,tmpS.order]
}

p2 <- barplot(layer.contributionsS,
        col=c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1"),
        xlab ="Non-spatial Model K-values",
        ylab = "layer contributions", 
        names.arg=paste0("K=",1:5))

p1+p2

##Making the admixture plots
###Spatial Models
load('outputs/conStruct/spk2_conStruct.results.Robj')
load('outputs/conStruct/spk2_data.block.Robj')
admix.props <- conStruct.results$chain_1$MAP$admix.proportions
match.layers.x.runs(admix.props1, admix.props)
make.structure.plot(admix.proportions = admix.props,
                    sample.order = order(data.block$coords[,1]),
                    sample.names = row.names(data.block$coords),
                    mar = c(4.5, 4, 2, 2))

load('outputs/conStruct/spk3_conStruct.results.Robj')
load('outputs/conStruct/spk3_data.block.Robj')
admix.props <- conStruct.results$chain_1$MAP$admix.proportions
match.layers.x.runs(admix.props, admix.props)
make.structure.plot(admix.proportions = admix.props,
                    sample.order = order(data.block$coords[,1]),
                    sample.names = row.names(data.block$coords),
                    mar = c(4.5, 4, 2, 2))

###Non-Spatial Models
load('outputs/conStruct/nspk4_conStruct.results.Robj')
load('outputs/conStruct/nspk4_data.block.Robj')
admix.props <- conStruct.results$chain_4$MAP$admix.proportions
match.layers.x.runs(admix.props, admix.props)
make.structure.plot(admix.proportions = admix.props,
                    sample.order = order(data.block$coords[,1]),
                    sample.names = row.names(data.block$coords),
                    mar = c(4.5, 4, 2, 2))

load('outputs/conStruct/nspk5_conStruct.results.Robj')
load('outputs/conStruct/nspk5_data.block.Robj')
admix.props <- conStruct.results$chain_1$MAP$admix.proportions
match.layers.x.runs(admix.props, admix.props)
make.structure.plot(admix.proportions = admix.props,
                    sample.order = order(data.block$coords[,1]),
                    sample.names = row.names(data.block$coords),
                    mar = c(4.5, 4, 2, 2))

#6 - Phylogenetic species delimitation models
#6.1 - Filtering 
glBFD <- gl.filter.callrate(glx, threshold = 1, method = 'loc')

#6.2 - Check to make sure you haven't removed EVERYTHING
nLoc(glBFD)
nInd(glBFD)
pc <- gl.pcoa(glBFD, nfactors = 5)
gl.pcoa.plot(pc, glBFD, labels = 'pop', xaxis = 1, yaxis = 2)

#6.3 - Select the randomly chosen indivs for BFD*
indlist <- glBFD$ind.names
write.csv(indlist, file='glBFD.csv')
drop_rand <- as.matrix(read.csv('glBFD.csv'))
glBFD <- gl.drop.ind(glBFD, ind.list = drop_rand, recalc = T)

nInd(glBFD)
nLoc(glBFD)
NA.posi(glBFD)

gl2snapp(glBFD, 'Ai.nex', outpath = 'data/derived', verbose = 2)

#2. Descriptive Statistics####
##COME BACK HERE AND SPLIT NTH BY STH TO LOOK INTO DIFFERENCES OF HEHO AND ALSO PAIRWISE FST
##FOR PAIRWISE FST HAVE STHN POPS BELOW THE DIAGONAL AND NTHN POPS ABOVE. NEED TO SPLIT BARRINGTON INTO NTH/STH POPS

#x.1 Splitting the dataset to Nth vs Sth
glNth <- gl.keep.pop(gl, pop.list = c('Dorrigo', 'Kerewong', 'NthBarrington'), recalc = T)
nInd(glNth)
nLoc(glNth)
levels(pop(glNth))
glNth$ind.names

glSth <- gl.drop.pop(gl, pop.list = c('Dorrigo', 'Kerewong', 'NthBarrington'), recalc = T)
nInd(glSth)
nLoc(glSth)
levels(pop(glSth))
glSth$ind.names

#x.2 Filtering out SNPs that are not of a good enough quality
glNth <- glNth %>% 
  gl.filter.reproducibility(., threshold = .95) %>% 
  gl.filter.monomorphs(.) %>% 
  gl.filter.secondaries(.) %>% 
  gl.filter.maf(., threshold = 0.02, verbose = 2) %>% 
  gl.filter.hwe(., alpha = .05, basis = 'any', bon = T, verbose = 2) %>% 
  gl.filter.hamming(., threshold = 0.2, rs = 5, pb = TRUE, verbose = 2) %>% 
  gl.filter.callrate(., threshold = .8, method = 'loc') %>% 
  gl.filter.callrate(., method = 'ind', threshold = .6, recalc = T) 

nInd(glNth)
nLoc(glNth)
levels(pop(glNth))
namesNth <- glNth$ind.names
write.csv(names, file = 'data/derived/namesNth.csv')
ind.callrate <- NA.posi(glNth)

Nthhierarchy <- read.csv('data/derived/heirarchy_Nth.csv')
glNth$strata <- Nthhierarchy
head(strata(glNth, ~ pop/gallery, combine = F))

glNthi <- gl2gi(glNth)
glNthi$strata <- Nthhierarchy
names(other(glNthi))
head(strata(glNthi, ~pop/gallery, combine = F))
glNthclone <- as.genclone(glNthi)

glSth <- glSth %>% 
  gl.filter.reproducibility(., threshold = .95) %>% 
  gl.filter.monomorphs(.) %>% 
  gl.filter.secondaries(.) %>% 
  gl.filter.maf(., threshold = 0.02, verbose = 2) %>% 
  gl.filter.hwe(., alpha = .05, basis = 'any', bon = T, verbose = 2) %>% 
  gl.filter.hamming(., threshold = 0.2, rs = 5, pb = TRUE, verbose = 2) %>% 
  gl.filter.callrate(., threshold = .8, method = 'loc') %>% 
  gl.filter.callrate(., method = 'ind', threshold = .6, recalc = T) 

nInd(glSth)
nLoc(glSth)
levels(pop(glSth))
namesSth <- glSth$ind.names
write.csv(namesSth, file = 'data/derived/namesSth.csv')
ind.callrate <- NA.posi(glSth)

Sthhierarchy <- read.csv('data/derived/hierarchy_Sth.csv')
glSth$strata <- Sthhierarchy
head(strata(glSth, ~ pop/gallery, combine = F))

glSthi <- gl2gi(glSth)
glSthi$strata <- Sthhierarchy
names(other(glSthi))
head(strata(glSthi, ~pop/gallery, combine = F))
glSthclone <- as.genclone(glSthi)

#2.1 Generate report of expected and observed heterozygosity
HeHoNth <- gl.report.heterozygosity(glNth, method = "pop", n.invariant = 0, boxplot = "adjusted", range = 1.5, cex.labels = 0.7,
                                 verbose = 2)
HeHoSth <- gl.report.heterozygosity(glSth, method = "pop", n.invariant = 0, boxplot = "adjusted", range = 1.5, cex.labels = 0.7,
                                    verbose = 2)
HeHo <- gl.report.heterozygosity(gl, method = "pop", n.invariant = 0, boxplot = "adjusted", range = 1.5, cex.labels = 0.7,
                                    verbose = 2)
par(mfrow = c(1,1))
#2.2 Calculation of Heterozygosity measures and Fis
DSNth <- glNthi %>% 
  basic.stats(., diploid = T)

DSSth <- glSthi %>% 
  basic.stats(., diploid = T)

#Custom Function to calculate SE of means of Ho, He, and Fis #####
SEmean <- function(x){
  x <- x[!is.na(x)]
  se <- sd(x)/sqrt(length(x))
  return(se)
}

#2.3 Generate a summary statistics table for loci, Heterozygosity, and inbreeding
##COME BACK HERE AND TEST FOR WAHLUND EFFECTS BY REGRESSING LOCI FIS BY FST. IF POSITIVE THEN WAHLUND
##Waples(2015)
NthPopStats <- select(HeHoNth, nInd:nLoc)
NthPopStats$Ho <- apply(DSNth$Ho, 2, mean, na.rm=T)
NthPopStats$HoSE <- (apply(DSNth$Ho, 2, SEmean))
NthPopStats$Hs <- (apply(DSNth$Hs, 2, mean, na.rm=T))
NthPopStats$HsSE <- apply(DSNth$Hs, 2, SEmean)
NthPopStats$Fis <- apply(DSNth$Fis, 2, mean, na.rm=T)
NthPopStats$FisSE <- apply(DSNth$Fis, 2, SEmean)

NthOverall <- as.data.frame(DSNth$overall) %>% 
  t()

SthPopStats <- select(HeHoSth, nInd:nLoc)
SthPopStats$Ho <- apply(DSSth$Ho, 2, mean, na.rm=T)
SthPopStats$HoSE <- (apply(DSSth$Ho, 2, SEmean))
SthPopStats$Hs <- (apply(DSSth$Hs, 2, mean, na.rm=T))
SthPopStats$HsSE <- apply(DSSth$Hs, 2, SEmean)
SthPopStats$Fis <- apply(DSSth$Fis, 2, mean, na.rm=T)
SthPopStats$FisSE <- apply(DSSth$Fis, 2, SEmean)

SthOverall <- as.data.frame(DSSth$overall) %>% 
  t()

write.csv(NthPopStats, file = 'outputs/Nth_summary_stats.csv')
write.csv(SthPopStats, file = 'outputs/Sth_summary_stats.csv')
write.csv(NthOverall, file = 'outputs/Nthglobal_stats.csv')
write.csv(SthOverall, file = 'outputs/Sthglobal_stats.csv')

#2.4 Pairwise Fst across populatons
Nthpwfst <- stamppFst(glNth, nboots = 100, percent = 95, nclusters = 6)
Sthpwfst <- stamppFst(glSth, nboots = 100, percent = 95, nclusters = 6)
Nthpwfst$Fsts
Sthpwfst$Fsts

Nthtemp <- Nthpwfst$Fsts
Nthtemp[upper.tri(Nthtemp)] <- t(Nthtemp)[upper.tri(Nthtemp)]
ord <- c('Dorrigo', 'Kerewong', 'NthBarrington')
NthDist <- Nthtemp[ord, ord]
NthDist

Sthtemp <- Sthpwfst$Fsts
Sthtemp[upper.tri(Sthtemp)] <- t(Sthtemp)[upper.tri(Sthtemp)]
ord <- c('SthBarrington', 'Olney', 'Cumberland', 'MtWilson', 'MaresForest', 'Termeil', 'Dampier', 'Broadwater')
SthDist <- Sthtemp[ord, ord]
SthDist

write.csv(NthDist, file = 'outputs/Nthdist.csv')
write.csv(SthDist, file = 'outputs/Sthdist.csv')

#2.5 AMOVA tests across populations
NthAMOVA <- poppr.amova(glNthclone, ~pop/gallery, cutoff = 0.6, nperm = 999)
NthAMOVA
Nthamova.test <- randtest(NthAMOVA)
Nthamova.test
plot(Nthamova.test)

SthAMOVA <- poppr.amova(glSthclone, ~pop/gallery, cutoff = 0.6, nperm = 999)
SthAMOVA
Sthamova.test <- randtest(SthAMOVA)
Sthamova.test
plot(Sthamova.test)

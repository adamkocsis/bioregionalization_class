# Paleobiogeographic analysis of PBDB data (last 10 Ma) 
# Adam T. Kocsis (Erlangen, 2020-06-17)
# CC-BY 4.0

library(icosa)
library(vegan)
library(igraph)
library(chronosphere)

# establish working directory
workdir <- "/mnt/sky/Dropbox/Teaching/FAU/Macroecology/material/src/bioregionalization"
setwd(workdir)

####################################################################################
# read in the data
# color scheme
load("data/allHex.RData")

# pbdb dataset prepared earlier
ceno6 <- readRDS("export/pbdb_species_49.rds")
# land polygons
land <- fetch("NaturalEarth")

source("R/methods/plots.R")

####################################################################################

# Geographic griding (icosa)

# create a grid
gr <- hexagrid(7, sp=TRUE)

	# visualize
	plot(gr)
	gridlabs(gr)

# locate the cells
ceno6$cells <- locate(gr, ceno6[, c("paleolng", "paleolat")])

# omit those entries where there is no species name or cells
ceno6use <- ceno6[!is.na(ceno6$cells) & !is.na(ceno6$trinomen),]

# contingency matrix
cont <- table(ceno6use$cells, ceno6use$trinomen)
samp <- cont[1:10,1:10]
samp


# incidence
cont[cont>1] <- 1

# Method 1. Compositional similarity
distmat <- vegdist(cont, method="jaccard")

# clustering
cluster <- hclust(distmat, "ward.D")

# plot this
par(mfrow=c(2,1))
plot(cluster)

# cutting the dendrogram-> membership vector
mem <- cutree(cluster, h=1.2)
abline(h=1.2, col="red")

biogeoplot(mem)

# try setting it to 
# clustering method set to "average"
# cutting height set to 0.99


####################################################################################
dev.off()
# Network- approach
# same contingency matrix 

# transform to a bipartite network
bipartite <- graph_from_biadjacency_matrix(cont)

# # you can cluster this directly
# infoBi <- cluster_infomap(bipartite)
# 
# # membership 
# info <- membership(infoBi)
# 
# # localities 
# binfoLoc <- info[names(info)%in%rownames(gr@faces)]
# 
# # with network anaysis
# biogeoplot(binfoLoc)

# project the graph (to look at only localities)
graph <- bipartite_projection(bipartite, which="false")
plot(graph)

# you can cluster this directly
infoGraph <- cluster_infomap(graph)

# membership 
info <- membership(infoGraph)

# with network anaysis
biogeoplot(info)

# outliers?

####################################################################################
# Tracing regions through time

# spatiotemporal cells
ceno6use$stc<- paste0("stg",ceno6use$stg,"_cell",  ceno6use$cells)

cont <- table(ceno6use$stc, ceno6use$trinomen)

# incidence
cont[cont>1] <- 1

# transform to a bipartite network
bipartite <- graph_from_biadjacency_matrix(cont)

# project the graph (to look at only localities)
graph <- bipartite_projection(bipartite, which="false")
plot(graph)

# you can cluster this directly
infoGraph <- cluster_infomap(graph)

# membership
info <- membership(infoGraph)

# display one time bine
bin <- 92

# get the subset corresponding to the focal bin
thisBin <- grep(paste0("stg", bin), names(info))
members <- info[thisBin]
names(members) <- gsub(paste0("stg", bin, "_cell"), "", names(members))
biogeoplot(members)


pdf("export/tracing.pdf", width=16, height=10)
for(bin in 92:94){
	thisBin <- grep(paste0("stg", bin), names(info))
	members <- info[thisBin]
	names(members) <- gsub(paste0("stg", bin, "_cell"), "", names(members))
	biogeoplot(members)
	mtext(text=paste0("Stage no. ", bin), side=3, line=1)
}
dev.off()






#######################################################
# Simplified code development:
# https://github.com/adamkocsis/obigeo

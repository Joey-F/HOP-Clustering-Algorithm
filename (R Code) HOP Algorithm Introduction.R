

#####
#####
##### Create some sample data #####


set.seed(3657)

dev <- .5

a <- cbind(rnorm(50, mean=1, sd=dev),rnorm(50, mean=1, sd=dev))
b <- cbind(rnorm(50, mean=4, sd=dev),rnorm(50, mean=1, sd=dev))
c <- cbind(rnorm(50, mean=1, sd=dev),rnorm(50, mean=4, sd=dev))
d <- cbind(rnorm(50, mean=4, sd=dev),rnorm(50, mean=4, sd=dev))

data <- rbind(a,b,c,d)



#####
#####
##### Run the HOP Algorithm #####

###
### (1) Calculate density
###

library(deldir)

voronoi <- deldir(x=data.frame(data))


# plot voronoi diagram for posterity

plot(data, type="n")
points(data, pch=20, col="red", cex=0.75)
plot(voronoi, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=1, col="black")


# density = inverse voronoi volume

density <- voronoi$summary$dir.area^(-1)


###
### (2) Calculate neighbors
###


delsgs <- voronoi$delsgs

# A data frame with 6 columns. The first 4 entries of each row are the coordinates 
# of the points joined by an edge of a Delaunay triangle, in the order (x1,y1,x2,y2). 
# The last two entries are the indices of the two points which are joined.


n <- nrow(data)

delaunay.TF <- matrix(rep(FALSE, n^2), nrow=n)

for (row in 1:nrow(delsgs)) {
    
    delaunay.TF[delsgs[row,5], delsgs[row,6]] <- TRUE   # cell [i,j]
    
    delaunay.TF[delsgs[row,6], delsgs[row,5]] <- TRUE   # cell [j,i] (symmetry)
    
}


# delaunay.TF[i,j] = TRUE if there exists a delaunay edge btwn vertices i,j



###
### (3) HOP to highest-density neighbor
###


HOP.onestep.forward <- function(point.density, edge.matrix) {

    # point.density is a scalar function applied to every vertex. HOP will move from low values to high values.

    # edge.matrix is NxN logical matrix with [i,j]=TRUE if there exists an edge between [i,j].

    HOPvec <- vector()

    for (row in 1:nrow(edge.matrix)) {

        tempvec <- edge.matrix[row,]        # iterate through edge matrix one row at a time

        index.TF <- which(tempvec==TRUE)

        greater.density <- index.TF[point.density[row]<point.density[index.TF]]     # constructs directed edges from low to high density

        newvec <- rep(FALSE, length(tempvec))

        newvec[greater.density]<- TRUE      # newvec only keeps edges going from low density -> high density

        edge.matrix[row,] <- newvec         # we only kept edge.matrix[i,j]=TRUE if density increases from i to j

        index.TF <- c(which(newvec==TRUE), row)

        HOPvec[row] <- index.TF[point.density[index.TF] == max(point.density[index.TF])]    # HOPvec = each point's 1-step HOP destination

    }

    output <- list(HOPvec, edge.matrix)

    names(output) <- c("HOPvec", "edge.matrix")

    return(output)

}


# HOP.onestep.forward <- function(point.density, edge.matrix) {
#     
#     # point.density is a scalar value applied to every vertex. HOP will move from low values to high values.
#     
#     # edge.matrix is NxN logical matrix with [i,j]=TRUE if there exists an edge between [i,j]. 
#     
#     HOPvec <- vector()
#     
#     onestep.neighbor <- function(tempvec) { # input = row of edge matrix; logical values indicating edge or no edge
#         
#         # tempvec <- edge.matrix[row,]
#         
#         index.TF <- which(tempvec==TRUE)  # index.TF = list of neighboring vertices
#         
#         greater.density <- index.TF[point.density[row]<point.density[index.TF]]     # compare each data point
#         
#         newvec <- rep(FALSE, length(tempvec))   
#         
#         newvec[greater.density]<- TRUE      # newvec only keeps edges going from low density -> high density
#         
#         edge.matrix[row,] <- newvec         # we only kept edge.matrix[i,j]=TRUE if density increases from i to j
#         
#         index.TF <- c(which(newvec==TRUE), row)
#         
#         HOPvec[row] <- index.TF[point.density[index.TF] == max(point.density[index.TF])]    # HOPvec = each point's 1-step HOP destination
#         
#     }
#     
#     output <- list(HOPvec, edge.matrix)     # vector of 1-step HOP destinations;  non-symmetric (directed) matrix of HOP edges
#     
#     names(output) <- c("HOPvec", "edge.matrix")
#     
#     return(output)
#     
# }


HOP.output <- HOP.onestep.forward(density, delaunay.TF)

HOPvec <- HOP.output$HOPvec

delaunay.TF <- HOP.output$edge.matrix



###
### (4) Chain HOP paths until all paths terminate at maxima
###


HOP.complete <- function(HOP.neighbor) {

    change <- 1
    
    n.iterations <- 0
    
    HOPvec.temp <- HOPvec

    while(change==1) {
        
        n.iterations <- n.iterations + 1
        
        HOPvec.new <- HOP.neighbor
        
        for (row in 1:length(HOP.neighbor)) {
            
            x <- HOP.neighbor[row]
            
            HOPvec.new[row] <- HOP.neighbor[x]
            
        }
        
        if (identical(HOP.neighbor, HOPvec.new)) {change <- 0}
        
        else {HOP.neighbor <- HOPvec.new}
        
    }
    
    output <- list(HOPvec.new, n.iterations)
    
    names(output) <- c("HOP.maxima", "n.iterations")
    
    return(output)

}


HOP.terminate <- HOP.complete(HOPvec)$HOP.maxima



plot(data, type="n")
points(data, pch=20, col=HOP.terminate, cex=.75)
plot(voronoi, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=1, col="black")

points(data[unique(HOP.terminate),], pch=1, col="red", cex=1.1)

for (max in unique(HOP.terminate)) {
    
    index <- which(HOP.terminate==max)
    
    for (item in index) {
        
        lines(x=data[c(item,max),], col="grey")
        
    }
    
}

        

#####
#####
##### Create graph for the HOP maxima #####


library(igraph)



distmat <- as.matrix(dist(data))


distmat[delaunay.TF==FALSE] <- 0    # only include edges that are (one-way) Delaunay-connected


graph <- graph_from_adjacency_matrix(distmat, mode="undirected", weighted=TRUE)  # connected graph of Delaunay edges


###
### Which clusters are connected via delaunay edges?
###

clusters <- make_clusters(graph, membership = HOP.terminate)    # assign each point to its HOP maximum

bond.edges <- crossing(clusters, graph)         # identify HOP clusters which share a Delaunay edge

bond.edges <- get.edgelist(graph)[bond.edges,]  # bond.edges is matrix with 2 columns listing vertex-vertex pairs on cluster boundaries

bond.edges <- apply(bond.edges, 2, as.numeric) 

bond.edges <- matrix(HOP.terminate[bond.edges],         # replace vertex indices with cluster labels, to show which clusters are joined
                     ncol=ncol(bond.edges), byrow=F)

bond.edges <- unique(bond.edges)                        # eliminate duplicate rows

bond.edges <- apply(bond.edges, 2, as.character) 




###
### What is the edge sequence connecting each pair of maxima? (even non-adjacent maxima, for now)
###

HOP.maxima <- sort(unique(HOP.terminate))       # ordered list of maxima indices


HOP.max.neighbors <- function(object) {         # generates all geodesics from a single maxima to all other maxima
    
    temp <- shortest_paths(graph=graph, from=object, to=HOP.maxima)
    
    return(temp$vpath)
    
}


paths <- lapply(HOP.maxima, HOP.max.neighbors)  # apply HOP.max.neighbors function over all starting maxima

names(paths) <- HOP.maxima  

# "paths" contains all paths from each maximum to all other maxima (list object)



###
### Geodesic lengths: the LONGEST edge on the shortest path
###


extract.max.edge <- function(HOP.path) {    # given a HOP path. identify the geodesic length
    
    all.edges <- E(graph, path=HOP.path)$weight
    
    max.edge <- max(all.edges)
    
    return(max.edge)
    
}


extract.geodesics <- function(list.of.HOP.paths) {    # apply the previous function over a list of HOP paths from a single HOP max
    
    sapply(list.of.HOP.paths, extract.max.edge)
    
}


geodesic.distance <- sapply(paths, extract.geodesics)       # apply the previous function over all possible starting HOP maxima

rownames(geodesic.distance) <- HOP.maxima
colnames(geodesic.distance) <- HOP.maxima



###
### Weighted connectivity matrix between HOP maxima 
###

temp <- matrix(rep(0, nrow(geodesic.distance)^2), nrow=nrow(geodesic.distance))
rownames(temp) <- rownames(geodesic.distance)
colnames(temp) <- colnames(geodesic.distance) 



temp[bond.edges] <- 1   # bond.edges specifies the [row,column] index of an existing edge (referenced by row/col name) 


temp <- pmax(temp, t(temp))     # force symmetry, reflecting the 1's


geodesic.neighbors <- geodesic.distance

geodesic.neighbors[temp==0] <- 0        # geodesic distance between neighboring HOP maxima



#####
#####
##### Chain together the HOP maxima #####


# Geodesics can be measured by their raw magnitude, or by their relative rank against all (sorted) geodesics.

# We generally choose rank over magnitude, but both options are presented here.


###
### Connect based on geodesic DISTANCE:
###

# geodesic.neighbors <- cbind(geodesic.neighbors, rep(max(geodesic.neighbors), nrow(geodesic.neighbors)))
# 
# n <- ncol(geodesic.neighbors)
# 
# colnames(geodesic.neighbors)[n] <- "data_max"
# 
# 
# connected.clusters <- matrix(rep(0, length(geodesic.distance)), nrow=nrow(geodesic.distance), dimnames=dimnames(geodesic.distance))
# 
# 
# 
# for (row in 1:nrow(geodesic.neighbors)) {
#     
#     temp <- sort(geodesic.neighbors[row,])
#     
#     temp <- diff(temp)
#     
#     connects <- names(temp)[1:(which(temp==max(temp))-1)]   # column names to connect with current max
#     
#     connected.clusters[row,connects] <- 1
#     
# }


###
### Connect based on geodesic RANK:
###

geodesic.ranks <- geodesic.neighbors

geodesic.ranks[geodesic.ranks==0 | geodesic.ranks==-Inf] <- NA  # insert NA for all non-neighboring maxima and self-loops

geodesic.ranks <- matrix(rank(geodesic.ranks, na.last = "keep", ties.method="min"),     # calculates geodesic ranks
               ncol=ncol(geodesic.ranks), byrow=TRUE)


# try things besides ties.method = min?



connected.clusters <- matrix(rep(0, length(geodesic.distance)), nrow=nrow(geodesic.distance), dimnames=dimnames(geodesic.distance))

# geodesic.ranks[is.na(geodesic.ranks)] <- 0        if kept, will measure 0 to n rather than 1 to n
                                                    # problem = distinguish between low-connected vs unconnected

colnames(geodesic.ranks) <- colnames(geodesic.distance)


# for (row in 1:nrow(geodesic.ranks)) {
#     
#     temp <- sort(geodesic.ranks[row,], na.last=FALSE)
#     
#     temp <- diff(temp)
#     
#     connects <- names(temp)[1:(which(temp==max(temp, na.rm=T))-1)]   # column names to connect with current max
#     
#     connected.clusters[row,connects] <- 1
#     
# }



geodesic.step.function <- function(row.index, rank.matrix) {    # applies to a single row of rank.matrix
    
    temp <- sort(rank.matrix[row.index,], na.last=TRUE)     # sorts columns according to rank; excludes NA
    
    temp <- diff(temp)
    
    connects <- names(temp)[1:(which(temp==max(temp, na.rm=T))-1)]   # column names to connect with current max
    
    newvec <- rep(0, ncol(rank.matrix))

    names(newvec) <- colnames(rank.matrix)

    newvec[connects] <- 1

    return(newvec)
    
}



row.counter <- 1:nrow(geodesic.ranks)

names(row.counter) <- colnames(geodesic.ranks)


connected.clusters <- sapply(X=row.counter, FUN=geodesic.step.function, geodesic.ranks) 

# iterate through row.counter, applying geodesic.step.function to each row of geodesic.ranks


connected.clusters <- pmax(connected.clusters, t(connected.clusters))   # only consider mutual connectivity

diag(connected.clusters) <- rep(1, ncol(connected.clusters))


# Now: need to (1) create graph and (2) connect all maxima which can be joined 

library(expm)

fully.connected <- (connected.clusters %^% 10)>0

fully.connected[fully.connected] <- 1



reassign.connected.maxima <- function(matrix.row) {
    
    newmax <- min(which(matrix.row==1))     # the column index corresponding to the first 1
    
    newmax <- names(matrix.row)[newmax]     # the actual maxima (column name)
    
    return(newmax)
}



max.reassignment <- apply(fully.connected, MARGIN=1, FUN=reassign.connected.maxima)


final.cluster <- rep(0, length(max.reassignment))
    
for (entry in 1:length(HOP.terminate)) {
    
    temp <- HOP.terminate[entry]
    
    temp <- as.character(temp)
    
    final.cluster[entry] <- max.reassignment[temp]
    
}

final.cluster <- as.numeric(final.cluster)


plot(data, type="n")
points(data, pch=20, col=final.cluster, cex=.75)

for (max in unique(final.cluster)) {
    
    index <- which(final.cluster==max)
    
    for (item in index) {
        
        lines(x=data[c(item,max),], col="grey")
        
    }
    
}



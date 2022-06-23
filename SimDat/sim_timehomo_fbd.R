# Sim fbd trees with FossilSim in R
# Orlando Schwery 10.6.2022

setwd("C:/Users/oschwery/Documents/SELU-Postdoc/FBD_PPS/SimDat/SimPseudoPPTrees/")
getwd()


library(FossilSim)
library(ape)
library(phangorn)

ntax = 117
lambda = 0.35
mu = 0.3
psi = 0.001
numsim = 100

fbd_trees <- sim.fbd.taxa(ntax, numsim, lambda, mu, psi, complete=FALSE)
#str(fbd_trees)
numtips <- c()
maxages <- c()
for (i in 1:length(fbd_trees)) {
    numtips <- c(numtips, Ntip(fbd_trees[[i]]))
    maxages <- c(maxages, max(node.depth.edgelength(fbd_trees[[i]])))
}
fosstips <- numtips-ntax
cbind(numtips, maxages, fosstips)

# par(mfrow=c(3,4))

# for (i in 1:length(fbd_trees)) {
#     plot(fbd_trees[[i]])
# }


# plot(fbd_trees[[1]], cex=0.5)
# nodelabels(cex=0.5)
# tiplabels(cex=0.5)
# axisPhylo(backward=FALSE)



treeinfo <- list()
rootages <- list()

for (i in 1:length(fbd_trees)) {
    fbd_trees[[i]]$node.label <- paste0("node", seq(1:(length(fbd_trees[[i]]$tip.label)-1)))

    select.tip.or.node <- function(element, tree) {
        ifelse(element < Ntip(tree)+1,
            tree$tip.label[element],
            tree$node.label[element-Ntip(tree)])
    }

    ## Making the edge table
    edge_table <- data.frame(
                    "parent" = fbd_trees[[i]]$edge[,1],
                    "par.name" = sapply(fbd_trees[[i]]$edge[,1],
                                        select.tip.or.node,
                                        tree = fbd_trees[[i]]),
                    "child" = fbd_trees[[i]]$edge[,2],
                    "chi.name" = sapply(fbd_trees[[i]]$edge[,2],
                                        select.tip.or.node,
                                        tree = fbd_trees[[i]])
                    )


    treedat <- cbind(edge_table, fbd_trees[[i]]$edge.length)

    treedat_ordered <- treedat[order(treedat[,3],decreasing=FALSE),]

    depths <- node.depth.edgelength(fbd_trees[[i]])
    ages <- cbind(treedat_ordered, depths[-(Ntip(fbd_trees[[i]])+1)])
    rootage <- max(ages[,6])
    ages <- cbind(ages, (round(ages[,6]-rootage, 7))*-1)

    n_sampanc <- table(ages[, 5])[1]

    status <- c(rep("extant", times=ntax), rep("fossil", times=length(fbd_trees[[i]]$tip.label)-n_sampanc-ntax), rep("sampanc", times=n_sampanc), rep("node", times=length(fbd_trees[[i]]$node.label)-1))

    ages <- cbind(ages, status)

    colnames(ages)[5:7] <- c("edge.length", "depth", "age")

    ages <- rbind(ages, c(NA, NA, Ntip(fbd_trees[[i]])+1, "node1", NA, 0, rootage, "root"))
    for (j in c(1,3,5,6,7)) {
        class(ages[, j]) <- "numeric"
    }

    treeinfo[[i]] <- ages

    rootages[[i]] <- rootage

}

fossils <- list()
for (i in 1:length(treeinfo)) {
    tempinfo <- treeinfo[[i]][treeinfo[[i]]$status == 'fossil',]
    tempinfo <- rbind(tempinfo, treeinfo[[i]][treeinfo[[i]]$status == 'sampanc',])
    tempinfo <- tempinfo[, c("chi.name", "age")]
    tempinfo$min_ma <- tempinfo$age-(rootages[[i]]/10)
    tempinfo$max_ma <- tempinfo$age+(rootages[[i]]/10)
    for (j in 1:length(tempinfo$min_ma)) {
        if (tempinfo$min_ma[j] < 0) {
            delta <- tempinfo$age[j]/2
            tempinfo$min_ma[j] <- tempinfo$age[j]-delta
            tempinfo$max_ma[j] <- tempinfo$age[j]+delta
            delta <- c()
        }
    }
    tempinfo <- tempinfo[, -2]
    colnames(tempinfo)[1] <- "specimen"
    fossils[[i]] <- tempinfo
    tempinfo <- c()
}

taxa <- fossils
for (i in 1:length(treeinfo)) {
    tempinfo <- treeinfo[[i]][treeinfo[[i]]$status == 'extant', ]
    tempinfo <- cbind(tempinfo[,"chi.name"], rep(0, times=nrow(tempinfo)), rep(0, times=nrow(tempinfo)))
    colnames(taxa[[i]]) <- c("taxon", "min", "max")
    colnames(tempinfo) <- c("taxon", "min", "max")
    taxa[[i]] <- rbind(taxa[[i]], tempinfo)
    tempinfo <- c()
}

ingroups <- list()
extanttrees <- list()
for (i in 1:length(treeinfo)) {
    extants <- treeinfo[[i]][treeinfo[[i]]$status == 'extant', ]
    extants <- extants[,"chi.name"]
    extanttrees[[i]] <- keep.tip(phy=fbd_trees[[i]], tip=extants)
    rootnode <- getMRCA(extanttrees[[i]], extants)
    firstnodes <- Children(extanttrees[[i]], rootnode)
    node_1 <- Descendants(extanttrees[[i]], firstnodes[1], type="tips")[[1]]
    node_2 <- Descendants(extanttrees[[i]], firstnodes[2], type="tips")[[1]]
    if (length(node_1) > length(node_2)) {
        ingroups[[i]] <- c(extanttrees[[i]]$tip.label[node_1])
    } else if (length(node_1) < length(node_2)) {
        ingroups[[i]] <- c(extanttrees[[i]]$tip.label[node_2])
    } else {
        ingroups[[i]] <- c(extanttrees[[i]]$tip.label[node_1])
    }

}


####### testplot
cols <- treeinfo[[1]]$status
for (i in 1:length(cols)) {
    if(cols[i] == "extant"){cols[i] <- "green"}
    if(cols[i] == "fossil"){cols[i] <- "red"}
    if(cols[i] == "sampanc"){cols[i] <- "blue"}
    if(cols[i] == "fossil"){cols[i] <- "red"}
    if(cols[i] == "sampanc"){cols[i] <- "blue"}
    if(cols[i] == "node") {cols[i] <- "gold"}
    if(cols[i] == "root") {cols[i] <- "orange"}
}
 plot(treeinfo[[1]]$age, col=cols)

for (i in 1:length(treeinfo)) {
    print(table(treeinfo[[i]]$status))
    print(nrow(fossils[[i]]))
}

#######
class(fbd_trees) <- "multiPhylo"
write.tree(fbd_trees, "FossSim.trees")

class(extanttrees) <- "multiPhylo"
write.tree(extanttrees, "Extant.trees")

for (i in 1:length(fbd_trees)) {
    write.tree(fbd_trees[[i]], paste("tree_", i, ".tre", sep=""))
}

for (i in 1:length(fbd_trees)) {
    write.tree(extanttrees[[i]], paste("extanttree_", i, ".tre", sep=""))
}


for (i in 1:length(fbd_trees)) {
    write.csv(treeinfo[[i]], paste("treeinfo_", i, ".csv", sep=""))
}


for (i in 1:length(fbd_trees)) {
    write.table(fossils[[i]], file=paste("fossils", i, ".tsv", sep=""), quote=FALSE, sep='\t', row.names=FALSE)
}

for (i in 1:length(fbd_trees)) {
    write.table(taxa[[i]], file=paste("taxa", i, ".tsv", sep=""), quote=FALSE, sep='\t', row.names=FALSE)
}


for (i in 1:length(fbd_trees)) {
    write.table(ingroups[[i]], file=paste("ingroups", i, ".tsv", sep=""), quote=FALSE, sep='\t', row.names=FALSE)
}

rootages

for (i in 1:length(treeinfo)) {
    print(table(treeinfo[[i]]$status))
}

for (i in 1:length(ingroups)) {
    print(length(ingroups[[i]]))
}

for (i in 1:length(treeinfo)) {
    print(i)
    print(aggregate(treeinfo[[i]]$age, list(treeinfo[[i]]$status), FUN=min))
    print(aggregate(treeinfo[[i]]$age, list(treeinfo[[i]]$status), FUN=mean))
    print(aggregate(treeinfo[[i]]$age, list(treeinfo[[i]]$status), FUN=max))
}

# par(mfrow=c(1,2))
# plot(fbd_trees[[1]])
# rangeplot.asymmetric(fbd_trees[[1]])

# par(mfrow=c(2,2))
# plot(fbd_trees[[1]])

# f = sim.fossils.poisson(rate=3, tree=fbd_trees[[1]])
# f
# plot(f, tree=fbd_trees[[1]])

# SAt = SAtree.from.fossils(tree = fbd_trees[[1]], fossils = f)
# #plot(SAt$tree, cex=0.0001)
# plot(SAt$fossils, tree=SAt$tree)


# SAt_pruned = prune.fossils(tree = SAt$tree)
# #plot(SAt_pruned, cex=0.0001)
# plot(SAt_pruned$fossils, tree=SAt_pruned, cex=0.0001)

# SAt_sampled = sampled.tree.from.combined(tree = SAt$tree)
# #plot(SAt_sampled, cex=0.0001)
# plot(SAt_sampled$fossils, tree=SAt_sampled, cex=0.0001)

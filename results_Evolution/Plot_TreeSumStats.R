# Plot distributions of Tree Sumstats
# Orlando Schwery 23. June 2022

library(scales)

fbd_dir <- "C:/Users/oschwery/Documents/SELU-Postdoc/FBD_PPS/results_Evolution/results_folders_drive/results_"
sim_dir <- "C:/Users/oschwery/Documents/SELU-Postdoc/FBD_PPS/SimDat/SimPseudoPPTrees/"

fbd_runs <- c("1e", "1g", "1i", "1k", "2e", "2g")
sim_runs <- c("time_homo/", "skyline/")

fbd_file <- "summarize_trees_PPS.csv"
sim_file <- "summarize_trees_FossSim.csv"

fbd_dists <- list()
for (i in 1:length(fbd_runs)) {
    fbd_dists[[i]] <- read.csv(paste(fbd_dir, fbd_runs[i], "/", fbd_file, sep=""), header=TRUE)
}

sim_dists <- list()
for (i in 1:length(sim_runs)) {
    sim_dists[[i]] <- read.csv(paste(sim_dir, sim_runs[i], "/", sim_file, sep=""), header=TRUE)
}

modelcols <- c("lightgray", "darkgray", "lightblue", "darkblue", "red", "green", "plum", "purple")
modelbordercols <- c("darkgray", "black", "darkblue", "black", "darkred", "darkgreen", "purple", "purple")


alph <- 0.3

plotComparePDFsSummStats <- function(sim_stats, emp_stats, statNo, modelcols, modelbordercols, modelsets) {

  plot(density(sim_stats[, statNo], na.rm=TRUE),
    xlim=c(min(c(sim_stats[, statNo], emp_stats[, statNo]), na.rm=TRUE),
      max(c(sim_stats[, statNo], emp_stats[, statNo]), na.rm=TRUE)),
      col=alpha(modelbordercols[modelsets[1]], alph*3),
      xlab=colnames(emp_stats[statNo]), ylab="Density", main="")

  polygon(density(sim_stats[, statNo], na.rm=TRUE), col=alpha(modelcols[modelsets[1]], alph), border=alpha(modelbordercols[modelsets[1]], alph))

  lines(density(emp_stats[, statNo], na.rm=TRUE), col=alpha(modelbordercols[modelsets[2]], alph*3))

  polygon(density(emp_stats[, statNo], na.rm=TRUE), col=alpha(modelcols[modelsets[2]], alph), border=alpha(modelbordercols[modelsets[2]], alph))

  #abline(v=emp_stats[, statNo], col=empcols[generating_model[1]], lwd=2.5)

  #abline(v=emp_stats[, statNo], col=empcols[generating_model[2]], lwd=2.5)

}

#plotComparePDFsSummStats(sim_dists[[1]], fbd_dists[[1]], 1, modelcols, modelbordercols, c(1,6))

# par(mfrow=c(3,4))
# for (i in 1:length(sim_dists)) {
#     for (j in 1:length(fbd_dists)) {
#         for (k in 1:ncol(sim_dists[[1]])) {
#             plotComparePDFsSummStats(sim_dists[[i]], fbd_dists[[j]], k, modelcols, modelbordercols, c(i, 7+j))
#         }
#     }
# }

pdf(file=paste(sim_dir, "TreeSumStats.pdf", sep=""), width=25, height=25)
par(mfrow=c(5,5))
for (i in 1) {
    for (j in c(1,2,5,6)) {
        for (k in 1:ncol(sim_dists[[1]])) {
            plotComparePDFsSummStats(sim_dists[[i]], fbd_dists[[j]], k, modelcols, modelbordercols, c(i, 2+j))
        }
        plot.new()
        plot.new()
    }
}
for (i in 2) {
    for (j in c(3,4)) {
        for (k in 1:ncol(sim_dists[[1]])) {
            plotComparePDFsSummStats(sim_dists[[i]], fbd_dists[[j]], k, modelcols, modelbordercols, c(i, 2+j))
        }
        plot.new()
        plot.new()
    }
}

dev.off()

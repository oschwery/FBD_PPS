

library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)

getwd()
setwd("C:/Users/oschwery/Documents/SELU-Postdoc/FBD_PPS/")

file <- "TimeHomoTrials/antsSA_"

sets <- 4

trials <- letters[1:14]

traces <- list()

for (i in 1:length(trials)) {
    traces[[i]] <- readTrace(paste(file, sets, trials[i], ".log", sep=""))
}

tracesMCMC <- list()

for (i in 1:length(traces)) {
    tracesMCMC[[i]] <- as.mcmc(traces[[i]][[1]])
}

# tracesin1 <- list()

# for (i in 1:length(traces)) {
#     tracesin1[[i]] <- traces[[i]][[1]]
# }

# parameters of interest
params_all <- c("Posterior", "Likelihood", "alpha_morpho", "diversification", "extinction_rate", "origin_time", "speciation_rate", "turnover")

params_overall <- c("Posterior", "Likelihood")

params_morpho <- c("alpha_morpho")

params_divrates <- c("diversification", "extinction_rate", "speciation_rate", "turnover")

params_age <- c("origin_time")

cols_all = c("purple", "blue", "darkgreen", "green", "gold", "orange", "red", "black")

scenarios = c("4a", "4b", "4c", "4d", "4e", "4f", "4g", "4h", "4i", "4j", "4k", "4l", "4m", "4n")
#plotTrace(traces[[1]], vars=params_overall)

#plotTrace(tracesin1, vars=params_overall)

####################
# Plot traces

pdf(file="TimeHomo_TrialTraces.pdf", width=10, height=30)
for (j in 1:length(params_all)) {
    mins <- c()
    maxs <- c()
    for (i in 1:length(tracesMCMC)) {
        mins <- c(mins, min(tracesMCMC[[i]][, params_all[j]]))
        maxs <- c(maxs, max(tracesMCMC[[i]][, params_all[j]]))
    }
    damins <- min(mins)
    damaxs <- max(maxs)

    par(mfrow=c(7, 2))
    for (i in 1:length(tracesMCMC)) {
        traceplot(tracesMCMC[[i]][,params_all[j]], ylim=c(damins, damaxs), ylab=params_all[j], col=cols_all[j], main=scenarios[i])
    }
}
dev.off()



prep_plots <- function(traces, daparams) {
    mins <- c()
    maxs <- c()
    for (i in 1:length(traces)) {
        mins <- c(mins, min(traces[[i]][[1]][, daparams]))
        maxs <- c(maxs, max(traces[[i]][[1]][, daparams]))
    }
    damins <- min(mins)
    damaxs <- max(maxs)
    daplotz <- list()
    for (i in 1:length(traces)) {
        daplotz[[i]] <- plotTrace(traces[[i]], vars=daparams)[[1]] +
        coord_cartesian(xlim=c(damins, damaxs))
    }
    return(daplotz)
}


exec_plots <- function(daplotz, ncols) {
    grid.arrange(daplotz[[1]], 
                daplotz[[2]], 
                daplotz[[3]], 
                daplotz[[4]], 
                daplotz[[5]], 
                daplotz[[6]], 
                daplotz[[7]], 
                daplotz[[8]], 
                daplotz[[9]], 
                daplotz[[10]], 
                daplotz[[11]], 
                daplotz[[12]], 
                daplotz[[13]], 
                daplotz[[14]], 
                ncol=ncols)
} 



overallplotz <- prep_plots(traces, params_overall)
morphoplotz <- prep_plots(traces, params_morpho)
divrateplotz <- prep_plots(traces, params_divrates)
ageplotz <- prep_plots(traces, params_age)

pdf(file="TimeHomo_TrialPosteriors.pdf", width=30, height=25)
    exec_plots(overallplotz, 3)
    exec_plots(morphoplotz, 3)
    exec_plots(divrateplotz, 3)
    exec_plots(ageplotz, 3)
dev.off()




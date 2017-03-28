#! /usr/bin/Rscript --vanilla

## Lightly edited version for human consumption.

## Strategy for testing: given a tree with an inferred mutation map,
## take the number of mutations at a given codon and simulate that
## many mutations under a null (neutral) model of mutation.  This
## becomes a test of the arrangement of mutations *given* that there
## are a certain number of them.  We then take these simulations and
## compute our LRT for those, in order to get at the null
## distribution.

## Outputs two .RData files, which can be directly loaded into
## R for further processing.  The key things defined in these
## files are:
## - null.dists, which provide the simulated null distributions of our statistic
##   for all the different numbers of mutations appearing for each codon
##   position
## - LLRs, which holds our computed statistics
## - branch.info, which is the input CSV file, in R form

## Helper routines defined in this file:
source("LRT_helpers_RealData.r")

library(optparse)
library(parallel)

## FIXME what will we be inputting and outputting and how?
option.list <-
  list(make_option("--numsim", default=1000, type="integer",
                   help=paste("Number of null simulations to perform",
                     "for each case (default 1000)")),
       make_option("--parallel", default=1, type="integer",
                   help=paste("Number of cores to use; if 1 no parallelization is done")),
        make_option("--seed", default=2233, type="integer",
                   help=paste("Random seed")))
parser <-
  OptionParser(usage="usage: %prog [--numsim <number of simulations>] [--parallel <number of cores>] [--seed <random seed>] <branch info CSV file> <output RData file>",
               option_list=option.list,
               add_help_option=TRUE)

opts <- parse_args(parser, args=commandArgs(TRUE),
                   positional_arguments=TRUE)

input.file <- opts$args[1]
output.file <- opts$args[2]
num.sim <- opts$options$numsim
num.cores <- opts$options$parallel
random.seed <- opts$options$seed
set.seed(random.seed)

branch.info <- read.csv(input.file, stringsAsFactors=FALSE)

## For each codon, compute the LLR; then simulate a bunch of mutation
## maps under neutrality to characterize the null distribution.
num.codons <- ncol(branch.info) - 4

## To save as much labour as possible on the computing of the
## null empirical distributions, we precompute the empirical distributions.
## They can be shared between codons with the same number of mutations.

null.dists <- list()

## Use sapply if we aren't using parallelism; otherwise use mclapply.
my.apply <- sapply
if (num.cores > 1)
  {
    my.apply <- function (x, y)
      {
        return(unlist(mclapply(x, y, mc.preschedule=FALSE, mc.cores=num.cores)))
      }
  }

branch.lengths <- branch.info$endtime - branch.info$starttime
for (codon.pos in 1:num.codons)
  {
    num.muts <- sum(branch.info[, codon.pos + 4] != 0)

    if ((as.character(num.muts) %in% names(null.dists)) || num.muts == 0)
      {
        next
      }
    
    cat("Simulating null distribution on the tree with ",
        num.muts, " mutations....\n", sep="")

    null.dists[[as.character(num.muts)]] <-
      my.apply(1:num.sim,
               function (sim.idx)
               {
                 branches.with.mut <- sample(1:nrow(branch.info), num.muts,
                                             prob=branch.lengths)
                 mut.branch.vector <- rep(0, nrow(branch.info))
                 ## Mark the branches with a mutation with a 4 (since 1, 2, and 3
                 ## denote mutations as found by ancestral reconstruction).
                 mut.branch.vector[branches.with.mut] <- 4
                 
                 curr.null.ll <- max.null.log.likelihood(branch.info$starttime,
                                                         branch.info$endtime,
                                                         mut.branch.vector)
                 
                 curr.alt.opt.time <- max.alt.log.likelihood(branch.info$starttime,
                                                             branch.info$endtime,
                                                             branch.info$timetotip,
                                                             mut.branch.vector)
                 curr.alt.ll <- curr.alt.opt.time$value
                 ## cat("Simulation ", sim.idx, " finished.\n", sep="")
                 
                 return(curr.alt.ll - curr.null.ll)
               })
  }

## Set up a checkpoint here so if it crashes later at least we didn't
## lose the null distributions.
save.image(file=paste(output.file, ".nulldists.RData", sep=""))

LLRs <- NULL
alt.opt.times <- NULL
for (codon.idx in 1:num.codons)
  {
    cat("Computing log-likelihood ratio for codon ", codon.idx, "....\n", sep="")
    num.muts <- sum(branch.info[,codon.idx+4] != 0)
    
    ## If there were no mutations, then we can stop.
    if (num.muts == 0)
      {
        cat("No mutations for this codon; skipping.\n")
        next
      }
    
    null.ll <- max.null.log.likelihood(branch.info$starttime,
                                       branch.info$endtime,
                                       branch.info[,codon.idx+4])
    
    ## This returns the results of an optimization.
    alt.opt.time <- max.alt.log.likelihood(branch.info$starttime,
                                           branch.info$endtime,
                                           branch.info$timetotip,
                                           branch.info[,codon.idx+4])
    ## The result now is a list with entries 'par' and 'value'; 'par'
    ## holds the maximizing switch times; 'value' holds the optimum
    ## result.
    alt.ll <- alt.opt.time$value
    LLR <- alt.ll - null.ll
    LLRs <- rbind(LLRs,
                  data.frame(null.ll=null.ll,
                             alt.ll=alt.ll,
                             LLR=LLR))
    alt.opt.times <- c(alt.opt.times, alt.opt.time)
                             
  }

## Save the output in a manner we can analyze.
save.image(file=output.file)

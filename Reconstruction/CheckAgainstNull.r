#! /usr/bin/Rscript --vanilla

## Lightly edited for public consumption.

## Strategy for testing: given a tree with an inferred mutation map,
## take the number of mutations at a given codon and simulate that
## many mutations under a null (neutral) model of mutation.  This
## becomes a test of the arrangement of mutations *given* that there
## are a certain number of them.  We then take these simulations and
## compute our LRT for those, in order to get at the null
## distribution.

## The input branch information CSV file contains one line per node that
## describes the parent edge of that node:
## 
## node: name of the node in the tree
## starttime: time (measured from some point in the past) that the branch starts
## endtime: time that the branch ends, i.e. the time of the node
## codon1, codon2, ... codonN: 0 if no mutation on that branch, non-0 if there is
##
## The node names of the leaf nodes are formatted as
## "[id]__type[type]_mut[0 or 1]_time[true time of node]"
## and the true times will be extracted using a regular expression.  Modify this
## as you see fit to suit your data.

source("LRT_helpers.r")

library(optparse)
library(parallel)

option.list <-
  list(make_option("--numsim", default=1000,
                   help=paste("Number of null simulations to perform",
                     "for each case (default 1000)")),
       make_option("--parallel", default=1,
                   help=paste("Number of cores to use; if 1 no parallelization is done")))
parser <-
  OptionParser(usage="usage: %prog [--numsim <number of simulations>] [--parallel <number of cores>] <branch info CSV file> <output RData file>",
               option_list=option.list,
               add_help_option=TRUE)

opts <- parse_args(parser, args=commandArgs(TRUE),
                   positional_arguments=TRUE)

input.file <- opts$args[1]
output.file <- opts$args[2]
num.sim <- opts$options$numsim
num.cores <- opts$options$parallel

branch.info <- read.csv(input.file, stringsAsFactors=FALSE)
min.time <- min(branch.info$starttime)
max.time <- max(branch.info$endtime)

## Our data had node names that were either of the format 
## 546__type2_mut0_time3000
## for leaf nodes, or
## Node163
## for internal nodes.  This code extracts the "true" time of the
## former nodes; modify accordingly for your data.  You could also
## modify it so that you simply tell it a value for real.max.time,
## as that's the only place this information is used.
actual.dates <-
  sapply(branch.info$node,
         function (x)
         {
           node.name.pattern <- ".*__type.*_mut.*_time(.*)"
           result <- NA
           if (grepl(node.name.pattern, x))
             {
               result <- as.numeric(sub(node.name.pattern, "\\1", x))
             }
           return(result)
         })

real.max.time <- max(actual.dates, na.rm=TRUE)

## For each codon, compute the LLR; then simulate a bunch of mutation
## maps under neutrality to characterize the null distribution.
num.codons <- ncol(branch.info) - 4

## Get the null distribution of the LLR for the selected codon.
null.emp.dist <- rep(NA, num.sim)
LLR <- NA

selected.codon.idx <- 101
cat("Computing log-likelihood ratio for selected codon....\n")
num.muts <- sum(branch.info[,selected.codon.idx+4] != 0)

## If there were no mutations, then we can stop.
if (num.muts != 0)
  { 
    null.ll <- max.null.log.likelihood(branch.info$starttime,
                                       branch.info$endtime,
                                       branch.info[,selected.codon.idx+4])
    
    # This returns the results of an optimization.
    alt.opt.time <- max.alt.log.likelihood(branch.info$starttime,
                                           branch.info$endtime,
                                           branch.info$timetotip,
                                           branch.info[,selected.codon.idx+4])
    ## The result now is a list with entries 'par' and 'value'; 'par'
    ## holds the maximizing switch times; 'value' holds the optimum
    ## result.
    alt.ll <- alt.opt.time$value
    LLR <- alt.ll - null.ll
    
    ## We observed num.muts mutations.  Simulate that many mutations under
    ## neutrality.
    branch.lengths <- branch.info$endtime - branch.info$starttime

    cat("Simulating null distribution....\n")

    my.apply <- sapply
    if (num.cores > 1)
      {
        my.apply <- function (x, y)
        {
          return(unlist(mclapply(x, y, mc.preschedule=FALSE, mc.cores=num.cores)))
        }
      }
    null.emp.dist <-
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
                 cat("Simulation ", sim.idx, " finished.\n", sep="")
                 
                 return(curr.alt.ll - curr.null.ll)
               })
  }

## Save the output in a manner we can analyze.
save.image(file=output.file)

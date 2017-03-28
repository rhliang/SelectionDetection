## FIXME this is the relevant code from GetNodeDepths.r, which includes
## some other stuff that doesn't get used in our analysis.

## Helper to extract the "true" numerical date from a node label.
get.actual.time <- function (x)
  {
    node.name.pattern <- ".*__type.*_mut.*_time(.*)"
    result <- NA
    if (grepl(node.name.pattern, x))
      {
        ## Replace "pt" with ".".
        time.with.pt <- sub(node.name.pattern, "\\1", x)
        time.str <- gsub("pt", ".", time.with.pt)
        result <- as.numeric(time.str)
      }
    return(result)
  }


## Function that extracts all the necessary information from a Newick
## tree to locate mutations in the tree.
get.node.depths <- function(my.tree)
  {
    ## my.tree$edge is a (# of nodes)x2 matrix describing edge:
    ## the first column is the ancestor index, and the second column
    ## is the descendent.  The indices run as:
    ## 1 to (number of tips): tips, as indexed by my.tree$tip.label
    ## (number of tips) + 1 to (number of tips) + (number of internal nodes):
    ##  nodes, as indexed by my.tree$node.label (plus number of tips).
    
    node.depths <-
      data.frame(label=c(my.tree$tip.label, my.tree$node.label),
                 first.ttt=c(rep(0, length(my.tree$tip.label)),
                   rep(NA, length(my.tree$node.label))),
                 second.ttt=c(rep(0, length(my.tree$tip.label)),
                   rep(NA, length(my.tree$node.label))),
                 depth=node.depth.edgelength(my.tree),
                 length=rep(NA, length(my.tree$tip.label) + length(my.tree$node.label)))
    node.depths$actual.time <- sapply(node.depths$label, get.actual.time)

    node.depths$nearest.tip.idx <-
      c(1:length(my.tree$tip.label), rep(NA, length(my.tree$node.label)))
    node.depths$orig.depth <- node.depths$depth
    
    node.depths$parent.idx <-
      rep(NA, length(my.tree$tip.label) + length(my.tree$node.label))

    node.depths$first.child <-
      rep(NA, length(my.tree$tip.label) + length(my.tree$node.label))
    node.depths$second.child <-
      rep(NA, length(my.tree$tip.label) + length(my.tree$node.label))
    
    ## Adjust the tips.
    tip.indices <- 1:length(my.tree$tip.label)
    unique.tip.times <- unique(node.depths$actual.time[tip.indices])
    
    for (utt in unique.tip.times)
      {
        ## Get the tips with this tip time and take their max depth.
        curr.group.indices <-
          intersect(tip.indices, which(node.depths$actual.time == utt))
        curr.group.depth <-
          max(node.depths$depth[curr.group.indices])
        node.depths$depth[curr.group.indices] <- curr.group.depth
        
        ## We also need to adjust the lengths of the edges.
        node.depths$length[curr.group.indices] <-
          sapply(curr.group.indices,
                 function (idx)
                 {
                   ## Find the edge in the tree that has this tip
                   ## as the child.
                   parent.edge <- which(my.tree$edge[,2] == idx)

                   unadjusted.length <- my.tree$edge.length[parent.edge]
                   offset <- node.depths$depth[idx] - node.depths$orig.depth[idx]
                   return(unadjusted.length + offset)
                 })

        node.depths$parent.idx[curr.group.indices] <-
          sapply(curr.group.indices,
                 function (idx)
                 {
                   parent.edge <- which(my.tree$edge[,2] == idx)

                   return(my.tree$edge[parent.edge, 1])
                 })
      }
    
    ## Sweep through all of the internal nodes, filling in the node depth
    ## if both of the children have node depths set already.
    changes.made <- TRUE
    unset.internal.nodes <- 1:length(my.tree$node.label) + length(my.tree$tip.label)
    while (changes.made && length(unset.internal.nodes) > 0)
      {
        changes.made <- FALSE
        new.unset.internal.nodes <- NULL
        for (node.idx in unset.internal.nodes)
          {
            desc.edge.indices <- which(my.tree$edge[,1] == node.idx)
            ## Sanity check:
            if (length(desc.edge.indices) != 2)
              {
                cat("Error: internal node",
                    my.tree$node.label[node.idx-length(my.tree$tip.label)],
                    "does not have two descendants!\n")
              }
            
            desc.1.idx <- my.tree$edge[desc.edge.indices[1],2]
            desc.2.idx <- my.tree$edge[desc.edge.indices[2],2]
            
            first.desc.ttt <- min(node.depths$first.ttt[desc.1.idx],
                                  node.depths$second.ttt[desc.1.idx])
                                  
            second.desc.ttt <- min(node.depths$first.ttt[desc.2.idx],
                                   node.depths$second.ttt[desc.2.idx])
            
            if (is.na(first.desc.ttt) || is.na(second.desc.ttt))
              {
                new.unset.internal.nodes <- c(new.unset.internal.nodes, node.idx)
                next
              }
            
            ## Now we know both depths are set.  Add the appropriate edge
            ## length to each and take the minimum.  For these nodes,
            ## we know that their length has been set.
            first.desc.length <- node.depths$length[desc.1.idx]
            second.desc.length <- node.depths$length[desc.2.idx]

            node.depths$first.child[node.idx] <- desc.1.idx
            node.depths$second.child[node.idx] <- desc.2.idx
            node.depths$first.ttt[node.idx] <- first.desc.ttt + first.desc.length
            node.depths$second.ttt[node.idx] <- second.desc.ttt + second.desc.length

            if (node.depths$first.ttt[node.idx] <= node.depths$second.ttt[node.idx])
              {
                node.depths$nearest.tip.idx[node.idx] <-
                  node.depths$nearest.tip.idx[desc.1.idx]
              } else {
                node.depths$nearest.tip.idx[node.idx] <-
                  node.depths$nearest.tip.idx[desc.2.idx]
              }
            
            ## Set this node's length.
            parent.edge <- which(my.tree$edge[,2] == node.idx)
            ## If this is the root, we skip it.
            if (length(parent.edge) == 1)
              {
                node.depths$length[node.idx] <- my.tree$edge.length[parent.edge]
                node.depths$parent.idx[node.idx] <- my.tree$edge[parent.edge, 1]
              } else if (length(parent.edge) > 1) {
                print("Warning: there is a node with more than one parent edge")
              }
            
            changes.made <- TRUE
          }
        unset.internal.nodes <- new.unset.internal.nodes
        
        ## print(node.depths)
      }
    return(node.depths)
  }


## Helper function that returns, given the depth of a node and the
## length of its parent branch, and the node's time to nearest tip,
## the lengths of the branch that are in the early zone, dead zone,
## and hot zone.  Return a vector with 3 elements being those lengths.
segment.lengths <- function(depth, length, time.to.tip, T.s, T.h)
  {
    if (is.na(length))
      return(c(0,0,0))
    
    start.time <- depth - length
    end.time <- depth
    ## The start of the region in the hot zone (factoring in that it
    ## may overlap the early zone) is:
    hot.zone.start <- max(end.time - (T.h - time.to.tip), T.s)
    early <- max(0, min(end.time, T.s) - start.time)
    dead <- max(0, min(end.time, hot.zone.start) - max(T.s, start.time))
    hot <- max(0, end.time - max(start.time, hot.zone.start))

    return(c(early, dead, hot))
  }


## A helper that takes a data frame of mutations (as read from
## changes.csv, filtered, and having a phy.time column added) and
## produces a data frame with columns (time, time.to.tip, node.idx)
## where node.idx is the index of the child node of the branch
## containing the mutation.
get.mutation.info <- function(node.depths, muts)
  { 
    ## Get branch info for each mutation.
    muts <- merge(muts, node.depths, by="label")
    ## Remove the mutations on the root node's parent branch.
    muts <- muts[!is.na(muts$length),]

    ## At this point, muts may be empty (all mutations were on
    ## the root branch).
    times.to.tip <- NULL
    if (nrow(muts) > 0)
      {
        times.to.tip <- sapply(1:nrow(muts),
                               function(idx)
                               {
                                 curr.mut <- muts[idx,]
                                 return(pmin(curr.mut$first.ttt, curr.mut$second.ttt) +
                                    curr.mut$depth -
                                        curr.mut$phy.time)
                               })
      }

    return(data.frame(depth=muts$phy.time, time.to.tip=times.to.tip,
                      node.idx=match(muts$label, node.depths$label)))
  }


## Generate null simulated data.  Returns a data frame of the form
## produced by get.mutation.info.
sim.null.data <- function(node.depths, num.muts, r.seed=NULL)
  {
    if (!is.null(r.seed))
      set.seed(r.seed)

    ## Sample from the branches, with replacement.  Remove the
    ## root from consideration.
    non.root.indices <- which(!is.na(node.depths$length))
    mut.branch.indices <-
      sample(non.root.indices,
             size=num.muts,
             replace=TRUE, prob=node.depths$length[non.root.indices])

    ## Go through each of these mutations and choose a precise
    ## location for it on the branch.  Then figure out how close it is
    ## to the root and how close it is to the tip.
    mut.times <- lapply(mut.branch.indices,
                        function(idx)
                        {
                          curr.branch <- node.depths[idx,]
                          ## Place the mutation on the branch.
                          mut.time <-
                            runif(1, min=curr.branch$depth-curr.branch$length,
                                  max=curr.branch$depth)
                          dist.to.tip <-
                            min(curr.branch$first.ttt, curr.branch$second.ttt) +
                              curr.branch$depth - mut.time

                          return(c(mut.time, dist.to.tip))
                        })

    mut.info <- data.frame(depth=sapply(mut.times, function(x) x[1]),
                           time.to.tip=sapply(mut.times, function(x) x[2]),
                           node.idx=mut.branch.indices)

    return(mut.info)
  }

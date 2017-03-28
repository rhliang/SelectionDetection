## R helpers that were pulled out of LikelihoodRatioTest.r and will
## be used in CheckAgainstNull.r.

## Function that computes the null likelihood for a given codon.
## muts is 0 for no mutation, 1, 2, or 3 for mutations.

## The null likelihood is given by, if B is the total branch length, n is the total
## number of mutations, b(i) is the branch that mutation i is on, and
## N ~ Poi(lambda * B),

## P(N = n) \prod_{i=1}^n (lambda b_i)/(lambda B)
## = exp(-lambda B) (lambda B)^n/n! \prod_{i=1}^n ...

## The product is independent of lambda, so we can simply choose
## lambda to maximize the P(N=n) term; the maximizer is given by
## lambda B = n.
max.null.log.likelihood <- function(start.times, end.times, muts)
  {
    total.bl <- sum(end.times - start.times)
    mut.branches <- which(muts != 0)
    ## n == total number of mutations
    n <- length(mut.branches)

    if (n == 0)
      {
        return(NA)
      }

    ## "Marked" branches are ones with mutations on them.
    marked.bl <- end.times[mut.branches] - start.times[mut.branches]

    log.likelihood <-
      n*log(n) - n - sum(log(1:n)) - n*log(total.bl) + sum(log(marked.bl))
    return(log.likelihood)
  }

## Helper function that returns, given a branch's start time, end
## time, and time to nearest tip, the lengths of the branch that are
## in the early zone, dead zone, and hot zone.  Return a vector with 3
## elements being those lengths.
segment.lengths <- function(start.time, end.time, time.to.tip, T.s, T.h)
  {
    ## The start of the region in the hot zone (factoring in that it
    ## may overlap the early zone) is:
    hot.zone.start <- max(end.time - (T.h - time.to.tip), T.s)
    early <- max(0, min(end.time, T.s) - start.time)
    dead <- max(0, min(end.time, hot.zone.start) - max(T.s, start.time))
    hot <- max(0, end.time - max(start.time, hot.zone.start))

    return(c(early, dead, hot))
  }


## Maximum likelihood under alternative hypothesis given values for
## T_s and T_h.  This requires an optimization over the parameters
## r_1, r_2, and r_3; constraints on the parameters make it actually
## just r_1 and r_2 that we need to optimize over.
## PRE: end.times >= start.times (entry-wise).

max.alt.log.likelihood.given.switches <-
  function(start.times, end.times, times.to.tip, muts, T.s, T.h)
  { 
    mut.branches <- which(muts != 0)
    ## As above, n == total number of mutations.
    n <- length(mut.branches)

    if (n == 0)
      {
        return(NA)
      }

    ## L.e, L.d, L.h = unadjusted branch lengths in each of the early,
    ## dead, and hot zones.

    ## This produces a list where each entry corresponds to a branch
    ## and looks like [pre, mid, post].
    
    branch.segments <-
      lapply(1:length(start.times),
             function (idx)
             {
               segment.lengths(start.times[idx], end.times[idx],
                               times.to.tip[idx],
                               T.s, T.h)
             })

    L.e <- sum(sapply(branch.segments, function (x) x[1]))
    L.d <- sum(sapply(branch.segments, function (x) x[2]))
    L.h <- sum(sapply(branch.segments, function (x) x[3]))

    ## DEBUGGING
    ## cat("Number of mutations: ", n, "\n", sep="")
    ## cat("L.e, L.d, L.h = ", L.e, ", ", L.d, ", ", L.h,
    ##     "\n", sep="")

    ## Values we will get via optimization.
    max.log.likelihood <- NULL
    opt.r1 <- NULL
    opt.r2 <- NULL
    opt.r3 <- NULL

    ## Note that this first case handles both L.d != 0 and L.d == 0.
    if (L.e != 0 && L.h != 0)
      {
        ## Peg the "adjusted total branch length" to n; this maximizes the
        ## Poisson part of the likelihood.  We now need to maximize the
        ## following.  The parameter r1.r2 is a vector [r1, r2].
        branch.logs <- function(r1.r2)
          {
            r1 <- r1.r2[1]
            r2 <- r1.r2[2]
            ## The value of r3 is a function of the stuff we have above.
            r3 <- (n - r1*L.e - r2*L.d)/L.h

            ## ## DEBUGGING
            ## cat("branch.logs called with r1 = ", r1, ", r2 = ", r2,
            ##     ", r3 = ", r3,
            ##     "\n", sep="")
            
            summands <-
              sapply(mut.branches,
                     function (idx)
                     {
                       l1 <- branch.segments[[idx]][1]
                       l2 <- branch.segments[[idx]][2]
                       l3 <- branch.segments[[idx]][3]
                       return(log(r1*l1 + r2*l2 + r3*l3))
                     })

            ## ## DEBUGGING
            ## print(summands)
            
            return(sum(summands))
          }
        
        ## On the (r1, r2)-plane (r1 is the x-axis, r2 is the y-axis), the
        ## admissible region is a triangle with corners at (0, 0), (n/L.e,
        ## 0), and (n/(L.e+L.d+L.h), n/(L.e+L.d+L.h)).  So, a comfortable point
        ## in the middle is:
        starting.value <- c(n/(L.e+L.d+L.h), n/(2*(L.e+L.d+L.h)))
        
        ## ## DEBUGGING
        ## cat("Starting value (inside triangular region) = (",
        ##     starting.value[1], ", ", starting.value[2], ")\n", sep="")
        ## cat("L.e, L.d, L.h = ", L.e, ", ", L.d, ", ", L.h, "\n", sep="")
        
        ## The constraints we need to satisfy are:
        ## r1 - r2 >= 0
        ## r2 >= 0
        ## r3 >= r2, or n - r1*L.e - r2*(L.d+L.h) >= 0.
        constr.coeffs <- matrix(
                                c(1, -1,
                                  0, 1,
                                  -L.e, -L.d-L.h),
                                nrow=3, byrow=TRUE)
        constr.consts <- c(0, 0, -n)

        opt.result <- constrOptim(starting.value, branch.logs, grad=NULL,
                                  ui=constr.coeffs, ci=constr.consts,
                                  control=list(fnscale=-1))
    
        max.log.likelihood <- -n - sum(log(1:n)) + opt.result$value
        opt.r1 <- opt.result$par[1]
        opt.r2 <- opt.result$par[2]
        opt.r3 <- (n - opt.r1*L.e - opt.r2*L.d)/L.h
        
      } else if (L.e == 0 && L.d != 0 && L.h != 0) {
        
        ## In this case, we're really only optimizing over one
        ## parameter (r2), because we are still pegging the "adjusted
        ## total branch length" to n, so r3 is just (n - r2*L.d)/L.h.
        branch.logs <- function(r2)
          {
            r3 <- (n - r2*L.d)/L.h

            ## This is as above, but obviously the length before T1 is
            ## 0.
            summands <-
              sapply(mut.branches,
                     function (idx)
                     {
                       l2 <- branch.segments[[idx]][2]
                       l3 <- branch.segments[[idx]][3]
                       return(log(r2*l2 + r3*l3))
                     })
            return(sum(summands))
          }

        ## In this regime, r2 must satisfy 0 <= r2 <= n/(L.d+L.h).
        opt.result <- optimize(branch.logs, interval=c(0, n/(L.d+L.h)),
                               maximum=TRUE)
        
        max.log.likelihood <- -n - sum(log(1:n)) + opt.result$objective
        opt.r1 <- NA
        opt.r2 <- opt.result$maximum
        opt.r3 <- (n - opt.r2*L.d)/L.h
        
      } else if (L.e != 0 && L.d != 0 && L.h == 0) {
        
        ## In this case, we're only optimizing over one
        ## parameter (r1), and r2 is just (n - r1*L.e)/L.d.
        branch.logs <- function(r1)
          {
            r2 <- (n - r1*L.e)/L.d

            ## This is as above, but obviously the length after T2 is
            ## 0.
            summands <-
              sapply(mut.branches,
                     function (idx)
                     {
                       l1 <- branch.segments[[idx]][1]
                       l2 <- branch.segments[[idx]][2]
                       return(log(r1*l1 + r2*l2))
                     })
            return(sum(summands))
          }

        ## In this regime, r1 must satisfy
        ## n/(L.e+L.d) <= r1 <= n/L.e.
        opt.result <- optimize(branch.logs, interval=c(n/(L.e+L.d), n/L.e),
                               maximum=TRUE)
        
        max.log.likelihood <- -n - sum(log(1:n)) + opt.result$objective
        opt.r1 <- opt.result$maximum
        opt.r2 <- (n - opt.r1*L.e)/L.d
        opt.r3 <- NA
        
      } else {
        ## Now we know that we only have none non-zero L.e, L.d, or L.h.
        ## The maximum is simply the same as that of the null
        ## hypothesis.
        max.log.likelihood <- max.null.log.likelihood(start.times, end.times,
                                                      muts)
        if (L.e != 0)
          {
            opt.r1 <- n/L.e
          } else {
            opt.r1 <- NA
          }
        if (L.d != 0)
          {
            opt.r2 <- n/L.d
          } else {
            opt.r2 <- NA
          }
        if (L.h != 0)
          {
            opt.r3 <- n/L.h
          } else {
            opt.r3 <- NA
          }
        
      }
    
    return(list(max.value=max.log.likelihood, r1=opt.r1, r2=opt.r2, r3=opt.r3))
  }

## Now, we'd want to optimize the above over T.s and T.h.
## We need:
## min.time <= T.s <= max.time
## T.h <= max.time - T.s
max.alt.log.likelihood <- function(start.times, end.times, times.to.tip, muts)
  {
    if (sum(muts != 0) == 0)
      {
        return(NA)
      }
    
    ## We need a wrapper for max.alt.log.likelihood.given.switches so
    ## it can be optimized over.
    ll.helper <- function(switch.times)
      {
        T.s <- switch.times[1]
        T.h <- switch.times[2]

        ## DEBUGGING
        ## cat("ll.helper called with parameters ", T.s, ", ", T.h, "\n", sep="")

        return(max.alt.log.likelihood.given.switches(start.times, end.times,
                                                     times.to.tip, muts,
                                                     T.s, T.h)$max.value)
      }

    ## Program in the above constraints: they translate to
    ## T.s - min.time >= 0
    ## -T.s + max.time >= 0
    ## -T.s - T.h + max.time >= 0.
    constr.coeffs <- rbind(c(1, 0),
                           c(-1, 0),
                           c(-1, -1))
    constr.consts <- c(min.time, -max.time, -max.time)

    ## Reasonable starting values:
    starting.value <-
      c(min.time + (max.time-min.time)/3, (max.time-min.time)/10)

    ## DEBUGGING
    ## cat("starting.value = (", starting.value[1], ", ",
    ##     starting.value[2], ")\n", sep="")

    opt.over.times <- constrOptim(starting.value, ll.helper, grad=NULL,
                                  ui=constr.coeffs, ci=constr.consts,
                                  control=list(fnscale=-1))

    return(opt.over.times)
  }

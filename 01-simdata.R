## This source code is licensed under the FreeBSD license (see "LICENSE" file)
## (c) 2013 Felix Schönbrodt

## ======================================================================
## This file runs the Monte-Carlo simulations on a multicore machine.
## The numbers of used processor cores is defined in nCore. It can also be set to 1.
## ======================================================================


setwd("~/cogsci/EM3/corridor of stability")
library(data.table)
library(ggplot2)
library(reshape2)

# functions needed for multicore support - the script can also be run on a single core if nCore is set to 1
library(pacman)
p_load(foreach, doMC)

nCore <- 1
registerDoMC(nCore)
print(paste("Registered multicores: ", getDoParWorkers()))

# Load Ruscio's function for inducing a correlation into a bivariate data set with empirical marginal distributions, as well as some helper functions
source("Ruscio_GenData.R")

## ======================================================================
## Simulation parameters
## ======================================================================

popsize <- 100000	# size of population
n.max <- 1000		# maximum sample size
n.min <- 20			# minimum sample size
B <- 10000			# number of bootstrapped trajectories
rs <- seq(.1, .7, by=.2)	# which correlations should be induced in the population?


## ======================================================================
## Load marginal distributions: Should be a list of length two with raw data. Each list element contains one vector of data which defines the marginal distribution. List elements can have differing lengths.
## ======================================================================

# Simulate Gaussian distributions
GAUSS <- list(rnorm(1000000), rnorm(1000000))

# Or, load empirical distributions, e.g. the data set provided by Ted Micceri:
# http://www.freewebs.com/tedstats/Files/Real_Data.zip

# TESTDIST is used in the functions below - root it to any list of desired distributions
TESTDIST <- GAUSS

## ======================================================================
## The functions for the simulation - no changes should be necessary below this point
## ======================================================================

# Helper: Transform correlation to Fisher's Z
r2Z <- function(r) {
	return(0.5 * log((1 + r)/(1 - r)))
}

# Helper: Recode  Fisher's to correlation
Z2r <- function(Z) {
	return((exp(2*Z)-1)/(exp(2*Z)+1))
}

# calculate all correlations when sample size increases from n.min to n.max
corEvol <- function(df, n.min=20, stepsize=1) {
	res <- matrix(NA, nrow=length(seq(n.min, nrow(df), by=stepsize)), ncol=2)
	for (i in seq(n.min, nrow(df), by=stepsize)) {
		res[i-n.min+1, 1:2] <- c(n=i, R=cor(df[1:i,1], df[1:i,2]))
	}
	colnames(res) <- c("n", "R") 	
	return(res)	
}

#' @param DIST A list of two vectors containing the marginal distributions
#' @param n.max Maximum sample size to be simulated
#' @param B Number of bootstrapped trajectories
#' @param rs A vector of correlations to be imposed on the bivariate data set
#' @param replace Sample with or without replacement?
simCorEvol <- function(DIST, n.min=20, n.max=1000, B=10, rs=seq(.2, .4, by=.2), replace=TRUE, nCore=1) {

	if (nCore > 1 & nCore > getDoParWorkers()) stop("nCore does not match the number of registered multicores!")

	# Replications are split up in nCore pieces to harvest multicores
	if (B %% nCore != 0) {
		B <- ceiling(B/nCore)*nCore
		print(paste("Warning: Replications is not dividible by", nCore, "(the number of used processor cores) - set replications to", B))
	}

	# outer loop for the distribution of replications across cores
	res1 <- foreach(batch=1:nCore, .combine=rbind) %dopar% {

		maxcount <- (n.max-n.min+1)*(B/nCore)*length(rs)
		res <- matrix(NA, ncol=6, nrow=maxcount)

		counter <- 1
		for (r in rs) {
				
			print(paste("\nComputing next population at ", Sys.time()))

			# Generate population with specified rho
			pop <- GenData(DIST, rho=r, N=popsize)

			for (rep in 1:(B/nCore)) {		# rep = replication	
				
				# draw a boostrap sample from the population
				sam <- pop[sample(1:nrow(pop), size=n.max, replace=replace), ]	
				f <- corEvol(sam)
				f <- cbind(f, rho = r)
				f <- cbind(f, repl=rep)
				f <- cbind(f, batch=batch)		# batch keeps the number of the multicore processor
				f <- cbind(f, unique=(batch-1)*(B/nCore)*length(rs) + counter)
				res[(((counter-1)*nrow(f))+1):(counter*nrow(f)), 1:6] <- f
				counter <- counter + 1
			}
		}
		
		colnames(res) <- c("n", "R", "rho", "repl", "batch", "unique")
		gc()
		return(res)
	}
	
	return(res1)
}

s1 <- system.time(sim1 <- simCorEvol(TESTDIST, n.max=n.max, B=B, rs=rs, replace=TRUE, nCore=nCore))
print(s1)



## ======================================================================
## For subsequent analyses: t0 keeps the raw data as a data.table
## ======================================================================

t0 <- data.table(sim1)
rm(sim1)
setkey(t0, rho, unique)

dir.create("results", showWarnings = FALSE)
save(t0, file="results/sim1.RData")
print("Simulation results saved to results/sim1.RData!")
print("Simulation finished.")

## This source code is licensed under the FreeBSD license (see "LICENSE" file)
## (c) 2013 Felix Schönbrodt

# ---------------------------------------------------------------------
# First: Load the simulation results. Maybe you have to adjust the path:

load("results/sim1.RData")


## ======================================================================
## Calculate the distribution of POS
## ======================================================================


# finds consecutive numbers in a vector, starting from the tail
findConsecutive <- function(x, consec=2) {
	if (length(x) < consec) return(NA)
	x.r <- x[length(x):1]
	x.r.diff <- x.r[1:(length(x)-1)] - x.r[2:length(x)]
	for (i in 1:(length(x.r.diff)-consec+2)) {
		if (all(x.r.diff[i:(i+consec-2)] == 1)) {return(x[length(x)-i+1])}
	}
	return(NA)
}


# findPOS returns the point of stability of a single trajectory
# @param n Vector of sample sizes
# @param R Vector of correlations
# @param rho The true correlation
# @param w Half-width of corridor of stability
# @param inarow How many breaks must take place in a row to qualify for a "real" break?
findPOS <- function(n, R, rho, w, inarow=1) {
	COS <- Z2r(r2Z(rho + c(-1, 1)*w))
	breaks <- which(R < COS[1] | R > COS[2])
	if (is.infinite(max(breaks))) {
		n.stable <- min(n)	# no break: POS = minimal sample size
	} else {
		if (inarow==1) {
				BREAK <- max(breaks)	# get first break (seen from the tail)
			} else {
				BREAK <- findConsecutive(breaks, inarow)
			}
		if (is.na(BREAK)) {
			n.stable <- min(n)	# no break? POS = minimal sample size
		} else {
			n.stable <- n[BREAK]
		}
	}
	return(data.frame(POS=n.stable, w=w))
}


# helper function: return quantiles of vector as data.frame
qu <- function(x, prob=c(.80, .90, .95)) {
	data.frame(POS=quantile(x, prob=prob, na.rm=TRUE), confidence=prob)
}

# ---------------------------------------------------------------------
# Find POS with 1 consecutive breaks

system.time(stab1.10 <- t0[, findPOS(n, R, rho, w=.10, inarow=1), by=c("unique", "rho")])
system.time(stab1.15 <- t0[, findPOS(n, R, rho, w=.15, inarow=1), by=c("unique", "rho")])
system.time(stab1.20 <- t0[, findPOS(n, R, rho, w=.20, inarow=1), by=c("unique", "rho")])

stab1 <- rbind(stab1.10, stab1.15, stab1.20)
stab1.summary <- stab1[, qu(POS), by=c("rho", "w")]
stab1.summary

# ---------------------------------------------------------------------
# Find POS with 2 consecutive breaks

system.time(stab2.10 <- t0[, findPOS(n, R, rho, w=.10, inarow=2), by=c("unique", "rho")])
system.time(stab2.15 <- t0[, findPOS(n, R, rho, w=.15, inarow=2), by=c("unique", "rho")])
system.time(stab2.20 <- t0[, findPOS(n, R, rho, w=.20, inarow=2), by=c("unique", "rho")])

stab2 <- rbind(stab2.10, stab2.15, stab2.20)

stab2.summary <- stab2[, qu(POS), by=c("rho", "w")]
stab2.summary


# ---------------------------------------------------------------------
# Find POS with 3 consecutive breaks

system.time(stab3.10 <- t0[, findPOS(n, R, rho, w=.10, inarow=3), by=c("unique", "rho")])
system.time(stab3.15 <- t0[, findPOS(n, R, rho, w=.15, inarow=3), by=c("unique", "rho")])
system.time(stab3.20 <- t0[, findPOS(n, R, rho, w=.20, inarow=3), by=c("unique", "rho")])

stab3 <- rbind(stab3.10, stab3.15, stab3.20)

stab3.summary <- stab3[, qu(POS), by=c("rho", "w")]
stab3.summary

# ---------------------------------------------------------------------
# Combine the POS distributions

stab.combined <- cbind(stab1.summary, POS2=stab2.summary[, POS], POS3=stab3.summary[, POS])
stab.combined$diff2 <- stab.combined$POS - stab.combined$POS2
stab.combined$diff3 <- stab.combined$POS - stab.combined$POS3

summary(stab.combined$diff2)
summary(stab.combined$diff3)

qplot(rho, POS, data=stab.combined, facets=~confidence, geom=c("point", "line"), group=w) + theme_bw()

## Plot the distribution of POS
ggplot(stab1, aes(x=POS)) + geom_density() + facet_grid(rho~w, scales="free") + theme_bw()


# ---------------------------------------------------------------------
# The final POS table

ft1 <- stab1.summary[, round(mean(POS)), by=c("rho", "w", "confidence")]
(ft2 <- dcast(ft1, rho ~ confidence + w, value.var="V1"))

## save results
save(stab1, file="results/stab1.RData")
save(stab2, file="results/stab2.RData")
save(stab3, file="results/stab3.RData")
write.csv(finaltable, file="results/finaltable.csv")




## ======================================================================
## Finally: Do some sanity checks ...
## ======================================================================


# for each n: compare theoretical CI and bootstrapped percentiles

RHO <- .6
PQ <- t0[J(RHO), qu(R, prob=c(.025, .975)), by="n"]

tCI <- data.frame(n=min(PQ$n):max(PQ$n))	# tQ = theoretical CI
tCI$CI.lower <- Z2r(r2Z(RHO) - qnorm(.975)*(1/sqrt(tCI$n-3)))
tCI$CI.upper <- Z2r(r2Z(RHO) + qnorm(.975)*(1/sqrt(tCI$n-3)))

ggplot(PQ, aes(x=n, y=POS, linetype=factor(confidence))) + geom_line() + theme_bw() + geom_hline(yintercept=RHO, color="blue") + geom_line(data=tCI, aes(x=n, y=CI.upper), color="blue", linetype="dotted") + geom_line(data=tCI, aes(x=n, y=CI.lower), color="blue", linetype="dotted")

## --> bootstrapped percentiles of Gaussian marginal distributions are identical to the theoretical quantiles
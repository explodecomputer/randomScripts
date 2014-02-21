n <- 1000

# A causes B

z <- rnorm(n)
A <- z + rnorm(n)
B <- A + rnorm(n)


Ahat <- lm(A ~ z)$fitted
anova(lm(B ~ Ahat))$P[1]

Bhat <- lm(B ~ z)$fitted
anova(lm(A ~ Bhat))



# B causes A

z <- rnorm(n)
B <- rnorm(n)
A <- z + B

anova(lm(B ~ A))

Ahat <- lm(A ~ z)$fitted
anova(lm(B ~ Ahat))

anova(lm(B ~ z))

Bhat <- lm(B ~ z)$fitted
anova(lm(A ~ Bhat))


# Function to infer causal direction

twoStageLS <- function(A, B, zA, zB)
{
	Ahat_zA <- lm(A ~ zA)$fitted
	Ahat_zB <- lm(A ~ zB)$fitted
	Bhat_zA <- lm(B ~ zA)$fitted
	Bhat_zB <- lm(B ~ zB)$fitted

	pvals <- array(0,4)
	nom <- c("A ~ Bhat : zA", "A ~ Bhat : zB", "B ~ Ahat : zA", "B ~ Ahat : zB")

	pvals[1] <- anova(lm(A ~ Bhat_zA))$P[1]
	pvals[2] <- anova(lm(A ~ Bhat_zB))$P[1]
	pvals[3] <- anova(lm(B ~ Ahat_zA))$P[1]
	pvals[4] <- anova(lm(B ~ Ahat_zB))$P[1]

	return(data.frame(nom = nom, pvals = pvals))
}

inferCausality <- function(dat, threshold = 0.01)
{
	sig <- dat$pvals < threshold
	AtoB <- c(T, F, T, T)
	BtoA <- c(T, T, F, T)
	conf <- c(T, F, F, T)

	if(all(sig == AtoB))
	{
		return("AtoB")
	} else if(all(sig == BtoA)) {
		return("BtoA")
	} else if(all(sig == conf)) {
		return("Confounder")
	} else {
		return("Unknown")
	}
}



# A causes B

zA <- rnorm(n)
zB <- rnorm(n)
A <- zA + rnorm(n)
B <- A + rnorm(n) + zB

inferCausality(twoStageLS(A, B, zA, zB))


# B causes A

zA <- rnorm(n)
zB <- rnorm(n)
B <- rnorm(n) + zB
A <- zA + rnorm(n) + B

inferCausality(twoStageLS(A, B, zA, zB))


# confounder causes A and B

zA <- rnorm(n)
zB <- rnorm(n)
confounder <- rnorm(n)
B <- rnorm(n) + zB + confounder
A <- zA + rnorm(n) + confounder

inferCausality(twoStageLS(A, B, zA, zB))


# pleiotropy of zB but A causes B

zA <- rnorm(n)
zB <- rnorm(n)
A <- zA + rnorm(n) + zB/10
B <- A + rnorm(n) + zB

inferCausality(twoStageLS(A, B, zA, zB))
 


# Testing
# A is methylation, B is trait
# Correlation 0.05, 0.1, ..., 0.5
# var(zA) = 0.05, 0.1, ..., 1.0
# var(zB) = 0.01, 0.03, ..., 0.21
# Pleiotropic 
# Confounding, AtoB, BtoA


dat <- expand.grid(
	n = c(100, 1000, 10000, 100000),
	rsq = seq(0, 0.1, by=0.05),
	varZA = seq(0, 1, by=0.05),
	varZB = seq(0, 0.25, by=0.05),
	threshold = c(0.05, 0.01, 0.001),
	causality = c("Confounder", "AtoB", "BtoA"),
	propCorrect = NA
)

test <- with(dat, varZA + rsq)

dim(dat)


#' Create a genotype in HWE
#'
#' @param n number of individuals
#' @param p allele frequency
#' @export
#' @return array of 0s,1s,2s of length n
makeGeno <- function(n, p)
{
	pr <- c(p^2, 2*(1-p)*p, (1-p)^2)
	x <- rbinom(n, 2, p)
	return(x)
}


#' Create a phenotype using arbitrary number of known causal inputs
#'
#' @param cors array of variances for each input
#' @param ... Inputs, each being a vector of the same length
#' @return array of length of inputs
#' @examples \dontrun{
#' g1 <- makeGeno(1000, 0.5)
#' g2 <- makeGeno(1000, 0.3)
#' p <- makePhen(cors=c(0.2, 0.1, 0.15, 0.4, 0.15), g1, g2, rnorm(1000), rnorm(1000))
#'}
makePhen <- function(cors, ...)
{
	li <- list(...)
	stopifnot(length(li) == (length(cors)))
	stopifnot(sum(cors) <= 1)
	cors <- c(cors, 1-sum(cors))

	inputs <- do.call(cbind, li)
	n <- nrow(inputs)

	inputs <- cbind(inputs, rnorm(n))
	l <- ncol(inputs)

	for(i in 1:l)
	{
		inputs[,i] <- (inputs[,i] - mean(inputs[,i])) / sd(inputs[,i]) * sqrt(cors[i])
	}

	y <- apply(inputs, 1, sum)
	return(y)
}


createAtoB <- function(n, varZA, varZB, rsq)
{
	zA <- makeGeno(n, 0.5)
	zB <- makeGeno(n, 0.5)
	A <- makePhen(cors=varZA, zA)
	B <- makePhen(cors=c(rsq, varZB), A, zB)
	return(list(stat = 1, zA = zA, zB = zB, A = A, B = B))
}

createBtoA <- function(n, varZA, varZB, rsq)
{
	zA <- makeGeno(n, 0.5)
	zB <- makeGeno(n, 0.5)
	B <- makePhen(cors=varZB, zB)
	A <- makePhen(cors=c(rsq, varZA), B, zA)
	return(list(stat = 1, zA = zA, zB = zB, A = A, B = B))
}

createConfounder <- function(n, varZA, varZB, rsq)
{
	zA <- makeGeno(n, 0.5)
	zB <- makeGeno(n, 0.5)
	confounder <- rnorm(n)
	A <- makePhen(cors=c(varZA, sqrt(rsq)), zA, confounder)
	B <- makePhen(cors=c(varZB, sqrt(rsq)), zB, confounder)
	return(list(stat = 1, zA = zA, zB = zB, A = A, B = B))
}

simulationFramework <- function(dat, row, nrep)
{
	result <- array("", nrep)
	func <- get(paste("create", dat$causality[row], sep=""))
	for(i in 1:nrep)
	{
		variables <- func(dat, row)
		reg <- twoStageLS(variables$A, variables$B, variables$zA, variables$zB)
		result[i] <- inferCausality(reg, dat$threshold[row])
	}
	dat$propCorrect[row] <- sum(result == dat$causality[row]) / nrep
	return(list(dat = dat[row], result = result))
}


arguments <- commandArgs(T)




runSimulation




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
	stopifnot(length(li) == (length(cors)-1))
	cors <- cors / sum(cors)

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


#================================================================#
#================================================================#


# Create a system

# Methylation causes expression
# Expression causes height
# Methylation, expression and height each caused by a SNP also

n <- 50000
g1 <- makeGeno(n, 0.3)
meth <- makePhen(c(0.9, 0.1), g1)

g2 <- makeGeno(n, 0.3)
expr <- makePhen(c(0.4, 0.4, 0.2), meth, g2)

g3 <- makeGeno(n, 0.3)
height <- makePhen(c(0.01, 0.01, 0.98), expr, g3)


summary(lm(height ~ g3 + expr))
summary(lm(height ~ g3 + g2))
summary(lm(expr ~ height))


mod <- lm(cbind(meth, expr, height) ~ cbind(g1, g2, g3))
mod2 <- manova(cbind(meth, expr, height) ~ cbind(g1, g2, g3))

summary(mod)
anova(mod)
summary(mod2)


g = rbinom(100,2,0.2)
y = matrix(rnorm(1000), nrow = 100)
summary(lm(g~y))
summary(manova(y~g))


mod1 <- anova(mod)
mod2 <- summary(mod)
mod2

summary(lm(meth ~ cbind(g1, g2, g3)))
summary(lm(expr ~ cbind(g1, g2, g3)))
summary(lm(height ~ cbind(g1, g2, g3)))


dat <- cbind(expr, meth, height, g1, g2, g3)

head(dat)

cor(dat)
solve(cov(dat))

heatmap(solve(cov(dat)))


library(car)
Anova(mod)
summary(Anova(mod))




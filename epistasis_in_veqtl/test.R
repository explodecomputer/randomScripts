library(Rcpp)
library(microbenchmark)
library(snpStats)

sourceCpp("source.cpp")

x <- sample(1:3, 100, rep=T)
y <- rnorm(100)
X <- matrix(sample(1:3, 2000, rep=T), 100, 20)

microbenchmark(
	scan8df(y, X, x),
	for(i in 1:20) anova(lm(y ~ as.factor(x)*as.factor(X[,i])))
)

microbenchmark(
	scan4df(y, X, x),
	for(i in 1:20) anova(lm(y ~ as.factor(x):as.factor(X[,i])))
)


dat <- read.plink(bed="/Users/ghemani/Documents/WORK/data/chinese_ra/dataSASfilter.Quality.hwe1e-3.autosome")

y <- rnorm(nrow(dat$fam))

test <- scan8df(y, dat$genotypes@.Data, dat$genotypes@.Data[,1])
dim(test)
head(test)


table(is.na(test))

table(test[,1] == 8)





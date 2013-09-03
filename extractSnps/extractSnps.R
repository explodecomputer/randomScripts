library(snpStats)


findChromosomes <- function(snp, plinkrt)
{
	cmd <- paste("grep ", snp, " ", plinkrt, ".bim | head -n 1 | cut -d \":\" -f 2 | cut -f 1", sep="")
	print(snp)
	chr <- system(cmd, intern=TRUE)
	return(chr)
}


findChromosomesAll <- function(snplist, plinkrt)
{
	chr <- snplist
	for(i in 1:length(snplist))
	{
		a <- findChromosomes(snplist[i], plinkrt)
		chr[i] <- ifelse(is.null(a), NA, a)
	}
	dat <- data.frame(snp = snplist, chr = as.numeric(chr))
	snpdat$snp <- as.character(snpdat$snp)
	dat <- dat[order(dat$chr), ]
	return(dat)
}


extractSnps <- function(snpnames, plinkrt)
{
	require(snpStats)
	rawdata <- read.plink(bed=plinkrt, select.snps=snpnames)
	return(rawdata)
}


extractSnpsAll <- function(snpdat, plinkrt)
{
	chr <- unique(snpdat$chr)
	l <- length(chr)
	dat <- NULL
	for(i in 1:l)
	{
		nom <- gsub("\\*", chr[i], plinkrt)
		snps <- snpdat$snp[snpdat$chr == chr[i]]
		cat(i, ":", chr[i], ":", snps, "\n")
		dat <- cbind(dat, extractSnps(snps, nom))
	}
	return(dat)
}


ar <- commandArgs(T)
snplistfile <- ar[1]
plinkrt <- ar[2]


snplist <- scan(snplistfile, what="character")
snpdat <- findChromosomesAll(snplist, plinkrt)
extractSnpsAll(snpdat, plinkrt)


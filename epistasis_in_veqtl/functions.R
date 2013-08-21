library(Rcpp)
library(snpStats)
sourceCpp("source.cpp")

load("data/residuals_all.RData")
geno <- read.plink(bed="data/clean_geno_final")
all(geno$fam$pedigree == rownames(resphen))

table(hits$SNP %in% geno$map$snp.name )

hits <- subset(read.table("data/variance_effects_only_cis5_plus_trans9.txt", header=T), TRAIT %in% probeinfo$PROBE_ID & SNP %in% geno$map$snp.name)

getProbeAndSnp <- function(hits, bim, probeinfo, i)
{
	a <- which(probeinfo$PROBE_ID == hits$TRAIT[i])
	b <- which(bim$snp.name == hits$SNP[i])
	return(c(a,b))
}

jid <- 2
index <- getProbeAndSnp(hits, geno$map, probeinfo, jid)
index

a <- Sys.time()
output <- scan4df(resphen[,index[1]], geno$genotypes@.Data, geno$genotypes@.Data[, index[2]])
Sys.time() - a

head(output)


makeOutput <- function(output, hits, bim, i)
{
	pval8 <- -log10(pf(output[,3], output[,1], output[,2], lower.tail=FALSE))
	pval4 <- -log10(pf(output[,4], 4, output[,2], lower.tail=FALSE))
	a <- data.frame(bim, vChr = hits$CHR[i], vSNP = hits$SNP[i], probe = hits$TRAIT[i], df1 = output[, 1], pval8 = pval8, pval4 = pval4)
	return(subset(a, df1 == 8))

}

a <- makeOutput(output, hits, geno$map, jid)
dim(a)
max(a$pval8)

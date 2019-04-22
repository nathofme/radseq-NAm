#### LFMM
library(foreach)
#install lfmm
devtools::install_github("bcm-uga/lfmm")
library(lfmm)
#install LEA to prepare input files
#source("https://bioconductor.org/biocLite.R")
#biocLite("LEA")
library(LEA)

# EUSTrad.env includes PCA scores from popenvPCA (see spatialstarlings.R)
# create input file for genotypes
vcf2lfmm("/Users/nataliehofmeister/Documents/starlingRAD/analysis/EUSTallSNPr8maf01p.vcf", force = TRUE)
vcf2lfmm("/Users/nataliehofmeister/Documents/starlingRAD/analysis/EUSTallSNPr8.vcf", force = TRUE)
vcf2lfmm("/Users/nataliehofmeister/Documents/starlingRAD/analysis/EUSTallSNPr95.vcf", force = TRUE)



# K is number of latent factors to model confounding effects like population structure, background variation
# K can be derived from STRUCTURE or snmf (as in LEA manual)
lfmm.PC1 = NULL
lfmm.PC1 = lfmm("EUSTallSNPr95.lfmm",
                "EUSTrad.PC1.env",
                K = 1:3,
                repetitions = 10,
                iterations = 10000,
                burnin = 5000,
                project = "continue")

lfmm.PC2 = NULL
lfmm.PC2 = lfmm("EUSTallSNPr8.lfmm",
                "EUSTrad.PC2.env",
                K = 1:3,
                repetitions = 10,
                iterations = 10000,
                burnin = 5000,
                project = "continue")

lfmm.PC3 = NULL
lfmm.PC3 = lfmm("EUSTallSNPr8.lfmm",
                "EUSTrad.PC3.env",
                K = 1:3,
                repetitions = 10,
                iterations = 10000,
                burnin = 5000,
                project = "continue")

lfmm.PC1 = NULL
lfmm.PC1 = lfmm("EUSTallSNPr8.lfmm",
                "EUSTrad.PC1.env",
                K = 1:3,
                repetitions = 30,
                iterations = 10000,
                burnin = 5000,
                project = "continue")


# CHECK NULL HYPOTHESIS
# combine z-scores from multiple runs 

# 1 # Record z-scores from the 5 runs in the zs matrix
# requires that you specify d if used more than one variable in .env
zs.K1 = z.scores(lfmm.PC1, K = 1)
zs.K2 = z.scores(lfmm.PC1, K = 2)
zs.K3 = z.scores(lfmm.PC1, K = 3)

# 2 #Combine z-scores using the median
zs.K1.median = apply(zs.K1, MARGIN = 1, median)

# 3 # compute genomic inflation factor
lambda.K1 = median(zs.K1.median^2)/qchisq(0.5, df = 1)
lambda.K1

# 4 # compute adjusted p-values from the combined z-scores
adj.p.values.K1 = pchisq(zs.K1.median^2/lambda.K1, df = 1, lower = FALSE)

# 5 # visualize distribution of p-values
par(mfrow = c(1,1))

quartz()
# histogram should be flat with a peak near zero 
hist(adj.p.values.K1, col = "lightblue")
# loci with significant associations on the right of this manhattan plot
plot(-log10(adj.p.values.K1), pch = 19, col = "blue", cex = .7)

library(car)
qqPlot(-log10(adj.p.values.K1))

# Benjamini-Hochberg to adjust p-values for multiple testing
## FDR control: Benjamini-Hochberg at level q
## L = number of loci
L = 15038
#fdr level q
q = 0.1
w.K1 = which(sort(adj.p.values.K1) < q * (1:L)/L)
candidates.K1.PC1 = order(adj.p.values.K1)[w.K1]
candidates.K1.PC1
# candidates.K1 stores the list of candidate loci

# estimate FDR and true positive rates
estimated.FDR.K1.PC1 = length(which(candidates.K1.PC1 <= 900))/length(candidates.K1.PC1)
estimated.TP.K1.PC1 = length(which(candidates.K1.PC1 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.K1.PC1, "True Positive:", estimated.TP.K1.PC1))
# FDR 0.051, TPR = 28.35


# run with K=2
zs.K2.median = apply(zs.K2, MARGIN = 1, median)
lambda.K2 = median(zs.K2.median^2)/qchisq(0.5, df = 1)
lambda.K2
adj.p.values.K2 = pchisq(zs.K2.median^2/lambda.K2, df = 1, lower = FALSE)
par(mfrow = c(2,1))
hist(adj.p.values.K2, col = "lightblue")
plot(-log10(adj.p.values.K2), pch = 19, col = "blue", cex = .7)
qqPlot(-log10(adj.p.values.K2))
w.K2 = which(sort(adj.p.values.K2) < q * (1:L)/L)
candidates.K2.PC1 = order(adj.p.values.K2)[w.K2]
estimated.FDR.K2.PC1 = length(which(candidates.K2.PC1 <= 900))/length(candidates.K2.PC1)
estimated.TP.K2.PC1 = length(which(candidates.K2.PC1 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.K2.PC1, "True Positive:", estimated.TP.K2.PC1))
# FDR 0.049, TPR = 25.53

# K=3
zs.K3.median = apply(zs.K3, MARGIN = 1, median)
lambda.K3 = median(zs.K3.median^2)/qchisq(0.5, df = 1)
lambda.K3
adj.p.values.K3 = pchisq(zs.K3.median^2/lambda.K3, df = 1, lower = FALSE)
par(mfrow = c(2,1))
hist(adj.p.values.K3, col = "lightblue")
plot(-log10(adj.p.values.K3), pch = 19, col = "blue", cex = .7)
qqPlot(-log10(adj.p.values.K3))
w.K3 = which(sort(adj.p.values.K3) < q * (1:L)/L)
candidates.K3.PC1 = order(adj.p.values.K3)[w.K3]
candidates.K3.PC1
estimated.FDR.K3.PC1 = length(which(candidates.K3.PC1 <= 900))/length(candidates.K3.PC1)
estimated.TP.K3.PC1 = length(which(candidates.K3.PC1 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.K3.PC1, "True Positive:", estimated.TP.K3.PC1))
# FDR 0.046, TPR 27.99

# export to use in manhattan plot
LFMMadjp <- cbind(-log10(adj.p.values.K1),-log10(adj.p.values.K2),-log10(adj.p.values.K3))
head(LFMMadjp)
write.csv(LFMMadjp,file="LFMMpvalues.PC1.csv")
LFMM.PC1.candidates <- Reduce(intersect, list(candidates.K1.PC1,candidates.K2.PC1,candidates.K3.PC1))
LFMM.PC1.candidates

head(LFMM.PC1.candidates)


# PC2
zs.PC2.K1 = z.scores(lfmm.PC2, K = 1)
zs.PC2.K2 = z.scores(lfmm.PC2, K = 2)
zs.PC2.K3 = z.scores(lfmm.PC2, K = 3)
zs.PC2.K1.median = apply(zs.PC2.K1, MARGIN = 1, median)

# 3 # compute genomic inflation factor
lambda.PC2.K1 = median(zs.PC2.K1.median^2)/qchisq(0.5, df = 1)
lambda.PC2.K1
adj.p.values.PC2.K1 = pchisq(zs.PC2.K1.median^2/lambda.PC2.K1, df = 1, lower = FALSE)
par(mfrow = c(2,1))

# histogram should be flat with a peak near zero 
hist(adj.p.values.PC2.K1, col = "lightblue")
plot(-log10(adj.p.values.PC2.K1), pch = 19, col = "blue", cex = .7)
qqPlot(-log10(adj.p.values.PC2.K1))
# Benjamini-Hochberg to adjust p-values for multiple testing
L = 15038
q = 0.01
w.PC2.K1 = which(sort(adj.p.values.PC2.K1) < q * (1:L)/L)
candidates.PC2.K1 = order(adj.p.values.PC2.K1)[w.PC2.K1]
candidates.PC2.K1

estimated.FDR.PC2.K1 = length(which(candidates.PC2.K1 <= 900))/length(candidates.PC2.K1)
estimated.TP.PC2.K1 = length(which(candidates.PC2.K1 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.PC2.K1, "True Positive:", estimated.TP.PC2.K1))


# run with K=2
zs.PC2.K2 = z.scores(lfmm.PC2, K = 2, d=1)
zs.PC2.K2.median = apply(zs.PC2.K2, MARGIN = 1, median)
lambda.PC2.K2 = median(zs.PC2.K2.median^2)/qchisq(0.5, df = 1)
lambda.PC2.K2
adj.p.values.PC2.K2 = pchisq(zs.PC2.K2.median^2/lambda.PC2.K2, df = 1, lower = FALSE)
par(mfrow = c(2,1))
hist(adj.p.values.PC2.K2, col = "lightblue")
plot(-log10(adj.p.values.PC2.K2), pch = 19, col = "blue", cex = .7)
w.PC2.K2 = which(sort(adj.p.values.PC2.K2) < q * (1:L)/L)
candidates.PC2.K2 = order(adj.p.values.PC2.K2)[w.PC2.K2]
estimated.FDR.PC2.K2 = length(which(candidates.PC2.K2 <= 900))/length(candidates.PC2.K2)
estimated.TP.PC2.K2 = length(which(candidates.PC2.K2 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.PC2.K2, "True Positive:", estimated.TP.PC2.K2))



# K=3
zs.PC2.K3 = z.scores(lfmm.PC2, K = 3, d=1)
zs.PC2.K3.median = apply(zs.PC2.K3, MARGIN = 1, median)
lambda.PC2.K3 = median(zs.PC2.K3.median^2)/qchisq(0.5, df = 1)
lambda.PC2.K3
adj.p.values.PC2.K3 = pchisq(zs.PC2.K3.median^2/lambda.PC2.K3, df = 1, lower = FALSE)
par(mfrow = c(2,1))
hist(adj.p.values.PC2.K3, col = "lightblue")
plot(-log10(adj.p.values.PC2.K3), pch = 19, col = "blue", cex = .7)
w.PC2.K3 = which(sort(adj.p.values.PC2.K3) < q * (1:L)/L)
candidates.PC2.K3 = order(adj.p.values.PC2.K3)[w.PC2.K3]
estimated.FDR.PC2.K3 = length(which(candidates.PC2.K3 <= 900))/length(candidates.PC2.K3)
estimated.TP.PC2.K3 = length(which(candidates.PC2.K3 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.PC2.K3, "True Positive:", estimated.TP.PC2.K3))

LFMM.PC2.candidates <- Reduce(intersect, list(candidates.PC2.K1,candidates.PC2.K2,candidates.PC2.K3))
LFMM.PC2.candidates

# PC3
zs.PC3.K1 = z.scores(lfmm.PC3, K = 1)
zs.PC3.K2 = z.scores(lfmm.PC3, K = 2)
zs.PC3.K3 = z.scores(lfmm.PC3, K = 3)
zs.PC3.K1.median = apply(zs.PC3.K1, MARGIN = 1, median)

# 3 # compute genomic inflation factor
lambda.PC3.K1 = median(zs.PC3.K1.median^2)/qchisq(0.5, df = 1)
lambda.PC3.K1
adj.p.values.PC3.K1 = pchisq(zs.PC3.K1.median^2/lambda.PC3.K1, df = 1, lower = FALSE)
par(mfrow = c(2,1))

# histogram should be flat with a peak near zero 
hist(adj.p.values.PC3.K1, col = "lightblue")
plot(-log10(adj.p.values.PC3.K1), pch = 19, col = "blue", cex = .7)

# Benjamini-Hochberg to adjust p-values for multiple testing
L = 15038
q = 0.01
w.PC3.K1 = which(sort(adj.p.values.PC3.K1) < q * (1:L)/L)
candidates.PC3.K1 = order(adj.p.values.PC3.K1)[w.PC3.K1]
candidates.PC3.K1

estimated.FDR.PC3.K1 = length(which(candidates.PC3.K1 <= 900))/length(candidates.PC3.K1)
estimated.TP.PC3.K1 = length(which(candidates.PC3.K1 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.PC3.K1, "True Positive:", estimated.TP.PC3.K1))


# run with K=2
zs.PC3.K2 = z.scores(lfmm.PC3, K = 2, d=1)
zs.PC3.K2.median = apply(zs.PC3.K2, MARGIN = 1, median)
lambda.PC3.K2 = median(zs.PC3.K2.median^2)/qchisq(0.5, df = 1)
lambda.PC3.K2
adj.p.values.PC3.K2 = pchisq(zs.PC3.K2.median^2/lambda.PC3.K2, df = 1, lower = FALSE)
par(mfrow = c(2,1))
hist(adj.p.values.PC3.K2, col = "lightblue")
plot(-log10(adj.p.values.PC3.K2), pch = 19, col = "blue", cex = .7)
w.PC3.K2 = which(sort(adj.p.values.PC3.K2) < q * (1:L)/L)
candidates.PC3.K2 = order(adj.p.values.PC3.K2)[w.PC3.K2]
estimated.FDR.PC3.K2 = length(which(candidates.PC3.K2 <= 900))/length(candidates.PC3.K2)
estimated.TP.PC3.K2 = length(which(candidates.PC3.K2 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.PC3.K2, "True Positive:", estimated.TP.PC3.K2))



# K=3
zs.PC3.K3.median = apply(zs.PC3.K3, MARGIN = 1, median)
lambda.PC3.K3 = median(zs.PC3.K3.median^2)/qchisq(0.5, df = 1)
lambda.PC3.K3
adj.p.values.PC3.K3 = pchisq(zs.PC3.K3.median^2/lambda.PC3.K3, df = 1, lower = FALSE)
par(mfrow = c(2,1))
hist(adj.p.values.PC3.K3, col = "lightblue")
plot(-log10(adj.p.values.PC3.K3), pch = 19, col = "blue", cex = .7)
w.PC2.K3 = which(sort(adj.p.values.PC2.K3) < q * (1:L)/L)
candidates.PC3.K3 = order(adj.p.values.PC2.K3)[w.PC2.K3]
estimated.FDR.PC3.K3 = length(which(candidates.PC3.K3 <= 900))/length(candidates.PC3.K3)
estimated.TP.PC3.K3 = length(which(candidates.PC3.K3 > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR.PC3.K3, "True Positive:", estimated.TP.PC3.K3))


LFMM.PC3.candidates <- Reduce(intersect, list(candidates.PC3.K1,candidates.PC3.K2,candidates.PC3.K3))
LFMM.PC3.candidates


# export to use in manhattan plot
LFMMadjp.PC2 <- cbind(-log10(adj.p.values.PC2.K1),-log10(adj.p.values.PC2.K2),-log10(adj.p.values.PC2.K3))
head(LFMMadjp.PC2)
write.csv(LFMMadjp.PC2,file="LFMMpvalues.PC2.csv")
LFMMadjp.PC3 <- cbind(-log10(adj.p.values.PC3.K1),-log10(adj.p.values.PC3.K2),-log10(adj.p.values.PC3.K3))
head(LFMMadjp.PC3)
write.csv(LFMMadjp.PC3,file="LFMMpvalues.PC3.csv")



# how many candidate SNPs
length(LFMM.PC1.candidates)
length(LFMM.PC2.candidates)
length(LFMM.PC3.candidates)

LFMM.all.candidates <- Reduce(intersect, list(LFMM.PC1.candidates,LFMM.PC2.candidates,LFMM.PC3.candidates))
length(LFMM.all.candidates)


# pull out all adjusted p-values
foo <- cbind(adj.p.values.K1,adj.p.values.K2,adj.p.values.K3,
             adj.p.values.PC2.K1,adj.p.values.PC2.K2,adj.p.values.PC2.K3,
             adj.p.values.PC3.K1,adj.p.values.PC3.K2,adj.p.values.PC3.K3)
SNP <- 1:15038
LFMM.adjpvalues <- cbind(SNP,foo)
LFMM.adjpvalues <- as.data.frame(LFMM.adjpvalues)
# pull out p-values for these candidates
LFMM.candidates <- subset(LFMM.adjpvalues, SNP %in% LFMM.all.candidates)

# apply -log10
LFMM.candidates[,2:9] <- lapply(LFMM.candidates[,2:9], function(x) -log10(x))
head(LFMM.candidates)

# keep candidate only if it's > 5 sd on any significance measure
mean(adj.p.values.K1) + (5*sd(adj.p.values.K1))
mean(adj.p.values.K2) + (5*sd(adj.p.values.K2))
mean(adj.p.values.K3) + (5*sd(adj.p.values.K3))

mean(adj.p.values.PC2.K1) + (5*sd(adj.p.values.PC2.K1))
mean(adj.p.values.PC2.K2) + (5*sd(adj.p.values.PC2.K2))
mean(adj.p.values.PC2.K3) + (5*sd(adj.p.values.PC2.K3))

mean(adj.p.values.PC3.K1) + (5*sd(adj.p.values.PC3.K1))
mean(adj.p.values.PC3.K2) + (5*sd(adj.p.values.PC3.K2))
mean(adj.p.values.PC3.K3) + (5*sd(adj.p.values.PC3.K3))

LFMM.candidates.5sd <- LFMM.candidates[which((LFMM.candidates$adj.p.values.K1 > 2.38 &
                                                LFMM.candidates$adj.p.values.K2 > 2.38 &
                                                LFMM.candidates$adj.p.values.K3 > 2.38 |
                                                LFMM.candidates$adj.p.values.PC2.K1 > 2.38 &
                                                LFMM.candidates$adj.p.values.PC2.K2 > 2.38 &
                                                LFMM.candidates$adj.p.values.PC2.K3 > 2.38 |
                                                LFMM.candidates$adj.p.values.PC3.K1 > 2.38 &
                                                LFMM.candidates$adj.p.values.PC3.K2 > 2.38 &
                                                LFMM.candidates$adj.p.values.PC3.K3 > 2.38 )),]
length(LFMM.candidates.5sd$SNP) # 1315

positionAF <- read.csv("PositionAF.csv",sep=",")
head(positionAF)

LFMM.candidates.5sd.pos <- merge(LFMM.candidates.5sd, positionAF, 
      by.x = "SNP", by.y = "SNP_ID")

head(LFMM.candidates.5sd.pos)
write.csv(LFMM.candidates.5sd.pos,"LFMM.candidates.csv")




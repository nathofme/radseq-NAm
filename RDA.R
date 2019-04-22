### RDA based on Brenna Forester's work
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html

library(vegan)
library(rda)
library(robust)
library(psych)

# load in genotypic response and environmental predictors
# gen is vcftools 012 output, combined .012 and .012.indv files
gen <- read.csv("EUSTallSNPr8.012.csv", header=FALSE, na.strings="-1")
sum(is.na(gen))
sum(!is.na(gen))

str(gen)

gen.nomiss <- read.csv("EUSTallSNPr8.nomiss.012.csv", header=FALSE, na.strings="-1")
sum(is.na(gen.nomiss))
str(gen.nomiss)

gen[1:20,1:20]

#calculate probabilities of the three genotypes for each column (SNP)
n <- nrow(gen)
p_0 <- apply(gen, 2, function(x){sum(x == 0, na.rm = T)/(n - sum(is.na(x)))})
p_1 <- apply(gen, 2, function(x){sum(x == 1, na.rm = T)/(n - sum(is.na(x)))})
p_2 <- apply(gen, 2, function(x){sum(x == 2, na.rm = T)/(n - sum(is.na(x)))})
p <- data.frame(p_0, p_1, p_2)
# make sure probabilities add up to 1
p <- t(apply(p, 1, function(x){x/sum(x)}))
#make a table for indices of missing genotypes
NA_indices <- which(is.na(gen), arr.ind = T)
#replace missing genotypes by sampling from (0,1,2) based on probabilities given in table p
for (i in 1:nrow(NA_indices)) {
  x <- NA_indices[i, ]
  data[x[1], x[2]] <- sample(c(0:2), 1, replace = T, prob = p[x[2], ])
}


str(gen)
write.csv(gen, "EUSTallSNPr8.012.imp.csv", row.names = F, quote = F)
gen.imp <- read.csv("EUSTallSNPr8.012.imp.csv",sep=",")
str(gen.imp)

#gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
#sum(is.na(gen.imp))
#gen.imp[1:30,1:10]
#gen.imp <- as.data.frame(gen.imp)
#gen.imp <- lapply(gen.imp, as.numeric(unlist(gen.imp[2:159,2:15038])))
#str(gen.imp)

# impute genotypes 
library(radiator)
library(grur)
library(ranger)
library(missRanger)
#EUSTallSNPr8.imp.iv <- tidy_genomic_data(data = "EUSTallSNPr8.vcf", vcf.metadata = TRUE)
#EUSTallSNPr8.imp.g <- grur_imputations(data = EUSTallSNPr8.imp.i, imputation.method = "max")
#str(EUSTallSNPr8.imp.i)
#tidy_genomic_data(EUSTallSNPr8.imp.g)
#genomic_converter(EUSTallSNPr8.imp)

#gen$Individual <- as.character(gen$Individual)
gen.imp[,1:15038] <- lapply(gen.imp[,1:15038], as.numeric)
gen.imp <- as.data.frame(gen.imp[2:159,])
str(gen.imp)
colnames(gen.imp) <- gsub("[a-zA-Z ]", "", colnames(gen.imp))
head(gen.imp)

# checking with no missing data
gen.nomiss[,1:2499] <- lapply(gen.nomiss[,1:2499], as.numeric)
gen.nomiss <- as.data.frame(gen.nomiss[2:159,])
colnames(gen.nomiss) <- gsub("[a-zA-Z ]", "", colnames(gen.nomiss))
head(gen.nomiss)
str(gen.nomiss)
sum(!is.na(gen.nomiss))

n.no <- nrow(gen.nomiss)
p.no_0 <- apply(gen.nomiss, 2, function(x){sum(x == 0, na.rm = T)/(n - sum(is.na(x)))})
p.no_1 <- apply(gen.nomiss, 2, function(x){sum(x == 1, na.rm = T)/(n - sum(is.na(x)))})
p.no_2 <- apply(gen.nomiss, 2, function(x){sum(x == 2, na.rm = T)/(n - sum(is.na(x)))})
p.no <- data.frame(p.no_0, p.no_1, p.no_2)
# make sure probabilities add up to 1
p.no <- t(apply(p.no, 1, function(x){x/sum(x)}))
#make a table for indices of missing genotypes
NA_indices <- which(is.na(gen.nomiss), arr.ind = T)
#replace missing genotypes by sampling from (0,1,2) based on probabilities given in table p
for (i in 1:nrow(NA_indices)) {
  x <- NA_indices[i, ]
  data[x[1], x[2]] <- sample(c(0:2), 1, replace = T, prob = p.no[x[2], ])
}


# load in environmental data
env <- read.csv("EUSTrad.env.lessvar.csv")
env.df <- as.data.frame(env)
str(env.df)
pred <- env.df[,2:7]
pred <- as.data.frame(lapply(pred, as.numeric))
# remove BIO4 since not running
pred <- pred[c(1,3,5,6)]
str(pred)
quartz()
pdf("FigureS2_EnvPairs.pdf",height=6,width=6)
pairs.panels(pred)
dev.off()


env.plot <- read.csv("EnvVarSites.csv",sep=",")
str(env.plot)

library(reshape2)
env.plot <- melt(env.plot,id="location")
library(scales)
scale.bio1 <- rescale(env.plot$BIO1, from=c(0, max(env.plot$BIO1)))
scale.bio7 <- rescale(env.plot$BIO7, from=c(0, max(env.plot$BIO7)))
scale.bio12 <- rescale(env.plot$BIO12, from=c(0, max(env.plot$BIO12)))
scale.bio16 <- rescale(env.plot$BIO16, from=c(0, max(env.plot$BIO16)))
scale.elev <- rescale(env.plot$elevation, from=c(0, max(env.plot$elevation)))

env.plot <- cbind(env.plot, scale.bio1, scale.bio12, scale.bio16, scale.bio7, scale.elev)

quartz()
pdf("Figure3.elev.pdf",width=4,height=3)
ggplot(env.plot, aes(x=location,y=scale.elev)) + 
  geom_bar(stat="identity",position="dodge",fill="black") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 90, hjust = 1, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

pdf("Figure3.BIO1.pdf",width=4,height=3)
ggplot(env.plot, aes(x=location,y=scale.bio1)) + 
  geom_bar(stat="identity",position="dodge",fill="#b2182b") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 1, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

pdf("Figure3.BIO7.pdf",width=4,height=3)
ggplot(env.plot, aes(x=location,y=scale.bio7)) + 
  geom_bar(stat="identity",position="dodge",fill="#d6604d") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 1, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

pdf("Figure3.BIO12.pdf",width=4,height=3)
ggplot(env.plot, aes(x=location,y=scale.bio12)) + 
  geom_bar(stat="identity",position="dodge",fill="#2166ac") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 1, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

pdf("Figure3.BIO16.pdf",width=4,height=3)
ggplot(env.plot, aes(x=location,y=scale.bio16)) + 
  geom_bar(stat="identity",position="dodge",fill="#92c5de") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 1, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()


#ggplot(env.plot, aes(x=location,y=scale.value,fill=variable)) + 
  geom_bar(stat="identity",position="stack") + theme_bw() +
  scale_fill_manual(values=c("#b2182b","#d6604d","#2166ac","#92c5de","black"))
dev.off()


# randomize predictor variables to compare fit
pred.shuffle <- pred 
pred.shuffle <- pred.shuffle[sample(nrow(pred.shuffle)),]
str(pred.shuffle)

# Step 1: ordinate genotypes
library(vegan)
eustrad.rda <- rda(gen.imp ~ BIO1 + BIO7 + BIO16 + elevation, data=pred, scale=T)
eustrad.rda
RsquareAdj(eustrad.rda)
# pretty uniform
screeplot(eustrad.rda,npcs=8)
library(vegan)

# test for significance
signif.full <- anova.cca(eustrad.rda, parallel=getOption("mc.cores"))
signif.full
# significance of each axis
signif.axis <- anova.cca(eustrad.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
# check variance inflation < 10
vif.cca(eustrad.rda)

# try with shuffled data
eustrad.rda.shuffle <- rda(gen.imp ~ BIO1 + BIO7 + BIO12 + BIO16 + elevation, data=pred.shuffle, scale=T)
eustrad.rda.shuffle
RsquareAdj(eustrad.rda.shuffle)
signif.full.shuffle <- anova.cca(eustrad.rda.shuffle, parallel=getOption("mc.cores"))
signif.full.shuffle

# check with no missing data 
eustrad.nomiss.rda <- rda(gen.nomiss ~ BIO1 + BIO7 + BIO16 + elevation, data=pred, scale=T)


# plot
plot(eustrad.rda, scaling=3)

eco <- env$location
eco
# color scheme broken down by regions
rdbu.plus.palette=colorRampPalette(c('#b2182b','#f7f7f7','#4393ac'))
bg <- c("#67001f","#b2182b","white","lightgoldenrodyellow","#4393ac","gray90","lightyellow3","#d6604d","gray60","gold3","gray30","#d6604d","#92c5de","black","gold","#2166ac","gray10")

### Plot individuals, all SNPs, and predictor vectors
quartz()
pdf("RDA_LocalAdaptation_novec.pdf")
plot(eustrad.rda, type="n",scaling=3)
points(eustrad.rda, display="species", pch=20, cex=0.5, col="gray60", scaling=3)           # the SNPs
points(eustrad.rda, display="sites", pch=21, cex=1, col="black", 
       scaling=3, bg=bg[eco]) # the individuals
text(eustrad.rda, scaling=3, display="bp", col="black", cex=1)                           # the predictors
legend("bottomleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

load.rda <- scores(eustrad.rda, choices=c(1:3), display="species")
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# check for normal distribution
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

# get number of candidates > 3 SDs from mean loading
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 
cand3 <- outliers(load.rda[,3],3)
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand # 191

# organize results into data frame with SNP, loading, correlation
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
head(cand)
cand$snp <- gsub("[a-zA-Z ]", "", cand$snp)
cand$snp <- as.numeric(cand$snp)
str(cand)
length(cand$snp[duplicated(cand$snp)]) # no duplicates across predictors

write.csv(cand,"RDAcand.3sd.csv")

# add correlations with predictors
foo <- matrix(nrow=(ncand), ncol=5)  # columns for # of predictors
colnames(foo) <- c("BIO1","BIO7","BIO12","BIO16","elevation")
head(foo)

str(cand)
str(gen.imp)
str(pred)

# here, want a correlation between each SNP and each predictor, stored in matrix foo
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

# SNPs detected on multiple axes:
cand$snp[duplicated(cand$snp)] # 0 duplicates

# Let's see which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,9] <- names(which.max(abs(bar[4:8]))) # gives the variable
  cand[i,10] <- max(abs(bar[4:8]))              # gives the correlation
}

colnames(cand)[9] <- "predictor"
colnames(cand)[10] <- "correlation"
cand$correlation <- abs(cand$correlation)

write.csv(cand,"RDAcand.3sd.csv", row.names=F)

table(cand$predictor) 
table(cand$axis)      


### pull out candidates that load > 5 sd from mean and correlation is > 2 sd from mean
mean(cand$correlation) # 0.126
sd(cand$correlation) # 0.060
# 3 sd from mean is 0.24
# 2 sd from mean is 0.18
mean(cand$loading) # 0.118
sd(cand$loading) # 0.048
# 4 sd from mean 0.31
# 3 sd from mean 0.262
# 2 sd from mean 0.214

# strongest correlations
cand.correl.load <- cand[which((cand$loading > 0.214) | (cand$correlation > 0.12)),]
length(cand.correl.load$snp)
cand.correl.load


cand.load <- cand[which((cand$loading > 0.214) | (cand$correlation > 0.06)),]
length(cand.load$snp)
cand.load

# add in AF, positions and scaffolds
positionAF <- read.csv("PositionAF.csv",sep=",")
head(positionAF)

RDA.cand.pos <- merge(cand, positionAF, by.x = "snp", by.y = "SNP_ID")
RDA.cand.pos
write.csv(cand.pos,"RDA.candidates.csv")

# shared in both LFMM and RDA
RDA.LFMM.cand <- Reduce(intersect, list(LFMM.candidates$SNP,RDA.cand.pos$snp))
length(RDA.LFMM.cand) # 14
unique(RDA.LFMM.cand)


### plot AF differences among candidates & genome-wide






### get FSTs for each candidate
# using fst measures from bayescan, saved in manhattan.csv 
manhattan <- read.table("manhattan.csv", header=TRUE, sep=",")
manhattansubset <- manhattan[complete.cases(manhattan),]
manhattan <- data.frame(manhattansubset)
head(manhattan)


sel <- cand$snp
env <- cand$predictor
env[env=="BIO1"] <- '#b2182b'
env[env=="BIO12"] <- '#4393c3'
env[env=="BIO16"] <- '#d1e5f0'
env[env=="BIO7"] <- '#f4a582'
env[env=="elevation"] <- '#f7f7f7'

# color by predictor:
rownames(eustrad.rda$CCA$v) <- gsub("[a-zA-Z ]", "", rownames(eustrad.rda$CCA$v))
rownames(eustrad.rda$CCA$v) <- as.numeric(rownames(eustrad.rda$CCA$v))
col.pred <- rownames(eustrad.rda$CCA$v)

for (i in 1:length(sel)) {
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("chr",col.pred)] <- '#f1eef6'
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a')

# axes 1 & 2
plot(eustrad.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(eustrad.rda, display="species", pch=20, cex=1, col="gray60", scaling=3)
points(eustrad.rda, display="species", pch=20, cex=0.5, col="black", scaling=3)
text(eustrad.rda, scaling = 3, display = "bp", col="black", cex=1.1)
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a')
legend("bottomright", legend = c("BIO1","BIO7","BIO12","BIO16","Elevation"), bty = "n", col="gray32", pch = 21, cex=1.3, pt.bg = bg)


### plot correlations
quartz()

snpcorrel <- cand$correlation
sel <- cand$snp
env <- cand$predictor
correlmat <- as.data.frame(cbind(sel,env,snpcorrel))
str(correlmat)

cand.forplot <- cand[order(cand$correlation,decreasing=TRUE),]

quartz()
pdf("RDA_SNPcorrelations_4.pdf")
barplot(cand.forplot$correlation,horiz=TRUE,col=env,width=0.5,xlab="Correlation",xlim=c(0,0.3))
legend("topright",legend=c("BIO1 (39)","BIO7 (42)","BIO12 (41)","BIO16 (33)","elevation (36)"), bty = "n", col="black", pch = 22, cex=1.3, pt.bg = c("#b2182b","#f4a582","#4393c3","#d1e5f0","#f7f7f7"))
dev.off()

env <- cand$predictor
env[env=="BIO1"] <- '#b2182b'
env[env=="BIO12"] <- '#4393c3'
env[env=="BIO16"] <- '#d1e5f0'
env[env=="BIO7"] <- '#f4a582'
env[env=="elevation"] <- '#f7f7f7'

# individuals' weighted average for each RDA
eustrad.rda$CCA$wa
eustrad.rda


for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[1:191]))) # gives the variable
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

# Capblancq method continuing (doesn't show individuals)
library(robust)
library(qvalue)
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

res_rdadapt<-rdadapt(eustrad.rda, 5)
which(res_rdadapt[,2] < 0.1)
# only 8 candidates with q < 0.1

quartz()
library(ggplot2)
pdf("RDA_OutlierSNPs.pdf",width=5,height=5)
ggplot() + 
  geom_point(aes(x=eustrad.rda$CCA$v[,1], y=eustrad.rda$CCA$v[,2]), col = "gray80") +
  geom_point(aes(x=eustrad.rda$CCA$v[which(res_rdadapt[,2] < 0.1),1],
                 y=eustrad.rda$CCA$v[which(res_rdadapt[,2] < 0.1),2]), col = "#b2182b") +
  geom_segment(aes(xend=eustrad.rda$CCA$biplot[,1]/10, yend=eustrad.rda$CCA$biplot[,2]/10, x=0, y=0),
               colour="black", size=0.5, linetype=1,
               arrow=arrow(length = unit(0.02, "npc"))) +
  xlab("RDA 1") + ylab("RDA 2") + 
  theme_classic() + theme(legend.position="none")
dev.off()

geom_text(aes(x=1.2*eustrad.rda$CCA$biplot[,1]/10, y=1.2*eustrad.rda$CCA$biplot[,2]/10,
              label = colnames(env[,2:6]))) +

  

plot(eustrad.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(eustrad.rda, display="species", pch=20, cex=1, col="gray60", scaling=3)
points(eustrad.rda, display="species", pch=20, cex=0.5, col="black", scaling=3)
text(eustrad.rda, scaling = 3, display = "bp", col="black", cex=1.1)
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a')
legend("bottomright", legend = c("BIO1","BIO7","BIO12","BIO16","Elevation"), bty = "n", col="gray32", pch = 21, cex=1.3, pt.bg = bg)

  
  
# Step 3: adjust p-values for FDR
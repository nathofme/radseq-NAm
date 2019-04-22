library(RColorBrewer)
library(SNPRelate)
library(gdsfmt)
library(scales)
library(ggplot2)

setwd("~/Documents/starlingRAD/analysis")

### TEST FOR MISSING DATA ###
# convert SNPs to binary format (gdsfmt or Genomic Data Structure data files) - accelerates computation
snpgdsVCF2GDS(vcf.fn="/Users/nataliehofmeister/Documents/starlingRAD/analysis/EUSTallSNPr8.nomiss.vcf.recode.vcf", out.fn="eustrad.nomiss.gds",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)
snpgdsSummary("eustrad.nomiss.gds")
genofile <- snpgdsOpen("eustrad.nomiss.gds")
# 6287 SNPs

# get number of loci for one SNP per locus
snpgdsVCF2GDS(vcf.fn="/Users/nataliehofmeister/Documents/starlingRAD/analysis/EUSToneSNPr8maf01p.vcf", out.fn="eust.gds.one",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)
snpgdsSummary("eust.gds.one")
# 3568 SNPs

snpgdsVCF2GDS(vcf.fn="/Users/nataliehofmeister/Documents/starlingRAD/analysis/EUSTallSNPr8_rarer.vcf", out.fn="eustrad.rarer.gds",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)
snpgdsSummary("eustrad.rarer.gds")
genofile.rare <- snpgdsOpen("eustrad.rarer.gds")

# how much data is missing? (fraction for each sample) - median 32%, range 18-56% 
miss <- snpgdsSampMissRate(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss
summary(miss)
#remove EUST_TX0202 since >50%

###################################### PCAs ####################################

# PCA with SNPRelate
pca <- snpgdsPCA(gdsobj = genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
pcatab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(pcatab)
pcatab
plot(pcatab$EV2, pcatab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

# add labels by population
# order in vcf (check cohort file)
# label particular populations

samplingsite <- c("NM","NM","NM","NM","NM","NM","NM","KS","KS","KS","KS","KS","KS","KS","KS","KS","MO","MO","MO","MO","MO","MO","MO","MO","MO","IL","IL","IL","IL","IL","IL","IL","NC","NC","NC","NC","NC","NC","NC","NC","NC","NC","NC","AZ","AZ","AZ","AZ","AZ","AZ","AZ","AZ","AZ","AZ","IA","IA","IA","IA","IA","IA","IA","IA","IA","ID","ID","ID","ID","ID","ID","ID","ID","NV","NV","NV","NV","NV","NV","NV","NV","NV","NV","WA","WA","WA","WA","WA","WA","WA","WA","WA","WA","CA","CA","CA","CA","CA","CA","CA","CA","CA","CA","CA","CO","CO","CO","CO","CO","CO","CO","CO","TX","TX","TX","TX","TX","TX","TX","TX","TX","NH","NH","NH","NH","NH","NH","NH","NH","NH","NH","NH","NY","NY","NY","NY","NY","NY","NY","NY","NY","NY","NE","NE","NE","NE","NE","NE","NE","NE","NE","NE","NE","WI","WI","WI","WI","WI","WI","WI","WI")
# alphabetical # population <- c("AZ","AZ","AZ","AZ","AZ","AZ","AZ","AZ","AZ","AZ","CA","CA","CA","CA","CA","CA","CA","CA","CA","CA","CA","CO","CO","CO","CO","CO","CO","CO","CO","IA","IA","IA","IA","IA","IA","IA","IA","IA","ID","ID","ID","ID","ID","ID","ID","ID","IL","IL","IL","IL","IL","IL","IL","KS","KS","KS","KS","KS","KS","KS","KS","KS","MO","MO","MO","MO","MO","MO","MO","MO","MO","NC","NC","NC","NC","NC","NC","NC","NC","NC","NC","NC","NE","NE","NE","NE","NE","NE","NE","NE","NE","NE","NE","NH","NH","NH","NH","NH","NH","NH","NH","NH","NH","NH","NM","NM","NM","NM","NM","NM","NM","NV","NV","NV","NV","NV","NV","NV","NV","NV","NV","NY","NY","NY","NY","NY","NY","NY","NY","NY","NY","TX","TX","TX","TX","TX","TX","TX","TX","TX","WA","WA","WA","WA","WA","WA","WA","WA","WA","WA","WI","WI","WI","WI","WI","WI","WI","WI")

pcatab2 <- cbind(pcatab,samplingsite)
pcatab2
write.csv(pcatab2)
# then copy .csv to analysis folder

### PCA plot with all populations 
colorlist<-colorRampPalette(colors)
colorlist<-viridis(17)
quartz()
pdf("PCAallpop.pdf")
plot(pcatab2$EV1, pcatab2$EV2, xlab="PC1 1.07%", ylab="PC2 1.03%",col="black",bg=colorlist,pch=21,cex=0.8)
legend("topright",legend=levels(pcatab2$samplingsite),pch=21,col=c("black"),pt.bg=colorlist,cex=0.8)
dev.off()

### PCA plot with regions

population<-c("west","west","west","west","west","west","west","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","west","west","west","west","west","west","west","west","west","west","east","east","east","east","east","east","east","east","east","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east")
region <- c("Southwest","Southwest","Southwest","Southwest","Southwest","Southwest","Southwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Southeast","Southeast","Southeast","Southeast","Southeast","Southeast","Southeast","Southeast","Southeast","Southeast","Southeast","Southwest","Southwest","Southwest","Southwest","Southwest","Southwest","Southwest","Southwest","Southwest","Southwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Northwest","Northwest","Northwest","Northwest","Northwest","Northwest","Northwest","Northwest","West","West","West","West","West","West","West","West","West","West","Northwest","Northwest","Northwest","Northwest","Northwest","Northwest","Northwest","Northwest","Northwest","Northwest","West","West","West","West","West","West","West","West","West","West","West","Mountain","Mountain","Mountain","Mountain","Mountain","Mountain","Mountain","Mountain","South","South","South","South","South","South","South","South","South","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Northeast","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest","Midwest")
pcatab3 <- cbind(pcatab,region,realpop)
pcatab3

# just east and west divisions
colors<-c("#4393c3","#f7f7f7")

quartz()
pdf("PCArealpop.pdf")
plot(pcatab3$EV1, pcatab3$EV2, xlab="PC1 1.07%", ylab="PC2 1.03%",col="black",bg=colors,pch=21,cex=2)
legend("topright",legend=levels(pcatab3$realpop),pch=21,col=c("black"),pt.bg=colors,cex=1)
dev.off()


###################################### ADEGENET ####################################
### includes dudi.pca, AMOVA, Mantel
library(adegenet)
library(ade4)
library(poppr)

# make genind object
genind<-read.structure("EUSToneSNPr8maf01p.stru",n.ind=158,n.loc=3570,onerowperind=FALSE,col.lab=1,col.pop=2,row.marknames=0,NA.char="-9",ask=FALSE,quiet=FALSE)
genind
# check number of loci if not reading properly
sum(!sapply(read.table("EUSToneSNPr8maf01p.stru", sep = "\t"), is.logical))
# 3570

genind.all<-read.structure("EUSTallSNPr8.stru",n.ind=158,n.loc=15040,onerowperind=FALSE,col.lab=1,col.pop=2,row.marknames=0,NA.char="-9",ask=FALSE,quiet=FALSE)
genind.all

sum(!sapply(read.table("EUSTallSNPr8.stru", sep = "\t"), is.logical))
# 15040

# modify genind object to update strata for AMOVA
pop(genind) # population info stored in this accessor
# create df to bind to genind with hierarchical levels needed
# vectors already created for PCAs above
individual<-c(as.factor(seq(1,158,1)))
individual
population<-c("west","west","west","west","west","west","west","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","west","west","west","west","west","west","west","west","west","west","east","east","east","east","east","east","east","east","east","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","west","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east","east")
strata_df<-as.data.frame(cbind(individual,samplingsite,population,region))
strata_df
head(strata_df)
strata(genind)<-strata_df
genind
#genclone<-as.genclone(genind)
#genclone

# add lat-long coordinates
coords <- read.csv("latlongxy.csv")
head(coords)
genind@other$xy <- coords
genind@other$xy

# check strata
genind
table(strata(genind,~population))
poppr(genind)

### Ready for analysis!

################################## more pop gen ############################# 

### manhattan plot
library(qqman)

manhattan <- read.table("manhattan.csv", header=TRUE, sep=",")
manhattansubset <- manhattan[complete.cases(manhattan),]
manhattan <- data.frame(manhattansubset)
head(manhattan)

snpsOfInterest <- c(1895,7208,10808,11103)

#for LFMM
snpsOfInterest <- c(1247,1246,14704,12657,10040,5382,10039,3104,5385,3105,4566,4567,4564,11768,4563,11771,11770,3603,14573,14574,3834,4565,3601,3835,7841,14575,1362,11769,10179,1361,3955,7211,11802,11798,11804,3953)

# Bayescan
manhattan(manhattan,chr="CHROM", bp="POS", snp="SNP", p="log10.PO.", logp=TRUE, ylim=c(0,2),ylab="log10p",xlab="Scaffold")


pdf("Manhattan_FST.pdf",height=3)
manhattan(manhattan,chr="CHROM", bp="POS", snp="SNP", p="fst", logp=FALSE, ylab="Weir and Cockerham Fst",xlab="Scaffold",ylim=c(0,0.04),cex=0.5,genomewideline=0.0079,suggestiveline=0.0126,highlight=snpsOfInterest)
dev.off()





#LFMM
library(scales)
quartz()
pdf("LFMM_LocalAdaptation_withlegend.pdf",width=12,height=3)
par(mar = c(5,5,2,5))
with(manhattan, plot(manhattan$SNP, manhattan$LFMM.K3.PC2, pch=16, axes=F, col=alpha("#4393c3", 0.8), cex=0.5,
     ylab="log10p",xlab="SNP", ylim=c(0,50)),xlim=c(0,16000))
axis(side=1,at=c(0,2000,4000,6000,8000,10000,12000,14000,16000))
axis(side=2)
par(new = T)
with(manhattan, plot(manhattan$SNP, manhattan$LFMM.K3.PC3, pch=16, axes=F, xlab=NA, ylab=NA, col=alpha("gray70",0.8), cex=0.5, ylim=c(0,50)))
abline(h=18, col="black")
par(new = T)
with(manhattan, plot(manhattan$SNP, manhattan$LFMM.PC1.K3, pch=16, ylab=NA, xlab=NA, axes=F, cex=0.5, col="black"))
legend("top",legend=c("PC1","PC2","PC3"),pch=16,col=c("black",alpha("#4393c3",0.8),alpha("gray70",0.8)),cex=1,bty="n",horiz=TRUE)
dev.off()

# plot LFMM and allele frequency
quartz()
pdf("SelectionSignif_v_AlleleFreq.pdf",height=4)
par(mar = c(5,5,2,5))
with(manhattan, plot(manhattan$SNP, manhattan$minorAF, type="l",lwd=0.2, axes=F, xlab=NA, ylab=NA, col="#92c5de", ylim=c(0,1), xlim=c(0,16000)))
par(new=T)
with(manhattan, plot(manhattan$SNP, manhattan$LFMM.K3.PC3, pch=16, cex=0.5, col="black", 
                     ylab="log10p",xlab="ScNP",
                     ylim=c(0,50), xlim=c(0,16000)))
abline(h=20, col="black")
par(new = T)
axis(side = 4)
mtext(side = 4, line = 3, 'Allele frequency')
#legend("topright",
#legend=c(expression(-log[10](italic(p))), "FST"),
#lty=c(1,0), pch=c(NA, 16), col=c("gray", "black"))
dev.off()



# plot bayescenv results
quartz()
pdf("Manhattan_Bayescenv.pdf",height=4)
par(mar = c(5,5,2,5))
with(manhattan, plot(manhattan$CHROM, manhattan$PEP_diff_g_alpha, pch=16, cex=0.5, col="black", 
             ylab="Difference in PEP",xlab="Scaffold",
             ylim=c(0,1), xlim=c(0,400)))
abline(h=0.9, col="black")
par(new = T)
with(manhattan, plot(manhattan$CHROM, manhattan$fst_bayescenv, pch=16, axes=F, xlab=NA, ylab=NA, col="#92c5de", cex=0.5, ylim=c(0,1), xlim=c(0,400)))
par(new = T)
with(manhattan, plot(manhattan$CHROM, manhattan$fst, type="l",lwd=2, axes=F, xlab=NA, ylab=NA, col="#2166ac", ylim=c(0,1), xlim=c(0,400)))
axis(side = 4)
mtext(side = 4, line = 3, 'FST')
#legend("topright",
       #legend=c(expression(-log[10](italic(p))), "FST"),
       #lty=c(1,0), pch=c(NA, 16), col=c("gray", "black"))
dev.off()

# compare LFMM and Bayescenv
quartz()
par(mar = c(5,5,2,5))
with(manhattan, plot(manhattan$CHROM, manhattan$LFMM.K3, pch=16, cex=1, col="gray", 
    ylab=expression(-log[10](italic(p))),xlab="Scaffold",
    ylim=c(0,25), xlim=c(0,400)))
par(new = T)
with(manhattan, plot(manhattan$CHROM, manhattan$PEP_diff_g_alpha, pch=16, axes=F, xlab=NA, ylab=NA, col="black", cex=0.3, ylim=c(0,1), xlim=c(0,400)))
axis(side = 4)
mtext(side = 4, line = 3, 'Difference in PEP')


######## summary stats beyond FST
library(mmod)
library(reshape2)
setPop(genind.all) <- ~population # Use ~population analyze by pop
diff_genind.all <- diff_stats(genind.all)
diff_genind.all

per.locus <- melt(diff_genind.all$per.locus, varnames = c("Locus", "Statistic"))
stats <- c("Hs", "Ht", "Gst", "Gprime_st", "D", "D")
glob <- data.frame(Statistic = stats, value = diff_genind.all$global)
head(per.locus)

ggplot(per.locus, aes(x = Statistic, y = value)) +
  geom_boxplot() +
  geom_point() +
  geom_point(size = rel(3), color = "red", data = glob) +
  theme_bw() +
  ggtitle("Estimates of population differentiation")


# observed vs expected heterozygosity
genind.smry <- summary(genind)
genind.smry
pdf("Heterozygosity.pdf")
plot(genind.smry$Hexp, genind.smry$Hobs, xlab="Expected heterozygosity",ylab="Observed heterozygosity",pch=16,cex=0.5)
abline(0, 1, col = "gray")
dev.off()
t.test(genind.smry$Hexp, genind.smry$Hobs, paired = TRUE, var.equal = TRUE)

# another PCA
genind.pca1 <- dudi.pca(genind, cent = FALSE, scale = FALSE, scannf = FALSE, nf=2)
barplot(genind.pca1$eig)
genind.pca1

s.label(genind.pca1$li)
s.kde2d(genind.pca1$li, add.p = TRUE, cpoint = 0)
add.scatter.eig(genind.pca1$eig, 2, 1, 2)

# plotting for figure
library(viridis)
library(factoextra)
dudi.pca.plot <- s.class(genind.pca1$li, fac=genind@strata$population, xax = 1, yax = 2,
     grid=FALSE, axesell=FALSE, col=viridis(17),clabel=1,addaxes=1)
scree.plot <- screeplot(genind.pca1,main="",npcs=25)


### Mantel test for IBD
library(vegan)
dist.geo <- dist(genind$other$xy)
dist.genet <- dist(genind.X)
mtest <- mantel(dist.genet, dist.geo)
mtest

plot(dist.geo,dist.genet,xlab="Spatial distance",ylab="Genetic distance")

mantel.correlog <- mantel.correlog(dist(genind.X),D.geo=dist.geo,r.type="pearson", nperm=999)
mantel.correlog
plot(mantel.correlog)

mantel.correlog.more <- mantel.correlog(dist(genind.X),D.geo=dist.geo,r.type="pearson", nperm=99999)
mantel.correlog.more
plot(mantel.correlog.more)

### AMOVA
# run AMOVA on genind (genclone is haploid...)
# check that "total" is 2N-1 where N is number of individuals
# real pop is just east-west based on EEMS break
amova.pop<-poppr.amova(genind,hier=~population,method="ade",nperm=1000)
amova.pop 
amova.pop.pegas<-poppr.amova(genind,hier=~population,method="pegas",nperm=1000)
amova.pop.pegas
amova.samples<-poppr.amova(genind,hier=~population,within=FALSE,method="ade",nperm=1000)
amova.samples

# test for significance
set.seed(1230)
amova.pop.signif<-randtest(amova.pop,nrepet=999)
amova.pop.signif
pdf("AMOVA.pop.signif.pdf")
plot(amova.pop.signif)
dev.off()

# rand test for variance within individuals
as.randtest(sim=rnorm(1000),obs=amova.pop.signif$obs[1],alter=c("less"))
# rand test for variance among individuals within population
as.randtest(sim=rnorm(1000),obs=amova.pop.signif$obs[2],alter=c("greater"))
# rand test for variance among populations
as.randtest(sim=rnorm(1000),obs=amova.pop.signif$obs[3],alter=c("less"))

# AMOVAs without within-individual variance
set.seed(20151219)
pegas_amova <- pegas::amova(dist.genet ~ population/individual, data = strata_df, nperm = 1000)
pegas_amova
adonis_amova <- adonis(dist.genet ~ population, data = strata_df, permutations = 1000)
adonis_amova

# shuffle population assignments to check if results sensitive
genind.shuffle <- genind
set.seed(9001)
head(strata(genind)[sample(nInd(genind)),-1])
strata(genind.shuffle) <- strata(genind.shuffle)[sample(nInd(genind)), -1]
strata(genind.shuffle2) <- strata(genind.shuffle2)[sample(nInd(genind)), -1]

amova.pop.shuffle<-poppr.amova(genind.shuffle,hier=~population,method="ade")
amova.pop.shuffle
amova.pop.shuffle2<-poppr.amova(genind.shuffle2,hier=~population,method="ade")
amova.pop.shuffle2
amova.pop.shuffle3<-poppr.amova(genind.shuffle3,hier=~population,method="ade")
amova.pop.shuffle3

set.seed(1734)
amova.pop.shuffle.signif<-randtest(amova.pop.shuffle,nrepet=999)
plot(amova.pop.shuffle.signif)
amova.pop.shuffle.signif
set.seed(8932)
amova.pop.shuffle2.signif<-randtest(amova.pop.shuffle2,nrepet=999)
plot(amova.pop.shuffle2.signif)
amova.pop.shuffle.signif
set.seed(3721)
amova.pop.shuffle3.signif<-randtest(amova.pop.shuffle3,nrepet=999)
plot(amova.pop.shuffle3.signif)
amova.pop.shuffle3.signif

# arbitrary regions: Midwest, Southeast, etc.
amova.region<-poppr.amova(genind,hier=~region,method="ade")
amova.region

# with region and population
amova.all<-poppr.amova(genind,hier=~region/population,method="ade")
amova.all

# test for significance
set.seed(1328)
amova.all.signif<-randtest(amova.all,nrepet=999)
plot(amova.all.signif)
amova.all.signif



################# GPhocs
library(HDInterval)
eustrad.mcmc <- read.csv("EUSTallSNPr95.mcmc.log",sep="\t")
head(eustrad.mcmc)

hdi(eustrad.mcmc$tau_root, credMass = 0.95)
hdi(eustrad.mcmc$theta_east, credMass = 0.95)
hdi(eustrad.mcmc$theta_west, credMass = 0.95)
hdi(eustrad.mcmc$theta_root, credMass = 0.95)
hdi(eustrad.mcmc$m_west..east, credMass = 0.95)
hdi(eustrad.mcmc$m_east..west, credMass = 0.95)



############################## conStruct ############################## 
library(conStruct)
data(conStruct.data)

# preparing input files
# allele frequency matrix
# no locus labels
conStruct.freqs <- structure2conStruct(infile = "EUSTallSNPr8maf01p.conStruct.str",
   start.loci = 3,
   onerowperind = TRUE,
   missing.datum = -9,
   outfile = "conStruct.freqs")
# looked at matrix to check that it worked

# geographic distance matrix 
library(fields)
conStruct.coord <- read.csv("latlong.conStruct.csv",stringsAsFactors=FALSE)
conStruct.coord <- as.matrix(conStruct.coord)

conStruct.geodist <- rdist.earth(conStruct.coord,miles=FALSE)
hist(conStruct.geodist)
conStruct.geodist[1:20,1:20]
# because kept several individuals in each pop, cells identical within blocks
# chose to run this way bc diff genetic distances even though same geographic

# make sure geodist is symmetric matrix (diagonal of 0s)
pmean <- function(x,y) (x+y)/2
conStruct.geodist[] <- pmean(conStruct.geodist, matrix(conStruct.geodist, nrow(conStruct.geodist), byrow=TRUE))
conStruct.geodist

library("rstan")
stan(conStruct(spatial = TRUE, 
               K = 3, 
               freqs = conStruct.freqs,
               geoDist = conStruct.geodist, 
               coords = conStruct.coord,
               prefix = "spK3",
               n.chains = 1, 
               n.iter = 3000), control = list(adapt_delta = 0.99))

conStruct(spatial = TRUE, 
          K = 2, 
          freqs = conStruct.freqs,
          geoDist = conStruct.geodist, 
          coords = conStruct.coord,
          prefix = "spK2-2000",
          n.chains = 3, 
          n.iter = 2000)

conStruct(spatial = TRUE, 
          K = 3, 
          freqs = conStruct.freqs,
          geoDist = conStruct.geodist, 
          coords = conStruct.coord,
          prefix = "spK3-5000",
          n.chains = 3, 
          n.iter = 5000)

#match.layers.x.runs(admix.mat1 = ,admix.mat2 = ,admix.mat1.order = NULL)

spK3.5000.conStruct.results <- load("spK3.5000.conStruct.results.Robj")
spK3.5000.data.block <- load("spK3.5000.data.block.Robj",verbose=TRUE)

str(spK3.5000.data.block)

library(maps)
quartz()
pdf("ConStruct_AdmixtureMap.pdf")
#maps::map(xlim = range(conStruct.coord[,1])+c(-4,4), ylim = range(conStruct.coord[,2])+c(-4,4), col="gray")
maps::map(xlim = c(-130,-60), ylim = c(25,50), col="gray")
make.admix.pie.plot(admix.proportions=conStruct.results$chain_3$MAP$admix.proportions, coords=conStruct.coord, add=TRUE)
dev.off()



library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue", suppressUpdates=T)
biocLite("SNPRelate", suppressUpdates=T)
install_github("green-striped-gecko/dartR")
library(dartR)
browseVignettes("dartR")

library(adegenet)
genlight <- read.PLINK()

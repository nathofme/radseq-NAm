#### Candidate Function & Variation
library(plyr)
# need HitTable.csv from BLAST and PantherGeneInfo.csv 
# updated candidates from RDA.R and LFMM.R with BED_NAME
# these files based on common nomenclature of fasta file output from bedtools getfasta

LFMM.candidates.final <- read.csv("LFMM.candidates.final.csv",sep=",") # different from LFMM.R output because includes BED_NAME
head(LFMM.candidates.final)

SNPfst <- read.csv("snp.fst.csv",sep=",")
LFMM.candidates.final <- merge(LFMM.candidates.final, SNPfst, 
                                 by.x = "SNP", by.y = "SNP")

LFMM.hit.table <- read.csv("LFMMCandidates-HitTable.csv",sep=",")
head(LFMM.hit.table)


# merge candidates and hit table to get gene IDs in same data frame
Suppl.Table.2.LFMM <- merge(LFMM.candidates.final, LFMM.hit.table, 
                                 by.x = "BED_NAME", by.y = "Window", all=TRUE)
head(Suppl.Table.2.LFMM)
# add in GO terms for genes
cand.GO <- read.csv("CandidateGeneGOs.csv",sep=",")
head(cand.GO)
write.csv(Suppl.Table.2.LFMM,"SupplementaryTable2_LFMM.csv")


RDA.candidates.final <- read.csv("RDA.candidates.final.csv",sep=",") # includes BED_NAME
head(RDA.candidates.final)
RDA.hit.table <- read.csv("RDACandidates-HitTable.csv",sep=",")
head(RDA.hit.table)
Suppl.Table.2.RDA <- merge(RDA.candidates.final, RDA.hit.table, 
                            by.x = "BED_NAME", by.y = "Window", all=TRUE)
write.csv(Suppl.Table.2.RDA,"SupplementaryTable2_RDA.csv")
head(Suppl.Table.2.RDA)

# check AF and FST of candidates
mean(LFMM.candidates.final$fst) # essentially the same as mean genome-wide FST
summary(LFMM.candidates.final$AF)
summary(positionAF$AF)

t.test(positionAF$AF,LFMM.candidates.final$AF,alternative="less") # decidedly not different



# see if certain GOs overrepresented
cand.GO.acc <- cand.GO[,1:2]
head(cand.GO.acc)
cand.geneID.standard <- read.csv("CandidateGeneIDSignificance.csv",sep=",")
head(cand.geneID.standard)

cand.revigo <- merge(cand.GO.acc, cand.geneID.standard, 
                           by.x = "Gene_name", by.y = "Gene.ID", all=TRUE)
head(cand.revigo)
write.csv(cand.revigo,"Candidates.REVIGO.csv")


# examining GO terms with PANTHER
PANTHER <- read.csv("PantherAnalysis.csv",sep=",")
head(PANTHER)

PANTHERuncorrected <- read.csv("PantherAnalysisUncorrected.csv",sep=",")
head(PANTHERuncorrected)
str(PANTHERuncorrected)
colnames(PANTHERuncorrected) <- c("GO_name","category","ref.no","observed",
                               "expected","Fold.Enrichment","over.under","FDR","adj.P.value","raw.P.value")

summary(PANTHER.topcand$adj.P.value)
PANTHER.topcand <- subset(PANTHERuncorrected, adj.P.value < 0.98)
colnames(PANTHER.topcand) <- c("GO_name","FDR","adj.P.value","ref.no","observed",
                                  "expected","Fold.Enrichment","over.under","raw.P.value","FDR2","adjP2")
#write.csv(PANTHER.topcand,"PantherAnalysisTopCand.csv")

str(PANTHER.topcand)

GOslim.acc <- read.csv("GOSlimAccession.csv",sep="\t")
head(GOslim.acc)
PANTHER.topcand.go <- join(PANTHER.topcand, GOslim.acc, 
                           by="GO_name", type="left",match="all")
length(PANTHER.topcand.go$GO_name)
head(PANTHER.topcand)

# plotting best-represented GO categories
quartz()
par(mar = c(5,5,2,5))
ggplot(PANTHER.topcand, aes(x=GO_name,y=Fold.Enrichment)) + 
  geom_bar(x=GO_name,y=Fold.Enrichment,fill=category) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





#PANTHER[,"adj.P.value"] <- p.adjust(PANTHER$raw.P.value, method = "fdr", n = length(PANTHER$raw.P.value))
# get GO accession numbers
#acc.go <- cand.GO[,2:3]
#head(acc.go)
#library(plyr)
#PANTHER.go <- plyr::join(PANTHER, acc.go,by="GO_name", type="left",match="first")
#head(PANTHER.go)
#write.csv(PANTHER.go,"PantherAnalysisGO.csv")


gene.sig <- read.csv("CandidateGeneIDSignificance.csv")
gene.go <- cand.GO[,1:2]
gene.GO.sig <- merge(cand.GO, gene.sig,
                     by.x = "Gene_name", by.y = "Gene.ID", all=TRUE)

str(gene.GO.sig)

gene.GO.sig <- subset(gene.GO.sig,select = -GO_definition)
write.csv(gene.GO.sig,"CandidateGeneSignificance.csv")

PANTHER.sig <- plyr::join(PANTHER, gene.GO.sig,by="GO_accession", type="left",match="first")
unique(PANTHER.sig)
head(PANTHER.sig)

#### Figure 3D
gene.sig.all <- read.csv("Table_in_Figure3_all.csv")
head(gene.sig)

gene.sig$Predictor <- as.factor(gene.sig$Predictor)

#### Testing for AF differences or NS with top candidates

snp.AFNS <- read.csv("snp.geneID.afns.csv",sep=",")
snp.AFNS <- snp.AFNS[,1:4]
head(snp.AFNS)
head(gene.sig)

gene.sig.AFNS <- plyr::join(gene.sig.all, snp.AFNS,by="Gene.ID", type="left",match="first")
head(gene.sig.AFNS)
summary(gene.sig.AFNS$AF)

quartz()
pdf("Figure3_AFcand.pdf",width=4.5,height=2.5)
ggplot(gene.sig.AFNS, aes(x=Significance,y=AF, col=Predictor)) + 
  geom_point(cex=2) +
  scale_color_manual(values=c("black","#92c5de","#b2182b","gray")) + 
  xlim(0.05,0.3) + theme_bw() + ylab("Allele frequency") +
  theme(legend.position = "none", legend.title = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1))
dev.off()


### GO terms in candidates

pdf("Figure3_GOcand.pdf",width=4.5,height=2.5)
ggplot(gene.sig, aes(x=Significance,y=Process, col=Predictor)) + 
  geom_jitter(cex=2,width=0,height=0.05) +
  scale_color_manual(values=c("black","#92c5de","#b2182b","gray")) + 
  xlim(0.12,0.28) + ylab(NULL) + xlab(NULL) + theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1))
dev.off()



library(GOplot)
head(gene.GO.sig)
length(gene.GO.sig$Gene_name)
Category <- rep("BP",3460)
str(Category)

str(EC$david)
GOplot.david <- cbind(Category,gene.GO.sig)
head(GOplot.david)
colnames(GOplot.david) <- c("Category","Genes","ID","adj_pval")
# need to adjust so GO IDs don't repeat and genes separated by commas

circ <- circle_dat(GOplot.david, GOplot.genelist)




# using GO Slim for PANTHER
# Significance column is standardized among RDA & LFMM where higher is more significant (>log10p and >r2)
signal <- filter(PANTHER.sig, grepl("signal",PANTHER$GO_name))
head(signal)
summary(signal$Fold.Enrichment)
summary(signal$adj.P.value)
summary(signal$Significance)



stimul <- filter(PANTHER, grepl("stimul",PANTHER$GO_name))
mean(stimul$Fold.enrichment)




# limma tests for overrepresentation of GO terms
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")
library(limma)
goana



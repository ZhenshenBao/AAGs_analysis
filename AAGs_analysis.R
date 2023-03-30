library(org.Hs.eg.db)
library(clusterProfiler)
library(limma)
library(igraph)

PANCANExp <- `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`
cancernames <- SurvivalTable$cancer.type.abbreviation
cancernames <- cancernames[!duplicated(cancernames)]
samnames <- SurvivalTable$sample
samsp <- strsplit(as.character(samnames),"-")
samname <- do.call(rbind,samsp)
samnames <- cbind(samname,samnames)
samnames <- as.matrix(samnames)

normalsams <- samnames[which(samnames[,4] == "11"),]
normalsams <- normalsams[,5]
normalsams <- gsub("-",".",normalsams)
normalsams <- intersect(normalsams,colnames(PANCANExp))
normalsams <- gsub("\\.","\\-",normalsams)
cancersams <- samnames[which(samnames[,4] != "11"),]
cancersams <- cancersams[,5]

Angiogenic_factors <- as.vector(t(Angiogenic_factors[,1]))
genesy <- bitr(PANCANExp$sample, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
PANCANExp[rownames(genesy),1] <- genesy[,2]

iAngiogenic_factors <- intersect(Angiogenic_factors,PANCANExp$sample)
Angiogenicexps <- PANCANExp[PANCANExp$sample %in% iAngiogenic_factors == TRUE,-1]
rownames(Angiogenicexps) <- iAngiogenic_factors
allEXP <- PANCANExp[!duplicated(PANCANExp$sample),]
rownames(allEXP) <- allEXP$sample
allEXP <- allEXP[,-1]

norAngexps <- Angiogenicexps
norAngexps[,2:ncol(Angiogenicexps)] <- (Angiogenicexps[,2:ncol(Angiogenicexps)]-t(normalavges))/t(normalsds)




######整体差异表达分析#########
#avgsums <- NULL

#Angiogenicexps <- as.data.frame(Angiogenicexps)
cancernames <- c("HNSC", "THCA", "BRCA", "UCEC", "LIHC", "CHOL", "COAD", "READ", "ESCA",
                 "STAD", "BLCA", "KICH", "KIRC", "KIRP", "LUAD", "LUSC")
DEG_AAGs <- NULL
for(j in 1:length(cancernames)){
  cancersurvival <- SurvivalTable[which(SurvivalTable$cancer.type.abbreviation == cancernames[j]), ]
  
  intersams <- intersect(cancersurvival$sample,cancersams)
  internsams <- intersect(cancersurvival$sample,normalsams)
  
  group <- c(rep("N",length(internsams)),rep("C",length(intersams)))
  design <- model.matrix(~0+factor(group))
  colnames(design)=levels(factor(group))
  rownames(design)=c(internsams,intersams)
  
  
  intersams <- gsub("-",".",intersams)
  intersams <- intersect(intersams,colnames(Angiogenicexps))
  intersams1 <- gsub("\\.","\\-",intersams)
  testsams <- c(internsams,intersams1)
  testdesign <- design[testsams,]
  contrast.matrix<-makeContrasts("C-N",levels=testdesign)
  testsams1 <- gsub("-",".",testsams)
  testdata <- allEXP[,testsams1]
  fit <- lmFit(testdata,testdesign)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  nrDEG <- nrDEG[iAngiogenic_factors,]
  nrDEG1 <-nrDEG
  nrDEG1[which(nrDEG1$adj.P.Val <= 0.05 & (nrDEG1$logFC <= -1)),1] <- -1
  nrDEG1[which(nrDEG1$adj.P.Val <= 0.05 & (nrDEG1$logFC >= 1)),1] <- 1
  nrDEG1[which(nrDEG1$adj.P.Val > 0.05 |  nrDEG1$logFC  > -1 & nrDEG1$logFC <  1),1]  <- 0
  nrDEG <- cbind(nrDEG, nrDEG1$logFC)
  #avgsum <- c(avgsum, logfcsum)
  DEG_AAGs <- cbind(DEG_AAGs,nrDEG1$logFC)
  write.csv(nrDEG,paste(cancernames[j],"_nrDEG.csv",sep=""))
}
colnames(DEG_AAGs) <- cancernames
rownames(DEG_AAGs) <- iAngiogenic_factors
write.csv(DEG_AAGs,"DEG_AAGs.csv")


###########T-SNE############
library(Rtsne)
library(ggplot2)
library(ggsci)
TS_Angiogenic_factors <- as.vector(t(P_Angiogenic_factors))
cintersams <- intersect(SurvivalTable$sample,cancersams)
cintersams1 <- gsub("-",".",cintersams)
cintersams1 <- intersect(cintersams1,colnames(allEXP))
cintersams <- gsub("\\.","\\-",cintersams1)
cancerEXP <- allEXP[,cintersams1]
allsurvival <- SurvivalTable[SurvivalTable$sample %in% cintersams == TRUE,]
allsurvivals <- NULL
for(i in 1: length(cancernames)){
  allsurvivals <- rbind(allsurvivals,allsurvival[allsurvival$cancer.type.abbreviation %in% cancernames[i] == TRUE,])
}
allsurvival <- allsurvivals
cintersams <- allsurvival$sample
cintersams1 <- gsub("-",".",cintersams)
cancertypes <- as.vector(allsurvival$cancer.type.abbreviation)
names(cancertypes) <- as.vector(allsurvival$sample)
cancertypes <- cancertypes[cintersams]
cancerEXP <- cancerEXP[,cintersams1]
cancerEXP <- cancerEXP[TS_Angiogenic_factors,]
cancerEXP <- t(cancerEXP)

set.seed(42)
cancerEXP1 <- cancerEXP[,1:21]
cancerEXP1[is.na(cancerEXP1)] <- 0
cancerEXP1 = apply(cancerEXP1,2,as.numeric)
tsne_out = Rtsne(cancerEXP1, dims = 2, pca = T, max_iter = 1000, theta = 0.4, perplexity = 40, verbose = F)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")
cancers <- cancertypes
ggplot(tsne_result,aes(tSNE1,tSNE2,color=cancers))  + scale_color_manual(values = c("#666633", "#CCCCCC", "#FFFF00", "#FF9900", "#FFCC00", "#99FFCC", 
                                                                                    "#33CCFF", "#3399FF", "#3366FF", "#33FFFF", "#FF33FF", "#993399",
                                                                                    "#FF6600", "#FF3333", "#336666", "#FFFFCC"))+ geom_point()+ theme(legend.key = element_blank())+ theme(panel.grid.major=element_line(colour=NA),
                                                                                                                                                                                           panel.background = element_rect(fill = "transparent",colour = NA),
                                                                                                                                                                                           plot.background = element_rect(fill = "transparent",colour = NA),
                                                                                                                                                                                           panel.grid.minor = element_blank())
#######血管生成分数########
cancernames <- c("HNSC", "THCA", "BRCA", "UCEC", "LIHC", "CHOL", "COAD", "READ", "ESCA",
                 "STAD", "BLCA", "KICH", "KIRC", "KIRP", "LUAD", "LUSC")
DEG_AAGs1 <- DEG_AAGs
DEG_AAGs1[is.na(DEG_AAGs1)] <- 0
classdata <- NULL
for(j in 1:length(cancernames)){
  cancersurvival <- SurvivalTable[which(SurvivalTable$cancer.type.abbreviation == cancernames[j]), ]
  intersams <- intersect(cancersurvival$sample,cancersams)
  intersams <- gsub("-",".",intersams)
  intersams <- intersect(colnames(allEXP),intersams)
  internsams <- intersect(cancersurvival$sample,normalsams)
  internsams <- gsub("-",".",internsams)
  internsams <- intersect(colnames(allEXP),internsams)
  cancerdata <- allEXP[,intersams]
  cancerdata <- apply(cancerdata,2,scale)
  rownames(cancerdata) <- rownames(allEXP)
  normaldata <- allEXP[,internsams]
  normaldata <- apply(normaldata,2,scale)
  rownames(normaldata) <- rownames(allEXP)
  diseasecode <- DEG_AAGs1[,cancernames[j]]
  diseasecode[which(diseasecode == 1)] <- 0
  diseasecode[which(diseasecode == -1)] <- 1
  cancerdata <- cancerdata[names(diseasecode),]
  normaldata <- normaldata[names(diseasecode),]
  cancerdata[is.na(cancerdata)] <- 0
  normaldata[is.na(normaldata)] <- 0
  cancerdata <- cancerdata*diseasecode
  normaldata <- normaldata*diseasecode
  cancerscore <- apply(cancerdata,2,sum)
  normalscore <- apply(normaldata,2,sum)
  cancerscore <- abs(cancerscore - mean(normalscore))
  normalscore <- abs(normalscore - mean(normalscore))
  alldata <- c(normalscore,cancerscore)
  types <- c(rep("Normal",length(normalscore)),rep("Tumor",length(cancerscore)))
  cancergroup <- rep(cancernames[j],length(alldata))
  groupdata <- cbind(alldata,types,cancergroup)
  classdata <- rbind(classdata,groupdata)
  write.csv(cancerscore,paste("cancerscore_",cancernames[j],".csv",sep=""))
  write.csv(normalscore,paste("normalscore_",cancernames[j],".csv",sep=""))
  #boxplot(normalscore,cancerscore,names = c("normal","cancer"),col = c(colors()[11],colors()[33]))
}
classdata <- as.data.frame(classdata)
colnames(classdata) <- c("score","group","Tumor")
classdata$Tumor = as.factor(classdata$Tumor)
classdata$group = as.factor(classdata$group)
classdata$score = as.numeric(classdata$score)
rownames(classdata) <- NULL
ggplot(classdata, aes(x=Tumor, y=score, fill=group)) +
  geom_boxplot()+
  scale_y_continuous()+
  scale_fill_manual(values=c("Tumor" = "red", "Normal" = "green"))+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())


###################################################
cancernames <- c("HNSC", "THCA", "BRCA", "UCEC", "LIHC", "CHOL", "COAD", "READ", "ESCA",
                 "STAD", "BLCA", "KICH", "KIRC", "KIRP", "LUAD", "LUSC")
DEG_AAGs1 <- DEG_AAGs
DEG_AAGs1[is.na(DEG_AAGs1)] <- 0
survivalscore <- NULL
for(j in 1:length(cancernames)){
  cancersurvival <- SurvivalTable[which(SurvivalTable$cancer.type.abbreviation == cancernames[9]), ]
  intersams <- intersect(cancersurvival$sample,cancersams)
  intersams <- gsub("-",".",intersams)
  intersams <- intersect(colnames(allEXP),intersams)
  internsams <- intersect(cancersurvival$sample,normalsams)
  internsams <- gsub("-",".",internsams)
  internsams <- intersect(colnames(allEXP),internsams)
  cancerdata <- allEXP[,intersams]
  cancerdata <- apply(cancerdata,2,scale)
  rownames(cancerdata) <- rownames(allEXP)
  normaldata <- allEXP[,internsams]
  normaldata <- apply(normaldata,2,scale)
  rownames(normaldata) <- rownames(allEXP)
  diseasecode <- DEG_AAGs[,cancernames[j]]
  diseasecode[which(diseasecode == 1)] <- 0
  diseasecode[which(diseasecode == -1)] <- 1
  cancerdata <- cancerdata[names(diseasecode),]
  normaldata <- normaldata[names(diseasecode),]
  cancerdata[is.na(cancerdata)] <- 0
  normaldata[is.na(normaldata)] <- 0
  normaldata <- apply(normaldata,1,mean)
  cancerdata <- (cancerdata-normaldata)*diseasecode
  cancerdata[is.na(cancerdata)] <- 0
  cancerdata[which(cancerdata < 0)] <- 0
  #normaldata <- normaldata*diseasecode
  cancerscore <- apply(cancerdata,2,sum)
  #normalscore <- apply(normaldata,2,sum)
  #meanns <- mean(normalscore)
  #cancerscore <- cancerscore-meanns
  
  internsams <- gsub("\\.","-",internsams)
  intersams <- gsub("\\.","-",intersams)
  
  #ns <- cancersurvival[cancersurvival$sample %in% internsams == TRUE,]$PFI.time
  cs <- cancersurvival[cancersurvival$sample %in% intersams == TRUE,]$DSS.time #* cancersurvival[cancersurvival$sample %in% intersams == TRUE,]$DSS
  ossts <- cbind(cancerscore,cs)
  ossts <- ossts[!is.na(ossts[,2]),]
  
  #ggplot(as.data.frame(ossts), aes(x=cs, y=as.numeric(cancerscore),fill=cs)) +
  #geom_boxplot()
  
  ossts <- ossts[which(ossts[,2] != 0),]
  survivalscore <- rbind(survivalscore,ossts)
  #write.csv(ossts,paste("cancersurvival_",cancernames[j],".csv",sep=""))
  
}
#cor(survivalscore[,1],survivalscore[,2])
survivalscore1 <- survivalscore
survivalscore1[which(survivalscore1[,2] <= 365*2), 2] <- 2
#survivalscore1[which(survivalscore1[,2] > 365 & survivalscore1[,2]<= 365*2), 2] <- 2
survivalscore1[which(survivalscore1[,2] > 365*2 & survivalscore1[,2]<= 365*3), 2] <-3
survivalscore1[which(survivalscore1[,2] > 365*3 & survivalscore1[,2]<= 365*4), 2] <-4
survivalscore1[which(survivalscore1[,2] > 365*4 & survivalscore1[,2]<= 365*5),2] <-5
survivalscore1[which(survivalscore1[,2] > 365*5), 2] <- 6

survivalscore1 <- as.data.frame(survivalscore1)
survivalscore1$cancerscore <- as.numeric(survivalscore1$cancerscore)
survivalscore1$cs <- as.numeric(survivalscore1$cs)
ggplot(survivalscore1,aes(x=cancerscore,y=cs))+
  geom_point(size=4,color="steelblue2",alpha=0.5)+
  geom_smooth(color="steelblue2",method="lm",se=FALSE)+
  theme_bw()+theme(plot.title=element_text(hjust=0.5,size=16),axis.title=element_text(size=14))


survivalscore1$cs <- as.factor(survivalscore1$cs)
ggbarplot(survivalscore1,x = "cs",y = "cancerscore",add = "mean_se",color = "cs")



####################################################
cancernames <- c("HNSC", "LIHC", "CHOL", "COAD", "READ", "ESCA",
                 "STAD", "KICH", "KIRC", "KIRP", "LUAD", "LUSC")
DEG_AAGs1 <- DEG_AAGs
DEG_AAGs1[is.na(DEG_AAGs1)] <- 0
Angscore <- NULL
for(j in 1:length(cancernames)){
  cancersurvival <- SurvivalTable[which(SurvivalTable$cancer.type.abbreviation == cancernames[j]), ]
  intersams <- intersect(cancersurvival$sample,cancersams)
  intersams <- gsub("-",".",intersams)
  intersams <- intersect(colnames(allEXP),intersams)
  internsams <- intersect(cancersurvival$sample,normalsams)
  internsams <- gsub("-",".",internsams)
  internsams <- intersect(colnames(allEXP),internsams)
  cancerdata <- allEXP[,intersams]
  cancerdata <- apply(cancerdata,2,scale)
  rownames(cancerdata) <- rownames(allEXP)
  normaldata <- allEXP[,internsams]
  normaldata <- apply(normaldata,2,scale)
  rownames(normaldata) <- rownames(allEXP)
  diseasecode <- DEG_AAGs[,cancernames[j]]
  diseasecode[which(diseasecode == 1)] <- 0
  diseasecode[which(diseasecode == -1)] <- 1
  cancerdata <- cancerdata[names(diseasecode),]
  normaldata <- normaldata[names(diseasecode),]
  cancerdata[is.na(cancerdata)] <- 0
  normaldata[is.na(normaldata)] <- 0
  normaldata <- apply(normaldata,1,mean)
  cancerdata <- abs(cancerdata-normaldata)*diseasecode
  cancerdata[is.na(cancerdata)] <- 0
  cancerdata[which(cancerdata < 0)] <- 0
  cancerscore <- apply(cancerdata,2,sum)
  write.csv(cancerscore,paste("Angscore_",cancernames[j],".csv",sep=""))
  #Angscore <- c(Angscore, cancerscore) 
}

##############
library(graphite)
#get kegg pathways
human_kegg <- pathways("hsapiens", "kegg")
keggps <- list()
for(i in 1:317){
  keggps[[i]] <- convertIdentifiers(human_kegg[[i]], "SYMBOL")
  cat("\n");
  cat(paste("Done  ",i,sep=""))
}

from <- as.vector(t(tissue_enrich_analysis[which(tissue_enrich_analysis[,2]==1),1]))
from <- as.vector(t(from))
agps <- list()
for(i in 1:length(from)){
  a <- from[i]
  gps <- NULL
  for(ii in 1:317){
    p1 <- keggps[[ii]]
    genes_p1<-nodes(p1)
    genes_p1 <- gsub("SYMBOL:","",genes_p1)
    a1 <- intersect(a,genes_p1)
    if(length(a1) != 0){
      gps <- c(gps, p1@id)
    }
  }
  agps[[i]] <- gps
}

to <- as.vector(t(tissue_enrich_analysis[which(tissue_enrich_analysis[,2]==2),1]))
to <- as.vector(t(to))
bgps <- list()
for(i in 1:length(to)){
  b <- to[i]
  gps <- NULL
  for(ii in 1:317){
    p1 <- keggps[[ii]]
    genes_p1<-nodes(p1)
    genes_p1 <- gsub("SYMBOL:","",genes_p1)
    b1 <- intersect(b,genes_p1)
    if(length(b1) != 0){
      gps <- c(gps, p1@id)
    }
  }
  bgps[[i]] <- gps
}



plotedges <- NULL
for(i in 1:length(from)){
  for(j in 1:length(to)){
    c <- intersect(agps[[i]],bgps[[j]])
    if(length(c)!=0){
      plotedges <- rbind(plotedges,c(from[i],to[j]))
    }
  }
  cat("\n");
  cat(paste("Done  ",i,sep=""))
}




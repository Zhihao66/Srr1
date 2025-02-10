rm(list=ls())
library(phytools)
library(dplyr)
library(parallel)
library("readr")
source("D:/6_1KFG/ALAPFILEOK/SFun220218.R")
library(DECIPHER)
setwd("D:/SRR1/ZhiaoQuery/")

# Open Reciproc Best Hit (RBH) Table: All CopciAB protein vs 189 species protein file
RBH<-read_csv("D:/SRR1/RBH_189/CopciAB_vs_189_RBH_MMseq0523.csv", col_names = T) 

# Open Differentiall Expressed gene table 
ZHQ<-read_tsv("Down_DEG_FC2_0408.txt")
colnames(ZHQ)<-"protID"

# ALL RBH Phylostratigraphy
ccirbh<-RBH[grep(".T0",RBH$Query),]
sz(ccirbh$Query)
ccirbh$protID<-sapply(str_split(ccirbh$Query,"\\."), "[", 1)
ccirbh$Qaz<-sapply(str_split(ccirbh$protID,"_"), "[", 2)
ccirbh<-dplyr::arrange(ccirbh,Query,Tspec,desc(bitScore),eValue)
szurt<-ccirbh[!duplicated(ccirbh[,c("Query",'Tspec')]),]  # we retain only the best RBH hits 
sz(szurt$Query)
sz(szurt$Target)
sz(szurt$Tspec)
save(list=ls(),file="FCesek_RBHi0409.saved")

#Open species tree
fa<-read.nexus("D:/SRR1/IQ1/SRR1_rooted.tre")
spt<-LadderizeTree(fa)

ir <- szurt %>% group_by(Query)
cldb<-group_split(ir)
length(cldb)

head(cldb)
x<-cldb[[6]]

Phylostrat<-function(x){
  tlbl<-which(spt$tip.label %in% unique(x$Tspec))
  #length(tlbl)
  if(length(tlbl)>1){
    mrca<-getMRCA(spt,tlbl)
  }else{
    mrca<-tlbl
  }
  return(c(unique(x$Query),mrca))
}

phs<-lapply(cldb,function(x) Phylostrat(x))
cls<-lfuzo(phs,1)


# Calculation of the enrichment of downregulated genes in each phylostratum------------------------------------------------------------
ZHQ$dataset<-"Down_DEG_FC2"
table(ZHQ$dataset)
ds<-"Down_DEG_FC2"
    
gy<-NULL
for(ds in unique(ZHQ$dataset)){
    cc<-ZHQ[which(ZHQ$dataset %in% ds),]
    
    clsZH<-bind_cols(cls,cc[match(cls$'1',cc$protID),"dataset"])
    colnames(clsZH)<-c("protID","mrca","dataset")
    
    # ALL
    aktdeg<-clsZH
    PTR<-data.frame(table(aktdeg$mrca))  
    sum(PTR$Freq)
    #as.character(PTR$Var1)
    PTT<-tibble(node=(length(spt$tip.label)+1):(length(spt$tip.label)+spt$Nnode))
    PTT<-rbind(PTT,178)
    PTT$db<-PTR[match(PTT$node,PTR$Var1),"Freq"]
    PTT[is.na(PTT$db),"db"]<-0
    
    unique(clsZH$dataset)
    aktdeg<-clsZH[!is.na(clsZH$dataset),]
    vPTR<-data.frame(table(aktdeg$mrca))
    sum(vPTR$Freq)
    PTT$akt<-vPTR[match(PTT$node,vPTR$Var1),"Freq"]
    PTT[is.na(PTT$akt),"akt"]<-0
    VP<-PTT[which(PTT$akt!=0|PTT$db!=0),]
    colSums(VP[,c(2,3)])
    
    File2<-paste0("nodelabels.pdf")
    pdf(file = File2, width=8,height=16, pointsize = 18, family = "sans", bg = "white")
    plot(spt,cex=0.5)
    nodelabels(node=VP$node,frame = "none",col="red",cex=0.8)
    dev.off()
    
    File2<-paste0("Phylostrat_scatter",ds,".pdf")
    pdf(file = File2, width=8,height=8, pointsize = 18, family = "sans", bg = "white")
      plot(VP$db,VP$akt,pch=16,col="white",main=ds)
      text(VP$db,VP$akt,label=VP$node,cex=0.45,col="green3")
    dev.off()
    
    VP$rate<-round(VP$akt/VP$db,3)
    write.table(VP, file=paste0('ZHQ_VP_',ds,'.csv'), quote=FALSE, sep=',', row.names = F, col.names = T)
    
    PTTm<-PTT[-nrow(PTT),]
    tblue<-t_col("blue",percent = 80)
    tred<-t_col("red",percent = 50)
    
      File2<-paste0("Phylostrat_",ds,".pdf")
      pdf(file = File2, width=12,height=20, pointsize = 18, family = "sans", bg = "white")
        plot.phylo(spt,cex=0.6, main=paste0(ds," proteins:",nrow(aktdeg)))
        nodelabels(pch=16,col=tblue,cex=(PTT$db)/500)
        nodelabels(pch=16,col=tred,cex=(PTT$akt)/10)
        nodelabels(paste0(PTTm$db,"/",PTTm$akt),col="black",cex=0.35, frame="none")
      dev.off()
    
      
      FETtab<-function(i){
          colSums(VP[,c(2,3)])
          ez<-VP[i,]
          others<-VP[-i,]
          oth<-colSums(others[c(2,3)])
          ABCD<-rbind(ez[,c(2,3)],oth)
          ABCD$db<-(ABCD$db-ABCD$akt)
          sum(colSums(ABCD))
          mx<-as.matrix(ABCD)
          aranyok<-round(mx[1,]/(mx[1,]+mx[2,]),3) 
          names(aranyok)<-c("background_IpA","DEG_IpA") #IpA = interesting in a certain node / sum of All (including the certain node)
          Fstat<-fisher.test(mx)
          PVAL<-Fstat$p.value
          names(PVAL)<-"pval"
          ODDS<-Fstat$estimate
          return(c(aranyok,PVAL,ODDS))
      }
      
      feter<-lapply(1:nrow(VP),function(h) FETtab(h))
      FETerr<-lfuzo(feter,1)
      FETerr$sig<-ifelse(FETerr$pval<0.05,"sig","non")
      colnames(VP)<-c("node","all","DEG","rate")
      DEGFET<-cbind(VP,FETerr)
      DEGFET$version<-ds
      write.table(DEGFET, file=paste0("DEGFET_",ds,SD,".csv"), quote=FALSE, sep=",", row.names=F, col.names=T)
      gy<-rbind(gy,DEGFET)
      save(list=ls(),file=paste0("SAVE_",ds,".saved"))
      
}
# output: Fisher's exact test results for each phylostratum
write.table(gy, file=paste0("DEGFET_ALL_",ds,SD,".csv"), quote=FALSE, sep=",", row.names=F, col.names=T)

sigek<-gy[which(gy$sig=="sig"),]
sigek$odddir<-ifelse(sigek$`odds ratio`>1,1,0)

# output: node names, number of all genes in the phylostratum and number of DEGs in the phylostratum 
PTT
write.table(PTT, file=paste0("PTT_",ds,SD,".csv"), quote=FALSE, sep=",", row.names=F, col.names=T)



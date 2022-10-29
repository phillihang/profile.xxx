library(limma)
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


gene=read.table("gene.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(data))
geneExp=data[sameGene,]

out=rbind(ID=colnames(geneExp),geneExp)
library("limma")
              
inputFile=".txt"                                             
fdrFilter=0.05                                                  
logFCfilter=1                                                   


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])      
treatNum=length(group[group==0])    
Type=c(rep(1,conNum), rep(2,treatNum))

outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	 }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)

heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffnecroptosisExp.txt",sep="\t",col.names=F,quote=F)


pdf(file="vol.pdf",height=5,width=5)
xMax=max(abs(as.numeric(as.vector(outTab$logFC))))
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="lightsalmon",cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="deepskyblue",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

library(pheatmap)
hmExp=data[as.vector(outDiff[,1]),]
hmExp=log2(hmExp+0.01)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=6,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         #scale = "row",
           color = colorRampPalette(c("deepskyblue", "white","lightsalmon"))(500),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=6,
         fontsize_col=10)
dev.off()

library(ggpubr)


rt=read.table("diffnecroptosisExp.txt",sep="\t",header=T,check.names=F,row.names=1)       #读取输入文件
rt=t(rt)
type=c(rep("N",conNum),rep("T",treatNum))

data=data.frame()
for(i in colnames(rt)){
  data=rbind(data,cbind(expression=log2(rt[,i]+1),gene=i,type))
}
write.table(data,file="diffnecroptosisExp1.txt",sep="\t",row.names=F,quote=F)

data=read.table("diffnecroptosisExp1.txt",sep="\t",header=T,check.names=F)      
p=ggboxplot(data, x="gene", y="expression", color = "type", 
     ylab="Gene expression",
     xlab="",
     palette = c("blue","orange") )
p=p+rotate_x_text(60)
pdf(file="boxplot.pdf",width=11,height=5)                          
p
dev.off()
install.packages("BiocManager", dependencies=T)
gene <- read.table
gene <- gene$V1
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
ego_ALL <- enrichGO(gene          = gene,
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Mm.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2)
ego_all <- data.frame(ego_ALL)
ego_all <- read.csv("GO.csv",row.names = 1,check.names = F)
library(ggplot2)
p <- ggplot(data=ego_all, aes(x= Description, y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = c("#8DA1CB", "#FD8D62")) + theme_bw() + 
  facet_grid(ONTOLOGY~.,scales = "free",space = "free")  +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")

genelist <- bitr(gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- genelist$ENTREZID          
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
kegg <- data.frame(ekegg)

#install.packages("survival")

library(survival)
pFilter=0.05                                                             

setwd                  

rt=read.table("geneexp.txt",header=T,sep="\t",check.names=F,row.names=1)  
rt$futime=rt$futime/365                                                    

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
	   if(sd(rt[,gene])<0.01){
	      next}
	   a=rt[,gene]<=median(rt[,gene])
	   diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	   pValue=1-pchisq(diff$chisq,df=1)
	   fit=survfit(Surv(futime, fustat) ~ a, data = rt)
	   cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	   coxSummary = summary(cox)
	   coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	   if((pValue<pFilter) & (coxP<pFilter)){
	         sigGenes=c(sigGenes,gene)
	       	 outTab=rbind(outTab,
	                      cbind(gene=gene,
	                            KM=pValue,
	                            B=coxSummary$coefficients[,"coef"],
	                            SE=coxSummary$coefficients[,"se(coef)"],
	                            HR=coxSummary$conf.int[,"exp(coef)"],
	                            HR.95L=coxSummary$conf.int[,"lower .95"],
	                            HR.95H=coxSummary$conf.int[,"upper .95"],
			                    pvalue=coxP) )
	  }
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)    
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)    
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)


rt <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))


pdf(file="forest.pdf", width = 6,height = 7)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()
#install.packages('survival')

library(survival)                                         
setwd("")      
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)


outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)


riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
    file="risk.txt",
    sep="\t",
    quote=F,
    row.names=F)
#install.packages("pheatmap")


summary(fit)
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)
setwd("")           
rt=read.table("Risk.txt",sep="\t",header=T,row.names=1,check.names=F)     
rt=rt[order(rt$riskScore),]                                   




riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])


line=rt[,"riskScore"]
line[line>10]=10


pdf(file="riskScore.pdf",width = 10,height = 3.5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("blue",lowLength),
     rep("orange",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("orange","blue"),cex=1.2)
dev.off()

color=as.vector(rt$fustat)
color[color==1]="orange"
color[color==0]="blue"
pdf(file="survStat.pdf",width = 10,height = 3.5)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("orange","blue"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

rt1=log2(rt[c(3:(ncol(rt)-2))]+0.01)
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap.pdf",width = 10,height = 3)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         color = colorRampPalette(c("blue", "white", "orange"))(50) )
dev.off()

data=read.table("risk.txt",header=T,sep="\t",check.names=F)
diff=survdiff(Surv(futime, fustat) ~risk,data = data)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = data)
pdf(file="survival.pdf",width=5.5,height=5)
plot(fit,
      lwd=2,
      col=c("red","blue"),
      xlab="Time (year)",
      ylab="Survival rate",
      main=paste("Survival curve (p=", pValue ,")",sep=""),
      mark.time=T)
legend("topright",
      c("High risk", "Low risk"),
     lwd=2,
      col=c("red","blue"))
dev.off()



inputFile="risk.txt"        
survFile="survival.pdf"        
rocFile="ROC.pdf"              



rt=read.table(inputFile,header=T,sep="\t")

diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=5,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table.height=.25)
pdf(file=survFile,onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()


ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(3,5,10),ROC=TRUE)
pdf(file=rocFile,width=5,height=5)
plot(ROC_rt,time=3,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=10,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
        c(paste0('AUC at 3 years: ',round(ROC_rt$AUC[1],3)),
          paste0('AUC at 5 years: ',round(ROC_rt$AUC[2],3)),
          paste0('AUC at 10 years: ',round(ROC_rt$AUC[3],3))),
        col=c("green",'blue','red'),lwd=2,bty = 'n')
dev.off()

#install.packages('survival')
#install.packages("survivalROC")

library(survival)
setwd("")                       
rt=read.table("indepInput.txt",header=T,sep="\t",check.names=F,row.names=1)


uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	 coxSummary = summary(cox)
	 uniTab=rbind(uniTab,
	              cbind(id=i,
	              HR=coxSummary$conf.int[,"exp(coef)"],
	              HR.95L=coxSummary$conf.int[,"lower .95"],
	              HR.95H=coxSummary$conf.int[,"upper .95"],
	              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	              )
}
write.table(uniTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)


multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)



bioForest=function(coxFile=null,forestFile=null){
		rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
		gene <- rownames(rt)
		hr <- sprintf("%.3f",rt$"HR")
		hrLow  <- sprintf("%.3f",rt$"HR.95L")
		hrHigh <- sprintf("%.3f",rt$"HR.95H")
		Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
		pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
		pdf(file=forestFile, width = 7,height = 4)
		n <- nrow(rt)
		nRow <- n+1
		ylim <- c(1,nRow)
		layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
		xlim = c(0,3)
		par(mar=c(4,2.5,2,1))
		plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
		text.cex=0.8
		text(0,n:1,gene,adj=0,cex=text.cex)
		text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
		text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
		xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
		plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
		arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
		abline(v=1,col="black",lty=2,lwd=2)
		boxcolor = ifelse(as.numeric(hr) > 1, 'orange', 'blue')
		points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
		axis(1)
		dev.off()
}



bioForest(coxFile="uniCox.txt",forestFile="uniForest.pdf")
bioForest(coxFile="multiCox.txt",forestFile="multiForest.pdf")




library(survivalROC)                 
rt=read.table("indepInput.txt",header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime/365
rocCol=rainbow(ncol(rt)-2)
aucText=c()

pdf(file="multiROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
  xlab="False positive rate", ylab="True positive rate",
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)


j=1
for(i in colnames(rt[,3:(ncol(rt)-1)])){
	roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
	j=j+1
	aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
	lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()

data <- read.table("clinical.txt",header = T,row.names = 1)
library(rms)
dd<-datadist(data)
options(datadist="dd")
options(na.action="na.delete")
summary(data$futime)
coxpbc<-cph(formula = Surv(futime,fustat) ~  age + stage + T + M + N + risk ,data=data,x=T,y=T,surv = T,na.action=na.delete)
print(coxpbc)
surv<-Survival(coxpbc) 
surv3<-function(x) surv(1095,x)
surv5<-function(x) surv(1825,x)
surv8<-function(x) surv(3650,x)


x<-nomogram(coxpbc,fun = list(surv3,surv5,surv8),lp=T,
            funlabel = c('3-year survival Probability','5-year survival Probability','10-year survival Probability'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
pdf("nomogram_classical.pdf",width = 12,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()
setwd("TCGA ESTIMATE")  

library(utils) 
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
expr <- read.table("BRCA_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

filterCommonGenes(input.f = "BRCA_fpkm_mRNA_01A.txt",   
                  output.f = "BRCA_fpkm_mRNA_01A.gct",   
                  id = "GeneSymbol")   
estimateScore("BRCA_fpkm_mRNA_01A.gct",   
              "BRCA_fpkm_mRNA_01A_estimate_score.txt",   
              platform="affymetrix")   

result <- read.table("BRCA_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))

rownames(result) <- colnames(expr)
write.table(result, file = "BRCA_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分

#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.pacakges('parallel')
#'       install.packages('preprocessCore')
#'       if preprocessCore is not available in the repositories you have selected, run the following:
#'           source("http://bioconductor.org/biocLite.R")
#'           biocLite("preprocessCore")
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#' This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
#' single-threaded in Windows.
#'
#' Usage:
#'       Navigate to directory containing R script
#'
#'   In R:
#'       source('CIBERSORT.R')
#'       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#'
#'       Options:
#'       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#'       ii) QN = Quantile normalization of input mixture (default = TRUE)
#'
#' Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
#' Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#' Core algorithm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
CoreAlg <- function(X, y){

  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)

  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#' do permutations
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #print(itor)

    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#' Main functions
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mixture_file heterogenous mixed expression
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @export
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){

  #read in data
  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- read.table(mixture_file, header=T, sep="\t", check.names=F)
  Y <- Y[!duplicated(Y[,1]),]
  rownames(Y)<-Y[,1]
  Y<-Y[,-1]
  X <- data.matrix(X)
  Y <- data.matrix(Y)

  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  # 矩阵X是LM22矩阵
  # 矩阵Y是待预测的矩阵
  # 然后对
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  while(itor <= mixtures){

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1

  }

  #save results
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}

setwd("cibersort")   
install.packages('e1071')
install.packages('parallel')
#install.packages("BiocManager")
BiocManager::install("preprocessCore", version = "3.13")
library(e1071)
library(parallel)
library(preprocessCore)
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'BRCA_fpkm_mRNA_01A.txt'  
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersort.Rdata")  
load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]  

mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) 
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, 
        names.arg = rep("",nrow(ciber.res)), 
        yaxt = "n", 
        ylab = "Relative percentage", 
        col = mycol) 
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), 
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, 
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   
setwd("cibersort_AY")
library(tidyverse)
a <- read.table("CIBERSORT-Results.txt", sep = "\t",row.names = 1,check.names = F,header = T)
a <- a[,1:22]
identical(rownames(a),rownames(group))
b <- group
class(b$group)
a$group <- b$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
#install.packages("ggsci")
#install.packages("tidyr")
#install.packages("ggpubr")

b <- gather(a,key=CIBERSORT,value = Proportion,-c(group,sample))

ggboxplot(b, x = "CIBERSORT", y = "Proportion",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "kruskal.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()

setwd("ssGSEA_AY")
BiocManager::install('GSVA')
library(tidyverse)
library(data.table)
library(GSVA)
cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
colnames(cellMarker)[2] <- "celltype"
type <- split(cellMarker,cellMarker$celltype)
cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})

save(cellMarker,file = "cellMarker_ssGSEA.Rdata")
load("immune_infiltration//cellMarker_ssGSEA.Rdata")

expr <- data.table::fread("LIHC_fpkm_mRNA_01A.txt",data.table = F)   
rownames(expr) <- expr[,1]   
expr <- expr[,-1]   
expr <- as.matrix(expr)   


gsva_data <- gsva(expr,cellMarker, method = "ssgsea")

a <- gsva_data %>% t() %>% as.data.frame()
identical(rownames(a),rownames(group))
a$group <- group$group
a <- a %>% rownames_to_column("sample")
write.table(a,"ssGSEA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=ssGSEA,value = Expression,-c(group,sample))

ggboxplot(b, x = "ssGSEA", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

library(GSVA)
library(tidyverse)
library(ggpubr)
library(ggsci)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

immunity <- read.table("immune_functions.txt", header = F)
colnames(immunity) <- c("CellType","Symbol")

immunity <- immunity %>% 
  split(., .$CellType) %>% 
  lapply(., function(x)(x$Symbol))
immunity <- lapply(immunity, unique)
tcga_expr <- read.table("BRCA_fpkm_mRNA_01A.txt", row.names = 1,check.names = F)
group <- read.table("risk.txt",header = T)
group$group <- ifelse(group$group == 1, "low risk","high risk")
tcga_expr1 <- tcga_expr[,group$id]

tcga_gsva <- as.data.frame(t(gsva(as.matrix(tcga_expr1), immunity, method = "ssgsea")))
tcga_gsva$group <- group$group
data <- reshape2::melt(tcga_gsva)
ggplot(data, aes(variable, value)) +
  geom_boxplot(aes(color = group, fill = group))+
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("#ED0000FF", "#00468bFF")) +
  stat_compare_means(aes(group = group), label = "p.signif") + theme_bw() + Seurat::RotatedAxis() +
  theme(title = element_text(size = 20,colour = 'black'),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.ticks.y = element_line(color='black', size=2, lineend = 10),
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.key.size=unit(.5,'cm'),
        legend.position = "top",
        legend.box = "",
        axis.text.y = element_text(size = 11,colour = 'black'),
        axis.text.x = element_text(size = 11,colour = 'black'),
        axis.title.x = element_text(size = 15,colour = 'black'),
        axis.title.y = element_text(size = 15,colour = 'black'),
        strip.text.x = element_text(size = 15,colour = 'black'))
#install.packages('survival')

library(survival)                                         
setwd("")       
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime/365


multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)


outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)


riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
    file="risk.txt",
    sep="\t",
    quote=F,
    row.names=F)

library(survival)
library(survminer)
library(timeROC)

inputFile="risk.txt"       
survFile="survival.pdf"        
rocFile="ROC.pdf"               
setwd("")      


rt=read.table(inputFile,header=T,sep="\t")


diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}

surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=5,
		           risk.table=TRUE,
		           legend.labs=c("High risk", "Low risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 1,
		           risk.table.title="",
		           palette=c("red", "blue"),
		           risk.table.height=.25)
pdf(file=survFile,onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()



ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file=rocFile,width=5,height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
        c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],3)),
          paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],3)),
          paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],3))),
        col=c("green",'blue','red'),lwd=2,bty = 'n')
dev.off()

library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("")         
files=grep("tsv",dir(),value=T)                                         
data = lapply(files,read.delim)                                         
names(data) = files

dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".tsv","",dataSet$.id)                           

gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,colour=pathway,group=pathway))+
  geom_line(size = 1.5) + scale_color_manual(values = gseaCol[1:nrow(dataSet)]) +   
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits = c(min(dataSet$RUNNING.ES - 0.02), max(dataSet$RUNNING.ES + 0.02))) +   
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) +   theme(legend.position = c(0,0),legend.justification = c(0,0)) + #legend注释的位值
  guides(colour = guide_legend(title = NULL)) + theme(legend.background = element_blank()) + theme(legend.key = element_blank())+theme(legend.key.size=unit(0.5,'cm'))
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "high expression<----------->low expression", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)

gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()

pdf('multipleGSEA.pdf',      
     width=7,                
     height=5.5)             
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="symbol.txt"      
riskFile="riskAll.txt"      
geneFile="gene.txt"      
setwd("")    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)


group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
data=avereps(data)
	

risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c(sameGene,"risk")]


sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
	if(sd(rt1[,i])<0.001){next}
	wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"risk"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "risk")
rt1=rt1[,sigGene]


rt1=melt(rt1,id.vars=c("risk"))
colnames(rt1)=c("risk","Gene","Expression")
	
group=levels(factor(rt1$risk))
rt1$risk=factor(rt1$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
	
boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Risk",
				  width=0.8,
				  palette = c("#00468bFF", "#ED0000FF") )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=risk),
	method="wilcox.test",
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
	
pdf(file="checkpoint.diff.pdf", width=12, height=7)
print(boxplot)
dev.off()


#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


library(limma)
library(ggpubr)
gene="CD274"              
expFile="symbol.txt"     
riskFile="riskAll.txt"      
setwd("")    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

data=data[gene,,drop=F]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


sameSample=intersect(row.names(risk), row.names(data))
risk=risk[sameSample, "risk",drop=F]
data=data[sameSample, , drop=F]
rt=cbind(risk, data)


rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#??????ͼ
pdf(file=paste0(gene, ".vioplot.pdf"), width=6, height=5)
ggviolin(rt, x="risk", y=gene, fill = "risk", 
         xlab="Risk", ylab=paste0(gene, " expression"),
         legend.title="Risk",
         palette=c("green", "red"),
         add = "boxplot", add.params = list(fill="white"))+ 
         stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

options(stringsAsFactors = FALSE)
setwd("")

library(readxl)
rt1 <- read_excel(path = "DTP_NCI60_ZSCORE.xlsx", skip = 7)
colnames(rt1) <- rt1[1,]
rt1 <- rt1[-1,-c(67,68)]

table(rt1$`FDA status`) 
rt1 <- rt1[rt1$`FDA status` %in% c("FDA approved", "Clinical trial"),]
rt1 <- rt1[,-c(1, 3:6)]
write.table(rt1, file = "drug.txt",sep = "\t",row.names = F,quote = F)

rt2 <- read_excel(path = "RNA__RNA_seq_composite_expression.xls", skip = 9)
colnames(rt2) <- rt2[1,]
rt2 <- rt2[-1,-c(2:6)]
write.table(rt2, file = "geneExp.txt",sep = "\t",row.names = F,quote = F)


rm(list = ls())
library(impute)
library(limma)

rt <- read.table("drug.txt",sep="\t",header=T,check.names=F)
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
drug <- rt[,2:ncol(rt)]
dimnames <- list(rownames(drug),colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)

mat <- impute.knn(data)
drug <- mat$data
drug <- avereps(drug)

exp <- read.table("geneExp.txt", sep="\t", header=T, row.names = 1, check.names=F)
dim(exp)
exp[1:4, 1:4]


gene <- read.table("gene.txt",sep="\t",header=F,check.names=F)
genelist <- as.vector(gene[,1])
genelist

genelist <- gsub(" ","",genelist)
genelist <- intersect(genelist,row.names(exp))
exp <- exp[genelist,]

outTab <- data.frame()
for(Gene in row.names(exp)){
  x <- as.numeric(exp[Gene,])
  for(Drug in row.names(drug)){
    y <- as.numeric(drug[Drug,])
    corT <- cor.test(x,y,method="pearson")
    cor <- corT$estimate
    pvalue <- corT$p.value
    if(pvalue < 0.01){
      outVector <- cbind(Gene,Drug,cor,pvalue)
      outTab <- rbind(outTab,outVector)
    }
  }
}
outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab, file="drugCor.txt", sep="\t", row.names=F, quote=F)

library(ggplot2)
library(ggpubr)

plotList_1 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab)
}
for(i in 1:corPlotNum){
  Gene <- outTab[i,1]
  Drug <- outTab[i,2]
  x <- as.numeric(exp[Gene,])
  y <- as.numeric(drug[Drug,])
  cor <- sprintf("%.03f",as.numeric(outTab[i,3]))
  pvalue=0
  if(as.numeric(outTab[i,4])<0.001){
    pvalue="p<0.001"
  }else{
    pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
  }
  df1 <- as.data.frame(cbind(x,y))
  p1=ggplot(data = df1, aes(x = x, y = y))+
    geom_point(size=1)+
    stat_smooth(method="lm",se=FALSE, formula=y~x)+
    labs(x="Expression",y="Z score",title = paste0(Gene,", ",Drug),subtitle = paste0("Cor=",cor,", ",pvalue))+
    theme(axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank())+
    theme_bw()
  plotList_1[[i]]=p1
}

plotList_2 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab)
}
for(i in 1:corPlotNum){
  Gene <- outTab[i,1]
  Drug <- outTab[i,2]
  x <- as.numeric(exp[Gene,])
  y <- as.numeric(drug[Drug,])
  df1 <- as.data.frame(cbind(x,y))
  colnames(df1)[2] <- "IC50"
  df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
  compaired <- list(c("low", "high"))
  p1 <- ggboxplot(df1,
                  x = "group", y = "IC50",
                  fill = "group", palette = c("#00AFBB", "#E7B800"),
                  add = "jitter", size = 0.5,
                  xlab = paste0("The_expression_of_", Gene),
                  ylab = paste0("IC50_of_", Drug)) +
    stat_compare_means(comparisons = compaired,
                       method = "wilcox.test",   
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")))
  plotList_2[[i]]=p1
}

nrow <- ceiling(sqrt(corPlotNum))
ncol <- ceiling(corPlotNum/nrow)
ggarrange(plotlist=plotList_1,nrow=nrow,ncol=ncol)
ggarrange(plotlist=plotList_2,nrow=nrow,ncol=ncol)


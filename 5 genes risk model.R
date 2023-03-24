#########################################################################################
rm(list=ls())
library(ggrisk)
library(survival)
library(rms)
library(pROC)
library(Seurat)
library(glmnet)

## The univariate Cox regression analysis was performed to determine the association between the cluster 10 marker genes

PDX_human_cell_tsne <- readRDS("PDX_human_cell_tsne.rds")
tcga_data <-readRDS("tcga_PAAD.rds")
scRNA=PDX_human_cell_tsne
logFCfilter=0.5
adjPvalFilter=0.05
scRNA.markers <- FindAllMarkers(object = scRNA,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)
sig.markers=scRNA.markers[(abs(as.numeric(as.vector(scRNA.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(scRNA.markers$p_val_adj))<adjPvalFilter),]
sig.markers10 <- scRNA.markers %>%dplyr::filter(cluster ==10)
ids=intersect(sig.markers10$gene,rownames(tcga_data))
tcga_dat_m=tcga_data[ids,]
tcga_clin=tcga_clin[,c(1:5)]
tcga_dat_m=cbind(tcga_clin,t(tcga_dat_m))

# Univariate Cox regression analysis
library(survival)
tcga.cox=c()
for (i in colnames(tcga_dat_m)[6:ncol(tcga_dat_m)]){
  dat=data.frame(time=tcga_dat_m$OS.time,
                 status=tcga_dat_m$OS,
                 gene=tcga_dat_m[,i])
  fmla <- as.formula("Surv(time, status) ~gene")
  cox <- survival::coxph(fmla, data = dat)
  re = c(summary(cox)[[7]][5], 
         summary(cox)[[7]][2], 
         summary(cox)[[8]][3], 
         summary(cox)[[8]][4])
  tcga.cox=rbind(tcga.cox,re)
}
head(tcga.cox)
rownames(tcga.cox)=colnames(tcga_dat_m)[6:ncol(tcga_dat_m)]
colnames(tcga.cox) = c("p.value", "HR", "Low 95%CI", "High 95%CI")
tcga.cox=as.data.frame(tcga.cox)
head(tcga.cox)

## LASSO-analysis
sig.names=rownames(tcga.cox[which(tcga.cox$p.value < 0.01),])
set.seed(2021)###do not change
fit1=glmnet(as.matrix(tcga_dat_m[,sig.names])
            ,cbind(time=tcga_dat_m$OS.time,
                   status=tcga_dat_m$OS)
            ,family="cox"
            ,nlambda=100
            , alpha=1) 
cv.fit<-cv.glmnet(as.matrix(tcga_dat_m[,sig.names])
                  ,cbind(time=tcga_dat_m$OS.time,
                         status=tcga_dat_m$OS)
                  ,family="cox"
                  ,nlambda=100
                  , alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min 
plot(fit1, xvar="lambda")
plot(cv.fit)
tcga_dat1 <- data.frame(time=tcga_dat_m$OS.time,
                        status=tcga_dat_m$OS,
                        tcga_dat_m[,names(sig.coef)])
fmla <- as.formula(paste0("Surv(time, status) ~"
                          ,paste0(names(sig.coef),collapse = '+')))
cox1 <- coxph(fmla, data =as.data.frame(tcga_dat1))
lan <- coef(cox1)
round(lan, 3)
genes <- names(cox1$coefficients)
paste0(round(lan, 3), '*', names(lan))

## The risk score formula constructed
# Calculated risk score
risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
risk.tcgaz <- mosaic::zscore(risk.tcga)
tcga_dat1$time=tcga_dat1$time/365
ROC_rt=timeROC::timeROC(T=tcga_dat1$time, 
                        delta=tcga_dat1$status,
                        marker=risk.tcgaz, cause=1,
                        weighting='marginal',
                        times=c(1,3,5), 
                        ROC=TRUE,iid = T)
p.dat=rbind(data.frame(V1=ROC_rt$FP[,1],V2=ROC_rt$TP[,1],Type=paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],2))),
            data.frame(V1=ROC_rt$FP[,2],V2=ROC_rt$TP[,2],Type=paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],2))),
            data.frame(V1=ROC_rt$FP[,3],V2=ROC_rt$TP[,3],Type=paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],2))))
p.dat=as.data.frame(p.dat)
p.dat$Type=as.factor(p.dat$Type)
roc_plot=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))+
  stat_smooth(aes(colour=Type),se = FALSE, size = 1)+
  theme_bw()+xlab('False positive fraction')+ylab('True positive fraction') +
  theme(axis.text.y=element_text(family="Times",face="plain"),
        axis.text.x=element_text(family="Times",face="plain"),
        axis.title.x=element_text(family="Times",face="plain"),
        axis.title.y=element_text(family="Times",face="plain"),
        plot.title=element_blank(),
        plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),
        legend.position=c(1,0),
        legend.justification=c(1,0),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.title = element_text(family="Times",face="plain"),
        legend.text = element_text(family="Times",face="plain"))

# Survival curve
dat_km<-data.frame(time=tcga_dat1$time,
                   status=tcga_dat1$status,
                   riskscore=risk.tcga,
                   risk=ifelse(risk.tcgaz>=0,'High','Low'))
rownames(dat_km)=rownames(tcga_dat1)
diff=survdiff(Surv(time, status) ~risk,data = dat_km)

pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(time, status) ~ risk, data = dat_km)
tcga_km=survminer::ggsurvplot(fit, 
                              data=dat_km,
                              conf.int=T,
                              pval=pValue,
                              pval.size=6,
                              legend.title="Risk",
                              legend.labs=c("High risk", "Low risk"),
                              xlab="Time(years)",
                              #break.time.by = 3,
                              palette=c("red", "blue"),
                              risk.table=TRUE,
                              risk.table.title="",
                              risk.table.col = "strata",
                              risk.table.height=.25)

# multivariate Cox regression analysis
ids=intersect(rownames(dat_km),tcga_clin$sample)
dat_km=dat_km[ids,]
match(rownames(dat_km),tcga_clin$sample)
tcga_clin=tcga_clin[,-c(1:5)]
colnames(tcga_clin)=c("age","grade","gender","size","stage")
rt=cbind(tcga_clin,dat_km)
multiCox=coxph(Surv(time, status) ~riskscore+stage+age+gender+grade+size, data = rt)
multiCoxSum=summary(multiCox)
outMultiTab=data.frame()
outMultiTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
outMultiTab=outMultiTab[,-2]
write.table(outMultiTab,file="TCGA.multiCox.txt",sep="\t",row.names=F,quote=F)

# Forest map
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  height=nrow(rt)/12.5+5
  pdf(file=forestFile, width = 7,height = height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  Ϣ
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
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}
bioForest(coxFile="TCGA.multiCox.txt",forestFile="forest.pdf",forestCol=c("red","green"))

# Chi-square test
library(ggstatsplot)
ggbarstats(rt, x = treatment.adj, y = risk)+theme_classic() 
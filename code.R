
library(sva)
library(limma)
library(glmnet)
library(Boruta)

rt1=read.table("GSE142530_fpkm.txt",sep="\t",header=T,check.names=F)    
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1=avereps(data1)

rt2=read.table("GSE28619_GPL570_sample_exp.txt",sep="\t",header=T,check.names=F)   
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=avereps(data2)

sameGene=intersect( row.names(data1),row.names(data2) )
data=cbind(data1[sameGene,],data2[sameGene,])

group1 <- read.table("GSE142530_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group2 <- read.table("GSE28619_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

tum.sam1 <- rownames(group1[which(group1$Group == 1),,drop = F])
tum.sam2 <- rownames(group2[which(group2$Group == 1),,drop = F])
nor.sam1 <- rownames(group1[which(group1$Group == 0),,drop = F])
nor.sam2 <- rownames(group2[which(group2$Group == 0),,drop = F])

tum.sam <- c(tum.sam1, tum.sam2)
nor.sam <- c(nor.sam1, nor.sam2)
data <- log2(data[,c(tum.sam1,nor.sam1,tum.sam2,nor.sam2)] + 1)     

write.table(data,file="merge.txt",sep="\t",quote=F,row.names = T,col.names = NA) 

rt=read.table("merge.txt",sep="\t",header=T,check.names=F, row.names = 1,stringsAsFactors = F) 
rt=as.matrix(rt)

batch <- data.frame(batch = rep(c("Group1","Group2"), times = c(nrow(group1),nrow(group2))))
modcombat = model.matrix(~1, data=batch)
outTab <- as.data.frame(ComBat(dat=as.matrix(data), batch=batch$batch, mod=modcombat)) 
outTab <- outTab[,c(nor.sam, tum.sam)]
write.table(outTab,file="normalize.txt",sep="\t",quote=F,row.names = T,col.names = NA)

library("FactoMineR")
library("factoextra")
pca.plot = function(data,col){
        
        df.pca <- PCA(t(data), graph = FALSE)
        fviz_pca_ind(df.pca,
                     geom.ind = "point",
                     col.ind = col ,
                     addEllipses = TRUE,
                     legend.title = "Groups"
        )
}
data<-read.table("merge.txt",sep = "\t",row.names = 1,header = T)
group<-read.table("group.txt",sep = "\t",row.names = 1,header = T)
a<-pca.plot(data,factor(group$Group))
pdf(file="PCA-raw.pdf",width=6,height=5)
a
dev.off()
data1<-read.table("normalize.txt",sep = "\t",row.names = 1,header = T)
b<-pca.plot(data1,factor(group$Group))
pdf(file="PCA-after.pdf",width=6,height=5)
b
dev.off()

logFoldChange=1                                     
adjustP=0.05
rt=read.table("normalize.txt",sep="\t",header=T,check.names=F)                         
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
rt=normalizeBetweenArrays(as.matrix(rt))

modType=c(rep("con",length(nor.sam)),rep("treat",length(tum.sam))) 
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F,row.names = T,col.names = NA)

diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & P.Value < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F,row.names = T,col.names = NA)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & P.Value  < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F,row.names = T,col.names = NA)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & P.Value  < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F,row.names = T,col.names = NA)

hmExp=rt[rownames(diffSig),]
diffExp=hmExp
write.table(diffExp,file="heatmap.txt",sep="\t",quote=F,col.names=NA,row.names = T)

tiff(file="vol.tiff",
     width = 12,            
     height =12,           
     units ="cm",
     compression="lzw",
     bg="white",
     res=600)
xMax=max(abs(allDiff$logFC))
yMax=max(-log10(allDiff$P.Value))
plot(allDiff$logFC, -log10(allDiff$P.Value), xlab="logFC",ylab="-log10(P.Value)",
     main="Volcano",xlim=c(-xMax,xMax),ylim=c(0,yMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(allDiff, P.Value<adjustP & logFC>logFoldChange)
points(diffSub$logFC, -log10(diffSub$P.Value), pch=20, col="palevioletred1",cex=0.8)
diffSub=subset(allDiff, P.Value<adjustP & logFC<(-logFoldChange))
points(diffSub$logFC, -log10(diffSub$P.Value), pch=20, col="dodgerblue",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()






library(limma)
library(pheatmap)
library(gplots)



plotdata=read.table("heatmap.txt",sep="\t",header=T,check.names=F,row.names = 1)
up=read.table("up.txt",header=F,stringsAsFactors = F)
down=read.table("down.txt",header=F,stringsAsFactors = F)

annCol <- data.frame(Tissue = rep(c("Control","Disease"),c(19,25)),
                     row.names = colnames(plotdata),
                     stringsAsFactors = F)

annRow <- data.frame(Direct = rep(c("Up","Down"),c(length(up$V1),length(down$V1))),
                     row.names = c(up$V1,down$V1),
                     stringsAsFactors = F)

annColors <- list("Tissue"=c("Control"="blue",
                             "Disease"="red"),
                  "Direct"=c("Up"="yellow",
                             "Down"="green"))

plotdata <- plotdata[c(up$V1,down$V1),]
plotdata <- t(scale(t(plotdata)))
plotdata[plotdata > 5] <- 5
plotdata[plotdata < -5] <- -5

pdf(file="heatmap_tn.pdf",width = 8,height = 8)
pheatmap(plotdata, 
         scale = "none",
         annotation_row=annRow, 
         annotation_col=annCol, 
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(20),
         fontsize_row=7,
         fontsize_col=5,
         fontsize=7,
         show_rownames = F,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = F)
dev.off()







library(tidyverse)
library(glmnet)
source('msvmRFE.R')   
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)

train<-read.csv("easy_input.csv",row.names = 1,as.is = F,header = T)
train[1:4,1:4]
x <- as.matrix(train[,-1])  
(y <- ifelse(train$group == "NR", 0,1))
x[1:4,1:4]

set.seed(94)
fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
plot(fit, xvar = "dev", label = TRUE)
cvfit = cv.glmnet(x, y, 
                  nfold=10, 
                  family = "binomial", type.measure = "class")
plot(cvfit)
cvfit$lambda.min 
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
write.csv(lasso_fea,"feature_lasso.csv")

set.seed(94)
predict <- predict(cvfit, newx = x[1:nrow(x),], s = "lambda.min", type = "class")
table(predict,y)
input <- train
svmRFE(input, k = 10, halve.above = 100) 

nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)

top.features = WriteFeatures(results, input, save=F)
head(top.features)
write.csv(top.features,"feature_svm.csv")

featsweep = lapply(1:5, FeatSweep.wrap, results, input)
featsweep
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

PlotErrors(errors, no.info=no.info)

Plotaccuracy(1-errors,no.info=no.info) 

which.min(errors)

(myoverlap <- intersect(lasso_fea, top.features[1:which.min(errors), "FeatureName"])) 

summary(lasso_fea%in%top.features[1:which.min(errors), "FeatureName"])

pdf("C_lasso_SVM_venn.pdf", width = 5, height = 3)
grid.newpage()
venn.plot <- venn.diagram(list(LASSO = lasso_fea, 
                               SVM_RFE = as.character(top.features[1:which.min(errors),"FeatureName"])), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("LASSO", "SVM_RFE"), 
                          main = "Overlap")
grid.draw(venn.plot)
dev.off()

library(ggplot2)
library(ggpubr)
pathway = read.csv("feature_svm.csv",header=TRUE,row.names=1,check.names = FALSE) 
p = ggplot(pathway,aes(reorder(FeatureName,AvgRank),AvgRank))
p=p + geom_point()  
p=p + geom_point(aes(size=AvgRank))
pbubble = p+ geom_point(aes(size=AvgRank,color=AvgRank))
pr = pbubble+scale_color_gradient(low="green",high = "red")
pr = pr+labs(color=expression(AvgRank),size="AvgRank",  
             x="FeatureName",y="AvgRank")
pr + theme_bw()+rotate_x_text(45)+scale_color_continuous(low='purple', high='green')+theme(axis.text.x=element_text(size=12))





library(randomForest)
set.seed(589)
exp <- read.table("ARGexp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- as.matrix(exp)
exp <- pmax(exp,0)*2
exp<-log2(exp+1)
group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
com_sam <- intersect(colnames(exp),rownames(group))
var <- apply(exp, 1, mad)
exp <- exp[var > quantile(var)[4],]
rownames(exp) <- gsub("-","_",rownames(exp))
dat <- cbind.data.frame(t(exp[,com_sam]),group[com_sam,,drop = F])
dat$Group <- ifelse(dat$Group == "control",0,1)
group <- factor(group$Group)
model_RF <- randomForest(Group ~ ., 
                         data = dat,
                         ntree = 1000,
                         nPerm = 50,
                         mtry = floor(sqrt(ncol(dat)-1)), 
                         proximity = T,
                         importance = T)

varImpPlot(model_RF)
imp <- as.data.frame(importance(model_RF, type = 1))
sel_var <- imp[order(imp$`%IncMSE`,decreasing = T),,drop = F]
imp.ori <- as.data.frame(importance(model_RF))
write.table(imp.ori,"importance of variables.txt",sep = "\t",row.names = T,col.names = NA,quote = F)





source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "normalize.txt", perm=100, QN=TRUE)






library(dplyr)
library(tidyr)
library(tibble)

res_cibersort<-read.table("CIBERSORT-Results.txt",sep = "\t",header = T,row.names = 1)
res_cibersort<-res_cibersort[,]
dd1 <- res_cibersort %>% 
        as.data.frame() %>% 
        rownames_to_column("sample") %>% 
        mutate(type=c(rep("Control",19),rep("Disease",25))) %>%                     
        pivot_longer(cols=2:23,
                     names_to= "celltype",
                     values_to = "Proportion")

library(ggplot2)
ggplot(dd1,aes(sample,Proportion,fill = celltype)) + 
        geom_bar(position = "fill",stat = "identity")+
        theme_bw()+
        guides(fill=guide_legend(ncol=1))+
        facet_wrap(~type, scales = "free_x", nrow = 1)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))






rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
library(corrplot)
pdf("corHeatmap.pdf",height=13,width=13)            
corrplot(corr=cor(rt),
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
dev.off()






library(vioplot)                                                
normal=19                                                       
tumor=25                                                       
rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)   
pdf("vioplot.pdf",height=8,width=15)            
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,61),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
text(seq(1,63,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)

for(i in 1:ncol(rt)){
        normalData=rt[1:normal,i]
        tumorData=rt[(normal+1):(normal+tumor),i]
        vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
        vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
        wilcoxTest=wilcox.test(normalData,tumorData)
        p=round(wilcoxTest$p.value,3)
        mx=max(c(normalData,tumorData))
        lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
        text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
}
dev.off()







library(ggplot2)
library(dplyr)
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
expr <- read.table("normalize.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("genelist.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = F)
ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
for (i in gene$V1) {
        message(paste0("analysis of ",i," starts..."))
        subexpr <- as.numeric(expr[i,])
        names(subexpr) <- colnames(expr)
        
        lsam <- names(subexpr[subexpr < median(subexpr)])
        hsam <- names(subexpr[subexpr >= median(subexpr)])

        dat <- as.numeric(expr[i,]); names(dat) <- colnames(expr)
        comsam <- intersect(names(dat), rownames(ciber))
        tmp1 <- dat[comsam]
        tmp2 <- ciber[comsam,]
        
        var <- colnames(ciber)
        data <- data.frame(var)
        for (j in 1:length(var)){
                test <- cor.test(as.numeric(tmp2[,j]),tmp1,method = "spearman")
                data[j,2] <- test$estimate                                            
                data[j,3] <- test$p.value
        }
        names(data) <- c("symbol","correlation","pvalue")
        data <- as.data.frame(na.omit(data))
        data %>% 
                filter(pvalue <0.05) %>%  
                ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
                geom_segment(aes(xend=0,yend=symbol)) +
                geom_point(aes(col=pvalue,size=abs(correlation))) +
                scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
                scale_size_continuous(range =c(2,8))  +
                theme_minimal() +
                ylab(NULL)
        ggsave(paste0("correlation between CIBERSORT and expression of ", i,".pdf"),width = 8,height = 6)
}






library(ggplot2)
library(stringr)
library(limma)
a_1 <- read.table("normalize.txt",header = T,row.names = 1,sep = "\t", quote = "",fill = T,check.names=F)
dim(a_1)
head(a_1[,1:3])
a_2 <- as.data.frame(t(a_1))
dim(a_2)
head(a_2[,1:3])

a_3 <- a_1
a_3$Id <- rownames(a_3)
dim(a_3)
head(a_3[,1:3])

b_1 <- read.table("Immunomodulator_and_chemokines(1).txt",header = T,sep = "\t", quote = "",fill = T)         
dim(b_1)
head(b_1)
b_2 <- b_1[b_1$type == "MHC",]
# b_2 <- b_1[b_1$type == "Immunoinhibitor",]
# b_2 <- b_1[b_1$type == "Immunostimulator",]
# b_2 <- b_1[b_1$type == "chemokine",]
# b_2 <- b_1[b_1$type == "receptor",]

dim(b_2)
head(b_2)
data1 <- dplyr::inner_join(b_2,a_3,by="Id")
dim(data1)
head(data1[,1:6])
data2 <- a_2[,c("MOXD1", "PDZK1IP1","SLC51B",data1$Id)]     
dim(data2)
head(data2[,1:5])
library(Hmisc)
CorMatrix <- function(cor,p) {
        ut <- upper.tri(cor) 
        data.frame(row = rownames(cor)[row(cor)[ut]] ,
                   column = rownames(cor)[col(cor)[ut]], 
                   cor =(cor)[ut], 
                   p = p[ut] )
}

res <- rcorr(as.matrix(data2))
result_1 <- CorMatrix(res$r, res$P)
head(result_1)

dim(result_1)    
result_2 <- result_1[result_1$row == "MOXD1" | result_1$row == "PDZK1IP1" | result_1$row == "SLC51B",]   
dim(result_2)
head(b_2)
b_2$column <- b_2$Id
head(b_2)
result_3 <- dplyr::inner_join(result_2,b_2,by="column")
dim(result_3)
result1 <- result_3[,1:4]
head(result1)
dim(result1)
result1$Regulation <- result1$cor
result1[,5][result1[,5] > 0] <- c("postive")
result1[,5][result1[,5] < 0] <- c("negative")
head(result1)
colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")

write.table(result1,file="corResult_MHC.xls",sep="\t",quote=F,col.names=T,row.names = F)  

a1 <- read.table("corResult_MHC.xls",header = T,sep = "\t", quote = "",fill = T)      
head(a1)

data2 <- a1

data2$pv <- -log10(data2$pvalue)
data2$type <- data2$cor
data2[,7][data2[,7] > 0 ] <- 2
data2[,7][data2[,7] < 0 ] <- 1
summary(data2)
data3 <- data2[order(data2$immuneGene,data2$cor),]
head(data3)
dim(data3)
data4 <- data3
dim(data4)
summary(data4)
p <- ggplot(data4,aes(x=gene,y=immuneGene)) +
        geom_point(aes(colour = cor, size=pv)) +
        labs(x="gene",y="MHC")    
p <- p + scale_colour_gradient2(low = "blue", high = "red", mid = "white",
                                midpoint = 0,  space = "Lab",
                                name="Pearson\nCorrelation")
p <- p + theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        theme(axis.text=element_text(size = 15)) +
        theme(axis.text.x=element_text(colour = "black",angle=0,hjust=0.5,size = 15)) +
        theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
        theme(axis.title =element_text(size = 20)) +
        theme(text = element_text(size = 15))

p

ggsave("MHC.pdf",height = 9,width = 7)








library(ggpubr)

pFilter=0.99

rt=read.table("ARGexp.txt",sep="\t",header=T,row.names=1,check.names=F)    
data=rt

Type=read.table("cluster.Immunity.txt",sep="\t",check.names=F,header=F,row.names = 1)
colnames(Type)=c("Subtype")

outTab=data.frame()
data=cbind(data,Type)
for(i in colnames(data[,1:(ncol(data)-1)])){
        rt1=data[,c(i,"Subtype")]
        colnames(rt1)=c("expression","Subtype")
        ksTest<-t.test(expression ~ Subtype, data = rt1)
        pValue=ksTest$p.value
        if(pValue<pFilter){
                outTab=rbind(outTab,cbind(rt1,gene=i))
                print(pValue)
        }
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)
data=read.table("data.txt",sep="\t",header=T,check.names=F)       
data$Subtype=factor(data$Subtype, levels=c("Normal","HIV"))

library(dplyr)
library(ggplot2)
data=read.table("data.txt",sep="\t",header=T,check.names=F)      
data$Subtype=factor(data$Subtype, levels=c("Control","Disease"))
data %>% 
        filter(Subtype %in% c("Control","Disease")) %>% 
        ggplot(aes(x= gene, y= expression, fill = Subtype, color = Subtype))+
        geom_boxplot()+
        scale_fill_manual(name= "Subtype", values = c("deepskyblue", "pink"))+
        scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3"))+
        theme_bw()+labs(x="", y="Fraction")+
        theme(axis.text.x = element_text(colour = "black",vjust = 1,size = 12, hjust = 1),legend.position="top")+theme(text = element_text(size = 15))+
        rotate_x_text(45)+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method = "t.test")







library(ggplot2)
library(stringr)
a_1 <- read.table("symbol.txt",header = T,row.names = 1,sep = "\t", quote = "",fill = T,check.names=F)
dim(a_1)
head(a_1[,1:3])
a_2 <- as.data.frame(t(a_1))
dim(a_2)
head(a_2[,1:3])

a_3 <- a_1
a_3$Id <- rownames(a_3)
dim(a_3)
head(a_3[,1:3])

b_1 <- read.table("Immunomodulator_and_chemokines.txt",header = T,sep = "\t", quote = "",fill = T)
dim(b_1)
head(b_1)
b_2 <- b_1[b_1$type == "chemokine",]
dim(b_2)
head(b_2)

data1 <- dplyr::inner_join(b_2,a_3,by="Id")
dim(data1)
head(data1[,1:6])
data2 <- a_2[,c("MOXD1","PDZK1IP1","SLC51B",data1$Id)]  
dim(data2)
head(data2[,1:5])
library(Hmisc)

CorMatrix <- function(cor,p) {
        ut <- upper.tri(cor) 
        data.frame(row = rownames(cor)[row(cor)[ut]] ,
                   column = rownames(cor)[col(cor)[ut]], 
                   cor =(cor)[ut], 
                   p = p[ut] )
}

res <- rcorr(as.matrix(data2),type = "pearson")
result_1 <- CorMatrix(res$r, res$P)
head(result_1)

dim(result_1)    
result_2 <- result_1[result_1$row == "MOXD1" | result_1$row == "PDZK1IP1"| result_1$row == "SLC51B"  ,]   
dim(result_2)
head(b_2)
b_2$column <- b_2$Id
head(b_2)
result_3 <- dplyr::inner_join(result_2,b_2,by="column")
dim(result_3)
result1 <- result_3[,1:4]
head(result1)
dim(result1)
result1$Regulation <- result1$cor
result1[,5][result1[,5] > 0] <- c("postive")
result1[,5][result1[,5] < 0] <- c("negative")
head(result1)
colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")

write.table(result1,file="corResult_Immunostimulator.xls",sep="\t",quote=F,col.names=T,row.names = F)

a1 <- read.table("corResult_Immunostimulator.xls",header = T,sep = "\t", quote = "",fill = T)
head(a1)

data2 <- a1
library(ggpubr)
data2$pvalue <- ifelse(data2$pvalue < 0.05,
                       ifelse(data2$pvalue < 0.01,"**","*"),
                       "")
data2$pvalue[1:20]
data2$type <- data2$cor
data2[,6][data2[,6] > 0 ] <- 2
data2[,6][data2[,6] < 0 ] <- 1
summary(data2)
data3 <- data2[order(data2$immuneGene,data2$cor),]
head(data3)
dim(data3)
data4 <- data3[data3$pvalue < 0.05,]
dim(data4)
summary(data4)
p <- ggplot(data4,aes(x=gene,y=immuneGene)) +
        geom_point(aes(colour = cor, size=pvalue)) +
        labs(x="",y="Disease regulatory genes")
p <- p + scale_colour_gradient2(low = "blue", high = "red", mid = "white",
                                midpoint = 0, limit = c(-1, 1), space = "Lab",
                                name="Pearson\nCorrelation")
p <- p + theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        theme(axis.text=element_text(size = 15)) +
        theme(axis.text.x=element_text(colour = "black",angle=0,hjust=0.5,size = 15)) +
        theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
        theme(axis.title =element_text(size = 20)) +
        theme(text = element_text(size = 15))

p+rotate_x_text(45)
ggsave("Immunostimulator.pdf")






library(ggplot2)
library(ggExtra)
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F,row.names = 1)
dat<-as.data.frame(t(rt))

corr_eqn <- function(x,y,digits=2) {
        test <- cor.test(x,y,type="pearson")
        paste(paste0("n = ",length(x)),
              paste0("r = ",round(test$estimate,digits),"(Pearson)"),
              paste0("p.value= ",round(test$p.value,digits)),
              sep = ", ")
}
gene<-as.numeric(dat$MOXD1) 
imucell<-dat$PSEN1  
corr_eqn(gene,imucell)

gg<-ggplot(dat, aes(x=gene, y=imucell)) + 
        geom_point(color = "black") + 
        geom_smooth(method="loess", se=F,color="blue") + 
      
        labs( 
                y="PSEN1",
                x="MOXD1", 
                title="Scatterplot")+ 
       
        labs(title = paste0(corr_eqn(gene,imucell)))+
        theme_bw()
gg
gg2 <- ggMarginal(gg, type="density")
gg2 <- ggMarginal(gg, type="density",xparams = list(fill ="orange"),
                  yparams = list(fill ="skyblue"))








workdir <- "C://"; setwd(workdir)
fig.path    <- file.path(workdir, "Figures")
if (!file.exists(fig.path)) { dir.create(fig.path) } 

library(pROC)

clin <- read.csv("临床症状.csv",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
meta <- read.csv("代谢物.csv",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
colnames(meta) <- gsub(",","_",colnames(meta)) 

samID <- intersect(rownames(clin),rownames(meta))
metaID <- colnames(meta)
clinID <- colnames(clin) 

auc.cutoff <- 0.5 
outTab <- NULL 
for (c in clinID) {
        for (m in metaID) {

                tmp <- data.frame(clin = clin[samID,c],
                                  meta = meta[samID,m],
                                  stringsAsFactors = F)

                tmp <- as.data.frame(na.omit(tmp))

                model_LR <- glm(clin ~ meta, 
                                data = tmp, 
                                family = "binomial")
                LR.sum <- summary(model_LR)

                pred_LR <- data.frame(prob = predict(model_LR, newdata = tmp,type="response"),
                                      group = tmp$clin,
                                      stringsAsFactors = F)

                LR.roc <- roc(pred_LR$group,pred_LR$prob)
                LR.roc.ci <- ci.auc(LR.roc) 
                outTab <- rbind.data.frame(outTab,
                                           data.frame(clin = c,
                                                      meta = m, 
                                                      auc = LR.roc.ci[2], 
                                                      low.ci = LR.roc.ci[1],
                                                      up.ci =LR.roc.ci[3], 
                                                      p.logit =  LR.sum$coefficients[2,4], 
                                                      stringsAsFactors = F),
                                           stringsAsFactors = F)

                if(auc(LR.roc) > auc.cutoff) {
                        pdf(file.path(fig.path,paste0(c," vs ",m," roc.pdf")),width = 6.5,height = 6.5)
                        LR.roc <- plot.roc(pred_LR$group,pred_LR$prob,ylim=c(0,1),xlim=c(1,0),
                                           smooth=F,
                                           ci=TRUE, 
                                           main=paste0("Performance of using ",m," to predict ",c),
                                           col="steelblue",
                                           lwd=2, 
                                           legacy.axes=T,
                                           print.auc = T)
                        invisible(dev.off())
                }
        }
}
write.csv(outTab,"outTab.csv",quote = F,row.names = F) 










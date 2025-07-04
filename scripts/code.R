rm(list = ls())

suppressMessages(library(devtools))
suppressMessages(library(customLayout))
suppressMessages(library(stringr))
suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(tidydr))
suppressMessages(library(openxlsx))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(fgsea))
suppressMessages(library(corrplot))
suppressMessages(library(colorspace))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(maftools))
suppressMessages(library(vegan))
suppressMessages(library(forcats))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(rstatix))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ggcor))
suppressMessages(library(ggstance))
suppressMessages(library(tidyverse))
suppressMessages(library(GOplot))
suppressMessages(library(caret))
suppressMessages(library(writexl))
suppressMessages(library(rcartocolor))
suppressMessages(library(ggcorrplot))
suppressMessages(library(psych))
suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))
suppressMessages(library(cols4all))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(scales))
suppressMessages(library(oncoPredict))
suppressMessages(library(gghalves))
suppressMessages(library(cowplot))
suppressMessages(library(IOBR))
suppressMessages(library(estimate))
suppressMessages(library(UpSetR))
suppressMessages(library(ggbiplot))
suppressMessages(library(ggsci))
options(stringsAsFactors = F)
########### TCGA 
tcga.exp=readMatrix('00_origin_datas/Preprocessed/tcga.exp.txt')
tcga_cli=readMatrix('00_origin_datas/Preprocessed/tcga.cli.txt')
identical(as.vector(tcga_cli$sampleID),as.vector(colnames(tcga.exp)))

genesets=read.table("00_origin_datas/GeneList.txt",sep="\t",header=F,as.is=T,quote="\"",fill=T,check.names = F,stringsAsFactors = F)
colnames(genesets)
#tcga.exp=tcga.exp[intersect(rownames(tcga.exp),rownames(icgc.t.exp)),]
#############ssgsea#####
tcga_exp=as.data.frame(tcga.exp)

ssgsea_list <- list()
ssgsea_list[['imm']] <- as.vector(unique(genesets$V1))
ssgsea_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga_exp,
                                                  genelist = ssgsea_list)
ssgsea_score <- as.data.frame(t(ssgsea_score))
ssgsea_score$"Group"=""
######### -01 
inds1=which(substr(rownames(ssgsea_score),14,15)=="01")
inds2=which(substr(rownames(ssgsea_score),14,15)=="11")

ssgsea_score$"Group"[inds1]="Case"
ssgsea_score$"Group"[inds2]="Control"
fig1a <- mg_violin(ssgsea_score[, c("Group", "imm")]
                   ,melt = T
                   ,xlab = ''
                   ,legend.pos = 'tl'
                   ,ylab = 'score')
fig1a
dev.off()
ggsave("PDFs/fig1a.pdf",fig1a, width = 8, height = 8)

######## #######

cox.pval=0.05
tcga_exp=tcga_exp[genesets$V1,inds1]
tcga_cli=tcga_cli[inds1,]
identical(as.vector(tcga_cli$sampleID),colnames(tcga_exp))
tcga.imm.cox=cox_batch(t(scale(t(as.matrix(tcga_exp))))
                       ,time = tcga_cli$OS.time/365
                       ,event = tcga_cli$OS)

table(tcga.imm.cox$p.value<0.05)
table(tcga.imm.cox$p.value<0.01)
table(tcga.imm.cox$p.value<0.001)
#FALSE  TRUE 
#1388   3157 
tcga.imm.cox=tcga.imm.cox[order(tcga.imm.cox$HR,decreasing = T),]
tcga.imm.cox.sig=tcga.imm.cox[which(tcga.imm.cox$p.value<cox.pval),]
nrow(tcga.imm.cox.sig)

write.table(tcga.imm.cox,'01_ConsensusClusterPlus/tcga.imm.cox.txt',sep = "\t",quote = F,row.names = T,col.names = T)
write.table(tcga.imm.cox.sig,'01_ConsensusClusterPlus/tcga.imm.cox.sig.txt',sep = "\t",quote = F,row.names = T,col.names = T)

##1.2 ##########
clusterAlg_name=c('hc','pam','km','kmdist')[3]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[1]
consen_gene=rownames(tcga.imm.cox.sig)
length(consen_gene)
tcga_consen_data=as.matrix(tcga_exp[intersect(consen_gene,rownames(tcga_exp)),])

tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   

#tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, median))
#tcga_consen_data=as.dist(1-cor(tcga_consen_data,method = 'spearman'))
#write.table(tcga_consen_data,'01_ConsensusClusterPlus/tcga_consen_data.txt',sep = "\t",quote = F,row.names = T,col.names = T)
#tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data, maxK = 10, reps = 500, pItem = 0.8, pFeature = 1, title = "01_ConsensusClusterPlus", clusterAlg = clusterAlg_name, distance = distance_name, plot = "pdf", writeTable = F, seed = 123456)#########
#save(tcga_clust_subtype,file='00_origin_datas/tcga.subtype0.05.RData')
load('00_origin_datas/tcga.subtype0.05.RData')#########


k=3
#subtype.cols1=rev(c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977", "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A", "#39737C", "#86B4A9", "#82853B", "#CCC94D"))
#subtype.cols=c("#A6CDE2","#F59899","#FCBF6E","#E11E26","#74C476","#86B4A9")
subtype.cols=c("#FF7F0F","#39737C","#32A251")
subtype.cols1=c("#F7A24F","#C6133B","#90162D","#93A5CB")
tcga.subtype <- data.frame( Samples=names(tcga_clust_subtype[[k]]$consensusClass),Subtype=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Subtype=paste0('C',tcga.subtype$Subtype)

table(tcga.subtype$Subtype)
colnames(tcga_cli)[1]='Samples'
tcga.subtype.cli=merge(tcga.subtype,tcga_cli,by='Samples')
write.table(tcga.subtype.cli[,1:2],'01_ConsensusClusterPlus/Subtype.txt',sep = "\t",quote = F,row.names = F,col.names = T)
fig1b=ggplotKMCox(data.frame(time = tcga.subtype.cli$OS.time/365
                             , event = tcga.subtype.cli$OS
                             , tcga.subtype.cli$Subtype)
                  ,add_text = '',show_confint = F,palette = subtype.cols)

fig1b
dev.off()
ggsave("PDFs/fig1b.pdf",fig1b, width = 8, height = 8)

##################PCA

tcga_exp_var=t(tcga_exp[,tcga.subtype.cli$Samples])

tcga_exp_var=tcga_exp_var[ , which(apply(tcga_exp_var, 2, var) != 0)]##
dim(tcga_exp_var)
cluster.pca <- prcomp(tcga_exp_var, scale=T)
cluster.pca.plot <- ggbiplot(cluster.pca, scale=1, groups = tcga.subtype.cli$Subtype,
                             ellipse = TRUE,ellipse.prob=0.3, circle = F,var.axes=F) +
  scale_color_manual(values = subtype.cols) + 
  theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top',
        panel.grid = element_blank(),text = element_text(family = 'Times')) +
  xlab('PCA1') + ylab('PCA2')+xlim(-3,3)+ylim(-3,3)
cluster.pca.plot
library(scatterplot3d)
cluster.pca1=as.data.frame(cluster.pca$x)
cluster.pca1$"Samples"=rownames(cluster.pca1)
cluster.pca2=merge(tcga.subtype.cli,cluster.pca1,by="Samples")
cluster.pca2=cluster.pca2[,c("Samples","Subtype","PC1","PC2","PC3")]
cluster.pca2$"color"=""
cluster.pca2$"color"[which(cluster.pca2$Subtype=="C1")]=subtype.cols[1]
cluster.pca2$"color"[which(cluster.pca2$Subtype=="C2")]=subtype.cols[2]
cluster.pca2$"color"[which(cluster.pca2$Subtype=="C3")]=subtype.cols[3]
dev.off()
pdf(file="PDFs/fig1c.pdf",width=8,height=8)

s3d <- scatterplot3d(cluster.pca2[,3:5],pch = 16
              ,color= cluster.pca2[,6],grid=TRUE#, box=FALSE
              ,angle=60, main= "3D PCA plot"
              ,cex.symbols= 0.6)
  legend(s3d$xyz.convert(40,60,30), title = "Subtype",
       xpd=TRUE,inset= -0.01,
       legend = unique(cluster.pca2$Subtype),
       col = subtype.cols , pch = 16)
dev.off()
##################### 

colnames(tcga.subtype.cli)
tcga.subtype.cli=tcga.subtype.cli[,-7]
head(tcga.subtype.cli)
tcga.subtype.cli.cmp=list()
for(i in c(3:7,11)){
  #group.color=subtype.cols[1:4]
  
  p=plotMutiBar_tmp(table(tcga.subtype.cli[,i],tcga.subtype.cli$Subtype)
                    ,fill.color = subtype.cols1
                    ,isAuto = F,showValue = F
                    ,legTitle=colnames(tcga.subtype.cli)[i])
  tcga.subtype.cli.cmp=c(tcga.subtype.cli.cmp,list(p$Bar))
}
length(tcga.subtype.cli.cmp)

fig1g=mg_merge_plot(tcga.subtype.cli.cmp,nrow = 2,ncol = 3
                    ,labels='D')
fig1g
dev.off()

savePDF('PDFs/fig1d.pdf',fig1g,height = 6,width=6)



#02##########################
inds1=which(substr(colnames(tcga.exp),14,15)=="01")
tcga.t.exp_use=tcga.exp[,inds1]
write.table(tcga.t.exp_use,file="00_origin_datas/tcga.t.exp_use.txt",sep="\t",quote=F,col.names = TRUE)
###TIDE
tcga.t.exp_TIDE=sweep(tcga.t.exp_use,1,apply(tcga.t.exp_use, 1, median))
write.table(tcga.t.exp_TIDE,'00_origin_datas/TIDE_exp.txt',sep = "\t",quote = F,row.names = T,col.names = T)

rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples
tcga.subtype.cli=tcga.subtype.cli[colnames(tcga.t.exp_use),]
identical(as.vector(tcga.subtype.cli$Samples),colnames(tcga.t.exp_use))

tcga.C1.geneList=getGeneFC_use(gene.exp=tcga.t.exp_use,group=tcga.subtype.cli$Subtype
                               ,ulab='C1',dlab = NULL)
tcga.C2.geneList=getGeneFC_use(gene.exp=tcga.t.exp_use,group=tcga.subtype.cli$Subtype
                               ,ulab='C2',dlab = NULL)
tcga.C3.geneList=getGeneFC_use(gene.exp=tcga.t.exp_use,group=tcga.subtype.cli$Subtype
                               ,ulab='C3',dlab = NULL)

kegmt<-read.gmt("00_origin_datas/h.all.v2023.2.Hs.entrez.gmt") #

tcga.C1.gsea<-GSEA(tcga.C1.geneList,TERM2GENE = kegmt,seed = 123456) #GSEA
tcga.C2.gsea<-GSEA(tcga.C2.geneList,TERM2GENE = kegmt,seed = 123456) #GSEA
tcga.C3.gsea<-GSEA(tcga.C3.geneList,TERM2GENE = kegmt,seed = 123456) #GSEA


tcga.C1.gsea.res=tcga.C1.gsea@result
tcga.C2.gsea.res=tcga.C2.gsea@result
tcga.C3.gsea.res=tcga.C3.gsea@result



mg_venn_plot(list(C1=rownames(tcga.C1.gsea.res)
                  , C2 = rownames(tcga.C2.gsea.res)
                  , C3 = rownames(tcga.C3.gsea.res)))

tcga.hallmark.union=Reduce(union,list(C1=rownames(tcga.C1.gsea.res)
                                      , C2 = rownames(tcga.C2.gsea.res)
                                      , C3 = rownames(tcga.C3.gsea.res)))


tcga.gsea.heatmap.dat=matrix(c(0),nrow =3,ncol = length(tcga.hallmark.union))
rownames(tcga.gsea.heatmap.dat)=c('C1', 'C2','C3')
colnames(tcga.gsea.heatmap.dat)=tcga.hallmark.union

tcga.gsea.heatmap.dat[1,match(rownames(tcga.C1.gsea.res),colnames(tcga.gsea.heatmap.dat))]=tcga.C1.gsea.res$NES
tcga.gsea.heatmap.dat[2,match(rownames(tcga.C2.gsea.res),colnames(tcga.gsea.heatmap.dat))]=tcga.C2.gsea.res$NES
tcga.gsea.heatmap.dat[3,match(rownames(tcga.C3.gsea.res),colnames(tcga.gsea.heatmap.dat))]=tcga.C3.gsea.res$NES


range(tcga.gsea.heatmap.dat)
tcga.gsea.heatmap.dat
colnames(tcga.gsea.heatmap.dat)=gsub("HALLMARK_","",colnames(tcga.gsea.heatmap.dat))

fig2a=Heatmap(as.matrix(t(scale(tcga.gsea.heatmap.dat)))
              , name = "NES"
              , rect_gp = gpar(col = "white")
              , cluster_rows = T
              , show_row_dend = F
              , cluster_columns = F
              , show_column_dend = F
              , show_column_names = T
              , show_row_names = T
              , row_names_gp = gpar(fontsize = 10)
              , row_names_side  = c("left")
              , col = circlize::colorRamp2(c(-3, 0, 3), c('navy', 'white', 'red'))
              , border = TRUE)

pdf(file="PDFs/fig2a.pdf",width=6,height=10)
plot(fig2a)
dev.off()


###################### 
dim(tcga.t.exp_use)

tcga.deg.c1=mg_limma_DEG_use(exp = tcga.t.exp_use,group=tcga.subtype.cli$Subtype,ulab = 'C1',dlab = NULL)
tcga.deg.c2=mg_limma_DEG_use(exp = tcga.t.exp_use,group=tcga.subtype.cli$Subtype,ulab = 'C2',dlab = NULL)
tcga.deg.c3=mg_limma_DEG_use(exp = tcga.t.exp_use,group=tcga.subtype.cli$Subtype,ulab = 'C3',dlab = NULL)


tcga.deg.c1.res=tcga.deg.c1$DEG
tcga.deg.c2.res=tcga.deg.c2$DEG
tcga.deg.c3.res=tcga.deg.c3$DEG

tcga.deg.c1.res=data.frame(tcga.deg.c1.res,gene=rownames(tcga.deg.c1.res),Group=rep("C1 vs Other",nrow(tcga.deg.c1.res)))
tcga.deg.c2.res=data.frame(tcga.deg.c2.res,gene=rownames(tcga.deg.c2.res),Group=rep("C2 vs Other",nrow(tcga.deg.c2.res)))
tcga.deg.c3.res=data.frame(tcga.deg.c3.res,gene=rownames(tcga.deg.c3.res),Group=rep("C3 vs Other",nrow(tcga.deg.c3.res)))

tcga.deg.res=rbind(tcga.deg.c1.res,tcga.deg.c2.res,tcga.deg.c3.res)

fig2b=mg_volcano_custom(diffData = tcga.deg.res
                        ,tile.col = as.character(subtype.cols)
                        ,log2FC.cutoff = log2(1.5),topGeneN = 5)
fig2b
dev.off()
pdf(file="PDFs/fig2b.pdf",width=8,height=8)
fig2b
dev.off()

tcga.deg.c1$Summary
tcga.deg.c2$Summary
tcga.deg.c3$Summary

get_DEG_p=function(df_deg,p.cutoff=0.05,logfc.cutoff=1){
  df.deg.res=df_deg$DEG
  df.deg.sig=df.deg.res[which(df.deg.res$P.Value<p.cutoff & abs(df.deg.res$logFC)>logfc.cutoff),]
}

tcga.deg.c1.sig=get_DEG(tcga.deg.c1,logfc.cutoff=log2(1.5),p.cutoff = 0.05)
tcga.deg.c2.sig=get_DEG(tcga.deg.c2,logfc.cutoff=log2(1.5),p.cutoff = 0.05)
tcga.deg.c3.sig=get_DEG(tcga.deg.c3,logfc.cutoff=log2(1.5),p.cutoff = 0.05)

deg.c1 = tcga.deg.c1$DEG
deg.c1$"gene"=rownames(deg.c1)
deg.c2 = tcga.deg.c2$DEG
deg.c2$"gene"=rownames(deg.c2)
deg.c3 = tcga.deg.c3$DEG
deg.c3$"gene"=rownames(deg.c3)

Deg_list <- list("deg.c1" = deg.c1[,c(7,1:6)], "deg.c2" = deg.c2[,c(7,1:6)],"deg.c3" = deg.c3[,c(7,1:6)])


write_xlsx(Deg_list, path = "02_Pathway/LIHC.deg.xlsx")






venn=list('C1 vs Other'=rownames(tcga.deg.c1.sig)
          ,'C2 vs Other'=rownames(tcga.deg.c2.sig)
          ,'C3 vs Other'=rownames(tcga.deg.c3.sig))

mg_venn_plot(venn,fill=subtype.cols)
dev.off()
tcga.deg.sig=Reduce(intersect ,venn)

length(tcga.deg.sig)#239
writeMatrix(tcga.deg.sig,'02_Pathway/tcga.deg.sig.txt',header=F)
pdf(file="PDFs/fig2c.pdf",width=6,height=6)
mg_venn_plot(venn,fill=subtype.cols)
dev.off()
####################RCircos################
#  ------------------------------------------------------
# 
tcga.deg.sig=read.table('02_Pathway/tcga.deg.sig.txt',header=F)
ann=readMatrix
gene=as.data.frame(ann)
gene<-gene[gene$gene_type=="protein_coding",]
pick_info<-c("seqnames","start","end","gene_name")
gtf_df_pick<-dplyr::distinct(gene[,pick_info]) 
gtf_sig <- gtf_df_pick[which(gtf_df_pick$gene_name%in%tcga.deg.sig$V1),]
dim(gtf_sig)
dim(tcga.deg.sig)
colnames(gtf_sig)=c("Chr","Start","End","Gene")
head(gtf_sig)
deg.c1_ann=deg.c1[tcga.deg.sig$V1,c(7,1)]
deg.c2_ann=deg.c2[tcga.deg.sig$V1,c(7,1)]
deg.c3_ann=deg.c3[tcga.deg.sig$V1,c(7,1)]

colnames(deg.c1_ann)=c("Gene","C1")
colnames(deg.c2_ann)=c("Gene","C2")
colnames(deg.c3_ann)=c("Gene","C3")
data_list <- list(deg.c1_ann, deg.c2_ann, deg.c3_ann,gtf_sig)

my_merge <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "Gene")
}

deg.ann=Reduce(my_merge, data_list) 
deg.ann=deg.ann[,c(5:7,1:4)]
writeMatrix(deg.ann,'02_Pathway/deg.sig.txt',header=T)

library(RCircos)
# 
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
# 
chr.exclude <- NULL
# 
tracks.inside <- 7
# 
tracks.outside <-0
# 
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)  
# 
RCircos.List.Plot.Parameters()
# 
RCircos.Set.Plot.Area()    
RCircos.Chromosome.Ideogram.Plot()
#####
params=RCircos.Get.Plot.Parameters()
#params$hist.color="orchid"###
params$scatter.color="lightblue"###
params$track.background="white"
RCircos.Reset.Plot.Parameters(params)
#####


side <- "in";
track.num <- 1;
RCircos.Gene.Connector.Plot(deg.ann, track.num, side);
name.col <- 4;
track.num <- 2;
RCircos.Gene.Name.Plot(deg.ann, name.col,track.num, side);

side <- "in";
params=RCircos.Get.Plot.Parameters()
params$line.color="orange"###
params$track.height=0.2
RCircos.Reset.Plot.Parameters(params)
RCircos.Line.Plot(deg.ann, 5,3, side)

params=RCircos.Get.Plot.Parameters()
params$line.color="cyan"###
RCircos.Reset.Plot.Parameters(params)
RCircos.Line.Plot(deg.ann, 6,4, side)

params=RCircos.Get.Plot.Parameters()
params$line.color="green"###
RCircos.Reset.Plot.Parameters(params)
RCircos.Line.Plot(deg.ann, 7,5, side)



####################RCircos##############
library(UpSetR)


venn=list('C1 vs Other'=rownames(tcga.deg.c1.sig)
          ,'C2 vs Other'=rownames(tcga.deg.c2.sig)
          ,'C3 vs Other'=rownames(tcga.deg.c3.sig))
upset_plot=upset(fromList(venn),
            sets = c("C1 vs Other", "C2 vs Other", "C3 vs Other"),
            main.bar.color = "black",
            sets.bar.color = subtype.cols)
dev.off()
pdf('PDFs/upset_plot.pdf',height = 8,width = 6)
upset_plot
dev.off()


###################
enrichment=mg_clusterProfiler(as.vector(tcga.deg.sig))
write.table(enrichment$Enrich_tab,file = '02_Pathway/enrichment.txt',sep = '\t',quote = F,row.names = T,col.names = T)

enrichment$KEGG@result$p.adjust=signif(enrichment$KEGG@result$p.adjust,digits = 3)
kegg_plot=barplot(enrichment$KEGG, showCategory=10)+theme(axis.text.y = element_text(size = 10)) 


#############
#enrich = enrichment$Enrich_tab 
#eKEGG <- enrich %>% 
#  filter(DB=="pathway_KEGG")%>%
#  filter(row_number() >= 1,row_number() <= 10) 
#eKEGG$DB="KEGG"
#eKEGG <- eKEGG[order(eKEGG$pValue),]
#kegg_plot<- ggplot(eKEGG,aes(x= enrichmentRatio ,y=reorder(description,enrichmentRatio) )) + 
#  geom_point(aes(size=size  ,color=-1*log10(pValue))) +
#  scale_color_gradient(low="#0072B5",high ="#BC3C29") + 
#  labs(color=expression(-log[10](pValue)),size="size", shape="DB",
#       x="EnrichmentRatio",y="",title="KEGG enrichment") + 
#  theme_bw()+scale_size_continuous(range=c(5,8))
#kegg_plot
#dev.off()
savePDF('PDFs/kegg_plot.pdf',kegg_plot,height = 4,width=8)
#########03######

#tcga.deg.sig=readMatrix('02_Pathway/string_interactions_short.txt default node.txt',header=F)
tcga.deg.sig=readMatrix('02_Pathway/tcga.deg.sig.txt',header=F)
colnames(tcga.deg.sig)="V1"
tcga.deg.sig=tcga.deg.sig$V1
length(tcga.deg.sig)
identical(as.vector(tcga.subtype.cli$Samples),colnames(tcga.t.exp_use))
tcga.cox=cox_batch(t(scale(t(tcga.t.exp_use[tcga.deg.sig,])))
                   ,time =  tcga.subtype.cli$OS.time/365
                   ,event =tcga.subtype.cli$OS)
dim(tcga.cox)

table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)
writeMatrix(tcga.cox,outpath = '03_Model/tcga.cox.txt')



p.cutoff=0.05
tcga.cox_use=tcga.cox
tcga.cox_use$coef=log(tcga.cox_use$HR)
tcga.cox_use$Gene=rownames(tcga.cox_use)
tcga.cox_use$type=rep('None',nrow(tcga.cox_use))
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef>0)]='Risk'
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef<0)]='Protective'
table(tcga.cox_use$type)

######### lasso
tcga.gene.sig=rownames(tcga.cox)[which(tcga.cox$p.value<p.cutoff)]
length(tcga.gene.sig)

table(tcga.cox_use$type)


tcga.cox_forVis=tcga.cox_use
tcga.cox_forVis=tcga.cox_forVis[which(tcga.cox_forVis$type %in% c('Risk','Protective')),]
tcga.cox_forVis$p.value=-log10(tcga.cox_forVis$p.value)
range(tcga.cox_forVis$p.value)


#################### LASSO
table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)

tcga.exp.sig=tcga.t.exp_use[tcga.gene.sig,]
tcga.exp.sig=t(tcga.exp.sig)
dim(tcga.exp.sig)


dim(tcga.exp.sig)
options(ggrepel.max.hnscerlaps = Inf)
tcga.subtype.cli$Samples=as.vector(tcga.subtype.cli$Samples)
identical(rownames(tcga.exp.sig),tcga.subtype.cli$Samples)
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples
tcga.lasso.res=mg_lasso_cox_use(t((t(tcga.exp.sig)))
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , nfolds = 10
                                , lambda.min = T
                                , figLabels=c('B','C'))
tcga.lasso.res$Genes

tcga.lasso.res$lambda

tcga.lasso.res$plot

tcga.exp.for.cox=tcga.t.exp_use[match(tcga.lasso.res$Genes,row.names(tcga.t.exp_use)),]
dim(tcga.exp.for.cox)
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)

lst.modl=createCoxModel_use((t(tcga.exp.for.cox))
                            , time = tcga.subtype.cli$OS.time/365
                            , event = tcga.subtype.cli$OS
                            , isStep = T)
lst.modl$Cox
lst.modl$Genes
lst.modl$fmla

lst.modl.Coef=lst.modl$Coef
names(lst.modl.Coef)=lst.modl$Genes
lst.modl.Coef

tcga.risk.score=lst.modl$Score
#tcga.risk.score=scale(tcga.risk.score)[,1]
tcga.risk.score=mosaic::zscore(tcga.risk.score)

range(tcga.risk.score)

lst.modl$Coef

gene.coef=data.frame(Gene=lst.modl$Genes,Coef=lst.modl$Coef)
gene.coef$Type=ifelse(lst.modl$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
table(gene.coef$Type)

fig3c=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  geom_text(aes(label=round(Coef,digits = 3)),color="black",hjust = "left")+
  
  coord_flip() +
  scale_fill_manual(values=pal_nejm(alpha = 0.9)(8)[c(1,4)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Lasso Cox coefficient") +
  theme_classic()+theme(legend.position = c(0,1))
# theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")


fig3AB=mg_plot_lasso_use(fit = tcga.lasso.res$Mode1
                         , cv_fit = tcga.lasso.res$Model2
                         , show_text = F
                         , figLabels = c('A', 'B'))
fig3AB
fig3abc=mg_merge_plot(fig3AB,fig3c,nrow = 1,ncol = 2,widths = c(2,1),labels = c('','C'))
#savePDF('PDFs/Fig7AB.pdf',fig7A,height = 4,width = 9)
savePDF('PDFs/fig3abc.pdf',fig3abc,height = 5,width = 15)


tcga.exp.forCox<- cbind(time=tcga.subtype.cli$OS.time/365,
                        status=tcga.subtype.cli$OS,
                        t(tcga.t.exp_use)[rownames(tcga.subtype.cli), lst.modl$Genes])



dim(tcga.exp.forCox)

fmla <- as.formula(paste0("Surv(time, status) ~",paste0(lst.modl$Genes,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga.exp.forCox))
fig3d=survminer::ggforest(cox,data=tcga.exp.forCox,noDigits = 3)
fig3abcd=mg_merge_plot(fig3abc,fig3d,nrow = 1,ncol = 2,widths = c(3,1),labels = c('','D'))
savePDF('PDFs/fig3abcd.pdf',fig3abcd,height = 5,width = 16)

fig3abd=mg_merge_plot(fig3AB,fig3d,nrow = 1,ncol = 2,widths = c(2,1),labels = c('','C'))

############### TCGA
#tcga.cutoff <- survminer::surv_cutpoint(data.frame(time=tcga.subtype.cli$OS.time/365,
 #   event=tcga.subtype.cli$OS,
 #  risk=tcga.risk.score), time = "time", event = "event",
#variables = c("risk"))
#tcga.cutoff=tcga.cutoff$cutpoint$cutpoint
tcga.cutoff=0
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)
risk.group.color=c("#E11E26","#1E78B4")
names(risk.group.color)=c('High','Low')
fig3ee=plotRiskScoreModel_use(riskScore = tcga.risk.score
                              ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                              , time = tcga.subtype.cli$OS.time/365
                              , event = tcga.subtype.cli$OS
                              , cutoff = tcga.cutoff
                              , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                              , pal = risk.group.color)
fig3e1=fig3ee[[2]]

fig3e1
dev.off()
fig3ee2=plotCoxModel_Batch_use(riskScore = tcga.risk.score
                               ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                               , time = tcga.subtype.cli$OS.time/365
                               , event = tcga.subtype.cli$OS
                               , cutoff = tcga.cutoff
                               , labs = c('High','Low')
                               , title = 'RiskType'
                               , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                               , pal = risk.group.color
                               , mks = c(1:5))
fig3e2=fig3ee2[[3]]
#pdf('PDFs/fig3e2.pdf',height = 6,width = 6)
fig3e2
dev.off()

tcga.group=ifelse(tcga.risk.score>tcga.cutoff,'High','Low')
tcga.group=data.frame(tcga.group)
colnames(tcga.group)='group'
table(tcga.group$group)
fig3e=mg_merge_plot(fig3e1,fig3e2,nrow = 1,ncol = 2,widths = c(1,2))
write.table(cbind(tcga.risk.score,tcga.group),file = '03_Model/tcga.group.txt',sep='\t',quote = F)

###########################
icgc.t.exp <- read.table("00_origin_datas/GEO_expression_GSE31210.txt",header = T,check.names = F,fill=T,sep = "\t")
icgc.t.cli <- read.table("00_origin_datas/GEO_cli_GSE31210.txt",header = T,check.names = F,fill=T,sep = "\t")
identical(colnames(icgc.t.exp),icgc.t.cli$Samples)
icgc.HCCDB18.t.exp=icgc.t.exp
match(lst.modl$Genes,row.names(icgc.HCCDB18.t.exp))
length(lst.modl$Genes)
icgc.HCCDB18.t.cli.os=icgc.t.cli
identical(icgc.HCCDB18.t.cli.os$Samples , colnames(icgc.HCCDB18.t.exp)) ####
icgc.HCCDB18.model.dat=icgc.HCCDB18.t.exp[match(lst.modl$Genes,row.names(icgc.HCCDB18.t.exp)),]

icgc.HCCDB18.risk.score=predictScoreByCoxModel(coxModel = lst.modl
                                               ,(t(icgc.HCCDB18.model.dat)))
#icgc.HCCDB18.risk.score=scale(icgc.HCCDB18.risk.score)
icgc.HCCDB18.risk.score=mosaic::zscore(icgc.HCCDB18.risk.score)

lst.modl$fmla
#lst.vd.mod1$fmla

#icgc.HCCDB18.cutoff <- survminer::surv_cutpoint(data.frame(time=as.numeric(icgc.HCCDB18.t.cli.os$OS1),
# event=icgc.HCCDB18.t.cli.os$os_status1,
# risk=icgc.HCCDB18.risk.score), time =                                                   "time", event = "event", variables = c("risk"))
#icgc.HCCDB18.cutoff=icgc.HCCDB18.cutoff$cutpoint$cutpoint
icgc.HCCDB18.cutoff=0

fig3gg1=plotRiskScoreModel_use(riskScore = icgc.HCCDB18.risk.score
                               , dat = t(icgc.HCCDB18.t.exp[intersect(lst.modl$Genes, row.names(icgc.HCCDB18.t.exp)),])
                               , time = icgc.HCCDB18.t.cli.os$OS.time/365
                               , event = icgc.HCCDB18.t.cli.os$OS
                               , cutoff = icgc.HCCDB18.cutoff
                               , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                               , pal = risk.group.color)
fig3g1=fig3gg1[[2]]
#pdf('PDFs/fig3g1.pdf',height = 6,width = 6)

fig3g1
dev.off()
fig3gg2=plotCoxModel_Batch_use(riskScore = icgc.HCCDB18.risk.score
                               , dat = t(icgc.HCCDB18.t.exp[intersect(lst.modl$Genes, row.names(icgc.HCCDB18.t.exp)),])
                               , time = as.numeric(icgc.HCCDB18.t.cli.os$OS.time/365) 
                               , event = as.numeric(icgc.HCCDB18.t.cli.os$OS)
                               , cutoff = icgc.HCCDB18.cutoff
                               , labs = c('High','Low')
                               , title = 'RiskType'
                               , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                               , pal = risk.group.color
                               , mks = c(1:5))
fig3g2=fig3gg2[[3]]
#pdf('PDFs/fig3g2.pdf',height = 6,width = 6)

fig3g2
dev.off()
icgc.HCCDB18.group=ifelse(icgc.HCCDB18.risk.score>icgc.HCCDB18.cutoff,'High','Low')
icgc.HCCDB18.group=data.frame(icgc.HCCDB18.group)
colnames(icgc.HCCDB18.group)='group'
table(icgc.HCCDB18.group)
fig3g=mg_merge_plot(fig3g1,fig3g2,nrow = 1,ncol = 2,widths = c(1,2))

write.table(cbind(icgc.HCCDB18.risk.score,icgc.HCCDB18.group),file = '03_Model//geo.group.txt',sep='\t',quote = F)


fig3eg=ggpubr::ggarrange(fig3e,fig3g, ncol = 1, nrow = 2)
fig3=ggpubr::ggarrange(fig3abd,fig3eg, ncol = 1, nrow = 2,heights = c(1,3))

savePDF('PDFs/fig3.pdf',fig3,height = 12,width = 16)

#######################

##################04#######
identical(rownames(as.data.frame(tcga.risk.score)),rownames(tcga.subtype.cli))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=tcga.risk.score)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')

Risktype.color=risk.group.color


##########
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)[c(3:7,11:13)]

table(tcga_cox_datas$pathologic_N)
tcga_cox_datas$pathologic_N[tcga_cox_datas$pathologic_N=='N0'|tcga_cox_datas$pathologic_N=='N1']<-'N0+N1'
tcga_cox_datas$pathologic_N[tcga_cox_datas$pathologic_N=='N2'|tcga_cox_datas$pathologic_N=='N3']<-'N2+N3'


table(tcga_cox_datas$pathologic_T)
tcga_cox_datas$pathologic_T[tcga_cox_datas$pathologic_T=='T1'|tcga_cox_datas$pathologic_T=='T2']<-'T1+T2'
tcga_cox_datas$pathologic_T[tcga_cox_datas$pathologic_T=='T3'|tcga_cox_datas$pathologic_T=='T4']<-'T3+T4'


table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='stge i'|tcga_cox_datas$Stage=='stge ii']<-'Stage i+ii'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='stge iii'|tcga_cox_datas$Stage=='stge iv']<-'Stage iii+iv'
table(tcga_cox_datas$Stage)


## 
univar_res<-unicox(vars=colnames(tcga_cox_datas)[c(3:7,11:13)],time = tcga_cox_datas$OS.time,event = tcga_cox_datas$OS,data=tcga_cox_datas)
univar_res
univar_res[which(univar_res$pvalue<0.05),]


## 
mutivar_res<-multicox(vars=rownames(univar_res[which(univar_res$pvalue<0.01),]),time = tcga_cox_datas$OS.time,event = tcga_cox_datas$OS,data=tcga_cox_datas,forest = F)
mutivar_res
mutivar_res[which(mutivar_res$pvalue<0.05),]


##
#dev.off()
mg_Forestplot(df_m = univar_res,outFile = 'PDFs/tcga.univar.forestplot.pdf',height = 4,width = 6)
mg_Forestplot(df_m = mutivar_res,outFile = 'PDFs/tcga.mutivar.forestplot.pdf',height = 4,width = 6)


###################nomo

dt=data.frame(RiskScore=tcga_cox_datas$Riskscore,
              pathologic_T=tcga_cox_datas$pathologic_T
            )

pdf('PDFs/fig4a.pdf', width = 12, height = 10)

nom.plot=mg_nomogram(clinical_riskscore=dt,
                     os = as.numeric(tcga_cox_datas$OS.time/365),
                     status = as.numeric(tcga_cox_datas$OS),
                     mks = c(1:5))


dev.off()

#pdf('PDFs/fig4b.pdf', width = 9, height = 6)
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1:5))
dev.off()
########################


####05 #############

#### ESTIMATE
#tcga.exp.estimate<-deconvo_estimate(eset=tcga.t.exp_use)
#save(tcga.exp.estimate,file='05_immune_drug/tcga.exp.estimate.RData')
#load('05_immune_drug/tcga.exp.estimate.RData')
#tcga.exp.estimate=get.IOBR.immu.format(tcga.exp.estimate)

### CIBERSORT
#tcga.exp.cibersort<-deconvo_cibersort(eset=tcga.t.exp_use,arrays=T)
#save(tcga.exp.cibersort,file='05_immune_drug/tcga.exp.cibersort.RData')
#load('05_immune_drug/tcga.exp.cibersort.RData')
#tcga.exp.cibersort=get.IOBR.immu.format(tcga.exp.cibersort)
############ TIMER 
tcga.exp.timer<-deconvo_timer(eset=as.matrix(tcga.t.exp_use),indications=rep('LUAD',ncol(tcga.t.exp_use)))
save(tcga.exp.timer,file='05_immune_drug/tcga.exp.timer.RData')
load('05_immune_drug/tcga.exp.timer.RData')
tcga.exp.timer=get.IOBR.immu.format(tcga.exp.timer)

############ MCP-counter 
tcga.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(tcga.t.exp_use))
save(tcga.exp.mcp,file='05_immune_drug/tcga.exp.mcp.RData')
load('05_immune_drug/tcga.exp.mcp.RData')
tcga.exp.mcp=get.IOBR.immu.format(tcga.exp.mcp)

#
#tcga.t.estimate=tcga.exp.estimate[rownames(tcga.risktype.cli),1:4]
#tcga.t.cibersort=tcga.exp.cibersort[rownames(tcga.risktype.cli),1:22]
tcga.t.timer=tcga.exp.timer[rownames(tcga.risktype.cli),]
tcga.t.mcp=tcga.exp.mcp[rownames(tcga.risktype.cli),]

#write.table(tcga.t.cibersort,'05_immune_drug/tcga.t.cibersort.txt',sep = "\t",quote = F,row.names = T,col.names = T)
#write.table(tcga.t.estimate,'05_immune_drug/tcga.t.estimate.txt',sep = "\t",quote = F,row.names = T,col.names = T)
write.table(tcga.t.mcp,'05_immune_drug/tcga.t.mcp.txt',sep = "\t",quote = F,row.names = T,col.names = T)
write.table(tcga.t.timer,'05_immune_drug/tcga.t.timer.txt',sep = "\t",quote = F,row.names = T,col.names = T)





fig5c=get_PlotMutiBoxplot(tcga.exp.mcp,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'MCP Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5c

fig5d=get_PlotMutiBoxplot(tcga.t.timer,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'TIMER Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5d

############ MCP-counter 
cr.mcp=psych::corr.test(x=tcga.risktype.cli[,'Riskscore'],
                    y=tcga.t.mcp[tcga.risktype.cli$Samples,]
                    ,method = 'spearman')
df.mcp=t(rbind(cr.mcp$r,cr.mcp$p))
colnames(df.mcp)=c('r','p.value')
df.mcp=data.frame(Riskscore='Riskscore',TIMER_count=rownames(df.mcp),df.mcp)
df.mcp
df.mcp <- df.mcp %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df.mcp)
#corrmat.color=colorRampPalette(c('#0072B5FF', '#BC3C29FF'))(100)
corrmat.color <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988","#BB4444" ))(5)

cor.mcp<-quickcor(tcga.t.mcp[tcga.risktype.cli$Samples,], cor.test = T,type = "lower") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "lower")) + 
  anno_link(df.mcp, mapping = aes(colour = col,
                              size = abs(r),
                              linetype = lty)) +
  scale_fill_gradientn(colours = corrmat.color,
                       values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25)) +
  remove_x_axis()+
  scale_size_area(max_size = 1) +
  scale_linetype_manual(values = c( "solid","dotted")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 3),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 1),
    linetype = "none")
savePDF('PDFs/cor.mcp.pdf',cor.mcp,height = 8,width = 8)

#cor.mcp
############ timer
cr.timer=psych::corr.test(x=tcga.risktype.cli[,'Riskscore'],
                          y=tcga.t.timer[tcga.risktype.cli$Samples,]
                          ,method = 'spearman')
df.timer=t(rbind(cr.timer$r,cr.timer$p))
colnames(df.timer)=c('r','p.value')
df.timer=data.frame(Riskscore='Riskscore',TIMER_count=rownames(df.timer),df.timer)
df.timer
df.timer <- df.timer %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df.timer)

cor.timer<-quickcor(tcga.t.timer[tcga.risktype.cli$Samples,], cor.test = T,type = "lower") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "lower")) + 
  anno_link(df.timer, mapping = aes(colour = col,
                                    size = abs(r),
                                    linetype = lty)) +
  scale_fill_gradientn(colours = corrmat.color,
                       values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25)) +
  remove_x_axis()+
  scale_size_area(max_size = 1) +
  scale_linetype_manual(values = c( "solid","dotted")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 3),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 1),
    linetype = "none")

savePDF('PDFs/cor.timer.pdf',cor.timer,height = 8,width = 8)





# ICGs####

tcga.icgs=immu_ICGs(tcga.t.exp_use)
colnames(tcga.icgs)
select_col=c('CTLA4','PDCD1',"TIGIT",'LGALS9','CD276','TNFSF4','TNFSF9')
icg.dat.RS=cbind(tcga.risktype.cli$Riskscore
                 ,t(tcga.exp)[tcga.risktype.cli$Samples,lst.modl$Genes]
                 ,tcga.icgs[tcga.risktype.cli$Samples,c(select_col)])
#c('CTLA4','PDCD1','PDCD1LG2',"TIGIT",'LGALS9','CD80','CD28','HAVCR2')
colnames(icg.dat.RS)[1]='Riskcsore'

icg_cor_res <- Hmisc::rcorr(as.matrix(icg.dat.RS),type = 'spearman')
icg_cor_res$P[is.na(icg_cor_res$P)] <- 0
icg_cor_res.p=icg_cor_res$P
icg_cor_res.p[1:5,1:5]
icg_cor_res.p<-ifelse(icg_cor_res.p<0.0001,'****',
                      ifelse(icg_cor_res.p<0.001,'***', 
                             ifelse(icg_cor_res.p<0.01,'**',
                                    ifelse(icg_cor_res.p<0.05,'*',''))))
ICG_plot <- pheatmap(icg_cor_res$r[-c(1:(length(lst.modl$Genes)+1)),c(lst.modl$Genes,'Riskcsore')],
         color = circlize::colorRamp2(c(-1,-0.5, 0,0.5, 1), c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988","#BB4444" )),
         main="Heatmap", # 
         display_numbers = icg_cor_res.p[-c(1:(length(lst.modl$Genes)+1)),c(lst.modl$Genes,'Riskcsore')], # 
         annotation_legend = FALSE ,
         cluster_cols = F, # 
         cluster_rows = F,
         show_rownames = T, #
         show_colnames = T,
         fontsize_row = 12, # 
         fontsize_col = 16)
dev.off()
pdf('PDFs/ICG_plot.pdf',height = 6,width = 6)
ICG_plot

dev.off()

#####TIDE######
tcga.tide<-read.csv('00_origin_datas/TIDE.csv',row.names = 1,stringsAsFactors = F)
tcga.tide=tcga.tide[rownames(tcga.risktype.cli),]
dim(tcga.tide)

tide.selected=c('Exclusion','Dysfunction','TIDE')

tcga.subtype.tide.p.all=c()
for(i in tide.selected){
  df_half=data.frame("RiskType"=tcga.risktype.cli$Risktype,"TIDE_score"=tcga.tide[rownames(tcga.risktype.cli),i])
  p <-  ggplot(data = df_half, aes(x = RiskType, y = TIDE_score, fill = RiskType)) +
    geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 1) +
    geom_point(aes(y = TIDE_score, color = RiskType), position = position_jitter(width = 0.15), size = 1, alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 1) +
    labs(y = i
         , x = 'RiskType') +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values = risk.group.color) +
    scale_colour_manual(values = risk.group.color) +theme_classic2()+
  ggpubr::stat_compare_means(comparisons = list(c(1,2)),method = 'wilcox.test',label= "p.signif")
  
  tcga.subtype.tide.p.all=c(tcga.subtype.tide.p.all,list(p))
}
length(tcga.subtype.tide.p.all)

TIDE_plot=mg_merge_plot(tcga.subtype.tide.p.all,nrow = 1,ncol = length(tide.selected),common.legend=T,legend =  "left")
TIDE_plot
savePDF('PDFs/TIDE_plot.pdf',TIDE_plot,height = 5,width = 6)

fig5ad=mg_merge_plot(fig5c,cor.mcp,fig5d,cor.timer,nrow = 2,ncol = 2,widths = c(1.5,1),labels = LETTERS[1:4])

savePDF('PDFs/fig5ad.pdf',fig5ad,height = 8,width = 9)
####
# 05######


#drug_exp=as.matrix(tcga.t.exp_use)

#GDSC2_Expr = readRDS(file=file.path(dir,'Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
#GDSC2_Res = readRDS(file = file.path(dir,"Training Data/GDSC2_Res.rds"))
#GDSC2_Res <- exp(GDSC2_Res)
#calcPhenotype(trainingExprData = as.matrix(GDSC2_Expr),
#              trainingPtype = as.matrix(GDSC2_Res),
#              testExprData = as.matrix(drug_exp),
#              batchCorrect = 'eb',  #   "eb" for ComBat
#              powerTransformPhenotype = TRUE,
#              removeLowVaryingGenes = 0.2,
#              minNumSamples = 10,
#              printOutput = TRUE,
#              removeLowVaringGenesFrom = 'rawData' )
tcga_durg_ic50_res=read.csv('05_immune_drug/calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
dim(tcga_durg_ic50_res)
tcga_durg_ic50_res[1:5,1:5]

IC50.mat=data.frame(Riskscore=tcga.risktype.cli$Riskscore,tcga_durg_ic50_res[tcga.risktype.cli$Samples,])

IC50_RS_cor <- corr.test(x =IC50.mat$Riskscore,
                         y = IC50.mat[,-1],
                         method = "spearman",adjust = "BH",ci = F)


IC50_RS_cor_res=data.frame(drugs=colnames( IC50.mat[,-1]))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
head(IC50_RS_cor_res)
table(IC50_RS_cor_res$p.adj<0.05,abs(IC50_RS_cor_res$cor)>0.3)
IC50_RS_cor_res=IC50_RS_cor_res[IC50_RS_cor_res$p.adj<0.05 & abs(IC50_RS_cor_res$cor)>0.3,]
IC50_RS_cor_res=IC50_RS_cor_res[order(IC50_RS_cor_res$cor),]
head(IC50_RS_cor_res)
IC50_plot <- ggplot(data=IC50_RS_cor_res,aes(x=cor,y=reorder(drugs,cor),
                                             color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_colour_gradient(low ='#ffc7c7' ,high = "#8785a2")+
  geom_segment(aes(yend=drugs,xend=0),size=1) +
  labs(x='spearman Correlation',y='')+theme_bw()+
  theme(text = element_text(family = 'Times'),legend.position = "bottom")

pdf('PDFs/IC50_RS_cor.pdf',height = 9,width = 4)
IC50_plot
dev.off()

######################
#######################################
########

tcga_tmb <- mg_getTCGATMBByCode('LUAD')
tcga_tmb <- as.data.frame(tcga_tmb)
head(tcga_tmb)
tcga_tmb$Sample <- paste0(tcga_tmb$Sample, '-01')
tcga.risktype.cli$Samples=substr(tcga.risktype.cli$Samples,0,15)
#colnames(tcga.risktype.cli)[1]="Sample"
tcga_tmb1 <- merge(tcga_tmb, tcga.risktype.cli[, c("Samples", "Risktype",'Riskscore')],
                  by.x = 'Sample', by.y = 'Samples')
head(tcga_tmb1)
tcga_tmb1$logTMB <- log2(tcga_tmb1$TMB + 1)


#tcga_tmb_plot <- wb_beeswarm_plot(tcga_tmb1[, c("Risktype", "logTMB")],
#                                  ylab = 'log2(Tumor mutation burden)',
#                                  show_compare = T,
#                                  xlab = "Risktype",
#                                  title = 'TCGA',
#                                  col = risk.group.color)

tcga_tmb_plot <- mg_violin(tcga_tmb1[, c("Risktype", "logTMB")]
                 ,melt = T
                 ,xlab = 'Risktype'
                 ,legend.pos = 'tl'
                 ,ylab = 'log2(Tumor mutation burden)')

tcga_tmb_plot
ggsave(plot = tcga_tmb_plot,
       filename = 'PDFs/tcga_tmb_plot.pdf',
       width = 5, height = 5)


##########
tcga.maf=getTCGAMAFByCode('LUAD')#
tcga.risktype.use=tcga.risktype.cli[,c('Samples','Risktype')]
table(tcga.risktype.use$Risktype)
colnames(tcga.risktype.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.use$Tumor_Sample_Barcode=substr(tcga.risktype.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.use.high=tcga.risktype.use[which(tcga.risktype.use$Risktype=='High'),]
tcga.risktype.use.low=tcga.risktype.use[which(tcga.risktype.use$Risktype=='Low'),]

write.table(tcga.risktype.use.high,file='06_risktype.mut/tcga.risktype.use.high.txt',row.names = F)
write.table(tcga.risktype.use.low,file='06_risktype.mut/tcga.risktype.use.low.txt',row.names = F)

tcga.maf.high=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.use.high$Tumor_Sample_Barcode))
tcga.maf.high<-read.maf(tcga.maf.high@data,isTCGA=T,clinicalData = '06_risktype.mut/tcga.risktype.use.high.txt')
tcga.maf.high@clinical.data

tcga.maf.low=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.use.low$Tumor_Sample_Barcode))
tcga.maf.low<-read.maf(tcga.maf.low@data,isTCGA=T,clinicalData = '06_risktype.mut/tcga.risktype.use.low.txt')
tcga.maf.low@clinical.data
dev.off()
#######
pdf('PDFs/tcga.maf.high.pdf',height = 6,width = 7,onefile = F)
oncoplot(maf = tcga.maf.high,top = 20,sortByAnnotation = T)
dev.off()
pdf('PDFs/tcga.maf.low.pdf',height = 6,width = 7,onefile = F)
oncoplot(maf = tcga.maf.low,top =20,sortByAnnotation = T)
dev.off()

######
tcga.geneList=getGeneFC(gene.exp=tcga.t.exp_use[,rownames(tcga.risktype.cli)], group=tcga.risktype.cli$Risktype,ulab='High',dlab='Low')
saveRDS(tcga.geneList,file = "06_risktype.mut/tcga.geneList.RDS")


tcga.geneList <- readRDS("06_risktype.mut/tcga.geneList.RDS")
library(clusterProfiler)
library(tidygraph)
library(cols4all)
tcga.gsea.KEGG=gseKEGG(
  tcga.geneList,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea")
saveRDS(tcga.gsea.KEGG,file = "06_risktype.mut/tcga.gsea.KEGG.RDS")
tcga.gsea.KEGG <- readRDS("06_risktype.mut/tcga.gsea.KEGG.RDS")

tcga.gsea.KEGG.res=tcga.gsea.KEGG@result
write.csv(tcga.gsea.KEGG.res,'06_risktype.mut/GSEA_res.csv',row.names = F)
table(tcga.gsea.KEGG.res$p.adjust<0.05 & tcga.gsea.KEGG.res$NES<0)
table(tcga.gsea.KEGG.res$p.adjust<0.05 & tcga.gsea.KEGG.res$NES>0)

ind1=tcga.gsea.KEGG.res %>% slice_max(n =6, order_by = NES)
#ind2=tcga.gsea.KEGG.res %>% slice_min(n =8, order_by = NES)  

esga_plot=enrichplot::gseaplot2(tcga.gsea.KEGG,ind1$ID, pvalue_table = F,color=c4a('brewer.paired',6) ,title ='KEGG enrichment')
#esga_plot2=enrichplot::gseaplot2(tcga.gsea.KEGG,ind2$ID, pvalue_table = T,color=c4a('brewer.paired',8),title ='KEGG enrichment in Low group')
pdf('PDFs/esga_plot.pdf', width = 6, height = 6)
esga_plot
dev.off()

save.image("project.Rdata")

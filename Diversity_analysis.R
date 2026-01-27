metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}
library("writexl")
library(phyloseq)
library(ggplot2)
library(ANCOMBC)
library(microbiomeSeq)
library(ggpubr)
library(ggsci)
library(cowplot)
library(factoextra)
library(FactoMineR)
library("vegan")
library("multcomp")
library("patchwork")
library("magrittr")
library("dplyr")
library(ggpubr)
mphlanin <- read.csv("input/merged.mfas.ra.txt", sep = "\t", strip.white = T, stringsAsFactors = F,row.names=1)
mphlanin_used_tmp<-data.frame(apply(mphlanin, 2, function(x) as.numeric(as.character(x))))
rownames(mphlanin_used_tmp) <- row.names(mphlanin)
# only Bacteria included
mphlanin_used_tmp2 <- mphlanin_used_tmp[!(grepl("k__Viruses", rownames(mphlanin_used_tmp),perl = TRUE)),]
mphlanin_used <-mphlanin_used_tmp2[grep("s__",rownames(mphlanin_used_tmp2),perl =TRUE,invert=FALSE),]
mphlanin_used <- mphlanin_used[grep("t__",rownames(mphlanin_used),perl =TRUE,invert=TRUE),]
metadata <- read.delim("input/metadatadf.txt", header=TRUE, sep = "\t")
metadatadf <- data.frame(metadata)
metadatadf<-data.frame(lapply(metadatadf, gsub, pattern = "-", replacement = ".", fixed = TRUE))

write_xlsx(metadatadf,"input/metadatadf.xlsx")
metadatadf_p <- read.csv("input/metadatadf_positive.txt", sep = "\t", strip.white = T, stringsAsFactors = F,row.names=1)
#SampleID should be corrected # corrected.
row.names(metadatadf) <- metadatadf$ID
library("phyloseq")
library("ggplot2")
sampledata <- sample_data(metadatadf)

mphlanin_used$Bacteria<-rownames(mphlanin_used)
write_xlsx(mphlanin_used,"input/mphlanin_used.xlsx")
mphlanin_used <- read.delim("input/mphlanin_used.txt", header=TRUE, sep = "\t")
rownames(mphlanin_used)<-mphlanin_used[,1]
mphlanin_used<-mphlanin_used[,-1]
identical(sort(rownames(sampledata)), sort(colnames(mphlanin_used)))
phyloseqin <- metaphlanToPhyloseq(mphlanin_used, metadat =sampledata)

#-----------
##Effector
#-----------
data <- as.data.frame(t(phyloseqin@otu_table))
all(row.names(data) ==phyloseqin@sam_data$ID)
res_pca <- prcomp(data, scale = T)
data_dist = as.matrix(vegdist(data, method = 'bray'))
res_adonis_all = adonis2(data_dist~class + age + sex + bmi,metadata,by = "term",permutations = 999)

tmp_class <- factor(phyloseqin@sam_data$class,levels = c("control","positive"))
fviz_pca_ind(res_pca, label="none",axes = c(1, 2), alpha.ind =1,habillage=tmp_class,invisible="quali",pointsize = 2.5, addEllipses = TRUE,ellipse.level=0.95,palette = color1)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.1))
res_adonis_all$item <- row.names(res_adonis_all)
tmp_adonis_result_dis_p <- res_adonis_all[!is.na(res_adonis_all$F),]
tmp_adonis_result_dis_p$item <- factor(tmp_adonis_result_dis_p$item,levels = tmp_adonis_result_dis_p$item)
tmp_adonis_result_dis_p$class <- ifelse(tmp_adonis_result_dis_p$`Pr(>F)`>0.05,"no_sig","sig")
tmp_adonis_result_dis_p <- tmp_adonis_result_dis_p%>%
  mutate(
    p_sig = symnum(`Pr(>F)`,
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,1),
      symbols = c("****", "***", "**", "*",""),
      na = "ns"  
    ) %>% as.character()   
  )
tmp_adonis_result_dis_p <- tmp_adonis_result_dis_p[order(-tmp_adonis_result_dis_p$R2), ]
tmp_adonis_result_dis_p$item <- factor(tmp_adonis_result_dis_p$item,levels = tmp_adonis_result_dis_p$item)
pa <- ggplot(tmp_adonis_result_dis_p)+
  geom_bar(aes(x = item, y = R2,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = R2+0.002,label=p_sig,size = 2,fontface = "bold",angle = 270,vjust =0))+
  theme_cowplot()+
  scale_fill_manual(values = color4)+
  labs(x = "Effector",y = "R2") +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none", #移除图例
        axis.title  = element_text(size = 14,face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text( hjust =0.5,vjust = 0.5,face = "bold"),
        axis.ticks.x = element_blank(),
       axis.text  = element_text( hjust =0.5,size =12 ,vjust = 0.5),axis.text.y =element_text( angle =90,size =12,vjust =1.5,hjust = 0.5))
ggsave(plot = pa,filename = "adonis_R2.pdf",width =6,height = 5)

data <- as.data.frame(t(phyloseqin@otu_table))
all(row.names(data) ==phyloseqin@sam_data$ID)
sample_ids_p <- rownames(sample_data(phyloseqin))[sample_data(phyloseqin)$class == "positive"]
data_group <- data[sample_ids_p, ]
data_p <- data_group[, apply(data_group, 2, var) != 0]
res_pca <- prcomp(data_p, scale = T)
data_dist = as.matrix(vegdist(data_p, method = 'bray'))
metadatadf_p <- metadata[metadata$class == "positive",]
res_adonis_all = adonis2(data_dist~parasite_load+ age + sex + bmi,metadatadf_p,by = "term",permutations = 999)
res_adonis_all$item <- row.names(res_adonis_all)
tmp_adonis_result_dis_p <- res_adonis_all[!is.na(res_adonis_all$F),]
tmp_adonis_result_dis_p$item <- factor(tmp_adonis_result_dis_p$item,levels = tmp_adonis_result_dis_p$item)
tmp_adonis_result_dis_p$class <- ifelse(tmp_adonis_result_dis_p$`Pr(>F)`>0.05,"no_sig","sig")
tmp_adonis_result_dis_p <- tmp_adonis_result_dis_p%>%
  mutate(
    p_sig = symnum(`Pr(>F)`,
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,1),
      symbols = c("****", "***", "**", "*","ns"),
      na = "ns"  
    ) %>% as.character()   
  )
tmp_adonis_result_dis_p <- tmp_adonis_result_dis_p[order(-tmp_adonis_result_dis_p$R2), ]
tmp_adonis_result_dis_p$item <- factor(tmp_adonis_result_dis_p$item,levels = tmp_adonis_result_dis_p$item)
pa <- ggplot(tmp_adonis_result_dis_p)+
  geom_bar(aes(x = item, y = R2,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = R2+0.002,label=`Pr(>F)`,size = 2,fontface = "bold",angle = 270,vjust =0))+
  theme_cowplot()+
  scale_fill_manual(values = color4)+
  labs(x = "Effector",y = "R2") +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none", 
        axis.title  = element_text(size = 14,face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text( hjust =0.5,vjust = 0.5,face = "bold"),
        axis.ticks.x = element_blank(),
       axis.text  = element_text( hjust =0.5,size =12 ,vjust = 0.5),axis.text.y =element_text( angle =90,size =12,vjust =1.5,hjust = 0.5))
ggsave(plot = pa,filename = "adonis_R2_p.pdf",width =6,height = 5)
##alpha  diversity
phyloseqin <- subset_samples(phyloseqin, class %in% c("control","positive"))
tmp_otu <- as.data.frame(phyloseqin@otu_table)
tmp_sample <- phyloseqin@sam_data$ID
tmp_type <- factor(phyloseqin@sam_data$class,levels = c("control","positive"))
tmp_load <- as.numeric(phyloseqin@sam_data$parasite_load)
tmp_shannon <- vegan::diversity(t(tmp_otu),index = "shannon")
tmp_simpson <-  vegan::diversity(t(tmp_otu),index = "simpson")
tmp_invsimpson <- 1/tmp_simpson
data_alpha_raw <- data.frame(sample=tmp_sample
                               ,shannon=tmp_shannon
                               ,simpson=tmp_simpson
                               ,invsimpson=tmp_invsimpson
                               ,type=tmp_type
                             ,load_p=tmp_load
)
library(ggpubr)
library(cowplot)
ggplot(data_alpha_raw,aes(x=type,y=shannon))+
  geom_boxplot(aes(fill = type),alpha=0.4,outlier.shape =NA)+
  geom_jitter(aes(color = type, shape = type), height = 0, width = 0.3, size = 3,alpha=0.8) + 
  theme_cowplot()+
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14), 
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size =15,face = "bold"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 20, unit = "pt"), size = 15)  # 调整标题字体大小
  ) +
  scale_color_manual(values =color1) +
  scale_fill_manual(values = color1)+
  labs(y ="Shannon index") +
  stat_compare_means(
    comparisons = list(c("control", "positive")),
    method = "wilcox.test", label = "p.signif"
  )
ggsave("alpha-shannon.pdf",width = 4,height = 5)
ggplot(data_alpha_raw,aes(x=type,y=simpson)) + 
  geom_boxplot(aes(fill = type),alpha=0.4,outlier.shape =NA)+
  geom_jitter(aes(color = type, shape = type), height = 0, width = 0.3, size = 3,alpha=0.8) + 
  theme_cowplot()+
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14), 
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size =15,face = "bold"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 20, unit = "pt"), size = 15)  # 调整标题字体大小
  ) +
  scale_color_manual(values =color1)+
  scale_fill_manual(values = color1)+
  labs(y="Simpson index")+
  stat_compare_means(
    comparisons = list(c("control","positive")),
    method = "wilcox.test", label = "p.signif")
ggsave("alpha-simpson.pdf",width = 4,height = 5)
ggplot(data_alpha_raw,aes(x=type,y=invsimpson)) + 
  geom_boxplot(aes(fill = type),alpha=0.4,outlier.shape =NA)+
  geom_jitter(aes(color = type, shape = type), height = 0, width = 0.3, size = 3,alpha=0.8) + 
  theme_cowplot()+
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14),  
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size =15,face = "bold"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 20, unit = "pt"), size = 15)  # 调整标题字体大小
  ) +
  scale_color_manual(values =color1)+
  scale_fill_manual(values = color1)+
  labs(y="Invsimpson index")+
  stat_compare_means(
    comparisons = list(c("control","positive")),
    method = "wilcox.test", label = "p.signif")
ggsave("alpha-invsimpson.pdf",width = 4,height = 5)
##-----------
##beta diversity
##-------------------
sample_data(phyloseqin)$ID<-row.names(sample_data(phyloseqin))
phyloseqin.pcoa <- transform_sample_counts(phyloseqin, function(x) x/sum(x) * 100)
PERMANOVA_metadata <- data.frame(sample_data(phyloseqin.pcoa))
PERMANOVA_dis<-phyloseq::distance(phyloseqin.pcoa, method="bray", weighted=FALSE)
adonis_result_dis<-adonis2(unname(PERMANOVA_dis)~ class, data = PERMANOVA_metadata, permutations = 999)
PCoA.bray.ord <- ordinate(phyloseqin.pcoa, "PCoA", "bray", weighted = FALSE)
#for Group
plot.data<-PCoA.bray.ord$vectors[,1:2]
plot.data<-cbind(plot.data,PERMANOVA_metadata)
PC1 = plot.data[,1]
PC2 = plot.data[,2]
plotdata.group <- data.frame(rownames(plot.data),plot.data[,1:2],plot.data$class)
colnames(plotdata.group) <-c("sample","PC1","PC2","class")
plotdata.group$class <- factor(plotdata.group$class,levels = c("control","positive"),labels = c("control","positive"))
yf <- plotdata.group
yd1 <- yf %>% group_by(class) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(class) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1
fit1 <- aov(PC1~class,data = plotdata.group)
tuk1<-glht(fit1,linfct=mcp(class="Tukey"))
res1 <- cld(tuk1,alpha=0.05)
fit2 <- aov(PC2~class,data = plotdata.group)
tuk2<-glht(fit2,linfct=mcp(class="Tukey"))
res2 <- cld(tuk2,alpha=0.05)
test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,yd1 = yd1$Max,yd2 = yd2$Max,class =c("control","positive"))
test$class <- factor(test$class,levels = c("control","positive"))

p1 <- ggplot(plotdata.group,aes(class,PC1)) +
  geom_boxplot(aes(fill = class),alpha=0.8) +
  #geom_text(data = test,aes(x = class,y = yd1,label = PC1),size = 7,color = "black",fontface = "bold") +
  coord_flip() +
  ylim(-0.65,0.65)+
  scale_fill_manual(values=color1) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=20,face = "bold"),
        axis.text.x=element_blank(),
        panel.grid=element_blank(),
        legend.position = "none")+
  stat_compare_means(method="wilcox.test",comparisons = list(c("control","positive")),label = "p.signf",size=5)

p3 <- ggplot(plotdata.group,aes(class,PC2)) +
  geom_boxplot(aes(fill = class),alpha=0.8) +
  #geom_text(data = test,aes(x =class,y = yd2,label = PC2),size = 7,color = "black",fontface = "bold") +
  ylim(-0.45,0.5)+
  scale_fill_manual(values=color1) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=20,angle = 45,vjust = 1,hjust = 1,face = "bold"),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none")+
  stat_compare_means(method="wilcox.test",comparisons = list(c("control","positive")),label = "p.signf",size=5)

p2<-ggplot(plotdata.group, aes(PC1, PC2,fill=class)) +
  geom_point(aes(PC1, PC2,color=class,shape=class),size=8,alpha=0.6)+#,pch = 21
  stat_ellipse(aes(color=class),geom = "polygon",level = 0.85,alpha=0.1)+
  scale_fill_manual(values=color1)+
  scale_color_manual(values = color1)+
  xlab(paste("PCoA1 ( ",round(PCoA.bray.ord$values$Relative_eig[1]*100,2),"%"," )",sep=""))+
  ylab(paste("PCoA2 (",round(PCoA.bray.ord$values$Relative_eig[2]*100,2),"%"," )",sep=""))+
  xlim(-0.65,0.65) +
  ylim(-0.45,0.5) +
  theme(text=element_text(size=30,face = "bold"))+
  geom_vline(aes(xintercept = 0),linetype="dashed")+
  geom_hline(aes(yintercept = 0),linetype="dashed")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(),
        axis.ticks.length = unit(0.5,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=20,vjust = 7,face = "bold"),
        axis.title.y=element_text(colour='black', size=20,vjust = -2,face = "bold"),
        axis.text=element_text(colour='black',size=20),legend.position="none")

p4 <- ggplot(plotdata.group, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",adonis_result_dis$Df[1],"\nR2 = ",round(adonis_result_dis$R2[1],4),"\np-value = ",adonis_result_dis$`Pr(>F)`[1],sep = "")),size =5,fontface ="bold") +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

p5 <-p1+p4+p2+p3+plot_layout(heights = c(1,4),widths = c(5,1),ncol = 2,nrow = 2)
p5
ggsave("beta_PERMANOVA.pdf",p5,height =10,width = 12)
library(vegan)
library(aPCoA)
library(factoextra)
library(microbiome)
library(microbiomeSeq)
library(pheatmap)
merge_samples_mean <- function(physeq, group){
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  merged <- merge_samples(physeq, group)
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}
library(phyloseq)
tmp_data_bar <- merge_samples_mean(phyloseqin,"class")
tmp_data_bar_count <- transform_sample_counts(tmp_data_bar, function(x) x/sum(x) * 100)
tmp_data_bar_count<-tax_glom(tmp_data_bar_count, 'Phylum', NArm = T)
tmp_data_bar_count@sam_data$class<- factor(tmp_data_bar_count@sam_data$class,levels = c("control","positive"))
color_bar<-sample(c("#51574a","#74a993","#447c69","#8e8c6d","#e4bf80","#e9d78e","#e2975d","#f19670","#e16552","#be5168","#c94a53","#cf89a9","#993767","#4e2472","#9163b6","#05387d","#5698c4","#9abf88","#328c97","#7c9fb0"),19)

plot_bar(tmp_data_bar_count, fill='Phylum')+
  geom_bar(stat="identity", position = 'fill',color ="transparent")+
  xlab("") +
  ylab("Relative Abundance") +
  theme_classic(base_size = 10) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle('Phylum') +
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1),
        legend.key.size = unit(10, "pt"))+
  guides(fill = guide_legend( ncol = 1, byrow = TRUE))+
  scale_fill_manual(values=c("#51574a","#74a993","#447c69","#8e8c6d","#e4bf80","#e9d78e","#e2975d","#f19670","#e16552","#be5168","#c94a53","#cf89a9","#993767","#4e2472","#9163b6","#05387d","#5698c4","#9abf88","#328c97","#7c9fb0"))+
  theme(legend.text = element_text(size = 8),text = element_text(face = "bold"))+
theme(axis.text.x = element_text(size = 12))
ggsave2("phylum.relative.abundance.pdf",height =6,width = 5,useDingbats=FALSE)
# 提取物种丰度矩阵
otu_matrix <- otu_table(tmp_data_bar_count)
# 提取样本数据
sample_data_df <- data.frame(sample_data(tmp_data_bar_count))
# 提取分类信息
tax_table_df <- data.frame(tax_table(tmp_data_bar_count))

# 将 OTU 矩阵转换为长格式数据框
library(reshape2)
otu_long <- melt(otu_matrix)
colnames(otu_long) <- c("OTU", "Sample", "Abundance")

# 合并分类信息和样本数据
otu_long <- merge(otu_long, tax_table_df, by.x = "OTU", by.y = 0)
otu_long <- merge(otu_long, sample_data_df, by.x = "Sample", by.y = 0)
otu_long$sample <- "sample"

# 绘制柱状图
ggplot(otu_long, aes(x=sample,y=Abundance,fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "transparent") +
  xlab("") +
  ylab("") +
  theme_classic(base_size = 10) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle('Phylum') +
  guides(fill = guide_legend(title = NULL)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  scale_fill_manual(values=c("#51574a","#74a993","#447c69","#8e8c6d","#e4bf80","#e9d78e","#e2975d","#f19670","#e16552","#be5168","#c94a53","#cf89a9","#993767","#4e2472","#9163b6","#05387d","#5698c4","#9abf88","#328c97","#7c9fb0"))
#TOP6 门水平
library(dplyr)
#install.packages("tidyverse")
library(tidyverse) 
mphlanin_used <- read.delim("input/mphlanin_used.txt", header=TRUE, sep = "\t")
rownames(mphlanin_used)<-mphlanin_used[,1]
mphlanin_used<-mphlanin_used[,-1]
mphlanin_used1<-rownames_to_column(mphlanin_used)
mphlanin_used1<-gsub(".*\\|p","",mphlanin_used1$rowname)
mphlanin_used1<-as.data.frame(mphlanin_used1)
mphlanin_used1<-gsub("\\|c_.*","",mphlanin_used1$mphlanin_used1)
mphlanin_used1<-as.data.frame(mphlanin_used1)                            
mphlanin_used1<-cbind(mphlanin_used1,mphlanin_used)                 
box_phylum<-aggregate(mphlanin_used1[, 2:41], list(mphlanin_used1=mphlanin_used1$mphlanin_used1),sum)
#修改第一列列名的操作：
colnames(box_phylum)[1] <- "Phylum"
box_phylum<-t(box_phylum)
col1<-box_phylum[1,]
colnames(box_phylum)<-col1
box_phylum<-box_phylum[-1,]
box_phylum<-as.data.frame(box_phylum)
box_phylum$sample<-rownames(box_phylum)
library("writexl")
#write_xlsx(box_phylum,"D:/桌面/R project/box_phylum.xlsx")
#桌面转到R的操作：
box_phylum_new<- read.delim("input/box_phylum_new.txt", header=TRUE, sep = "\t",na.strings="")
#install.packages("reshape2")
library(reshape2)
box_phylum_new<-box_phylum_new[,-8]
box_phylum_new<-melt(box_phylum_new,id.vars = "class")
colnames(box_phylum_new)[3] <- "Abundance"
colnames(box_phylum_new)[2] <- "Species"
library(dplyr)
library(rstatix)
library(ggpubr)
ggplot(box_phylum_new, aes(x=class, y= Abundance,fill=class))+
  geom_boxplot(alpha=0.3)+
  geom_jitter(aes(color = class, shape = class), height = 0, width = 0.3, size = 3,alpha=0.6) + 
  xlab("")+
  facet_wrap(vars(Species),scales="free_y")+
  scale_fill_manual(values = color2)+
  scale_color_manual(values = color2)+
  theme_bw()+stat_compare_means(method = "wilcox.test",comparisons = list(c("control", "positive")),label = "p.signif",show.legend = F)+
  theme(legend.position="none",text = element_text(face = "bold"))
ggsave2("top6.phylum.pdf",height =8,width = 8,useDingbats=FALSE)

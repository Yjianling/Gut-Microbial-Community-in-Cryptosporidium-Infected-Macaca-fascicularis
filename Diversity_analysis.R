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
Pielou<-evenness(phyloseqin,index = "pielou", zeroes = TRUE, detection = 0)
data_alpha_raw <- data.frame(sample=tmp_sample
                               ,shannon=tmp_shannon
                               ,simpson=tmp_simpson
                               ,invsimpson=tmp_invsimpson
                             ,Pielou=Pielou
                               ,type=tmp_type
                             ,load_p=tmp_load
)
data_alpha <- data.frame(sample=tmp_sample
                               ,shannon=tmp_shannon
                               ,simpson=tmp_simpson
                               ,invsimpson=tmp_invsimpson
                             ,pielou=Pielou
                               ,type=tmp_type)
library(ggpubr)
library(cowplot)

pielou <- ggplot(data_alpha,aes(x=type,y=pielou))+
  geom_boxplot(aes(fill = type),alpha=0.4,outlier.shape =NA)+
  geom_jitter(aes(color = type, shape = type), height = 0, width = 0.3, size = 3,alpha=0.8) + 
  theme_cowplot()+
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14),  # 调整X轴标签字体大小
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14),  # 调整Y轴标签字体大小
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    text = element_text(size =15,face = "bold"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 20, unit = "pt"), size = 15)  # 调整标题字体大小
  ) +
  scale_color_manual(values =color1) +
  scale_fill_manual(values = color1)+
  labs(y ="Pielou index") +
  stat_compare_means(
    comparisons = list(c("control", "positive")),
    method = "wilcox.test", label = "p.signif"
  )
ggsave("alpha-Pielou.pdf",width = 4,height = 5)
shannon <- ggplot(data_alpha_raw,aes(x=type,y=shannon))+
  geom_boxplot(aes(fill = type),alpha=0.4,outlier.shape =NA)+
  geom_jitter(aes(color = type, shape = type), height = 0, width = 0.3, size = 3,alpha=0.8) + 
  theme_cowplot()+
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14),  # 调整X轴标签字体大小
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14),  # 调整Y轴标签字体大小
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
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
simpson <- ggplot(data_alpha_raw,aes(x=type,y=simpson)) + 
  geom_boxplot(aes(fill = type),alpha=0.4,outlier.shape =NA)+
  geom_jitter(aes(color = type, shape = type), height = 0, width = 0.3, size = 3,alpha=0.8) + 
  theme_cowplot()+
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14),  # 调整X轴标签字体大小
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14),  # 调整Y轴标签字体大小
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
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
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14),  # 调整X轴标签字体大小
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14),  # 调整Y轴标签字体大小
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
pielou <- pielou + ggtitle("Pielou")
shannon <- shannon + ggtitle("Shannon")
simpson <- simpson + ggtitle("Simpson")
figure.alpha_diversity<-ggarrange(shannon,pielou,simpson,ncol = 3)
pdf("Figure.alpha_diversity.pdf",height=5,width=8,useDingbats=FALSE)
figure.alpha_diversity
dev.off()
##-----------
##感染隐孢子虫与未感染的beta diversity
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
#top 10 phylum & genus
#食蟹猴感染与未感染门水平相对丰度堆叠柱状图
#library(mecodev)
library("microeco")
library(file2meco)
library(cowplot)
library(ggplot2)
phyloseqin
#mycol1 <-c("#1D2FA1","#5061D0","#22755A","#099B8C","#328c97","#AF5D33","#F48248","#F4A076","#e16552","#c04a53","#995168")
mycol1 <- c("#A3B18A","#5B8E4D","#93B7BE","#3C6E71","#E0A458","#C57B57","#6D4C3D","#DAD7CD","#BCB8B1","#8E9AAF","#2E4A62")
meco_physeq1 <- phyloseq2meco(phyloseqin)
Phylum <- trans_abund$new(dataset = meco_physeq1, taxrank = "Phylum", ntaxa = 10, groupmean ="class")
Phylum.plot <- Phylum$plot_bar(others_color = rev(mycol1), legend_text_italic = FALSE)+
  theme_classic() + 
  ylim(c(0,101)) +ylab('Mean Relative Abundance Ratio (%)')+
  theme_cowplot()+ 
  theme(legend.key=element_blank())+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.y = element_text(hjust = 0.5, face = "bold",size = 12))

pdf("Figure.bar.group.phylum.pdf",height=5,width=5)
Phylum.plot 
dev.off()

Genus.10 <- trans_abund$new(dataset = meco_physeq1, taxrank = "Genus", ntaxa = 10, groupmean = "class")
Genus.10.plot <- Genus.10$plot_bar(others_color = rev(mycol1), legend_text_italic = FALSE)+
  theme_classic() +
  ylim(c(0,101)) +
  ylab('Mean Relative Abundance Ratio (%)')+
  theme_cowplot()+ 
  theme(legend.key=element_blank())+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(hjust = 0.5, face = "bold",size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"))

pdf("Figure.bar.group.genus.top10.pdf",height=5,width=5)
Genus.10.plot 
dev.off()
#统计
tax_df <- as.data.frame(tax_table(phyloseqin))

n_distinct(tax_df$Phylum)
n_distinct(tax_df$Genus)
n_distinct(tax_df$Species)
#门属表
df_g <- Genus.10$data_abund
write_xlsx(df_g,"df_g.xlsx")
df_p <- Phylum$data_abund
write_xlsx(df_p,"df_p.xlsx")
#core genus
library("microbiome")
library("reshape2")
library("ggalluvial")
library(ggpubr)
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
core_color <- c("#1D2FA1","#5061D0","#22755A","#099B8C","#AF5D33","#F48248","#F4A076","#EA8D0D","#F4AC48","#F4BF76")
core.phyloseqin<-tax_glom(phyloseqin, 'Genus', NArm = F)#509
core.phyloseqin.genus <- transform_sample_counts(core.phyloseqin, function(x) x/sum(x) * 100)
core.phyloseqin.genus.t<- subset_taxa(core.phyloseqin.genus, !grepl("^GGB", Genus)&!grepl("_unclassified|_GGB|_sp", Genus))#172
dim(otu_table(core.phyloseqin.genus.t))
core.phyloseq.genus<-core(core.phyloseqin.genus.t,prevalence = 0.50,detection = 1)
ps.merge.core <- merge_samples_mean(core.phyloseq.genus, "class")
tax.table.core<-as.data.frame(tax_table(ps.merge.core))
core.phyloseq.genus.df<-as.data.frame(otu_table(ps.merge.core))
names.core<-row.names(core.phyloseq.genus.df)
core.phyloseq.genus.df<-cbind(names.core,core.phyloseq.genus.df)
sampleData<-reshape2::melt(core.phyloseq.genus.df,id=c(colnames(core.phyloseq.genus.df)[1]))
newColNames<-c("OTU","Group","Abundance")
colnames(sampleData)<-newColNames
sampleData$OTU<-as.factor(sampleData$OTU)
sampleData$Genus <-tax.table.core$Genus
sampleData$Group<-factor(sampleData$Group,levels = c("control","positive"),labels =  c("control","positive"))

ggplot(sampleData, aes(x = Group, y=Abundance,stratum = Genus, alluvium = Genus,fill = Genus, label = Genus))+ylab('Mean Relative Abundance')+xlab('')+
  scale_fill_manual(values = core_color)+
  geom_flow(stat = "alluvium", lode.guidance = "frontback",color = "white") +
  geom_stratum(width = 1/3,linetype=1,size=0.5,alpha =0.8,color='white') +theme_bw()+
   theme(text = element_text(face = "bold"),
    panel.background = element_blank(),         # 去除背景
    panel.grid.major = element_blank(),         # 去除主网格线
    panel.grid.minor = element_blank(),         # 去除次网格线
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    strip.background = element_blank()
    # legend.position = "none"  # 如需要隐藏图例可启用
  )
ggsave2("Figure.core.pdf",height=4,width=5)
library()
data<-psmelt(core.phyloseq.genus)
df<-data %>% dplyr::select(!(4:9))
df<-df%>%dplyr::select(!(5:9))
df<-df%>%dplyr::select(!(1:2))
ggplot(df, aes(x = class, y = Abundance, fill = class,colour = class)) +
  geom_boxplot(alpha =0.3,width = 0.5,outlier.shape  = NA ) +
  geom_point(position = position_jitter(width = 0.1,height =0.1),size = 2.5,alpha = 0.5 
  ) + 
  facet_wrap(~ Genus, ncol =4, scales = "free_y") + # 关键：按Genus分面，y轴自由缩放 
  scale_color_manual(values =c("#101c5e","#5e101c"))+
    scale_fill_manual(values =c("#101c5e","#5e101c"))+
  labs(y = "Relative Abundance") +
  theme_bw() + # 使用白色背景主题
  theme(legend.position = "none", # 隐藏图例（因x轴已标注）
        strip.text = element_text(face = "bold", size = 11), # 分面标题加粗
        axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black",face = "bold", size = 11))+ # 调整x轴标签 
  stat_compare_means(method = "wilcox.test", # 使用非参数检验
                     comparisons = list(c("control", "positive")), 
                     label = "p.signif",
                     hide.ns =F, size=4#,label.y = max(df$Abundance) * 1.005
                     )

ggsave2("core2.pdf",height =7,width =8,useDingbats=FALSE)

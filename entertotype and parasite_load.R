#entertotype
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x<=0.0000001,pseudocount,x))
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
pam.medoids=function(x,k) {
  require(cluster)
  medoids = as.vector(pam(as.dist(x), k, diss=TRUE)$medoids)
  return(medoids)
}
library(ade4)
library(aPCoA)
library(cluster)
library(clusterSim)
library(vegan)
#library(ggstatsplot)
library(vcd)
set.seed(12345)
tmp_meta <-metadatadf
#phyloseqin_genus<-tax_glom(phyloseqin, 'Genus', NArm = F)
library(tibble)
tmp_otu<-mphlanin_used
row.names(tmp_otu) <- gsub(".*\\|g__", "", row.names(tmp_otu))
tmp_otu<-rownames_to_column(tmp_otu)
library(tidyr)
tmp_otu<- separate(tmp_otu, rowname, into = c("Genus", "age"), sep = "\\|s_")
tmp_otu<-tmp_otu[,-2]
tmp_otu<- aggregate(tmp_otu[, 2:41], by = list(tmp_otu$Genus), FUN = sum)
rownames(tmp_otu)<-tmp_otu[,1]
tmp_otu<-tmp_otu[,-1]
identical(tmp_meta$ID,colnames(tmp_otu))
tmp_otu_genus<-tmp_otu
tmp_otu_genus<-tmp_otu_genus[-c(92:388),]
noise_removal <- function(dataframe, percent=0.2, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  ncol(Matrix) -> ncols
  bigones <- apply(Matrix,1, function(x){ length(x[x > 0]) > percent*ncols})
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}
tmp_otu_genus<- noise_removal(tmp_otu_genus,percent = 0.2)

tmp_jsd = dist.JSD(tmp_otu_genus)#计算JSD矩阵,JSD值越小，表示两个样本的分布越相似。
tmp_nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    tmp_nclusters[k]=0 
  } else {
    tmp_data_clusters=pam.clustering(tmp_jsd, k)
    tmp_nclusters[k]=index.G1(t(tmp_otu_genus),tmp_data_clusters,  d = tmp_jsd,
                              centrotypes = "medoids")
  }
}
tmp_ncluster_d <- data.frame(n = 1:20,ch = tmp_nclusters)
library(ggplot2)
ggplot(tmp_ncluster_d,aes(x = n,y = ch))+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(1,20,1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x = "k clusters", y = "Calinski-Harabasz index", title = "Optimal number of clusters")
ggsave2("cluster_number.pdf",width = 6,height = 5)
#将聚类结果整合到元数据
k_best = which(tmp_nclusters == max(tmp_nclusters), arr.ind = TRUE)
tmp_cluster=pam.clustering(tmp_jsd, k = k_best)
tmp_medoids=pam.medoids(tmp_jsd, k = k_best)
tmp_silhouette=mean(silhouette(tmp_cluster, tmp_jsd)[,3])
cat(tmp_silhouette)
tmp_meta$enterotype <- tmp_cluster
tmp_meta$enterotype <- ifelse(tmp_meta$enterotype == 1,"ET1","ET2")
tmp_jsd_data <- as.data.frame(as.matrix(tmp_jsd))
identical(row.names(tmp_jsd_data),tmp_meta$ID)
tmp_jsd_data$class <- tmp_meta$enterotype

tmp_pca=dudi.pca(as.data.frame(t(tmp_otu_genus)), scannf=F, nf=10)
tmp_bet=bca(tmp_pca, fac=as.factor(tmp_cluster), scannf=F, nf=2) 
tmp_bet_res <- as.data.frame(t(tmp_bet$tab))
tmp_bet_res$taxon = row.names(tmp_bet_res)

tmp_et1 <- tmp_bet_res[which(tmp_bet_res[,1] == max(tmp_bet_res[,1])),3]
tmp_et2 <- tmp_bet_res[which(tmp_bet_res[,2] == max(tmp_bet_res[,2])),3]
res_et <- c(tmp_et1,tmp_et2)

tmp_jsd_sub <- as.dist(as.matrix(tmp_jsd_data))
tmp_meta_sub <- tmp_meta

tmp_apcoa_res <- aPCoA(tmp_jsd_sub ~ age,tmp_meta_sub,maincov = enterotype)
tmp_apcoa_matrix <- data.frame(tmp_apcoa_res$plotMatrix)
identical(rownames(tmp_apcoa_matrix),tmp_meta$ID)
library(dplyr)
tmp_res1 <- tmp_meta_sub %>% mutate(PC1 = tmp_apcoa_matrix$X1,PC2 = tmp_apcoa_matrix$X2) %>% 
  dplyr::select(PC1,PC2,enterotype,class,ID)
tmp_res1$enterotype <- as.character(tmp_res1$enterotype)
tmp_adonis_res <- adonis2(tmp_jsd_sub~enterotype+class, tmp_meta_sub)
tmp_adonis_res
#color0 <- c("#0e8fc2","#ff0044")
ggplot(data=tmp_res1 ,aes(x=PC1,y=PC2,color=enterotype))+
  geom_point(aes(shape=class),size=3.5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PCoA1"),
       y=c("PCoA2")#,title = "PAM cluester"
  )+
  stat_ellipse(data=tmp_res1,
               geom = "polygon",
               aes(fill=enterotype,colour = enterotype),
               alpha=0.3)+
  scale_fill_manual(values = color6)+
  scale_color_manual(values = color6)+
  theme(text = element_text(face = "bold"),
    panel.background = element_blank(),         
    panel.grid.major = element_blank(),         
    panel.grid.minor = element_blank(),         
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), 
    strip.background = element_blank()
    # legend.position = "none"  
  )
ggsave("et.pcoa.pdf")
#箱线图
plot_l <- list()
plot_l <- list()
for(i in 1:length(res_et)){
  tmp_d <- as.data.frame(t(tmp_otu_genus))
  tmp_d$ID = row.names(tmp_d)
  tmp_d <- tmp_d[,c("ID",res_et[i])]
  colnames(tmp_d) <- c("ID","value")
  tmp_plot <- inner_join(tmp_res1,tmp_d,by = "ID")
  p1 <- ggplot(tmp_plot, aes(x=enterotype, y =value))+
    geom_boxplot(aes(color=enterotype,fill = enterotype),alpha=0.4,outlier.shape = NA)+
    geom_jitter(aes(color=enterotype),width = 0.2,height = 0,size = 2.5)+
    stat_compare_means(method="wilcox.test",
                       comparisons = list(c("ET1","ET2")),
                       label = "p.signif"
    )+
    theme_bw()+
    theme(legend.position="none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank()
    )+
    scale_color_manual(values = color6)+
    scale_fill_manual(values = color6)+
    labs(y = "Relative Abundance",title = res_et[i])+
    theme(text = element_text(face = "bold",size = 11),
    panel.background = element_blank(),         # 去除背景
    panel.grid.major = element_blank(),         # 去除主网格线
    panel.grid.minor = element_blank(),         # 去除次网格线
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), 
    plot.title  = element_text(hjust = 0.5),
    strip.background = element_blank()
    # legend.position = "none"  # 如需要隐藏图例可启用
  )
  #facet_wrap(.~species)
  plot_l[[i]] <- p1
}
plot_grid(plot_l[[1]],plot_l[[2]],nrow = 1)
ggsave("ET_Relative_abdance.pdf",height =5,width = 6)

df_plot <- tmp_meta%>%
  dplyr::group_by(class, enterotype) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(prop = n / sum(n))
desired_order <- c("control","positive")#"CVCV","CVVR","veax","VRVR"
df_plot$class <- factor(df_plot$class, levels = desired_order )
ggplot(df_plot, aes(x = class, y = prop, fill = enterotype)) +
  geom_bar(stat = "identity",alpha = 0.65) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = color6) +
  labs(x = "group", y = "ET proportion", fill = "Enterotype") +
  theme_bw()+
  theme(text = element_text(face = "bold",color = "black", size = 15),
    panel.background = element_blank(),         # 去除背景
    panel.grid.major = element_blank(),         # 去除主网格线
    panel.grid.minor = element_blank(),         # 去除次网格线
    axis.text= element_text(face = "bold",color = "black"),  # 属名斜体
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), 
    plot.title  = element_text(hjust = 0.5),
    strip.background = element_blank()
    # legend.position = "none"  # 如需要隐藏图例可启用
  )
ggsave("ET_group_402k2.pdf",height = 5,width =5)


tmp_dist_data <- as.data.frame(as.matrix(tmp_jsd))##data.dist
identical(row.names(tmp_dist_data),tmp_meta$ID)
tmp_dist_data$class <- tmp_meta$enterotype

tmp_pca=dudi.pca(as.data.frame(t(tmp_otu_genus)), scannf=F, nf=10)#执行主成分分析
tmp_bet=bca(tmp_pca, fac=as.factor(tmp_cluster), scannf=F, nf=2) #执行置换多元方差分析data.cluster
tmp_bet_res <- as.data.frame(t(tmp_bet$tab))
tmp_bet_res$taxon = row.names(tmp_bet_res)
et_driver_used <-tmp_bet_res
colnames(et_driver_used) <- c("ET1","ET2","taxon")#,"ET3","ET5","ET4"
et_color <- color6
library(tidyverse)
df_long <- et_driver_used %>%
  pivot_longer(
    cols = starts_with("ET"),
    names_to = "Enterotype",
    values_to = "Loading"
  )
top10_each_ET <- df_long %>%
  group_by(Enterotype) %>%
  slice_max(order_by = Loading, n = 10, with_ties = FALSE) %>%
  ungroup()
library(tidytext)
top10_each_ET <- top10_each_ET %>%
  mutate(taxon_reorder = reorder_within(taxon, Loading, Enterotype))
top10_each_ET <- top10_each_ET %>%
  mutate(taxon_reorder = reorder_within(taxon, Loading, Enterotype))
ggplot(
  top10_each_ET,
  aes( x =taxon_reorder,y = Loading,color = Enterotype,fill = Enterotype)) +
  geom_col(width = 0.7,alpha = 0.8) +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~ Enterotype, scales = "free_y") +
  scale_fill_manual(values = et_color ) +
  scale_color_manual(values = et_color ) +
  labs(
    x = NULL,
    y = "BCA loading",
    title = "Dominant genera driving each enterotype"
  ) +
  theme_bw() +
  theme(
     plot.title = element_text(hjust = 0.5),
    text = element_text(size = 18,face = "bold",color = "black"),
    axis.text.y = element_text(face = "italic",color = "black"),  # 属名斜体
 strip.text = element_text(face = "bold"),
    legend.position = "none" )
ggsave("Figure.top10.ET.4020.002&pre50.pdf",height =8,width =15)
  
top10_et1 <- et_driver_used%>%
  arrange(desc(et_driver_used$ET1)) %>%
  slice_head(n=10)%>%
 dplyr::select(-2)
top10_et2 <- et_driver_used%>%
  arrange(desc(et_driver_used$ET2)) %>%
  slice_head(n=10)%>%
 dplyr::select(-1)
Figure.ET1<-ggplot(top10_et1, aes(x = reorder(taxon, ET1), y = ET1)) +
  geom_col(fill =et_color[1]) +
  coord_flip() +
  labs(x = "", y = "") +
  ggtitle(label = "Top 10 drivers of ET1")+
  theme_classic()+
  theme(
     plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    axis.text.y = element_text(face = "italic")  # 属名斜体
  )

ggsave("Figure.top10.ET1.pdf", Figure.ET1 , height=5, width = 7 )
Figure.ET2<-ggplot(top10_et2, aes(x = reorder(taxon, ET2), y = ET2)) +
  geom_col(fill =et_color[2]) +
  coord_flip() +
  labs(x = "", y = "") +
  ggtitle(label = "Top 10 drivers of ET2")+
  theme_classic()+
  theme(
     plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    axis.text.y = element_text(face = "italic")  # 属名斜体
  )
ggsave("Figure.top10.ET2group.pdf",Figure.ET2 , height=5, width = 7 )

tmp_data_ET <- cbind(tmp_meta,data_alpha_raw)
shannon <-ggplot(tmp_data_ET,aes(x=enterotype,y=shannon)) +
  geom_boxplot(aes(color=enterotype),outlier.shape = NA)+
  geom_jitter(aes(color=enterotype),size = 3,height=0,width = 0.2,alpha = 0.7)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(),panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    text = element_text(size =15,face = "bold"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 20, unit = "pt"), size = 15) , # 调整标题字体大小
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(values =et_color)+
  labs(y="Shannon value")+
  stat_compare_means(
    comparisons =list(c("ET1","ET2")),
     method = "wilcox.test", label = "p.signif",hide.ns = F )
simpson <-ggplot(tmp_data_ET,aes(x=enterotype,y=simpson)) +
  geom_boxplot(aes(color=enterotype),outlier.shape = NA)+
  geom_jitter(aes(color=enterotype),size = 3,height=0,width = 0.2,alpha = 0.7)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(),panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    text = element_text(size =15,face = "bold"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 20, unit = "pt"), size = 15) , # 调整标题字体大小
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(values =et_color)+
  labs(y="Pielou value")+
  stat_compare_means(
    comparisons =list(c("ET1","ET2")),
     method = "wilcox.test", label = "p.signif",hide.ns = F )
pielou <-ggplot(tmp_data_ET,aes(x=enterotype,y=pielou)) +
  geom_boxplot(aes(color=enterotype),outlier.shape = NA)+
  geom_jitter(aes(color=enterotype),size = 3,height=0,width = 0.2,alpha = 0.7)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(),panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    text = element_text(size =15,face = "bold"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 20, unit = "pt"), size = 15) , # 调整标题字体大小
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(values =et_color)+
  labs(y="Simpson value")+
  stat_compare_means(
    comparisons =list(c("ET1","ET2")),
     method = "wilcox.test", label = "p.signif",hide.ns = F )
  #facet_wrap(.~key,scales = "free")
pielou2 <- pielou + ggtitle("Pielou")
shannon2 <- shannon + ggtitle("Shannon")
simpson2 <- simpson + ggtitle("Simpson")
figure.alpha_ET<-ggarrange(shannon2,pielou2,simpson2,ncol = 3)
pdf("Figure.alpha_ET.pdf",height=5,width=8,useDingbats=FALSE)
figure.alpha_ET
dev.off()
data_ET <- tmp_data_ET[c("ID","shannon","simpson","pielou","class","enterotype")]
ggplot(data_ET%>%filter(enterotype == "ET1"),aes(x=class,y=shannon)) +
  geom_boxplot(aes(color=class),outlier.shape = NA)+
  geom_jitter(aes(color=class),size = 3,height=0,width = 0.2,alpha = 0.7)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(),panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    text = element_text(size =15,face = "bold"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 20, unit = "pt"), size = 15) , # 调整标题字体大小
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(values =et_color)+
  labs(y="Shannon value")+
  stat_compare_means(
    comparisons =list(c("control","positive")),
     method = "wilcox.test", label = "p.signif",hide.ns = F )
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
# 封装一个函数（强烈推荐）
plot_alpha <- function(data, index_name, ylab_name){
  ggplot(data, aes(x = class, y = .data[[index_name]])) +
    geom_boxplot(aes(color = class), outlier.shape = NA) +
    geom_jitter(aes(color = class),
                size = 3, width = 0.2, alpha = 0.7) +
    theme_cowplot() +
    theme(
      axis.title.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      text = element_text(size = 15, face = "bold"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    scale_color_manual(values = et_color) +
    labs(y = ylab_name) +
    stat_compare_means(
      comparisons = list(c("control", "positive")),
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE
    )
}

data_ET2 <- data_ET %>% filter(enterotype == "ET2")

p1 <- plot_alpha(data_ET2, "shannon", "Shannon index")
p2 <- plot_alpha(data_ET2, "pielou", "Pielou index")
p3 <- plot_alpha(data_ET2, "simpson", "Simpson index")

p_all <- plot_grid(p1, p2, p3, ncol = 3)
title <- ggdraw() + 
  draw_label(paste0("Enterotype 2"),
             fontface = 'bold', size = 16)
plot_grid(title, p_all, ncol = 1, rel_heights = c(0.1, 1))
ggsave("ET2_alpha_diversity.pdf", width = 8, height = 4)

tmp_meta <- tmp_meta %>% mutate(parasite_load = as.numeric(parasite_load))
tmp_meta$log_parasite <- log1p(tmp_meta$parasite_load)

ggplot(tmp_meta %>% filter(class == "positive"),aes(x = enterotype, y =parasite_load,color = enterotype))+
  geom_boxplot(aes(y=parasite_load,fill = enterotype),outlier.shape = NA,alpha = 0.4)+
  geom_jitter(size = 3,width = 0.2,height = 0)+
  stat_compare_means(comparisons = list(c("ET1","ET2")),
                     method = "t.test", label = "p.signif")+
  xlab("Enterotype")+
  ylab("Parasite_load")+
  theme_bw()+
  theme(text = element_text(face = "bold",size = 11),
    panel.background = element_blank(),         # 去除背景
    panel.grid.major = element_blank(),         # 去除主网格线
    panel.grid.minor = element_blank(),         # 去除次网格线
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), 
    plot.title  = element_text(hjust = 0.5),
    strip.background = element_blank()
    # legend.position = "none"  # 如需要隐藏图例可启用
  )+
  scale_color_manual(values = color6)+
  scale_fill_manual(values = color6)
ggsave("ET_log_parasite_load.pdf",height =5,width =4)

#parasite_load & alpha diversity
data_alpha_load <- data.frame(data_alpha_raw[c(1:20),])
data_alpha_load$load_p <- as.numeric(data_alpha_load$load_p)
data_alpha_load$load_p <- log1p(data_alpha_load$load_p)
cor.test(data_alpha_load$shannon, data_alpha_load$load_p, method = "spearman")
cor.test(data_alpha_load$simpson, data_alpha_load$load_p, method = "spearman")
cor.test(data_alpha_load$pielou, data_alpha_load$load_p, method = "spearman")
library(ggpmisc)
ggplot(data_alpha_load, aes(x =load_p, y =shannon)) +
  geom_point( colour = "#143f69",size = 3, alpha = 0.8,na.rm = T) +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, level = 0.99, fill = "#bbb",color="#741119",na.rm = T) + 
 geom_text(aes(label = "rho = -0.06015038 ; p-value = 0.8016"),x =4.5, y =4.8)+
  #scale_color_manual(values = c("#669933","#663399","#339966","#996699")) +#,
  labs( x = "log (parasite_load + 1)") +
  #ylim(0.7,1.1)+
  theme_bw(base_size = 15)+
  theme(text = element_text(face = "bold"),
    panel.background = element_blank(),         # 去除背景
    panel.grid.major = element_blank(),         # 去除主网格线
    panel.grid.minor = element_blank(),         # 去除次网格线
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    strip.background = element_blank()
    # legend.position = "none"  # 如需要隐藏图例可启用
  )
ggsave("Shannon与载虫量log2.pdf",width =5,height =4)
ggplot(data_alpha_load, aes(x =load_p, y =simpson)) +
  geom_point( colour = "#143f69",size = 3, alpha = 0.8,na.rm = T) +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, level = 0.95, fill = "#bbb",color="#741119",na.rm = T) +
  geom_text(aes(label = "rho = 0.007518797 ; p-value = 0.9772"),x =4.5,y =0.99)+
  labs( x = "log (parasite load + 1)") +
  theme_bw(base_size = 15)+
  theme(text = element_text(face = "bold"),
    panel.background = element_blank(),         # 去除背景
    panel.grid.major = element_blank(),         # 去除主网格线
    panel.grid.minor = element_blank(),         # 去除次网格线
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    strip.background = element_blank()
    # legend.position = "none"  # 如需要隐藏图例可启用
  )
ggsave("simpson与载虫量log2.pdf",width =5,height =4)
ggplot(data_alpha_load, aes(x =load_p, y =pielou)) +
  geom_point( colour = "#143f69",size = 3, alpha = 0.8,na.rm = T) +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, level = 0.99, fill = "#bbb",color="#741119",na.rm = T) + 
 geom_text(aes(label = "rho = -0.1458647 ; p-value = 0.538"),x =4.5, y =4.8)+
  #scale_color_manual(values = c("#669933","#663399","#339966","#996699")) +#,
  labs( x = "log (parasite load + 1)",y="Pielou") +
  #ylim(0.7,1.1)+
  theme_bw(base_size = 15)+
  theme(text = element_text(face = "bold"),
    panel.background = element_blank(),         # 去除背景
    panel.grid.major = element_blank(),         # 去除主网格线
    panel.grid.minor = element_blank(),         # 去除次网格线
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), # 显示外边框
    strip.background = element_blank()
  )
ggsave("pielou与载虫量log10.pdf",width =5,height =4)

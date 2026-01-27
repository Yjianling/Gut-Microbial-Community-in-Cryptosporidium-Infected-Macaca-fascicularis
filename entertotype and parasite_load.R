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

tmp_jsd = dist.JSD(tmp_otu_genus)
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
    panel.background = element_blank(),        
    panel.grid.major = element_blank(),        
    panel.grid.minor = element_blank(),      
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), 
    plot.title  = element_text(hjust = 0.5),
    strip.background = element_blank()
    # legend.position = "none" 
  )
  #facet_wrap(.~species)
  plot_l[[i]] <- p1
}
plot_grid(plot_l[[1]],plot_l[[2]],nrow = 1)
ggsave("ET_Relative_abdance.pdf",height =5,width = 6)

ggstatsplot::ggbarstats(
  data = tmp_res1,
  x = enterotype,
  y = class,
  title = "Mfas"
)+coord_flip()+scale_fill_manual(
  values = color1
)

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
    panel.background = element_blank(),         
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),         
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), 
    plot.title  = element_text(hjust = 0.5),
    strip.background = element_blank()
    # legend.position = "none"
  )+
  scale_color_manual(values = color6)+
  scale_fill_manual(values = color6)
ggsave("ET_log_parasite_load.pdf",height =5,width =4)

data_alpha_load <- data.frame(data_alpha_raw[c(1:20),])
data_alpha_load$load_p <- log1p(data_alpha_load$load_p)
cor.test(data_alpha_load$shannon, data_alpha_load$load_p, method = "spearman")
cor.test(data_alpha_load$simpson, data_alpha_load$load_p, method = "spearman")
library(ggpmisc)
ggplot(data_alpha_load, aes(x =load_p, y =shannon)) +
  geom_point( colour = "#143f69",size = 3, alpha = 0.8,na.rm = T) +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, level = 0.99, fill = "#bbb",color="#741119",na.rm = T) + 
 geom_text(aes(label = "rho = -0.06015038 ; p-value = 0.8016"),x =4.5, y =4.8)+
  #scale_color_manual(values = c("#669933","#663399","#339966","#996699")) +
  labs( x = "log (parasite load + 1)") +
  #ylim(0.7,1.1)+
  theme_bw(base_size = 15)+
  theme(text = element_text(face = "bold"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),     
    panel.grid.minor = element_blank(),      
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5), 
    strip.background = element_blank()
    # legend.position = "none" 
  )
ggsave("Shannon与载虫量log2.pdf",width =5,height =4)
ggplot(data_alpha_load, aes(x =load_p, y =simpson)) +
  geom_point( colour = "#143f69",size = 3, alpha = 0.8,na.rm = T) +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, level = 0.95, fill = "#bbb",color="#741119",na.rm = T) +
  geom_text(aes(label = "rho = 0.007518797 ; p-value = 0.9772"),x =4.5,y =0.99)+
  labs( x = "log (parasite load + 1)") +
  #ylim(0.7,1.1)+
  theme_bw(base_size = 15)+
  theme(text = element_text(face = "bold"),
    panel.background = element_blank(),    
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),       
    panel.border = element_rect(color = "black", fill = NA,linewidth =1.5),
    strip.background = element_blank()
    # legend.position = "none" 
  )
ggsave("simpson与载虫量log2.pdf",width =5,height =4)

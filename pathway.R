###pathway
library(phyloseq)
library(Maaslin2)
library(factoextra)
library(vegan)
library(aPCoA)
library(microeco)
library(file2meco)
library(magrittr)
library(aplot)
#install.packages("xlsx")
library(readr)
library(dplyr)
library(openxlsx)
tmp_et <- tmp_meta %>% filter(class%in% c("positive","control"))
tmp_et <- as.data.frame(tmp_et)
row.names(tmp_et) <- tmp_et$ID
tmp_et <-tmp_et[,c(1,7)]
data_pathway <- read_tsv("input/SAMPLE_pathabundance_relab.tsv")
#write.xlsx(data_pathway,"D:/桌面/R project/data_pathway.xlsx")
data_pathway <- read.csv("input/data_pathway.txt", sep = "\t", strip.white = T, stringsAsFactors = F,row.names=1)
colnames(data_pathway)[1] <- "Pathway"
data_pathway <- as.data.frame(data_pathway)
colnames(data_pathway) <- gsub("_merge_Abundance","",colnames(data_pathway))
write_tsv(data_pathway, path = "input/data_pathway.tsv")
meco_pathway_metacyc <- humann2meco(feature_table = "input/data_pathway.tsv", 
                                    db = "MetaCyc", 
                                    sample_table = tmp_et)#将代谢途径丰度表转换成 MECO 格式
#数据清洗
meco_pathway_metacyc$tidy_dataset()
meco_pathway_metacyc$cal_abund(rel = T)
meco_pathway_metacyc$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Superclass2 %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$pathway %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Species %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Genus %<>% .[!grepl("unclass", rownames(.)), ]

meco_pathway_metacyc_used <- trans_abund$new(meco_pathway_metacyc, taxrank = 
                      "Superclass2", ntaxa = 10, use_percentage = FALSE)#转化TOP10的分类单元
meco_pathway_metacyc_used$ylabname <- "Relative Abundance"
#每个样本都有的条形图
meco_pathway_metacyc_used$plot_bar(facet ="class", bar_type = "notfull")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))
tmp_plot_data <- meco_pathway_metacyc_used$data_abund %>% arrange(desc(all_mean_abund))#根据all_mean_abund对数据进行排序
tmp_plot_data$Taxonomy <- factor(tmp_plot_data$Taxonomy,levels = rev(unique(tmp_plot_data$Taxonomy)))#将分类学信息（Taxonomy）转换为因子，并且逆序排序
tmp_plot_pathway <- unique(tmp_plot_data$Taxonomy)[1:15]#选择排序后的前15个分类学组
#两个分组的条形图                                   
ggplot(tmp_plot_data %>% filter(Taxonomy %in% tmp_plot_pathway),
       aes(x=class,y=Abundance,fill = Taxonomy)) +
  stat_summary(fun=median, geom="bar" ,width = 0.6,position = "stack")+
  theme_bw()+
  ylab("Relative Abundance")+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = rev(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                   "#A6761D","#666666","#8DCEBB","#ECAF80","#BAB7D9","#F394C4",
                   "#B2D28E","#F2D480","#D2BA8E")))
ggsave("Taxonomy_Relative_Abbundance_class.pdf")
tmp_data_stat <- data.frame()
summary_stat <- function(data,value_col = 1,class_col = 2){
  zero_count <- function(x){
    c = length(which(x == 0) )
    return(c)
  }
  
  data <- as.data.frame(data)
  median <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=median)
  mean <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=mean)
  sd <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=sd)
  len <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=length)
  zero <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=zero_count)
  
  data_stat <- data.frame(mean,median=median$x, sd=sd$x, len=len$x,zero=zero$x)
  colnames(data_stat) = c("Group","Mean","Median", "Sd", "Count","ZeroCount")
  data_stat$Se <- data_stat$Sd/sqrt(data_stat$Count)
  data_stat$q25 = 0
  data_stat$q75 = 0
  for(i_t in 1:nrow(data_stat)){
    tmp_class = data_stat[i_t,1]
    tmp_data <- data[data[,class_col] == tmp_class,value_col]
    tmp_q25 <- quantile(tmp_data,probs = 0.25) 
    tmp_q75 <- quantile(tmp_data,probs = 0.75)
    data_stat[i_t,8] = tmp_q25
    data_stat[i_t,9] = tmp_q75
  }
  return(data_stat)
}
for(i in 1:length(tmp_plot_pathway)){
  tmp_data <- tmp_plot_data %>% filter(Taxonomy == tmp_plot_pathway[i])
  tmp <- summary_stat(tmp_data,3,5) %>% mutate(Taxonomy = tmp_plot_pathway[i]) 
  tmp_data_stat <- rbind(tmp_data_stat,tmp)
}
library(tidyr)
tmp_data_stat_short <- tmp_data_stat %>% 
  dplyr::select(Group,Mean,Taxonomy) %>% 
  spread(Group,Mean) %>%
  mutate(Var=positive-control)
tmp_data_stat_short$Taxonomy <-  factor(tmp_data_stat_short$Taxonomy,levels = rev((tmp_data_stat_short$Taxonomy)))
#能看出两组差异的柱形图
ggplot(tmp_data_stat_short)+
  geom_bar(aes(x=Taxonomy,y = Var,fill =Taxonomy),
           stat = "identity",position = "identity")+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = "none",text = element_text(face = "bold"))+
  scale_fill_manual(
    values = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
               "#A6761D","#666666","#8DCEBB","#ECAF80","#BAB7D9","#F394C4",
               "#B2D28E","#F2D480","#D2BA8E")
  )+
  coord_flip()+
  labs(y="",x="CPM variation")
ggsave("CPM_vaeiation.pdf")
#lefse图
meco_pathway_metacyc_diff <- trans_diff$new(meco_pathway_metacyc,
                                            method =  "lefse",group = "class",
                                            taxa_level = "pathway",
                                            alpha = 0.05,p_adjust_method = "fdr"
                                            #reference = c("Type,Health"),fixed_effects = c("Type","Sex")
)
color7 <- c("#5E6A8B","#b80422")
meco_pathway_metacyc_diff$plot_diff_bar(use_number = 1:30)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values =color5)+
  scale_fill_manual(values =color5)
ggsave("pathway_metacyc_diff_lefse.pdf")
#PCOA
library(ape)
meco_pathway_metacyc$cal_betadiv(method = "bray")
meco_pathway_metacyc_distance <- as.dist(meco_pathway_metacyc$beta_diversity$bray)
meco_pathway_metacyc_beta <- trans_beta$new(dataset = meco_pathway_metacyc, 
                                            group = "class", measure = "bray")

meco_pathway_metacyc_beta$cal_ordination(ordination = "PCoA")
meco_pathway_metacyc_beta_adonis <- adonis2(meco_pathway_metacyc_distance~class,tmp_meta)

meco_pathway_metacyc_beta_tmp <- meco_pathway_metacyc_beta$res_ordination$scores
meco_pathway_metacyc_beta_tmp2 <- trans_env$new(dataset = meco_pathway_metacyc, 
                                                add_data = meco_pathway_metacyc_beta_tmp[, 1:2])
# 'KW_dunn' for non-parametric test
tmp_list = combine_list(meco_pathway_metacyc_beta_tmp$class)
color_1 <-"#0e8fc2"
color_2<-"#ff0044"
ggplot()+geom_point(data = meco_pathway_metacyc_beta_tmp,
             mapping = aes(x=PCo1,y=PCo2,color=class,shape=class),
             size=3.5,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.92,0.90),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PCoA1"),y=c("PCoA2"),title ="adonis R2:0.10962,p=0.002")+
  stat_ellipse(data=meco_pathway_metacyc_beta_tmp,geom = "polygon",
               aes(x=PCo1,y=PCo2,fill=class),
               color =NA,alpha=0.2)+
  scale_fill_manual(values = color5)+
  scale_color_manual(values = color5)
ggsave("pathway_metacyc_beta_PCOA.pdf")

#Randomforest
# install.packages("remotes")
# remotes::install_github("jeffreyevans/rfUtilities")
# install.packages("rfUtilities")
dir.create("output/randomforest_yjl",showWarnings = F)
out_dir <- "output/randomforest_yjl/"

rm(list=ls())
library(tidyverse)
library(magrittr)
library(phyloseq)
library(ggpubr)
library(cowplot)
library(ggrepel)
install.packages("ggpmisc")
library(ggpmisc)
install.packages("ggview")
library(ggview)
library(randomForest)
install.packages("rfUtilities")
library(rfUtilities)
library(rfPermute)
library(caret)
library(pROC)
source("common_custom_function.R")
set.seed(10001)
noise_removal <- function(data, percent=0.1,low = 0, method = "pre_cut",dim = 1){
  matrix <- data
  tmp_col_len <- ncol(matrix)
  if(method == "max_cut"){
    bigones <- apply(matrix,dim,function(x){max(x) > low})
    print(c(method,low))
  }
  if(method == "pre_cut"){
    bigones <- apply(matrix,dim,function(x){length(x[x > low]) >= percent*tmp_col_len})
    print(c(method,percent,low))
  }
  if(method == "mean_cut"){
    bigones <- apply(matrix,dim,function(x){sum(x,na.rm = T)/sum(!is.na(x)) > low})
    print(c(method,low))
  }
  if(method == "na_cut"){
    bigones <- apply(matrix,dim,function(x){sum(!is.na(x)) >= percent*tmp_col_len})
    print(c(method,percent))
  }
  matrix_return <- matrix[bigones,]
  return(matrix_return)
}
combine_list <- function(x,y = "none"){
  x = as.character(x)
  tmp_c=unique(x)
  tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
    dplyr::mutate(tmp = str_c(V1,",",V2)) %>% 
    dplyr::select(tmp)
  if(y != "none"){
    tmp_a <- as.data.frame(tmp_a[grepl(y,tmp_a[,1]),])
    colnames(tmp_a) <- "tmp"
  }
  t_list <- (str_split(tmp_a$tmp,pattern = ","))
  return(t_list)
}
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
zero_one_norm <- function(x){
  max = max(x)
  min = min(x)
  t = (x-min)/(max-min)
  return(t)
}
plot_specify_sp <- function(otu,meta,sp,color_m,sample_col = 1,group_col = 8){
  
  if(! sp %in% row.names(otu)){
    c <- print(paste(c("sp is not exist:",sp),collapse = " "))
    return(c)
  }
  
  meta <- meta[,c(sample_col,group_col)]
  colnames(meta) <- c("sample","group")
  tmp_data <- as.data.frame(t(otu[which(row.names(otu)==sp),meta$sample]))
  tmp_data$sample <- row.names(tmp_data)
  colnames(tmp_data) <- c("value","sample")
  tmp_data_p <- inner_join(meta,tmp_data,by = "sample")
  tmp_list <- combine_list(tmp_data_p$group)
  #tmp_data_p$value2 <- log10(tmp_data_p$value + 1e-6)
  res_l <- list()
  p2 <- ggplot(tmp_data)+
    geom_bar(aes(x=sample,y=value),
             stat="identity", position="stack")+
    theme_cowplot()+
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
    )+
    labs(title = sp)
  p <- ggplot(tmp_data_p,aes(x=group,y=value2,color = group))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(width = 0.2,height = 0)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="",title = sp)+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color_m)
  res_l[[1]] <- tmp_data_p
  res_l[[2]] <- p
  res_l[[3]] <- p2
  return(res_l)
}
plot_correalation <- function(d,col_1 = 1,col_2 = 2, meth = "pearson"){
  dd <- d[,c(col_1,col_2)]
  colnames(dd) <- c("g1","g2")
  res_l <- list()
  tmp_formula <- y ~ x
  p <- ggplot(dd,aes(x=g1,y=g2)) + 
    stat_smooth(method='lm',formula = tmp_formula,color = "navy")+
    geom_point(color = "gray40",size=3)+
    stat_cor(cor.coef.name = "rho",method = meth,
             label.y.npc = 1,label.x.npc = 0.35)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
    )
  t <- cor.test(dd$g1,dd$g2)
  m <- lm(g2~g1,dd)
  res_l[[1]] <- p
  res_l[[2]] <- t
  res_l[[3]] <- m
  return(res_l)
}
list_pair_test <- function(t_data,t_list,value_col,group_col,method = "np" ){
  t_res <- list()
  t_data <- t_data[,c(group_col,value_col)]
  colnames(t_data) <- c("group","value")
  for(i in 1:length(tmp_list)){
    i = 1
    t_g <- tmp_list[[i]]
    t_data1 <- t_data[t_data$group == t_g[1],]
    t_data2 <- t_data[t_data$group == t_g[2],]
    if(method == "np"){
      t_res[[i]] <- wilcox.test(t_data1$value,t_data2$value)
    }
    if(method == "p"){
      t_res[[i]] <- t.test(t_data1$value,t_data2$value)
    }
  }
  return(t)
}



# load("phyloseqin.Rdata")
# load("sampledata.Rdata")
tmp_otu <- as.data.frame(otu_table(phyloseqin))
tmp_otu <- noise_removal(tmp_otu,percent=0.2,low = 0.01, method = "pre_cut")
tmp_otu <- sweep(tmp_otu,2,colSums(tmp_otu),"/") * 100
tmp_meta <- sampledata %>% dplyr::select(ID,class)
group <- unique(tmp_meta$class)
tmp_data <- as.data.frame(t(tmp_otu))
tmp_data$ID <- row.names(tmp_data)
tmp_data <- inner_join(tmp_data,tmp_meta,by="ID") 
row.names(tmp_data) <- tmp_data$ID
tmp_data <- tmp_data %>% select(!ID)

tmp_train_use <- sample(nrow(tmp_data), nrow(tmp_data)*0.7)
tmp_data_train <- tmp_data[tmp_train_use,]

tmp_errrate <- c(1)
for(i in 1:ncol(tmp_data_train)-1){
  tmp_model <- randomForest(x=tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                            y=factor(tmp_data_train$class,levels = group[c(1,2)]),
                            ntree=1000, 
                            importance=TRUE, 
                            proximity=TRUE, 
                            mtry=i,
                            na.action=na.omit)
  tmp_err<-mean(tmp_model$err.rate)
  tmp_errrate[i] <- mean(tmp_err)
}
which.min(tmp_errrate)
model_classify <- randomForest(x=tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                               y=factor(tmp_data_train$class,levels = group[c(1,2)]),
                               ntree=1000, 
                               importance=TRUE, 
                               proximity=TRUE, 
                               replace = TRUE,
                               mtry=which.min(tmp_errrate),
                               na.action=na.omit
)
plot(model_classify)
model_classify2 <- randomForest(x=tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                                y=factor(tmp_data_train$class,levels = group[c(1,2)]),
                                ntree=1000, 
                                importance=TRUE, 
                                proximity =TRUE, 
                                replace = TRUE,
                                mtry=which.min(tmp_errrate),
                                na.action=na.omit,
                                oob_score=TRUE
)
model_imp <- importance(model_classify2)
model_imp <- data.frame(predictors = rownames(model_imp), model_imp)
model_imp_sort <- arrange(model_imp, desc(MeanDecreaseAccuracy))
#model_imp_sort_20 <- model_imp_sort[1:45,]
#model_imp_sort_20$predictors <- gsub("s__","",model_imp_sort_20$predictors)
model_classify_pm <- rfPermute(x=tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                               y=factor(tmp_data_train$class,levels = group[c(1,2)]),
                               ntree=1000, 
                               importance=TRUE, 
                               proximity =TRUE, 
                               replace = TRUE,
                               mtry=which.min(tmp_errrate),
                               na.action=na.omit,
                               oob_score=TRUE
)
model_imp_pvalue <- as.data.frame(model_classify_pm$pval[,,2]) %>% select(MeanDecreaseAccuracy)
model_imp_pvalue$predictors <- row.names(model_imp_pvalue)
colnames(model_imp_pvalue) <- c("Pvalue","predictors")
model_imp_sort <- inner_join(model_imp_sort,model_imp_pvalue,by = "predictors")
model_imp_sort <- model_imp_sort %>% mutate(
  sig = case_when(
    Pvalue < 0.01 ~ "**",
    Pvalue < 0.05 & Pvalue >= 0.01 ~ "*",
    TRUE ~ "ns"
  )
)
write_csv(model_imp_sort,"model_imp.csv")

model_imp_sort_top <- model_imp_sort[model_imp_sort$MeanDecreaseAccuracy>2,]
model_imp_sort_top$predictors <- factor(model_imp_sort_top$predictors,levels = rev(model_imp_sort_top$predictors))

tmp_data_train_rfcv <- list()
for(i in 1:5){
  tmp_data_train_rfcv[[i]] <- randomForest::rfcv(tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                                                 factor(tmp_data_train$class,levels = group[c(1,2)]),
                                                 cv.fold = 10, 
                                                 step = 1.2,
                                                 ntree=1000,
                                                 na.action=na.omit, 
                                                 oob_score=TRUE
  )
  print(i)
}

  
tmp_data_train_rfcv_cv <- data.frame(sapply(tmp_data_train_rfcv, '[[', 'error.cv'))
tmp_data_train_rfcv_cv$otu_num <- row.names(tmp_data_train_rfcv_cv)
tmp_data_train_rfcv_cv <- gather(tmp_data_train_rfcv_cv,key = "sp_num",value = "cv",!otu_num)
tmp_data_train_rfcv_cv_stat <- summary_stat(tmp_data_train_rfcv_cv,value_col = 3,class_col = 1) %>% arrange(Mean)
tmp_data_train_rfcv_cv_stat$Group <- as.numeric(tmp_data_train_rfcv_cv_stat$Group)
ggplot(tmp_data_train_rfcv_cv_stat)+
  geom_point(aes(x = Group, y = Mean))+
  geom_line(aes(x = Group, y = Mean),linetype = "dashed")+
  geom_errorbar(aes(x = Group, y = Mean,ymin = Mean-Se,ymax = Mean+Se))+
  scale_x_continuous(limits = c(0,100),
                     breaks=tmp_data_train_rfcv_cv_stat$Group
  )+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  labs(x = "Sp num", y = "Mean Cv")
ggsave("randomforest.error_cv.pdf")
library(ggplot2)
library(dplyr)

# 筛选出前23个 Pvalue 较小的条目
model_imp_top13 <- model_imp_sort_top %>%
  arrange(Pvalue) %>%  # 按 Pvalue 升序排列（Pvalue 越小越显著）
  head(13) 
ggplot(model_imp_top13) +
  geom_bar(aes(x = predictors, y = MeanDecreaseAccuracy),stat = "identity", fill = "grey80") +
  geom_text(aes(x = predictors, y = MeanDecreaseAccuracy+0.25,label = sig),angle = 270)+
  coord_flip() +
  labs(title= "The important species",x = "Species")+
  theme(strip.background = element_blank())+
  theme_bw()
ggsave(paste(out_dir,"cross.rf_imp13.pdf",sep = '/'),width = 8,height = 8)

plot_l1 <- list()
for(i in 1:nrow(model_imp_sort_top)){
  tmp_t <- as.data.frame(t(tmp_otu))
  tmp_t$ID <- row.names(tmp_t)
  tmp_d <- inner_join(tmp_meta,tmp_t,by ="ID") %>% dplyr::select(model_imp_sort_top[i,1],ID,class)
  
  tmp_c=unique(tmp_d$class)
  tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
    mutate(tmp = str_c(V1,",",V2)) %>% 
    select(tmp)
  tmp_list <- (str_split(tmp_a$tmp,pattern = ","))
  
  colnames(tmp_d) <- c("value","sample","groups")
  p <- ggplot(tmp_d,aes(x=groups,y=value,color = groups))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(width = 0.2,height = 0)+
    scale_colour_manual(values = color1)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="Relative abundance",title = model_imp_sort_top[i,1])+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")
  plot_l1[[i]] <- p
}
plot_grid(plot_l1[[1]],plot_l1[[2]],plot_l1[[3]],plot_l1[[4]],
          plot_l1[[5]],plot_l1[[6]],plot_l1[[7]],plot_l1[[8]],
          plot_l1[[9]],plot_l1[[10]],plot_l1[[11]],plot_l1[[12]],
          plot_l1[[13]],plot_l1[[14]],plot_l1[[15]],plot_l1[[16]]#,
          #plot_l1[[17]],plot_l1[[18]],plot_l1[[19]],plot_l1[[20]]
          ,nrow = 4)
ggsave("cross.rf_imp_sp16ALL.pdf",width = 20,height = 16)

model_imp_best <- model_imp_sort[model_imp_sort$MeanDecreaseAccuracy>2 & model_imp_sort$Pvalue< 0.05,]
#model_imp_best <- model_imp_best[1:28,] 
#tmp_data_train_f <- tmp_data_train[,c(model_imp_best$predictors,"groups")] 



tmp_data_train<-  tmp_data[ row.names(tmp_data) %in% row.names(tmp_data_train),]
tmp_data_train_f <- tmp_data[row.names(tmp_data) %in% row.names(tmp_data_train_f),c(model_imp_best$predictors,"class")]
tmp_data_pre <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train),]
tmp_data_pre_f <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train_f),c(model_imp_best$predictors,"class")]

model_train_2 <- randomForest(x = tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],
                              y = factor(tmp_data_train_f$class,levels = group[c(1,2)]),
                              ntree=1000,
                              importance=TRUE,
                              proximity=TRUE,
                              mtry = 20
)
train_predict <- predict(model_train_2,tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],type="response")
test_predict <- predict(model_train_2,tmp_data_pre_f[,1:(ncol(tmp_data_pre_f)-1)],type="response")

confusionMatrix(train_predict, factor(tmp_data_train_f$class,levels = group[c(1,2)]))
confusionMatrix(test_predict, factor(tmp_data_pre_f$class,levels = group[c(1,2)]))

tmp_mds_data <- MDSplot(model_train_2, factor(tmp_data_train_f$class),k=4,pch=30)
tmp_mds_data <- as.data.frame(tmp_mds_data$points) %>% 
  mutate(ID  = row.names(.)) %>% 
  inner_join(.,tmp_meta,by="ID")

ggplot(tmp_mds_data,aes(`Dim 1`,`Dim 2`,color=class))+
  geom_point(size = 4)+
  theme_bw()+
  scale_color_manual(values = color1)+
  theme(
    panel.grid = element_blank()
  )+
  labs(title = "RandomForest MDS plot")
ggsave("cross.rf_mds.pdf",width = 6,height = 4)
train_predict2 <- predict(model_train_2,tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],type="vote")
test_predict2 <- predict(model_train_2,tmp_data_pre_f[,1:(ncol(tmp_data_pre_f)-1)],type="vote")

test_roc <- roc(factor(tmp_data_pre_f$class,levels = group[c(1,2)]), 
                test_predict2[,1],  
                plot=T,
                ci.method="bootstrap",
                ci =TRUE,
                conf.level =0.95
)
train_roc <- roc(factor(tmp_data_train_f$class,levels = group[c(1,2)]), 
                 train_predict2[,1],  
                 plot=T,
                 ci.method="bootstrap",
                 ci =TRUE,
                 conf.level =0.95
)
library(pROC)
ci.auc(test_roc, method="bootstrap", boot.n=2000)
ci.auc(train_roc, method="bootstrap", boot.n=2000)
roc_with_ci <- function(obj,color_p) {
  ciobj <- ci.se(obj, specificities = seq(0, 1, l = 20))
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  dat.ci$lower[dat.ci$lower==1] <- 0.99                     
  ggroc(obj, colour =color_p) +theme_minimal() +
    geom_abline(
      slope = 1,
      intercept = 1,
      linetype = "dashed",
      alpha = 0.7,
      color = "grey"
    ) + coord_equal() +
    geom_ribbon(
      data = dat.ci,
      aes(x = x, ymin = lower, ymax = upper),
      fill = color_p,
      alpha = 0.1
    ) + ggtitle(capture.output(obj$auc))
}

roc_with_ci(test_roc,"firebrick3")
ggsave("cross.rf_roc_pre.pdf")
roc_with_ci(train_roc,"navy")
ggsave("cross.rf_roc_train.pdf")

save(model_train_2,
     tmp_data_train_f,tmp_data_pre_f,
     tmp_data,tmp_meta,
     file = "rf.Rdata")
save(mphlanin, file = "mphlanin.RData")
save(metadatadf, file = "metadatadf.monkey.RData")


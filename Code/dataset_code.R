#Name: dataset_code.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Concordance index (C-index) analysis for Nick's clusters of genes
#Load needed packages----
library(ggplot2)
library(survival)
library(survminer)
library(TCGAbiolinks)
library(tidyverse)

#Loading the cancer specific data frame----
#For TCGA-COAD (colorectal cancer)
load("Data/TCGA-Data/coad_df.RData",
     verbose = TRUE)

#For TCGA-READ (rectal cancer)
load("Data/TCGA-Data/read_df.RData",
     verbose = TRUE)

#For TCGA-LUSC (lung cancer)
load("Data/TCGA-Data/lusc_df.RData",
     verbose = TRUE)

#For TCGA-GBM (brain cancer)
load("Data/TCGA-Data/gbm_df.RData",
     verbose = TRUE)

#For TCGA-BRCA (breast cancer)
load("Data/TCGA-Data/brca_df.RData",
     verbose = TRUE)


#Loading the cox_model_fitter() function from its file----
source("Code/cox_model.R")

#Loading Nick's data and removing any data frames that have no genes----
clust_int_em_nod <- readRDS(file = "Data/Clusters_Integration_EMandC1_NoDup.rds")
clust_int_em_nod[[3]] <- NULL
clust_int_em_nod[[4]] <- NULL

clust_int_em_and_c1 <- readRDS(file = "Data/Clusters_Integration_EMandC1.rds")
clust_int_em_and_c1[[3]] <- NULL
clust_int_em_and_c1[[4]] <- NULL

clust_int_nod <- readRDS(file = "Data/Clusters_Integration_NoDup.rds")
clust_int_nod[[3]] <- NULL

#For loop to get all of the gene names from each data frame and submit them
#to the cox_model_fitter() function to get concordance index----
counter <- 1
outputs <- list()
my_cindicies <- c()

for (df in clust_int_em_nod){
  current_df <- df
  if(dim(current_df)[1] != 0){
    cox_model <- cox_model_fitter(my.seed = 1, cox.df = cox_df,
                                  gene.num = length(rownames(current_df)),
                                  cox.predictors = current_df$gene,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = FALSE,
                                  save.regular.cox.genes = FALSE,
                                  my.filename = paste0("Data/Outputs/COAD/",counter,"_clust_int_nod_coad.csv"))
    
    outputs[[counter]] <- cox_model
    perc_done <- paste0(round((100*counter)/length(clust_int_em_nod)), "% done with this dataset")
    print(perc_done)
    counter <- counter + 1
    
    #Storing all of the c-index values in a vector that we can use later to build the plot
    c_finder <-cox_model$CV$index[1]
    current_c <- cox_model$CV$cvm[c_finder]
    current_c <- round(current_c, digits = 4)
    my_cindicies <- c(my_cindicies, current_c)
    
  }else{
    print("This input doesn't have any genes and the cox model is skipping to the next data frame")
    perc_done <- paste0(round((100*counter)/length(clust_int_em_nod)), "% done with this dataset")
    print(perc_done)
    counter <- counter + 1
  }
  
  #Getting the top performing data frame 
  top_cindex <-max(my_cindicies)
  top_index <- which(my_cindicies==top_cindex)
  
}

print(paste0("The index of the top performing dataset is: ", top_index))
print(paste0("It's concordance index is: ",my_cindicies[top_index]))



#Now doing it for just the top Concordance index result----
cox_models <- list()
my_cindicies <- c()
counter <- 1

for (x in clust_int_nod[top_index]) {
  current_df <- x
  current_cox <- cox_model_fitter(my.seed = 1,
                                  cox.df = cox_df,
                                  gene.num = length(rownames(current_df)),
                                  cox.predictors = current_df,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = paste0("Data/Outputs/READ/",counter,"_clust_int_nod_read_top.csv"))

  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1

  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
}

cox_models$`1`$CV
c_index_df <- data.frame(c_index=my_cindicies, active_coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
write.csv(c_index_df, file = paste0("Data/Outputs/READ/",counter,"clust_int_nod_read_c_index_df.csv"))



#For plotting the coefficients----
coef_df <- data.frame(coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
coef_df_sub <- filter(coef_df, abs(coefs)>0.0001)


colnames(coef_df_sub)[2] <- "Gene"
coef_plot <- ggplot(data = coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("Nick Top READ Cox Coefficients")+
  ylab("Coefficients")+
  xlab("Gene")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


coef_plot

#Saving the finished plot as a .SVG----
svg("Data/Outputs/READ/read_coef_plot_clust_int_nod.svg")
coef_plot
dev.off()

#Saving the finished plot as a .PNG----
png("Data/Outputs/READ/read_coef_plot_clust_int_nod.png")
coef_plot
dev.off()

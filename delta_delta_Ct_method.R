library(writexl)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(plotrix)

#Set working directory and upload data
setwd("C:/Users/EUCHO/Desktop/3T3L-1 Ph and 9P exposure/Time course qPCR/e2 data")
dat <- read.csv("C:/Users/EUCHO/Desktop/3T3L-1 Ph and 9P exposure/Time course qPCR/e2 raw/syed time course 24h phe actin cyp1a1 paary fabp4 e2 -  Quantification Cq Results_0.csv") %>%
  select(Well, Cq) %>%
  filter(Cq != "NaN")

#Upload plate layout
plate_map <- read_xlsx("C:/Users/EUCHO/Desktop/3T3L-1 Ph and 9P exposure/Time course qPCR/plate map/plate_24h_phe.xlsx") %>%
  filter(Gene!= 0 & Gene!="NA")
  #plate number refers to the exposure plate

##########Experimental info#########
#sampling time
time <- "24h"
#housekeeping gene
hk_gene <- "Actin"  
#vehicle solvent
vc <- "DMSO"
#experiment number
exp <- "e2"

#######################################################################
#Add plate info to raw data
dat <- dplyr::left_join(dat, plate_map) %>%
  filter(Cq != "NaN"& Gene!= "NA") %>%
  filter(Chemical != "MI")

colnames(dat)[colnames(dat)=="Plate"] <- "Exp_plate"

#Generate a list of genes of interest (GOI)
genes <- data.frame(genes=unique(dat$Gene)) %>%
  filter(genes!=hk_gene)

genes_str <- as.character(genes$genes)

#Generate a list of chemicals
Chemical <- c(unique(dat$Chemical))
  
#Number of datasets/experimental plates included in the overall input data
plate <- data.frame(plate=unique(plate_map$Plate))
num_plate <- nrow(plate)

exp_plate <- as.character(plate$plate)

#Generate csv file for annotated raw data
write.csv(dat, paste0(time, "_dat_", paste0(exp_plate, collapse = "_"),"_", exp,".csv"))

colnames(dat)[colnames(dat)=="Exp_plate"] <- "Plate"

#Average Cq and calculate standard error in Cq
stderr_Cq <- dat %>%
  subset(select = -c(Well)) %>%
  group_by(Plate, Chemical, Gene, Concentration) %>%
  mutate(Std_dev = sd(Cq), Std_error = std.error(Cq)) %>%
  subset(select= -c(Cq))

stderr_Cq <- stderr_Cq[!duplicated(stderr_Cq), ]

avg_Cq <- dat %>%
  subset(select = -c(Well)) %>%
  group_by(Plate, Chemical, Gene, Concentration) %>%
  summarise_at(vars(Cq), list(Avg_Cq = mean)) %>%
  dplyr::left_join(stderr_Cq)

write.csv(avg_Cq, paste0(time,"_avg_Cq_",paste0(exp_plate, collapse = "_"),"_",exp,".csv"))

#Average Cq for the housekeeping gene
hk_Cq <- avg_Cq %>%
  filter(Gene == hk_gene)

#write.csv(hk_Cq, paste0(time,"_hk_Cq_",exp,".csv"))

#Average Cq for GOI
GOI_Cq <- avg_Cq %>%
  filter(Gene != hk_gene)

#write.csv(GOI_Cq, paste0(time,"_GOI_Cq_",exp,".csv"))

#####Calculate delta and standard error on delta#####
hk_Cq_chem <- data.frame(hk_Cq["Chemical"])
hk_Cq_conc <- data.frame(hk_Cq["Concentration"])

n1 <- nrow(hk_Cq_chem)
np <- nrow(plate)
np2 <- nrow(plate)

#Add new columns: delta and Std_error_dt
GOI_delta <- GOI_Cq %>%
  group_by(Chemical, Concentration) %>%
  mutate(delta=0, Std_error_dt=0)

#Calculate delta

while(n1!=0){
  n2 <- nrow(genes)
  
  a <- as.numeric(hk_Cq[hk_Cq$Plate==plate[np2,] & hk_Cq$Chemical==hk_Cq_chem[n1,] & hk_Cq$Concentration==hk_Cq_conc[n1,], "Avg_Cq"])
  GOI_delta[GOI_delta$Plate==plate[np2,] & GOI_delta$Chemical==hk_Cq_chem[n1,] & GOI_delta$Concentration==hk_Cq_conc[n1,],"delta"] <- GOI_delta[GOI_delta$Plate==plate[np2,] & GOI_delta$Chemical==hk_Cq_chem[n1,] & GOI_delta$Concentration==hk_Cq_conc[n1,],"Avg_Cq"] - a

#Calculate standard error on delta
  while(n2!=0){
    if(num_plate!=1){
      while(num_plate!=0){
        b <- as.numeric(GOI_delta[GOI_delta$Chemical==hk_Cq_chem[n1,] & GOI_delta$Concentration==hk_Cq_conc[n1,] & GOI_delta$Gene==genes[n2,] & GOI_delta$Plate==plate[num_plate,], "Std_error"])
        c <- as.numeric(hk_Cq[hk_Cq$Chemical==hk_Cq_chem[n1,] & hk_Cq$Concentration==hk_Cq_conc[n1,] & hk_Cq$Plate==plate[num_plate,], "Std_error"])
        
        GOI_delta[GOI_delta$Chemical==hk_Cq_chem[n1,] & GOI_delta$Concentration==hk_Cq_conc[n1,] & GOI_delta$Gene==genes[n2,] & GOI_delta$Plate==plate[num_plate,],"Std_error_dt"] <- sqrt(b^2+c^2)
        
        num_plate <- num_plate-1
      }
#if there is more than one dataset/plate included in the overall data  
    }else{
        
    b <- as.numeric(GOI_delta[GOI_delta$Chemical==hk_Cq_chem[n1,] & GOI_delta$Concentration==hk_Cq_conc[n1,] & GOI_delta$Gene==genes[n2,], "Std_error"])
    c <- as.numeric(hk_Cq[hk_Cq$Chemical==hk_Cq_chem[n1,] & hk_Cq$Concentration==hk_Cq_conc[n1,], "Std_error"])
  
    GOI_delta[GOI_delta$Chemical==hk_Cq_chem[n1,] & GOI_delta$Concentration==hk_Cq_conc[n1,] & GOI_delta$Gene==genes[n2,],"Std_error_dt"] <- sqrt(b^2+c^2)
    }
    n2 <- n2-1
    num_plate <- nrow(plate)
  }
  
  if(np>=1){
    n1 <- n1-1
    
    if(n1==0 & np2!=0){
      np2 <- np2-1
      n1 <- nrow(hk_Cq_chem)
      if(np2==0){
        n1 <- 0
      }
    }
  }else{
    n1 <- n1-1
  }
}

write.csv(GOI_delta, paste0(time,"_GOI_delta_", paste0(exp_plate, collapse = "_"),"_",exp,".csv"))

treat_delta <- GOI_delta %>%
  filter(Chemical!=vc)

ctrl_delta <- GOI_delta %>%
  filter(Chemical==vc) %>%
  subset(select=-c(Concentration)) %>%
  mutate(delta_delta=delta-delta, Std_error_dtdt=sqrt(2*Std_error_dt^2), Concentration=vc)

delta_delta <- treat_delta %>%
  group_by(Gene) %>%
  mutate(delta_delta=0, Std_error_dtdt=0)

#Calculate delta-delta
n <- nrow(genes)
num_plate <- nrow(plate)
n2 <- nrow(plate)

while(n>0 & n2>0){
  
  a <- lapply(ctrl_delta[ctrl_delta$Plate==plate[n2,] & ctrl_delta$Gene==genes[n,],"delta"], as.numeric)
  delta_delta[delta_delta$Plate==plate[n2,] & delta_delta$Gene==genes[n,],"delta_delta"] <- delta_delta[delta_delta$Plate==plate[n2,] & delta_delta$Gene==genes[n,],"delta"] - a

  if(nrow(genes)>1){
    n <- n-1
    n2 <- n2-1
    if(n2==0){
      n2 <- nrow(plate)
    }
    if(n==0){
      n2 <- 0
    }
  }else{ 
    n2 <- n2-1
    if(n2==0){
      n <- 0
      }
  }
}

#Calculate standard error on delta-delta
n <- nrow(genes)
n1 <- nrow(hk_Cq_chem)
num_plate <- nrow(plate)
n2 <- nrow(plate)

if(num_plate==1){
    n3 <- nrow(GOI_delta)
  }

if(n2>1){
      
  while(num_plate>0){
    b <- as.numeric(delta_delta[delta_delta$Chemical==hk_Cq_chem[n1,] & delta_delta$Concentration==hk_Cq_conc[n1,] & delta_delta$Gene==genes[n,] & delta_delta$Plate==plate[num_plate,], "Std_error_dt"])
    c <- as.numeric(ctrl_delta[ctrl_delta$Gene==genes[n,] & ctrl_delta$Plate==plate[num_plate,], "Std_error_dt"])
      
    delta_delta[delta_delta$Chemical==hk_Cq_chem[n1,] & delta_delta$Concentration==hk_Cq_conc[n1,] & delta_delta$Gene==genes[n,] & delta_delta$Plate==plate[num_plate,],"Std_error_dtdt"] <- sqrt(b^2+c^2)
        
    if(nrow(genes)>1){
      n1 <- n1-1
      if(n1==0){
        n <- n-1
        n1 <- nrow(hk_Cq_chem)
      }
      if(n==0&num_plate!=0){
        n <- nrow(genes)
        n1 <- nrow(hk_Cq_chem)
        num_plate <- num_plate-1
      }
    }else{
      n1 <- n1-1
      if(n1==0&num_plate!=0){
        num_plate <- num_plate-1
        n1 <- nrow(hk_Cq_chem)
      }
    }
  }
  
  }else{
  
    while(n1!=0){  
    b <- as.numeric(GOI_delta[GOI_delta$Chemical==hk_Cq_chem[n1,] & GOI_delta$Concentration==hk_Cq_conc[n1,] & GOI_delta$Gene==genes[n,], "Std_error_dt"])
    c <- as.numeric(ctrl_delta[ctrl_delta$Gene==genes[n,], "Std_error_dt"])
    
    delta_delta[delta_delta$Chemical==hk_Cq_chem[n1,] & delta_delta$Concentration==hk_Cq_conc[n1,] & delta_delta$Gene==genes[n,],"Std_error_dtdt"] <- sqrt(b^2+c^2)
    
    n1 <- n1-1
    n3 <- n3-1
  
    if(n1==0){
      n <- n-1
    }
    if(n1==0 & n3!=0){
      n1 <- nrow(hk_Cq_chem)
    }
    if(n3==0){
      n1 <- 0
    }
  }
}

#Calculate fold change
ctrl_delta <- ctrl_delta %>%
  subset(select=-c(Concentration)) %>%
  mutate(fold_change=2^-delta_delta, fc_stderr=0, Concentration="0")

ctrl_delta_neg <- ctrl_delta %>%
  mutate(fold_change=-1)
  
delta_delta <- delta_delta %>%
  mutate(fold_change=2^-delta_delta) %>%
  mutate(fc_stderr=2^-Std_error_dtdt)

#Generate summary Excel file and graph
write.csv(delta_delta, paste0(time,"_delta_delta_", paste0(exp_plate, collapse = "_"),"_",exp, ".csv"))

#Negative reciprocals of fold changes <1
delta_delta_1 <- delta_delta %>%
  filter(fold_change<1) %>%
  mutate(fold_change = -1/fold_change)

delta_delta <- delta_delta %>%
  filter(fold_change>1)%>%
  rbind(ctrl_delta, ctrl_delta_neg, delta_delta_1)

#Exclude certain chemicals e.g. BaP
delta_delta_nBaP <- delta_delta %>%
  filter(fold_change>1)%>%
  rbind(ctrl_delta, ctrl_delta_neg, delta_delta_1) %>%
  filter(Chemical!="0h") %>% filter (Chemical!="BaP")

# Fold change bar graphs with st error bars
ggsave(filename = paste0(time, "_", paste0(genes_str, collapse = "_"), "_",paste0(exp_plate, collapse = "_"),"_", exp, ".pdf"),
       ggplot(data = delta_delta, mapping = aes(x = Concentration, y = fold_change, fill = Chemical)) + 
         geom_bar(stat = "identity") + 
         geom_errorbar(aes(ymin=fold_change-fc_stderr, ymax=fold_change+fc_stderr), width=.2,
                       position=position_dodge(.9)) +
         labs(title = paste0(time, " ", paste0(genes_str, collapse = " "), " expression levels"),
              x = "Concentration", 
              y = "Fold change relative to vehicle ctrl",
              fill = "Chemical") +
         theme(text = element_text(size = 10),
               axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1)) +
         facet_grid(Plate~Gene),
       width = 25, height = 15, units = "cm")

#Exclude certain chemicals e.g. BaP
ggsave(filename = paste0(time, "_", paste0(genes_str, collapse = "_"), "_",paste0(exp_plate, collapse = "_"),"_", exp, "_noBAP.pdf"),
       ggplot(data = delta_delta_nBaP, mapping = aes(x = Concentration, y = fold_change, fill = Chemical)) + 
         geom_bar(stat = "identity") + 
         geom_errorbar(aes(ymin=fold_change-fc_stderr, ymax=fold_change+fc_stderr), width=.2,
                 position = position_dodge(.9)) +
         labs(title = paste0(time, " ", paste0(genes_str, collapse = " "), " expression levels"),
              x = "Concentration", 
              y = "Fold change relative to vehicle ctrl",
              fill = "Chemical") +
         theme(text = element_text(size = 10),
               axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1)) +
         facet_grid(Plate~Gene),
       width = 25, height = 15, units = "cm")


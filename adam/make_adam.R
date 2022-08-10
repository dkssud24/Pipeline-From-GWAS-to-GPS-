setwd(getwd())

library(bigstatsr)
library(bigsnpr)
library(dplyr)
library(data.table)
library(magrittr)
library(R.utils)
library(stringr)


hwe_list <- fread("test_mode.hwe")

hwe_list2 <- list()
hwe_list2 <- str_split(hwe_list[,6],by="/")
hwe_list3 <- do.call("rbind",hwe_list2)
hwe_list4 <- as.numeric(hwe_list3)
dim(hwe_list4) <- c(1141242,3)
names(hwe_list4) <- c("N1","N2","N3")
hwe_list5 <- cbind(hwe_list,hwe_list4)


# Count maximum N
hwe_list5$N_max <- apply(hwe_list5[,10:12],1,max)
hwe_list5 <- hwe_list5 %>% mutate(count = ifelse(N_max == N1,1,ifelse(N_max==N1,1,0)))

hwe_list6 <- hwe_list5 %>% mutate(geno1 = ifelse(count_min==0,A2,ifelse(count_min==1,A2,A1)),geno2= ifelse(count_min==0,A2,A1))
geno6 <- hwe_list6[,c("geno1","geno2")]
geno7 <- c(t(geno6))
geno8 <- cbind(ped,t(data.frame(geno7)))


fwrite(geno8,'adam.ped',quote=F,sep='\t',row.names=F,na="NA",col.names=F)

#!/usr/bin/env Rscript

argv = commandArgs(trailingOnly=TRUE)

action = argv[1]
env = argv[2]
traitID = argv[3]
name = argv[4]

library(data.table)
library(dplyr)

a <- fread("root_discovery.txt")
a <- data.frame(a)

b <- fread(env)
b <- data.frame(b)

int <- function(x){
        qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

c <- select(b,c('eid',traitID))

if(action == "--asthma"){
        names(c)[1] <- "asthma_eid"
        d <- left_join(a,c,by="asthma_eid")
        d <- d  %>% mutate(trans = int(get(traitID)))
        d <- na.omit(d)
        e <- d[,c(1,1,3:18)]
        names(e)[2] <- "FID"
        fwrite(e,paste(name,"_dis_raw.txt",sep=""),quote=F,sep='\t',row.names=F,na="NA")
        f <- e[,which(colnames(e)==traitID)]
        print(summary(f))
        print(sd(f))

        pheno <- e[,c(1,1,18)]
        names(pheno)[2] <- "IID"
        cov <- e[,c(1:16)]
        names(cov)[2] <- "IID"

        fwrite(pheno,paste(name,'_dis.pheno',sep=""),quote=F,sep='\t',row.names=F,na="NA")
        fwrite(cov,paste(name,'_dis.cov',sep=""),quote=F,sep='\t',row.names=F,na="NA")
}else if(action == "--ob"){
        names(c)[1] <- "FID"
        d <- left_join(a,c,by="FID")
        d <- d %>% mutate(trans = int(get(traitID)))
        d <- na.omit(d)
        e <- d[,c(1,1,3:18)]
        names(e)[2] <- "FID"
        fwrite(e,paste(name,'_dis_raw.txt',sep=""),quote=F,sep='\t',row.names=F,na="NA")
        
        f <- e[,which(colnames(e)==traitID)]
        print(summary(f))
        print(sd(f))
        
        pheno <- e[,c(1,1,18)]
        names(pheno)[2] <- "IID"
        cov <- e[,c(1:16)]
        names(cov)[2] <- "IID"
        fwrite(pheno,paste(name,'_dis.pheno',sep=""),quote=F,sep='\t',row.names=F,na="NA")
        fwrite(cov,paste(name,'_dis.cov',sep=""),quote=F,sep='\t',row.names=F,na="NA")
}
else{print("option error")}


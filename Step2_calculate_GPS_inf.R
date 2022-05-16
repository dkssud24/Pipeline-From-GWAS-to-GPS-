setwd(getwd())

library(bigstatsr)
library(bigsnpr)
library(dplyr)
library(data.table)
library(magrittr)
library(R.utils)

## Create correlation matrix

snp_readBed("QCref2.bed")
> info <- snp_attach("QCref2.rds")
> str(info)
##  List of 3
##   $ genotypes:Reference class 'FBM.code256' [package "bigstatsr"] with 16 fields
##    ..$ extptr      :<externalptr> 
##    ..$ extptr_rw   :<externalptr> 
##    ..$ nrow        : int 404
##    ..$ ncol        : int 1141242
##    ..$ type        : Named int 1
##    .. ..- attr(*, "names")= chr "unsigned char"
##    ..$ backingfile : chr "/home/adam/BiO/enju07/9_md/0_root/ref_pannel_bedbimfam/QCref2.bk"
##    ..$ is_read_only: logi FALSE
##    ..$ address     :<externalptr> 
##    ..$ address_rw  :<externalptr> 
##    ..$ bk          : chr "/home/adam/BiO/enju07/9_md/0_root/ref_pannel_bedbimfam/QCref2.bk"
##    ..$ rds         : chr "/home/adam/BiO/enju07/9_md/0_root/ref_pannel_bedbimfam/QCref2.rds"
##    ..$ is_saved    : logi TRUE
##    ..$ type_chr    : chr "unsigned char"
##    ..$ type_size   : int 1
##    ..$ file_size   : num 4.61e+08
##    ..$ code256     : num [1:256] 0 1 2 NA NA NA NA NA NA NA ...
##    ..and 26 methods, of which 12 are  possibly relevant:
##    ..  add_columns, as.FBM, bm, bm.desc, check_dimensions,
##    ..  check_write_permissions, copy#envRefClass, initialize, initialize#FBM,
##    ..  save, show#envRefClass, show#FBM
##   $ fam      :'data.frame':	404 obs. of  6 variables:
##    ..$ family.ID  : chr [1:404] "HG00096" "HG00097" "HG00099" "HG00100" ...
##    ..$ sample.ID  : chr [1:404] "HG00096" "HG00097" "HG00099" "HG00100" ...
##    ..$ paternal.ID: int [1:404] 0 0 0 0 0 0 0 0 0 0 ...
##    ..$ maternal.ID: int [1:404] 0 0 0 0 0 0 0 0 0 0 ...
##    ..$ sex        : int [1:404] 0 0 0 0 0 0 0 0 0 0 ...
##    ..$ affection  : int [1:404] -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 ...
##   $ map      :'data.frame':	1141242 obs. of  6 variables:
##    ..$ chromosome  : int [1:1141242] 1 1 1 1 1 1 1 1 1 1 ...
##    ..$ marker.ID   : chr [1:1141242] "rs3115860" "rs3131969" "rs2286139" "rs12562034" ...
##    ..$ genetic.dist: num [1:1141242] 0.489 0.49 0.493 0.496 0.501 ...
##    ..$ physical.pos: int [1:1141242] 753405 754182 761732 768448 779322 785989 798959 838555 846808 853954 ...
##    ..$ allele1     : chr [1:1141242] "C" "A" "C" "A" ...
##    ..$ allele2     : chr [1:1141242] "A" "G" "T" "G" ...
##   - attr(*, "class")= chr "bigSNP"

sumstats2 <- sumstats[,c(1,3,2,5,4,9,10,12,8)]
## Classes ??data.table?? and 'data.frame':	1352168 obs. of  9 variables:
##  $ #CHROM: int  1 1 1 1 1 1 1 1 1 1 ...
##  $ ID    : chr  "rs6650104" "rs11510103" "rs12565286" "rs3094315" ...
##  $ POS   : int  564477 567753 721290 752566 752721 753405 754182 760912 761147 761732 ...
##  $ ALT   : chr  "G" "G" "C" "G" ...
##  $ REF   : chr  "A" "A" "G" "A" ...
##  $ BETA  : num  NA NA -0.00646 -0.00545 -0.00581 ...
##  $ SE    : num  NA NA 0.00677 0.00332 0.00331 ...
##  $ P     : num  NA NA 0.34 0.1007 0.0794 ...
##  $ OBS_CT: int  174080 174083 172285 173675 174083 173019 172909 172751 172700 166497 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
names(sumstats2) <- c("chr","rsid","pos","a0","a1","beta","beta_se","p","n_eff")

info_snp <- snp_match(sumstats2, map)
1,352,168 variants to be matched.
0 ambiguous SNPs have been removed.
1,141,242 variants have been matched; 0 were flipped and 22,078 were reversed.


NCORES <- nb_cores()
tmp <- tempfile(tmpdir = "/home/adam/BiO/enju07/9_md/0_root/ref_pannel_bedbimfam/tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
corr <- NULL
ld <- NULL
fam.order <- NULL
CHR <- map$chr
chr <- map$chr
POS <- map$pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
corr0 <- NULL
CHR <- info$map$chromosome
POS <- info$map$physical.pos
G2 <- snp_fastImputeSimple(G)

for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
            G2,
            ind.col = ind.chr2,
            ncores = NCORES,
            infos.pos = POS2[ind.chr2],
            size = 3 / 1000
        )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc( ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)




snp_readBed("match_mer_vali.bed")
obj.bigSNP <- snp_attach("match_mer_vali.rds")
G_ht <- obj.bigSNP$genotypes
CHR_ht <- obj.bigSNP$map$chromosome
POS_ht <- obj.bigSNP$map$physical.pos
G2_ht <- snp_fastImputeSimple(G_ht)     //30min 정도 걸림

map_ht <- obj.bigSNP$map[,-c(2:3)]
names(map_ht) <- c("chr", "pos", "a0", "a1")
info_snp_ht <- snp_match(sumstats2, map_ht)


beta_inf <- snp_ldpred2_inf(corr,df_beta,h2 = h2_est)
ind.test <- 1:nrow(G2)
pred_inf <- big_prodVec( G2, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)

fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
fam.order2 <- fam.order[,c(1,2)]
score <- cbind(fam.order2,pred_inf)

phenotype <- fread("../pheno/Standing_height_vali_raw.txt")
pheno <- left_join(phenotype,score,by=c("FID","IID"))

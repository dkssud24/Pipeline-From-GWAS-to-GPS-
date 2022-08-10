#!/usr/bin/perl

use warnings;
use strict;

my $string;
my $string2;
for( my $i = 6; $i <= 22; $i = $i+1){
        $string = $i."_v2";
        $string2 = $i."_v2_s488288";

`plink --bed /BiO/00_original_data/2_2019_UK_biobank/download_genotypes/ukb_cal_chr$string.bed --fam /BiO/00_original_data/2_2019_UK_biobank/download_genotypes/ukb48422_cal_chr$string2.fam --bim /BiO/00_original_data/2_2019_UK_biobank/download_genotypes/ukb_snp_chr$string.bim --keep ../final/base_with_eve.fam --make-bed --out ukb_geno_adam_set_chr$i`}

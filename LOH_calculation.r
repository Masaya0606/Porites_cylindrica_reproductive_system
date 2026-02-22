#bcftools_LOH_vcf.sh CloneGroup2 and 7 このファイルがbcftools_LOH.shで計算している
#/Users/moritamasaya/NGS_analyses/Porites_LOHにあるbcftools_LOH_vcf.shで計算したものである。
#/Users/moritamasaya/NGS_analyses/Porites_LOHにあるbcftools_LOH_vcf.shで計算したものである。
#NIGのサーバーで計算しており、MAF0.05 geno0.1 HWE 10^6でフィルタリングしたVCFファイルを利用した。SNP数はおよそ6800

#CloneGourp7
df <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_7_M4_MH2_F2_GT.txt", header = FALSE,
                 col.names = c("CHROM","POS","M1","M2","M3","M4","MH1","MH2","F1","F2"))

# もし missing がありそうなら、先に除外
non_missing <- (df$M1 !="./.") &(df$M2 !="./.") & (df$M3 !="./.") & (df$M4 !="./.") & (df$MH1 != "./.") & (df$MH2 != "./.") & (df$F1 != "./.") & (df$F2 != "./.")
# Mがヘテロの座位
M_het <- (df$M1 == "0/1") | (df$M2 == "0/1")| (df$M3 == "0/1")| (df$M4 == "0/1")&non_missing
# MHがヘテロの座位
MH_het <- (df$MH1 == "0/1") | (df$MH2 == "0/1")&non_missing

# Fが全ホモの座位
F_hom <- (df$F1 %in% c("0/0","1/1")) & (df$F2 %in% c("0/0","1/1"))&non_missing
# Mが全ホモの座位
M_hom <- (df$M1 %in% c("0/0","1/1")) & (df$M2 %in% c("0/0","1/1"))& (df$M3 %in% c("0/0","1/1"))& (df$M4 %in% c("0/0","1/1"))&non_missing


# LOH発生
LOH_site <- MH_het & F_hom
LOH_site_M<-MH_het & M_hom

# LOH 率
N_total      <- sum(non_missing)          # 比較に使えたSNP総数
N_MH_het     <- sum(MH_het)               # 分母：MH側ヘテロ座位数
N_M_het      <-sum(M_het)
N_LOH_MHF        <- sum(LOH_site)             # 分子：LOH座位数
N_LOH_MF        <- sum(LOH_site_M)

LOH_rate_MH_F <- N_LOH_MHF / N_MH_het
LOH_rate_M_F <- N_LOH_MF / N_MH_het

N_total
N_MH_het
N_M_het
N_LOH_MHF
N_LOH_MF
LOH_rate_MH_F
LOH_rate_M_F

#pairwiseLOH
# LOH率
calc_loh_pair <- function(mh_gt, f_gt) {
  non_missing <- (mh_gt != "./.") & (f_gt != "./.")
  mh_het <- (mh_gt == "0/1") & non_missing
  f_hom <- (f_gt %in% c("0/0","1/1")) & non_missing
  
  denom <- sum(mh_het)
  num   <- sum(mh_het & f_hom)
  
  c(LOH_rate = ifelse(denom > 0, num / denom, NA),
    N_het_parent = denom,
    N_LOH = num)
}

loh_MH1_F1 <- calc_loh_pair(df$MH1, df$F1)
loh_MH2_F1 <- calc_loh_pair(df$MH2, df$F1)
loh_MH1_F2 <- calc_loh_pair(df$MH1, df$F2)
loh_MH2_F2 <- calc_loh_pair(df$MH2, df$F2)
loh_MH1_M1 <- calc_loh_pair(df$MH1, df$M1)
loh_MH1_M2 <- calc_loh_pair(df$MH1, df$M2)
loh_MH1_M3 <- calc_loh_pair(df$MH1, df$M3)
loh_MH1_M4 <- calc_loh_pair(df$MH1, df$M4)
loh_MH2_M1 <- calc_loh_pair(df$MH2, df$M1)
loh_MH2_M2 <- calc_loh_pair(df$MH2, df$M2)
loh_MH2_M3 <- calc_loh_pair(df$MH2, df$M3)
loh_MH2_M4 <- calc_loh_pair(df$MH2, df$M4)
loh_M1_F1 <- calc_loh_pair(df$M1, df$F1)
loh_M2_F1 <- calc_loh_pair(df$M2, df$F1)
loh_M3_F1 <- calc_loh_pair(df$M3, df$F1)
loh_M4_F1 <- calc_loh_pair(df$M4, df$F1)
loh_M1_F2 <- calc_loh_pair(df$M1, df$F2)
loh_M2_F2 <- calc_loh_pair(df$M2, df$F2)
loh_M3_F2 <- calc_loh_pair(df$M3, df$F2)
loh_M4_F2 <- calc_loh_pair(df$M4, df$F2)

loh_MH1_F1
loh_MH2_F1
loh_MH1_F2
loh_MH2_F2
loh_MH1_M1
loh_MH1_M2
loh_MH1_M3
loh_MH1_M4
loh_MH2_M1
loh_MH2_M2
loh_MH2_M3
loh_MH2_M4
loh_M1_F1
loh_M2_F1
loh_M3_F1
loh_M4_F1
loh_M1_F2
loh_M2_F2
loh_M3_F2
loh_M4_F2

#reclaculating
df7 <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_7_M4_MH2_F2_GT.txt",
                  header = FALSE,
                  col.names = c("CHROM","POS","M1","M2","M3","M4","MH1","MH2","F1","F2"))
calc_loh_pair <- function(mh_gt, f_gt) {
  non_missing <- (mh_gt != "./.") & (f_gt != "./.")
  mh_het <- (mh_gt == "0/1") & non_missing
  f_hom <- (f_gt %in% c("0/0","1/1")) & non_missing
  
  denom <- sum(mh_het)
  num   <- sum(mh_het & f_hom)
  
  c(LOH_rate = ifelse(denom > 0, num / denom, NA),
    N_het_parent = denom,
    N_LOH = num)
}

make_loh_row <- function(parent_name, offspring_name,
                         parent_gt, offspring_gt,
                         clone_group = "G2") {
  res <- calc_loh_pair(parent_gt, offspring_gt)
  data.frame(
    CloneGroup    = clone_group,
    Parent        = parent_name,
    Offspring     = offspring_name,
    OffspringType = ifelse(substr(offspring_name, 1, 1) == "MH", "F", "M"),
    LOH_rate      = res["LOH_rate"],
    N_het_parent  = res["N_het_parent"],
    N_LOH         = res["N_LOH"],
    row.names = NULL
  )
}
library(dplyr)
loh_G7 <- rbind(
  # MH → F
  make_loh_row("MH1", "F1", df7$MH1, df7$F1, "Group7"),
  make_loh_row("MH2", "F1", df7$MH2, df7$F1, "Group7"),
  make_loh_row("MH1", "F2", df7$MH1, df7$F2, "Group7"),
  make_loh_row("MH2", "F2", df7$MH2, df7$F2, "Group7"),
  
  # MH → M1–M4
  make_loh_row("MH1", "M1", df7$MH1, df7$M1, "Group7"),
  make_loh_row("MH1", "M2", df7$MH1, df7$M2, "Group7"),
  make_loh_row("MH1", "M3", df7$MH1, df7$M3, "Group7"),
  make_loh_row("MH1", "M4", df7$MH1, df7$M4, "Group7"),
  make_loh_row("MH2", "M1", df7$MH2, df7$M1, "Group7"),
  make_loh_row("MH2", "M2", df7$MH2, df7$M2, "Group7"),
  make_loh_row("MH2", "M3", df7$MH2, df7$M3, "Group7"),
  make_loh_row("MH2", "M4", df7$MH2, df7$M4, "Group7"),
  # F → M1–M4
  make_loh_row("F1", "M1", df7$MH1, df7$M1, "Group7"),
  make_loh_row("F1", "M2", df7$MH1, df7$M2, "Group7"),
  make_loh_row("F1", "M3", df7$MH1, df7$M3, "Group7"),
  make_loh_row("F1", "M4", df7$MH1, df7$M4, "Group7"),
  make_loh_row("F2", "M1", df7$MH2, df7$M1, "Group7"),
  make_loh_row("F2", "M2", df7$MH2, df7$M2, "Group7"),
  make_loh_row("F2", "M3", df7$MH2, df7$M3, "Group7"),
  make_loh_row("F2", "M4", df7$MH2, df7$M4, "Group7")
)

loh_G7



#CloneGroup2
df <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_2_M2_MH3_F1_GT.txt", header = FALSE,
                 col.names = c("CHROM","POS","M1","M2","MH1","MH2","MH3","F1"))

# removing loci if missing
non_missing <- (df$M1 !="./.") &(df$M2 !="./.") & (df$MH1 != "./.") & (df$MH2 != "./.") & (df$MH3 != "./.") & (df$F1 != "./.")

# heterozygous loci in M
M_het <- (df$M1 == "0/1") | (df$M2 == "0/1") & non_missing

# heterozygous loci in MH
MH_het <- (df$MH1 == "0/1") | (df$MH2 == "0/1")|(df$MH3 == "0/1") & non_missing

# heterozygous loci in F
F_het <- (df$F1 == "0/1")  & non_missing

# homozygous loci in F
F_hom <- (df$F1 %in% c("0/0","1/1")) & non_missing
M_hom <- (df$M1 %in% c("0/0","1/1")) & (df$M2 %in% c("0/0","1/1")) & non_missing
MH_hom <- (df$MH1 %in% c("0/0","1/1")) & (df$MH2 %in% c("0/0","1/1")) & (df$MH3 %in% c("0/0","1/1")) & non_missing

# LOH発生
LOH_site_MH_F <- MH_het & F_hom
LOH_site_MH_M <- MH_het & M_hom
LOH_site_F_M <- F_het & M_hom
LOH_site_F_MH <- F_het & MH_hom
N_total      <- sum(non_missing)          # 比較に使えたSNP総数
N_MH_het     <- sum(MH_het)               # 分母：MH側ヘテロ座位数
N_M_het      <-sum(M_het)
N_F_het      <-sum(F_het)
N_LOH_MH_F        <- sum(LOH_site_MH_F)             # 分子：LOH座位数
N_LOH_MH_M        <- sum(LOH_site_MH_M)             # 分子：LOH座位数
N_LOH_F_M        <- sum(LOH_site_F_M)             # 分子：LOH座位数
N_LOH_F_MH        <- sum(LOH_site_F_MH)             # 分子：LOH座位数
LOH_rate_MH_F     <- N_LOH_MH_F / N_MH_het
LOH_rate_MH_M     <- N_LOH_MH_M / N_MH_het
LOH_rate_F_M     <- N_LOH_F_M / N_F_het
LOH_rate_F_MH     <- N_LOH_F_MH / N_F_het

N_total
N_MH_het
N_M_het
N_F_het
N_LOH_MH_F
N_LOH_MH_M
N_LOH_F_M
N_LOH_F_MH
LOH_rate_MH_F
LOH_rate_MH_M
LOH_rate_F_M
LOH_rate_F_MH

# LOH rate
calc_loh_pair <- function(mh_gt, f_gt) {
  non_missing <- (mh_gt != "./.") & (f_gt != "./.")
  mh_het <- (mh_gt == "0/1") & non_missing
  f_hom <- (f_gt %in% c("0/0","1/1")) & non_missing
  
  denom <- sum(mh_het)
  num   <- sum(mh_het & f_hom)
  
  c(LOH_rate = ifelse(denom > 0, num / denom, NA),
    N_het_parent = denom,
    N_LOH = num)
}

"M1","M2","MH1","MH2","MH3","F1"
Group2_loh_MH1_F1 <- calc_loh_pair(df$MH1, df$F1)
Group2_loh_MH2_F1 <- calc_loh_pair(df$MH2, df$F1)
Group2_loh_MH3_F1 <- calc_loh_pair(df$MH3, df$F1)
Group2_loh_MH1_M1 <- calc_loh_pair(df$MH1, df$M1)
Group2_loh_MH2_M1 <- calc_loh_pair(df$MH2, df$M1)
Group2_loh_MH3_M1 <- calc_loh_pair(df$MH3, df$M1)
Group2_loh_MH1_M1 <- calc_loh_pair(df$MH1, df$M2)
Group2_loh_MH2_M1 <- calc_loh_pair(df$MH2, df$M2)
Group2_loh_MH3_M1 <- calc_loh_pair(df$MH3, df$M2)
Group2_loh_F1_M1 <- calc_loh_pair(df$F1, df$M1)
Group2_loh_F1_M2 <- calc_loh_pair(df$F1, df$M2)
Group2_loh_F1_MH1 <- calc_loh_pair(df$F1, df$MH1)
Group2_loh_F1_MH2 <- calc_loh_pair(df$F1, df$MH2)
Group2_loh_F1_MH3 <- calc_loh_pair(df$F1, df$MH3)
Group2_loh_M1_F1 <- calc_loh_pair(df$M1, df$F1)
Group2_loh_M2_F1 <- calc_loh_pair(df$M2, df$F1)
Group2_loh_M1_MH1 <- calc_loh_pair(df$M1, df$MH1)
Group2_loh_M1_MH2 <- calc_loh_pair(df$M1, df$MH2)
Group2_loh_M1_MH3 <- calc_loh_pair(df$M1, df$MH3)
Group2_loh_M2_MH1 <- calc_loh_pair(df$M2, df$MH1)
Group2_loh_M2_MH2 <- calc_loh_pair(df$M2, df$MH2)
Group2_loh_M2_MH3 <- calc_loh_pair(df$M2, df$MH3)

Group2_loh_MH1_F1
Group2_loh_MH2_F1
Group2_loh_MH3_F1
Group2_loh_MH1_M1
Group2_loh_MH2_M1
Group2_loh_MH3_M1
Group2_loh_MH1_M1
Group2_loh_MH2_M1
Group2_loh_MH3_M1
Group2_loh_F1_M1
Group2_loh_F1_M2
Group2_loh_F1_MH1
Group2_loh_F1_MH2
Group2_loh_F1_MH3
Group2_loh_M1_F1
Group2_loh_M2_F1
Group2_loh_M1_MH1
Group2_loh_M1_MH2
Group2_loh_M1_MH3
Group2_loh_M2_MH1
Group2_loh_M2_MH2
Group2_loh_M2_MH3

df <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_2_M2_MH3_F1_GT.txt",
                 header = FALSE,
                 col.names = c("CHROM","POS","M1","M2","MH1","MH2","MH3","F1"))

# Group 2 の LOH 結果をまとめる
loh_G2 <- rbind(
  # MH → F
  make_loh_row("MH1", "F1", df$MH1, df$F1, "Group2"),
  make_loh_row("MH2", "F1", df$MH2, df$F1, "Group2"),
  make_loh_row("MH3", "F1", df$MH3, df$F1, "Group2"),
  
  # MH → M1
  make_loh_row("MH1", "M1", df$MH1, df$M1, "Group2"),
  make_loh_row("MH2", "M1", df$MH2, df$M1, "Group2"),
  make_loh_row("MH3", "M1", df$MH3, df$M1, "Group2"),
  
  # MH → M2
  make_loh_row("MH1", "M2", df$MH1, df$M2, "Group2"),
  make_loh_row("MH2", "M2", df$MH2, df$M2, "Group2"),
  make_loh_row("MH3", "M2", df$MH3, df$M2, "Group2")
)

loh_G2

loh_all <- rbind(loh_G2, loh_G7)

ggplot(loh_all, aes(x = OffspringType, y = LOH_rate,
                    color = CloneGroup)) +
  geom_jitter(width = 0.1, size = 4, alpha = 0.9) +
  geom_hline(yintercept = 0.05, 
             linetype = "dashed", 
             color = "black", 
             linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(
    x = "Offspring type",
    y = "LOH rate (MH → offspring)",
    color = "Clone group"
  )

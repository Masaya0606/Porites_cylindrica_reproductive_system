#CloneGroup2
#/Users/moritamasaya/NGS_analyses/Porites_LOHにあるbcftools_LOH_vcf.shで計算したものである。
#NIGのサーバーで計算しており、MAF0.05 geno0.1 HWE 10^6でフィルタリングしたVCFファイルを利用した。SNP数はおよそ6800

df <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_2_M2_MH3_F1_GT.txt", header = FALSE,
                 col.names = c("CHROM","POS","M1","M2","MH1","MH2","MH3","F1")) 
samples<-c("M1","M2","MH1","MH2","MH3","F1")
calc_Hobs <- function(gt_vec) {
  non_missing <- gt_vec != "./."
  het <- gt_vec == "0/1"
  sum(het & non_missing) / sum(non_missing)
}

Hobs_vec <- sapply(df[samples], calc_Hobs)
Hobs_vec

Hobs_df_group2 <- data.frame(
  Sample = names(Hobs_vec),
  Hobs   = as.numeric(Hobs_vec),
  stringsAsFactors = FALSE
)

Hobs_df_group2

Hobs_df_group2$Group <- NA

Hobs_df_group2$Group[Hobs_df_group2$Sample %in% c("M1","M2")]      <- "M"
Hobs_df_group2$Group[Hobs_df_group2$Sample %in% c("MH1","MH2","MH3")] <- "MH"
Hobs_df_group2$Group[Hobs_df_group2$Sample %in% c("F1")]           <- "F"

Hobs_df_group2$Group <- factor(Hobs_df_group2$Group, levels = c("M","MH","F"))
Hobs_df_group2

library(ggplot2)

ggplot(Hobs_df_group2, aes(x = Group, y = Hobs)) +
  geom_jitter(width = 0.1, height = 0, size = 2) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Observed heterozygosity (Hobs)") +
  xlab("Reproductive type")


#CloneGroup7
df <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_7_M4_MH2_F2_GT.txt", header = FALSE,
                 col.names = c("CHROM","POS","M1","M2","M3","M4","MH1","MH2","F1","F2"))

samples <- c("M1","M2","M3","M4","MH1","MH2","F1","F2")

calc_Hobs <- function(gt_vec) {
  non_missing <- gt_vec != "./."
  het <- gt_vec == "0/1"
  sum(het & non_missing) / sum(non_missing)
}

Hobs_vec <- sapply(df[samples], calc_Hobs)
Hobs_vec

Hobs_df_group7 <- data.frame(
  Sample = names(Hobs_vec),
  Hobs   = as.numeric(Hobs_vec),
  stringsAsFactors = FALSE
)

Hobs_df_group7

Hobs_df_group7$Group <- NA

Hobs_df_group7$Group[Hobs_df_group7$Sample %in% c("M1","M2","M3","M4")]      <- "M"
Hobs_df_group7$Group[Hobs_df_group7$Sample %in% c("MH1","MH2")] <- "MH"
Hobs_df_group7$Group[Hobs_df_group7$Sample %in% c("F1", "F2")]           <- "F"

Hobs_df_group7$Group <- factor(Hobs_df_group7$Group, levels = c("M","MH","F"))
Hobs_df_group7


ggplot(Hobs_df_group7, aes(x = Group, y = Hobs)) +
  geom_jitter(width = 0.1, height = 0, size = 2) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Observed heterozygosity (Hobs)") +
  xlab("Reproductive type")

#group2と7を一緒に
library(dplyr)

# --- 1: Clone group 名を追加 ---
Hobs_df_group2$CloneGroup <- "group2"
Hobs_df_group7$CloneGroup <- "group7"

# --- 2: データ結合 ---
Hobs_all <- bind_rows(Hobs_df_group2, Hobs_df_group7)

# --- 3: プロット ---
ggplot(Hobs_all, aes(x = Group, y = Hobs, color = CloneGroup)) +
  geom_jitter(width = 0.1, size = 4, alpha = 0.8) +
  scale_y_continuous(limits = c(0,0.5)) +
  theme_classic() +
  labs(
    x = "Reproductive category",
    y = "Observed heterozygosity (Hobs)",
    color = "Clone group"
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 14)
  )

#glmで検定
#group2
Hobs_group2_glm<-glm(Hobs~Group, data=Hobs_df_group2)
library(multcomp)
mult_comp_Hobs_group2_glm <- glht(Hobs_group2_glm, linfct = mcp(Group = "Tukey"))
summary(mult_comp_Hobs_group2_glm)
#Ns
#念の為
Hobs_group2_glm_null<-glm(Hobs~1, data=Hobs_df_group2)
anova(Hobs_group2_glm, Hobs_group2_glm_null)

#group7
Hobs_group7_glm<-glm(Hobs~Group, data=Hobs_df_group7)
mult_comp_Hobs_group7_glm <- glht(Hobs_group7_glm, linfct = mcp(Group = "Tukey"))
summary(mult_comp_Hobs_group7_glm)
#Ns
#念の為
Hobs_group7_glm_null<-glm(Hobs~1, data=Hobs_df_group7)
anova(Hobs_group7_glm, Hobs_group7_glm_null)


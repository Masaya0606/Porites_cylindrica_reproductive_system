#Observed heterozygosity (Hobs) was calculated for each individual as the proportion of heterozygous sites among non-missing SNPs, 
#and compared among reproductive categories using generalized linear models with Tukey’s post hoc tests.

# -----------------------------
# Clone Group 2
# -----------------------------
# These data were generated using bcftools_LOH_vcf.sh
# located in /Users/moritamasaya/NGS_analyses/Porites_LOH.
# Analyses were conducted on the NIG supercomputer.
# The input VCF file was filtered with MAF = 0.05, geno = 0.1,
# and HWE = 10^6. The total number of SNPs is approximately 6,800.

df <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_2_M2_MH3_F1_GT.txt",
                 header = FALSE,
                 col.names = c("CHROM","POS","M1","M2","MH1","MH2","MH3","F1"))

samples <- c("M1","M2","MH1","MH2","MH3","F1")

# Function to calculate observed heterozygosity (Hobs)
calc_Hobs <- function(gt_vec) {
  non_missing <- gt_vec != "./."
  het <- gt_vec == "0/1"
  sum(het & non_missing) / sum(non_missing)
}

# Calculate Hobs for each individual
Hobs_vec <- sapply(df[samples], calc_Hobs)
Hobs_vec

# Convert to data frame
Hobs_df_group2 <- data.frame(
  Sample = names(Hobs_vec),
  Hobs   = as.numeric(Hobs_vec),
  stringsAsFactors = FALSE
)

Hobs_df_group2

# Assign reproductive category
Hobs_df_group2$Group <- NA
Hobs_df_group2$Group[Hobs_df_group2$Sample %in% c("M1","M2")] <- "M"
Hobs_df_group2$Group[Hobs_df_group2$Sample %in% c("MH1","MH2","MH3")] <- "MH"
Hobs_df_group2$Group[Hobs_df_group2$Sample %in% c("F1")] <- "F"

Hobs_df_group2$Group <- factor(Hobs_df_group2$Group,
                               levels = c("M","MH","F"))
Hobs_df_group2

library(ggplot2)

# Plot Hobs for Clone Group 2
ggplot(Hobs_df_group2, aes(x = Group, y = Hobs)) +
  geom_jitter(width = 0.1, height = 0, size = 2) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Observed heterozygosity (Hobs)") +
  xlab("Reproductive type")


# -----------------------------
# Clone Group 7
# -----------------------------

df <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_7_M4_MH2_F2_GT.txt",
                 header = FALSE,
                 col.names = c("CHROM","POS","M1","M2","M3","M4","MH1","MH2","F1","F2"))

samples <- c("M1","M2","M3","M4","MH1","MH2","F1","F2")

# Calculate Hobs
Hobs_vec <- sapply(df[samples], calc_Hobs)
Hobs_vec

Hobs_df_group7 <- data.frame(
  Sample = names(Hobs_vec),
  Hobs   = as.numeric(Hobs_vec),
  stringsAsFactors = FALSE
)

Hobs_df_group7

# Assign reproductive category
Hobs_df_group7$Group <- NA
Hobs_df_group7$Group[Hobs_df_group7$Sample %in% c("M1","M2","M3","M4")] <- "M"
Hobs_df_group7$Group[Hobs_df_group7$Sample %in% c("MH1","MH2")] <- "MH"
Hobs_df_group7$Group[Hobs_df_group7$Sample %in% c("F1","F2")] <- "F"

Hobs_df_group7$Group <- factor(Hobs_df_group7$Group,
                               levels = c("M","MH","F"))
Hobs_df_group7

# Plot Hobs for Clone Group 7
ggplot(Hobs_df_group7, aes(x = Group, y = Hobs)) +
  geom_jitter(width = 0.1, height = 0, size = 2) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Observed heterozygosity (Hobs)") +
  xlab("Reproductive type")


# -----------------------------
# Combine Clone Groups 2 and 7
# -----------------------------

library(dplyr)

# Add clone group identifiers
Hobs_df_group2$CloneGroup <- "group2"
Hobs_df_group7$CloneGroup <- "group7"

# Merge datasets
Hobs_all <- bind_rows(Hobs_df_group2, Hobs_df_group7)

# Plot combined data
ggplot(Hobs_all, aes(x = Group, y = Hobs, color = CloneGroup)) +
  geom_jitter(width = 0.1, size = 4, alpha = 0.8) +
  scale_y_continuous(limits = c(0, 0.5)) +
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


# -----------------------------
# Statistical tests using GLM
# -----------------------------

# Clone Group 2
Hobs_group2_glm <- glm(Hobs ~ Group, data = Hobs_df_group2)

library(multcomp)
mult_comp_Hobs_group2_glm <- glht(Hobs_group2_glm,
                                 linfct = mcp(Group = "Tukey"))
summary(mult_comp_Hobs_group2_glm)
# Not significant

# Null model comparison (for confirmation)
Hobs_group2_glm_null <- glm(Hobs ~ 1, data = Hobs_df_group2)
anova(Hobs_group2_glm, Hobs_group2_glm_null)


# Clone Group 7
Hobs_group7_glm <- glm(Hobs ~ Group, data = Hobs_df_group7)

mult_comp_Hobs_group7_glm <- glht(Hobs_group7_glm,
                                 linfct = mcp(Group = "Tukey"))
summary(mult_comp_Hobs_group7_glm)
# Not significant

# Null model comparison (for confirmation)
Hobs_group7_glm_null <- glm(Hobs ~ 1, data = Hobs_df_group7)
anova(Hobs_group7_glm, Hobs_group7_glm_null)

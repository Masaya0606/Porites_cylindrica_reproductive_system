# Clone groups 2 and 7
# Files for this analysis "Clone_group_7_M4_MH2_F2_GT.txt" and "Clone_group_7_M4_MH2_F2_GT.txt" were generated using bcftools_LOH.sh.
# The calculations were performed using bcftools_LOH_vcf.sh
# Analyses were conducted on the NIG supercomputer.
# The input VCF files were filtered with MAF = 0.05, geno = 0.1,
# and HWE = 10^6. The total number of SNPs is approximately 6,800.

# -----------------------------
# Clone Group 7
# -----------------------------

df <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_7_M4_MH2_F2_GT.txt",
                 header = FALSE,
                 col.names = c("CHROM","POS","M1","M2","M3","M4","MH1","MH2","F1","F2"))

# Remove sites with missing genotypes in any individual
non_missing <- (df$M1 != "./.") & (df$M2 != "./.") & (df$M3 != "./.") &
               (df$M4 != "./.") & (df$MH1 != "./.") & (df$MH2 != "./.") &
               (df$F1 != "./.") & (df$F2 != "./.")

# Sites heterozygous in at least one M individual
M_het <- ((df$M1 == "0/1") | (df$M2 == "0/1") |
          (df$M3 == "0/1") | (df$M4 == "0/1")) & non_missing

# Sites heterozygous in at least one MH individual
MH_het <- ((df$MH1 == "0/1") | (df$MH2 == "0/1")) & non_missing

# Sites homozygous in all F individuals
F_hom <- (df$F1 %in% c("0/0","1/1")) &
         (df$F2 %in% c("0/0","1/1")) & non_missing

# Sites homozygous in all M individuals
M_hom <- (df$M1 %in% c("0/0","1/1")) &
         (df$M2 %in% c("0/0","1/1")) &
         (df$M3 %in% c("0/0","1/1")) &
         (df$M4 %in% c("0/0","1/1")) & non_missing

# LOH sites
LOH_site_F <- MH_het & F_hom
LOH_site_M <- MH_het & M_hom

# Summary statistics
N_total   <- sum(non_missing)   # Total number of comparable SNPs
N_MH_het  <- sum(MH_het)        # Number of heterozygous sites in MH (denominator)
N_M_het   <- sum(M_het)

N_LOH_MH_F <- sum(LOH_site_F)   # LOH sites from MH to F
N_LOH_MH_M <- sum(LOH_site_M)   # LOH sites from MH to M

LOH_rate_MH_F <- N_LOH_MH_F / N_MH_het
LOH_rate_MH_M <- N_LOH_MH_M / N_MH_het

# -----------------------------
# Pairwise LOH calculation
# -----------------------------

# Function to calculate LOH rate between a parent and an offspring
calc_loh_pair <- function(parent_gt, offspring_gt) {
  non_missing <- (parent_gt != "./.") & (offspring_gt != "./.")
  parent_het  <- (parent_gt == "0/1") & non_missing
  offspring_hom <- (offspring_gt %in% c("0/0","1/1")) & non_missing

  denom <- sum(parent_het)
  num   <- sum(parent_het & offspring_hom)

  c(LOH_rate = ifelse(denom > 0, num / denom, NA),
    N_het_parent = denom,
    N_LOH = num)
}

# -----------------------------
# Recalculate and summarize Clone Group 7
# -----------------------------

df7 <- read.table("~/NGS_analyses/Porites_LOH/Clone_group_7_M4_MH2_F2_GT.txt",
                  header = FALSE,
                  col.names = c("CHROM","POS","M1","M2","M3","M4","MH1","MH2","F1","F2"))

# Function to format LOH results as a data frame
make_loh_row <- function(parent_name, offspring_name,
                         parent_gt, offspring_gt,
                         clone_group) {

  res <- calc_loh_pair(parent_gt, offspring_gt)

  OffspringType <- dplyr::case_when(
    startsWith(offspring_name, "F")  ~ "F",
    startsWith(offspring_name, "M")  ~ "M",
    startsWith(offspring_name, "MH") ~ "MH",
    TRUE ~ NA_character_
  )

  data.frame(
    CloneGroup    = clone_group,
    Parent        = parent_name,
    Offspring     = offspring_name,
    OffspringType = OffspringType,
    LOH_rate      = as.numeric(res["LOH_rate"]),
    N_het_parent  = as.numeric(res["N_het_parent"]),
    N_LOH         = as.numeric(res["N_LOH"])
  )
}

library(dplyr)

# Compile LOH results for Clone Group 7
loh_G7 <- rbind(
  # MH → F
  make_loh_row("MH1", "F1", df7$MH1, df7$F1, "Group7"),
  make_loh_row("MH2", "F1", df7$MH2, df7$F1, "Group7"),
  make_loh_row("MH1", "F2", df7$MH1, df7$F2, "Group7"),
  make_loh_row("MH2", "F2", df7$MH2, df7$F2, "Group7"),

  # MH → M
  make_loh_row("MH1", "M1", df7$MH1, df7$M1, "Group7"),
  make_loh_row("MH1", "M2", df7$MH1, df7$M2, "Group7"),
  make_loh_row("MH1", "M3", df7$MH1, df7$M3, "Group7"),
  make_loh_row("MH1", "M4", df7$MH1, df7$M4, "Group7"),
  make_loh_row("MH2", "M1", df7$MH2, df7$M1, "Group7"),
  make_loh_row("MH2", "M2", df7$MH2, df7$M2, "Group7"),
  make_loh_row("MH2", "M3", df7$MH2, df7$M3, "Group7"),
  make_loh_row("MH2", "M4", df7$MH2, df7$M4, "Group7")
)

# -----------------------------
# Combine Clone Groups and plot
# -----------------------------

loh_all <- rbind(loh_G2, loh_G7)

ggplot(loh_all, aes(x = OffspringType, y = LOH_rate, color = CloneGroup)) +
  geom_jitter(width = 0.1, size = 4, alpha = 0.9) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "black", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(
    x = "Offspring type",
    y = "LOH rate (MH → offspring)",
    color = "Clone group"
  )

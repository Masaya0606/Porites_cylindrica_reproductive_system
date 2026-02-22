# Load data
library(readxl)
SexChange <- read_excel(
  "~/path/Sex_of_Pcyl_update2601.xlsx",
  sheet = "Sex_of_P.cylindrica (2)"
)

library(dplyr)

# Filter out rows lacking hermaphroditism rate data (2024)
SexChange_filter <- SexChange %>%
  filter(!is.na(`2024_hermaphroditism_rate`))

# Define sex transition categories between 2024 and July 2025 histology
SexChange_filter <- SexChange_filter %>%
  mutate(
    transition = case_when(
      `2024` == "M"  & `2025_July_histo` == "MH" ~ "M → MH",
      `2024` == "MH" & `2025_July_histo` == "M"  ~ "MH → M",
      `2024` == "M"  & `2025_July_histo` == "M"  ~ "M → M",
      `2024` == "MH" & `2025_July_histo` == "MH" ~ "MH → MH",
      TRUE ~ NA_character_
    )
  )

table(SexChange_filter$transition)

library(tidyr)

# Convert hermaphroditism rates into long format (2024 vs 2025)
SexChange_long <- SexChange_filter %>%
  pivot_longer(
    cols = c(`2024_hermaphroditism_rate`, `2025_hermaphroditism_rate`),
    names_to = "Year",
    values_to = "Herm_rate"
  ) %>%
  mutate(Year = ifelse(Year == "2024_hermaphroditism_rate", "2024", "2025"))

# (Optional) Exclude FH if needed (note: current filter expression is likely incorrect)
SexChange_long_exFH <- SexChange_long %>%
  filter(`2024` == !is.na("FH"))

library(ggplot2)

# Plot changes in hermaphroditism rate across years for each colony,
# colored by the 2024 sex category
ggplot(
  SexChange_long,
  aes(
    x = Year, y = Herm_rate,
    group = Colony_ID,   # individual ID column
    color = `2024`       # color by sex category in 2024
  )
) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  theme_classic()

# Create a long-format dataset with sex category mapped to the corresponding year
df_long <- SexChange_filter %>%
  transmute(
    Colony_ID,
    Herm_2024 = `2024_hermaphroditism_rate`,
    Herm_2025 = `2025_hermaphroditism_rate`,
    Sex_2024  = `2024`,
    Sex_2025  = `2025_July_histo`
  ) %>%
  pivot_longer(
    cols = c(Herm_2024, Herm_2025),
    names_to = "Year",
    values_to = "Herm_rate"
  ) %>%
  mutate(
    Year = ifelse(Year == "Herm_2024", "2024", "2025"),
    Sex  = ifelse(Year == "2024", Sex_2024, Sex_2025)
  )

# Exclude FH category
df_long_exFH <- df_long %>%
  filter(Sex != "FH")

# Plot (including FH)
ggplot(df_long, aes(x = Year, y = Herm_rate, group = Colony_ID)) +
  geom_line(alpha = 0.35) +
  geom_point(aes(color = Sex), size = 2.5) +
  theme_classic() +
  labs(x = NULL, y = "Hermaphroditism rate")

# Plot (excluding FH)
ggplot(df_long_exFH, aes(x = Year, y = Herm_rate, group = Colony_ID)) +
  geom_line(alpha = 0.35) +
  geom_point(aes(color = Sex), size = 2.5) +
  theme_classic() +
  labs(x = NULL, y = "Hermaphroditism rate")

# Calculate the number of male polyps (total minus hermaphrodite polyps)
SexChange_long$`2024_number_of_male_polyp` <- c(
  SexChange_long$`2024_total` - SexChange_long$`2024_number_of_hermaphrodite_polyp`
)
SexChange_long$`2025_number_of_male_polyp` <- c(
  SexChange_long$`2025_total` - SexChange_long$`2025_number_of_hermaphrodite_polyp`
)

SexChange_long <- SexChange_long %>%
  mutate(
    `2024_number_of_male_polyp` = `2024_total` - `2024_number_of_hermaphrodite_polyp`,
    `2025_number_of_male_polyp` = `2025_total` - `2025_number_of_hermaphrodite_polyp`
  )

# Define transition groups and whether sex changed between years
SexChange_long <- SexChange_long %>%
  mutate(
    transition_group = case_when(
      `2024` == "MH" & `2025_July_histo` == "M"  ~ "MH_to_M",
      `2024` == "M"  & `2025_July_histo` == "MH" ~ "M_to_MH",
      `2024` == "MH" & `2025_July_histo` == "MH" ~ "MH_to_MH",
      `2024` == "M"  & `2025_July_histo` == "M"  ~ "M_to_M",
      TRUE ~ NA_character_
    ),
    change_binary = ifelse(`2024` == `2025_July_histo`,
                           "No_change", "Sex_change")
  )

# Compute hermaphroditism rates for 2024 and 2025, and their changes
SexChange_long <- SexChange_long %>%
  mutate(
    herm_rate_2024 = `2024_number_of_hermaphrodite_polyp` / `2024_total`,
    herm_rate_2025 = `2025_number_of_hermaphrodite_polyp` / `2025_total`,
    delta_herm = herm_rate_2025 - herm_rate_2024,
    abs_delta_herm = abs(delta_herm)
  )

library(ggplot2)

# Remove rows with missing transition labels
SexChange_long_filter <- SexChange_long %>%
  filter(!is.na(transition))

# Visualize changes in hermaphroditism rate by transition group
ggplot(SexChange_long_filter,
       aes(x = transition_group, y = abs_delta_herm, group = Colony_ID)) +
  geom_line(alpha = 0.35) +
  geom_point(aes(color = transition_group), size = 2.5) +
  theme_classic() +
  labs(x = NULL, y = "Changes in hermaphroditism rate")

library(glmmTMB)

# Fit a Tweedie GLMM: absolute change in hermaphroditism rate ~ transition group,
# with Genet as a random intercept
glmm_tweedie <- glmmTMB(
  abs_delta_herm ~ transition_group + (1 | Genet),
  data = SexChange_long_filter,
  family = tweedie(link = "log")
)

# Null model (random effect only)
glmm_tweedie_null <- glmmTMB(
  abs_delta_herm ~ 1 + (1 | Genet),
  data = SexChange_long_filter,
  family = tweedie(link = "log")
)

# Likelihood-ratio test and model summary
anova(glmm_tweedie, glmm_tweedie_null)
summary(glmm_tweedie)

# Post hoc comparisons (Tukey-adjusted)
em_glmm_model <- emmeans(glmm_tweedie, ~ transition_group)
pairwise_glmm_model <- contrast(em_glmm_model, method = "pairwise", adjust = "tukey")
summary(pairwise_glmm_model)

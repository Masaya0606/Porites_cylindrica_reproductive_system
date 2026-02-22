# -----------------------------
# Data import
# -----------------------------
library(readxl)

clone_sex_diameter <- read_excel(
  "~/path/ColonyID_Sex_Color_GPS_Density_Clones.xlsx"
)

clone_sex_diameter_2412 <- read_excel(
  "~/path/PoritesInohaRowData2412_colony_size.xlsx",
  sheet = "Sheet1"
)

# Load required libraries
library(ggplot2)
library(dplyr)

# -----------------------------
# Boxplot of colony size by sex (dataset 1)
# -----------------------------
ggplot(
  clone_sex_diameter,
  aes(
    x = factor(sex, levels = c("M", "F", "MH", "FH")),
    y = MeanDiameter,
    fill = sex
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Size", y = "Mean diameter (cm)") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

# -----------------------------
# Boxplot of colony size by sex (dataset 2: 2412)
# -----------------------------
ggplot(
  clone_sex_diameter_2412,
  aes(
    x = factor(Sex, levels = c("M", "F", "MH", "FH")),
    y = Mean_radius,
    fill = Sex
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Size", y = "Mean diameter (cm)") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

library(lmerTest)
library(dplyr)

# -----------------------------
# Exclude FH individuals
# -----------------------------
filtered_clone_sex_diameter <- clone_sex_diameter %>%
  filter(sex != "FH")

filtered_clone_sex_diameter_2412 <- clone_sex_diameter_2412 %>%
  filter(Sex != "FH")

# -----------------------------
# Linear mixed-effects models
# -----------------------------
# Refit the model (structure unchanged)
Sex_diameter_lmer <- lmer(
  MeanDiameter ~ sex + (1 | color),
  data = filtered_clone_sex_diameter
)

# Residual diagnostics
res <- resid(Sex_diameter_lmer)
qqnorm(res)
qqline(res)
hist(res)

Sex_diameter_lmer_2412 <- lmer(
  Mean_radius ~ Sex + (1 | Color),
  data = filtered_clone_sex_diameter_2412
)

Sex_diameter_lmer_2412_Null <- lmer(
  Mean_radius ~ 1 + (1 | Color),
  data = filtered_clone_sex_diameter_2412
)

anova(Sex_diameter_lmer_2412, Sex_diameter_lmer_2412_Null)

res <- resid(Sex_diameter_lmer_2412)
qqnorm(res)
qqline(res)
hist(res)

# Inspect fixed effects
summary(Sex_diameter_lmer)
summary(Sex_diameter_lmer_2412)

# -----------------------------
# Multiple comparisons (Tukey)
# -----------------------------
em_sex_size <- emmeans(Sex_diameter_lmer, ~ sex)
pairwise_comparisons_sex_size <-
  contrast(em_sex_size, method = "pairwise", adjust = "tukey")
summary(pairwise_comparisons_sex_size)

em_sex_size_2412 <- emmeans(Sex_diameter_lmer_2412, ~ Sex)
pairwise_comparisons_sex_size_2412 <-
  contrast(em_sex_size_2412, method = "pairwise", adjust = "tukey")
summary(pairwise_comparisons_sex_size_2412)

# -----------------------------
# Summarize data at the clone level
# -----------------------------
# Aggregate by CloneID to obtain a single mean diameter
# and the dominant sex (most frequent category)
clone_summary <- clone_sex_diameter %>%
  filter(!is.na(CloneID) & !is.na(sex)) %>%
  group_by(CloneID) %>%
  summarise(
    MeanDiameter = mean(MeanDiameter, na.rm = TRUE),
    MajorSex = names(sort(table(sex), decreasing = TRUE))[1]
  )

# Plot histogram of clone sizes by dominant sex
ggplot(clone_summary, aes(x = MeanDiameter, fill = MajorSex)) +
  geom_histogram(binwidth = 5, color = "black", alpha = 0.8) +
  labs(
    title = "Distribution of Clone Sizes by Sex",
    x = "Mean Diameter",
    y = "Number of Clones",
    fill = "Sex"
  ) +
  theme_minimal()

# -----------------------------
# Histogram of colony sizes per CloneID
# -----------------------------
filtered_data <- clone_sex_diameter %>%
  filter(!is.na(CloneID), !is.na(sex), !is.na(MeanDiameter))

ggplot(filtered_data, aes(x = MeanDiameter, fill = sex)) +
  geom_histogram(binwidth = 5, color = "black", alpha = 0.8) +
  facet_wrap(~ CloneID) +
  labs(
    title = "Histogram of Mean Diameter by Sex per CloneID",
    x = "Mean Diameter",
    y = "Count",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

library(ggplot2)
library(dplyr)

# Use only rows without missing values
filtered_data <- clone_sex_diameter %>%
  filter(!is.na(Clone), !is.na(sex), !is.na(MeanDiameter))

# Faceted histogram by clone category, colored by sex
ggplot(filtered_data, aes(x = MeanDiameter, fill = sex)) +
  geom_histogram(
    binwidth = 5,
    color = "black",
    alpha = 0.8,
    position = "stack"
  ) +
  facet_wrap(~ Clone) +
  labs(
    title = "Distribution of Mean Diameter by Clone Category and Sex",
    x = "Mean Diameter",
    y = "Count",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# -----------------------------
# Spatial structure of clones (GPS)
# -----------------------------
# Load data (CSV file specified as needed)
Density_Clone <- read_excel(
  "~/path/ColonyID_Sex_Color_GPS_Density_Clones.xlsx"
)

# Remove rows with missing values
df_clean <- Density_Clone %>%
  filter(
    !is.na(CloneID),
    !is.na(sex),
    !is.na(MeanDiameter),
    !is.na(Latitude),
    !is.na(Longitude)
  )

# Plot clone structure by GPS location
ggplot(df_clean, aes(x = Longitude, y = Latitude)) +
  geom_path(
    aes(group = CloneID, color = CloneID),
    linewidth = 0.6,
    alpha = 0.5
  ) +
  geom_point(
    aes(color = sex, size = MeanDiameter),
    alpha = 0.8
  ) +
  scale_size(range = c(2, 10)) +
  labs(
    title = "Clone Structure by GPS Location",
    x = "Longitude",
    y = "Latitude",
    size = "Mean Diameter",
    color = "Sex / CloneID"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

# -----------------------------
# Spatial plots faceted by CloneID
# -----------------------------
ggplot(df_clean, aes(x = Longitude, y = Latitude)) +
  geom_path(
    aes(group = CloneID),
    linewidth = 0.6,
    alpha = 0.5,
    color = "gray50"
  ) +
  geom_point(
    aes(color = sex, size = MeanDiameter),
    alpha = 0.9
  ) +
  facet_wrap(~ CloneID) +
  scale_size(range = c(2, 10)) +
  labs(
    title = "Clone Structure per CloneID",
    x = "Longitude",
    y = "Latitude",
    size = "Mean Diameter",
    color = "Sex"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold")
  )

# -----------------------------
# Spatial plots excluding non_Clone
# -----------------------------
df_clean <- Density_Clone %>%
  filter(
    !is.na(CloneID),
    CloneID != "non_Clone",
    !is.na(sex),
    !is.na(MeanDiameter),
    !is.na(Latitude),
    !is.na(Longitude)
  )

ggplot(df_clean, aes(x = Longitude, y = Latitude)) +
  geom_path(
    aes(group = CloneID),
    linewidth = 0.6,
    alpha = 0.5,
    color = "gray50"
  ) +
  geom_point(
    aes(color = sex, size = MeanDiameter),
    alpha = 0.9
  ) +
  facet_wrap(~ CloneID) +
  scale_size(range = c(2, 10)) +
  labs(
    title = "Clone Structures (excluding non_Clone)",
    x = "Longitude",
    y = "Latitude",
    size = "Mean Diameter",
    color = "Sex"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold")
  )

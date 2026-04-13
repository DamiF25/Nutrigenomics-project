# if required, install packages
install.packages("dplyr")
install.packages("survival")
install.packages("survminer")

# Load required packages
library(dplyr)
library(survival)
library(survminer)

# Set working directory

# Load data
raw_data <- read.csv("Pooled snail mortality.csv")

# Transform to long format: one row per death event
long_data <- raw_data %>%
  pivot_longer(cols = starts_with("Week"), names_to = "Week", values_to = "Deaths") %>%
  mutate(Week = as.numeric(str_remove(Week, "Week"))) %>%
  uncount(weights = Deaths) %>%
  mutate(status = 1)  # 1 = death event

# Add censored observations (snails that survived all 8 weeks)
censored_data <- raw_data %>%
  mutate(Survivors = Total_snail - Total_death) %>%
  select(Diet, Survivors) %>%
  uncount(weights = Survivors) %>%
  mutate(Week = 8, status = 0)  # 0 = censored

# Combine death and censored data
km_data <- bind_rows(long_data, censored_data)

# Fit Kaplan-Meier model
km_fit <- survfit(Surv(Week, status) ~ Diet, data = km_data)
names(km_fit$strata) <- gsub("Diet=", "", names(km_fit$strata))

# Perform pairwise log-rank test with BH adjustment
pairwise_survdiff(Surv(Week, status) ~ Diet, data = km_data)

# Calculate % survival
survival_percent <- km_data %>%
  group_by(Diet) %>%
  summarise(
    Total = n(),
    Survived = sum(status == 0),
    Percent_Survival = round((Survived / Total) * 100, 1)
  ) %>%
  arrange(desc(Percent_Survival))
View(survival_percent)

# Plot survival curves

# Manually assign colour palettes to diets 
diet_colors <- c(
  "DL" = "#1b9e77",
  "FF" = "#FF2A00", 
  "FL" = "#7385BA",
  "SG" = "#6A00FF",
  "WL" = "#A6D690",
  "WL+FF" = "#e6ab02",
  "WL+SG" = "#874F39"
)

# Build the plot first
plot <- ggsurvplot(
  km_fit,
  data = km_data,
  conf.int = FALSE,
  pval = TRUE,
  risk.table = FALSE,
  xlab = "Week",
  ylab = "Survival probability",
  legend.title = "Diet",
  legend = "right"
)

# Modify the ggplot object
plot$plot <- plot$plot +
  scale_color_manual(values = diet_colors) +
  scale_fill_manual(values = diet_colors) +
  theme_classic(base_size = 14) + 
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 1.0))

# Increase p-value font size
plot$pval <- TRUE
plot$pval.coord <- c(1, 0.1)  # Optional: reposition p-value
plot$pval.size <- 20           # Increase font size here

# Display the plot
print(plot)


ggsave("snail_survival_plot.png", plot = plot$plot, width = 8, height = 5.5, dpi = 300)




# Install necessary packages
install.packages("ggplot2")
install.packages("ggpubr") # for comparison significance on graph
install.packages("dplyr") 
install.packages("car") #For Levene's test
install.packages("multcompView") #To extract Tukey's test result
install.packages("ggstatsplot")

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(multcompView)
library(ggstatsplot)

# Set working directory

# Load data
data <- read.csv("Pooled R data on embryo per egg mass.csv")

# Shapiro-Wilk test for normality of EM size per diet
shapiro_test <- shapiro.test(data$Embryo_per_EM)

# Check results (if p<0.05, data is NOT normally distributed)
print(shapiro_test)

# Perform a Levene's test (if p<0.05, variances are NOT equal)
leveneTest(Embryo_per_EM ~ Diet, data = data)

# If data are normally distributed, and variances are equal, we can use ANOVA
anova_result <- aov(Embryo_per_EM ~ Diet, data = data)

anova_result

# Perform Tukey's post hoc test (all groups vs each other)
tukey_result <- TukeyHSD(anova_result, conf.level=.95)

# Print the results
print(tukey_result)

# Plot Tukey's results
plot(tukey_result, las = 1, col = "blue")

# Convert p-values to stars
p_to_stars <- function(p) {
  ifelse(p <= 0.001, "***",
         ifelse(p <= 0.01, "**",
                ifelse(p <= 0.05, "*", "")))
}

# Convert Tukey's results to a data frame
tukey_df <- as.data.frame(tukey_result$Diet)

# Add group1 and group2 columns by splitting the row names
tukey_annotations <- data.frame(
  group1 = sapply(strsplit(rownames(tukey_df), "-"), `[`, 1),
  group2 = sapply(strsplit(rownames(tukey_df), "-"), `[`, 2),
  p.adj = tukey_df$`p adj`
) %>%
  mutate(label = p_to_stars(p.adj)) %>%  # Convert to stars
  filter(label != "")                    # Remove "ns"

print(tukey_annotations)

# Calculate mean and standard error for each Diet group (FOR ERROR BAR PLOTTING)
data_summary <- data %>%
  group_by(Diet) %>%
  summarise(
    mean_embryo_per_egg_mass = mean(Embryo_per_EM),
    se_embryo_per_egg_mass = sd(Embryo_per_EM) / sqrt(n())  # Standard error
  )
View(data_summary)

# Define colors
diet_colors <- c("DL" = "#54BAA4", "FF" = "#DB6A24", "SG" = "#D97BE0",
                 "WL+FF" = "#FFD000", "WL+SG" = "#D4BA7F")

# Bar plot with error bars
ggplot(data_summary, aes(x = Diet, y = mean_embryo_per_egg_mass)) +
  geom_bar(stat = "identity", fill = diet_colors, linewidth = 1.5, width = 0.6) +
  geom_errorbar(aes(ymin = mean_embryo_per_egg_mass - se_embryo_per_egg_mass, 
                    ymax = mean_embryo_per_egg_mass + se_embryo_per_egg_mass), 
                width = 0.3, size = 0.5, color = "black") +  # Error bars
  stat_pvalue_manual(
    data = tukey_annotations,
    label = "label", size = 10,
    y.position = seq(max(data$Embryo_per_EM) + 0.2, by = 0.2, length.out = nrow(tukey_annotations)),
    step.increase = 0.15
  ) +
  labs(x = "Diet", y = "Average embryo per egg mass") +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme(legend.position = "None") +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.line = element_line(linewidth = 1.0))         # Thicker axis lines)

# Save the final figure
ggsave("Average embryo per egg mass.png", width = 6, height = 5.5, dpi = 300)

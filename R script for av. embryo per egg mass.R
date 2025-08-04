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

# Convert Tukey's results to a data frame
tukey_df <- as.data.frame(tukey_result$Diet)

# Add group1 and group2 columns by splitting the row names
tukey_annotations <- data.frame(
  group1 = sapply(strsplit(rownames(tukey_df), "-"), `[`, 1),
  group2 = sapply(strsplit(rownames(tukey_df), "-"), `[`, 2),
  p.adj = tukey_df$`p adj`
)

# Filter significant comparisons
significant_comparisons <- tukey_annotations %>%
  filter(p.adj < 0.05)

print(significant_comparisons)

# Round p-values to 3 decimal places
significant_comparisons$p.adj <- round(significant_comparisons$p.adj, 4)

# Calculate mean and standard error for each Diet group (FOR ERROR BAR PLOTTING)
data_summary <- data %>%
  group_by(Diet) %>%
  summarise(
    mean_embryo_per_egg_mass = mean(Embryo_per_EM),
    se_embryo_per_egg_mass = sd(Embryo_per_EM) / sqrt(n())  # Standard error
  )
View(data_summary)

# Bar plot with error bars
ggplot(data_summary, aes(x = Diet, y = mean_embryo_per_egg_mass)) +
  geom_bar(stat = "identity", fill = NA, color = "#F58427", linewidth = 1.5, width = 0.6) +
  geom_errorbar(aes(ymin = mean_embryo_per_egg_mass - se_embryo_per_egg_mass, 
                    ymax = mean_embryo_per_egg_mass + se_embryo_per_egg_mass), 
                width = 0.3, size = 0.5, color = "#F58427") +  # Error bars
  stat_pvalue_manual(
    data = significant_comparisons,
    label = "p.adj", size = 5.5,
    y.position = seq(max(data$Embryo_per_EM) + 0.2, by = 0.2, length.out = nrow(significant_comparisons)),
    step.increase = 0.15
  ) +
  labs(x = "Diet", y = "Average embryo per egg mass") +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme(legend.position = "None") +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.line = element_line(size = 1.0))         # Thicker axis lines)

# Save the final figure
ggsave("Average embryo per egg mass.tif", width = 6, height = 5.5, dpi = 300)

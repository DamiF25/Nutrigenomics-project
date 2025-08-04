# Install necessary packages
install.packages("ggplot2")
install.packages("ggpubr") # for comparison significance on graph
install.packages("dplyr") 
install.packages("car")
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

# Load the data
data <- read.csv("Pooled R data on egg mass output.csv")

# Shapiro-Wilk test for normality of eggs per snail
shapiro_test <- shapiro.test(data$Eggs)

# Check results (if p<0.05, data is NOT normally distributed)
print(shapiro_test)

# Perform a Levene's test (if p<0.05, variances are NOT equal)
leveneTest(Eggs ~ Diet, data = data)

# If data are normally distributed, and variances are equal, we can use ANOVA
anova_result <- aov(Eggs ~ Diet, data = data)
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
significant_comparisons$p.adj <- round(significant_comparisons$p.adj, 3)

# Calculate mean and standard error for each Diet group
data_summary <- data %>%
  group_by(Diet) %>%
  summarise(
    mean_eggs = mean(Eggs),
    se_eggs = sd(Eggs) / sqrt(n())  # Standard error
  )
View(data_summary)

# Box plot with significant comparisons
boxplot <- ggplot(data, aes(x = Diet, y = Eggs)) +
  geom_boxplot(aes(color = Diet), outlier.shape = NA, size = 1.5, alpha = 1) +  # Box outline
  geom_jitter(aes(color = Diet), position = position_jitter(width = 0.2), size = 1.5) +  # Add jitter points
    stat_summary(fun = mean, geom = "point", shape = 18, color = "black", size = 3) +
  stat_pvalue_manual(
    data = significant_comparisons,
    label = "p.adj", size = 5.5,
    y.position = seq(max(data$Eggs) + 8, by = 1, length.out = nrow(significant_comparisons)),
    step.increase = 0.08
  ) +
  labs(x = "Diet", y = "Average egg mass count") +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size =20))

# Print the plot
print(boxplot)

# Save the final figure
ggsave("Average egg mass count.png", width = 7, height = 6.5, dpi = 300)


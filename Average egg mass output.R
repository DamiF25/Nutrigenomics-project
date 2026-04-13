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

# Perform ANOVA
anova_result <- aov(Eggs ~ Diet, data = data)

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

# Add group1 and group2 columns by splitting the row names (also add stars)
tukey_annotations <- data.frame(
  group1 = sapply(strsplit(rownames(tukey_df), "-"), `[`, 1),
  group2 = sapply(strsplit(rownames(tukey_df), "-"), `[`, 2),
  p.adj = tukey_df$`p adj`
) %>%
  mutate(label = p_to_stars(p.adj)) %>%  # Convert to stars
  filter(label != "")                    # Remove "ns"

print(tukey_annotations)

# Calculate mean and standard error for each Diet group
data_summary <- data %>%
  group_by(Diet) %>%
  summarise(
    mean_eggs = mean(Eggs),
    se_eggs = sd(Eggs) / sqrt(n())  # Standard error
  )
View(data_summary)

# Define colors
diet_colors <- c("DL" = "#54BAA4", "FF" = "#DB6A24", "FL" = "#7385BA", "SG" = "#D97BE0",
                 "WL" = "#A6D690", "WL+FF" = "#FFD000", "WL+SG" = "#D4BA7F")

# Box plot with significant comparisons
boxplot <- ggplot(data, aes(x = Diet, y = Eggs)) +
  geom_boxplot(aes(color = Diet, fill = Diet), outlier.shape = NA, size = 1.5, alpha = 0.5) + 
  geom_jitter(aes(color = Diet), position = position_jitter(width = 0.2), size = 1.5) +  # Add jitter points
    stat_summary(fun = mean, geom = "point", shape = 18, color = "black", size = 3) +
  stat_pvalue_manual(
    data = tukey_annotations,
    label = "label", size = 10,
    y.position = seq(max(data$Eggs) + 8, by = 1, length.out = nrow(tukey_annotations)),
    step.increase = 0.08
  ) +
  scale_color_manual(values = diet_colors) +   # manual colors for outlines and jitter
  scale_fill_manual(values = diet_colors) +    # manual fill colors for boxplots
  labs(x = "Diet", y = "Average egg mass count") +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size =20),
        panel.border = element_rect(size = 2.0, fill = NA))

# Print the plot
print(boxplot)

# Save the final figure
ggsave("Average egg mass count_2.png", width = 7, height = 6.5, dpi = 300)

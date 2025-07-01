# Load necessary libraries
install.packages("ggplot2")
install.packages("dplyr") #For %>% function
install.packages("tidyr") #For the 'separate' function
install.packages("FSA") #For Dunn's test 
install.packages("car") #For Levene's test
install.packages("rstatix") #For stat_pvalue_manual
install.packages("ggpubr") #For stat_pvalue_manual (sometimes loaded via rstatix)

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(FSA) 
library(car)
library(rstatix)     
library(ggpubr)     


# Load and filter your data
snail_data <- read.csv("Pooled R data on final snail growth.csv")
week_8_data <- snail_data %>% filter(Week == 8)
week_8_data$Diet <- factor(week_8_data$Diet)  # Ensure Diet is a factor

# Check levels
str(week_8_data$Diet)
levels(week_8_data$Diet)

# Perform Shapiro_Wilk normality test for each Diet at Week 8
normality_test_results <- week_8_data %>%
  group_by(Diet) %>%
  summarize(p_value = shapiro.test(Size)$p.value)

# View the results
print(normality_test_results)

# Perform Shapiro-Wilk normality test on the entire Week 8 dataset
shapiro_test_result <- shapiro.test(week_8_data$Size)

# Display the result
shapiro_test_result

# Use Levene's test to check if the variance of Size is equal across the different Diet groups
# Perform Levene's Test to check homogeneity of variances, grouping by Diet ONLY
levene_test_result <- leveneTest(Size ~ Diet, data = week_8_data)

# View the result
print(levene_test_result)

#NON NORMALLY DISTRIBUTED DATA - MULTIPLE COMPARISONS - KRUSKAL-WALLIS
kruskal.test(Size ~ Diet, data = week_8_data)

# Save test result in a variable
kruskal_result <- kruskal.test(Size ~ Diet, data = week_8_data)
print(kruskal_result)

# Perform Dunn's test with Holm correction
dunn_results <- dunnTest(Size ~ Diet, data = week_8_data, method = "holm")$res

# Filter only significant comparisons and assign increasing y positions manually
significant_comparisons <- dunn_results %>%
  filter(P.adj <= 0.05) %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  mutate(
    p.adj = P.adj,
    label = paste0("p = ", signif(p.adj, 2))
  )

# Convert group names to numeric positions based on factor levels
group_levels <- levels(week_8_data$Diet)

# Sort comparisons to apply ascending y positions with spacing
significant_comparisons <- significant_comparisons %>%
  mutate(
    x1 = as.numeric(factor(group1, levels = group_levels)),
    x2 = as.numeric(factor(group2, levels = group_levels)),
    x_mid = (x1 + x2) / 2  # midpoint used for ordering
  ) %>%
  arrange(x_mid) %>%  # order left to right
  mutate(
    y.position = seq(
      from = max(week_8_data$Size, na.rm = TRUE) + 1,
      by = 1.7, # increase vertical space between significance comparison bars
      length.out = n()))

# Plot
ggplot(week_8_data, aes(x = Diet, y = Size, fill = Diet)) +
  geom_violin(trim = FALSE, alpha = 1, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 21, fill = "red", color = "black", size = 3) +
  # Add significance bars and labels
  geom_segment(data = significant_comparisons,
               aes(x = x1, xend = x2, y = y.position + 1.5, yend = y.position + 1.5),
               inherit.aes = FALSE) +
  geom_text(data = significant_comparisons,
            aes(x = (x1 + x2) / 2, y = y.position + 2.1, label = label),
            size = 6.0,
            inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Diet", y = "Shell diameter (mm)") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 22),
    legend.position = "none")

# Save the final figure
ggsave("Final growth_Dunn_sig.png", width = 9, height = 7, dpi = 300)

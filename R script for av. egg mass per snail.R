# Install necessary packages
install.packages("dplyr") 
install.packages("ggplot2")
install.packages("ggpubr") # for comparison significance on graph
install.packages("FSA") #For Dunn's test
install.packages("rstatix") #For stat_pvalue_manual

install.packages("multcompView") #To extract Tukey's test result
install.packages("ggstatsplot")

# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)    # For stat_pvalue_manual or stat_compare_means
library(FSA)        # For dunnTest
library(rstatix)    # For p-value adjustment (optional)

# Set working directory

# Load data
data <- read.csv("Pooed R data on egg mass per snail.csv")

# To estimate number of egg-laying snails per Diet
data %>%
  filter(Egg_mass > 0) %>%
  group_by(Diet) %>%
  summarise(estimated_egg_laying_snails = sum(Snail), .groups = "drop")

# To get total egg masses per diet
data %>%
  group_by(Diet) %>%
  summarise(total_egg_masses = sum(Egg_mass), .groups = "drop")

# To combine estimated number of egg-laying snails per Diet and total egg masses per diet in one table
summary_table <- data %>%
  group_by(Diet) %>%
  summarise(
    total_snails = sum(Snail),
    total_egg_masses = sum(Egg_mass),
    estimated_egg_laying_snails = sum(Snail[Egg_mass > 0]),
    avg_egg_mass_per_snail = total_egg_masses / total_snails,
    avg_egg_mass_per_egg_layer = ifelse(estimated_egg_laying_snails > 0,
                                        total_egg_masses / estimated_egg_laying_snails,
                                        NA),
    .groups = "drop"
  )

View(summary_table)

write.csv(summary_table, "EM_per_snail_summary_table.csv", row.names = FALSE)

# To visualise data

# Calculate avg egg mass per egg-laying snail per replicate (row)
data2 <- data %>%
  mutate(
    avg_egg_mass_per_egg_layer = ifelse(Egg_mass > 0, Egg_mass / Snail, NA_real_)
  ) %>%
  filter(!is.na(avg_egg_mass_per_egg_layer))  # only egg-laying snails

# With non-parametric data assumption
# Kruskal-Wallis test
kruskal_res <- data2 %>% kruskal_test(avg_egg_mass_per_egg_layer ~ Diet)
print(kruskal_res)

  # Dunn test with adjusted p-values
dunn_res <- data2 %>% dunn_test(avg_egg_mass_per_egg_layer ~ Diet, p.adjust.method = "bonferroni")
print(dunn_res)

# Order Diet factor alphabetically in the main data
data2 <- data2 %>%
  mutate(Diet = factor(Diet, levels = sort(unique(Diet))))

# Prepare Dunn’s test pairwise result
pval_table <- dunn_res %>%
  filter(p.adj < 0.05) %>%
  mutate(
    p = round(p.adj, 3),
    y.position = c(9.0, 9.7, 10.7),  # Adjust heights if needed
    group1 = as.character(group1),
    group2 = as.character(group2)
  ) %>%
  select(group1, group2, p, y.position)  # Include 'p' not 'p.adj'

# Plot
ggplot(data2, aes(x = Diet, y = avg_egg_mass_per_egg_layer)) + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 3, color = "orange") +
  stat_summary(fun = mean, geom = "point", size = 5, color = "black") +
  stat_pvalue_manual(pval_table, label = "p", tip.length = 0.02, size = 5.5) +
  labs(y = "Average egg mass per egg-laying snail", x = "Diet") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 20),
    legend.position = "none"
  )

# Save the final figure
ggsave("Average egg mass per egg-laying snail.tiff", width = 6.5, height = 6.5, dpi = 300)
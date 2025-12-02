library(tidyverse)
library(dplyr)
library(ggplot2)
library(vegan)
library(viridis)
library(RColorBrewer)
library(patchwork)

library(showtext)
showtext_auto()  # enables custom fonts in ggplot


# Read tviridis# Read the data
data <- read.csv("working_data_avg1.csv", stringsAsFactors = FALSE)



# Filter for CTL and ND treatments only
data_filtered <- data %>%
  filter(NutDisp %in% c("CTL", "ND")) %>%
  filter(MeanValue > 0)  # Remove zero abundance ASVs



# Create a community matrix (samples x ASVs)
# Rows = samples, Columns = ASVs, Values = abundance
community_matrix <- data_filtered %>%
  select(Sample_base, ASV, MeanValue) %>%
  pivot_wider(names_from = ASV, values_from = MeanValue, values_fill = 0) %>%
  column_to_rownames("Sample_base")

#cat("Community matrix dimensions:", dim(community_matrix), "\n")

# Calculate diversity indices
diversity_results <- data.frame(
  Sample_base = rownames(community_matrix),
  Shannon = diversity(community_matrix, index = "shannon"),
  Simpson = diversity(community_matrix, index = "simpson"),
  stringsAsFactors = FALSE
)

# Add treatment information
sample_metadata <- unique(data_filtered[,c('Sample_base', 'NutDisp', 'Time')])

diversity_results <- diversity_results %>%
  left_join(sample_metadata, by = "Sample_base")

# Display summary statistics

summary_stats <- diversity_results %>%
  group_by(NutDisp) %>%
  summarise(
    n_samples = n(),
    Shannon_mean = round(mean(Shannon, na.rm = TRUE), 3),
    Shannon_sd = round(sd(Shannon, na.rm = TRUE), 3),
    Simpson_mean = round(mean(Simpson, na.rm = TRUE), 3),
    Simpson_sd = round(sd(Simpson, na.rm = TRUE), 3),
    .groups = 'drop'
  )
diversity_results$NutDisp <- factor(diversity_results$NutDisp,
                                    levels = c("CTL", "ND"),
                                    labels = c("Control", "Nutri+Disp"))


print(summary_stats)

# Create violin plots
# Shannon diversity
p_shannon <- ggplot(diversity_results, aes(x = NutDisp, y = Shannon, fill = NutDisp)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  scale_fill_brewer(palette="Purples")+
  #scale_fill_manual(values = c("CTL" = "#E69F00", "ND" = "#56B4E9")) +
  labs(
    #title = "Shannon Diversity",
    #subtitle = "Comparison between CTL and ND treatments",
    x = "Treatment",
    y = "Shannon Diversity Index",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    #plot.title = element_text(family= "Arial", size = 40, face = "bold", hjust = 0.5),
    #plot.subtitle = element_text(family= "Arial", size = 30, hjust = 0.5),
    axis.title = element_text(size = 45, face= "bold"),
    axis.text = element_text( size = 45),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank()
  )

# Simpson diversity
p_simpson <- ggplot(diversity_results, aes(x = NutDisp, y = Simpson, fill = NutDisp)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  
  scale_fill_brewer(palette="Purples")+
  labs(
    #title = "Simpson Diversity Index",
    #subtitle = "Comparison between CTL and ND treatments", 
    x = "Treatment",
    y = "Simpson Diversity Index",
    fill = "Treatment"
  ) +
  #coord_cartesian(ylim = c(0, max(alpha_div$Shannon) + 1)) +
  theme_minimal() +
  theme(
    #plot.title = element_text(family = "Arial", size = 40, face = "bold", hjust = 0.5),
    #plot.subtitle = element_text(family = "Arial", size = 30, hjust = 0.5),
    axis.title = element_text(size = 45, face= "bold"),
    axis.text = element_text(size = 45),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank()
  )

# Display the plots
print(p_shannon)
print(p_simpson)



ggsave("Shanon.png", plot = p_shannon, 
       width = 10, height = 5, dpi = 300)
ggsave("Simpson.png", plot = p_simpson, 
       width = 10, height = 5, dpi = 300)


# Combined plot showing both diversity indices
div_long <- diversity_results %>%
  select(Sample_base,NutDisp, Shannon, Simpson) %>%
  pivot_longer(cols = c(Shannon, Simpson), names_to = "Diversity_Index", values_to = "Value")

p_combined2 <- ggplot(div_long, aes(x = NutDisp, y = Value, fill = NutDisp)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_grid(Diversity_Index ~ NutDisp, scales = "free") +
  scale_fill_brewer(palette="Purples")+
  labs(
    x = "6 weeks incubations", 
    y = "Diversity Index Value"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size= 45, face= "bold"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 45, face = "bold"),
    legend.position = "none")
    #legend.text = element_text(family= "Arial", size= 14),
    #legend.title = element_text(family= "Arial", size=14)
    #plot.title = element_text(size = 14, face = "bold")
  

print(p_combined2)

ggsave("Combined.png", plot = p_combined2, 
       width = 10, height = 5, dpi = 300)



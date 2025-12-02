library(tidyverse)

# Read the raw concentration data
bdloss <- read.csv("E:/Thami-Uni/Oil Microcosm Data/BDLOSS.csv")
colnames(bdloss) <- gsub("^X([0-9])", "\\1", colnames(bdloss))

# Clean up - remove empty rows
bdloss <- bdloss %>% 
  filter(!is.na(Treatment_type) & Treatment_type != "")

# Clean up Treatment_type (remove extra spaces)
bdloss$Treatment_type <- trimws(bdloss$Treatment_type)

# Check treatment types and counts
print("Treatment types in data:")
print(table(bdloss$Treatment_type))

# Store compound names (excluding Treatment_type)
compound_names <- colnames(bdloss)[-1]

# Convert to Long Format

bdloss_long <- bdloss %>%
  pivot_longer(cols = -Treatment_type,
               names_to = "Compound",
               values_to = "Concentration")


# 3. Calculate Mean Concentrations by Treatment (Baselines)


mean_by_treatment <- bdloss_long %>%
  group_by(Treatment_type, Compound) %>%
  summarise(mean_conc = mean(Concentration, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = Treatment_type, values_from = mean_conc)

# Rename columns to remove spaces
colnames(mean_by_treatment) <- gsub(" ", "_", colnames(mean_by_treatment))

print("Column names after cleaning:")
print(colnames(mean_by_treatment))

# Get mean Raw Diesel and mean UV weathered as baselines
raw_means <- mean_by_treatment %>% select(Compound, Raw_Diesel)
uv_means <- mean_by_treatment %>% select(Compound, UV_weathered)
                                         
#Calculate Percentage Loss for Each Replicate


pct_loss_data <- bdloss_long %>%
  left_join(raw_means, by = "Compound") %>%
  left_join(uv_means, by = "Compound") %>%
  mutate(
   
    Pct_Loss = case_when(
      Treatment_type == "UV weathered" ~ ((Raw_Diesel - Concentration) / Raw_Diesel) * 100,
      Treatment_type == "CTRL" ~ ((UV_weathered - Concentration) / UV_weathered) * 100,
      Treatment_type == "ND" ~ ((UV_weathered - Concentration) / UV_weathered) * 100,
      TRUE ~ NA_real_
    ),
    # Create comparison labels
    Comparison = case_when(
      Treatment_type == "UV weathered" ~ "Raw to UV",
      Treatment_type == "CTRL" ~ "UV to Control",
      Treatment_type == "ND" ~ "UV to Nutri+Disp",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Comparison))  # Remove Raw Diesel rows


# Calculate Summary Statistics (Mean, SD, SE)


pct_loss_summary <- pct_loss_data %>%
  group_by(Compound, Comparison) %>%
  summarise(
    mean_loss = mean(Pct_Loss, na.rm = TRUE),
    sd_loss = sd(Pct_Loss, na.rm = TRUE),
    n = n(),
    se_loss = sd_loss / sqrt(n),
    .groups = "drop"
  )

# Set factor order for compounds (original order from data)
pct_loss_summary$Compound <- factor(pct_loss_summary$Compound, 
                                    levels = rev(compound_names))

# Set factor order for Comparison
pct_loss_summary$Comparison <- factor(pct_loss_summary$Comparison,
                                      levels = c("Raw to UV", 
                                                 "UV to Control", 
                                                 "UV to Nutri+Disp"))

# Filter out compounds with all zeros or NaN (no variation possible)
pct_loss_summary <- pct_loss_summary %>%
  filter(!is.nan(mean_loss) & !is.infinite(mean_loss))

# View the summary data
print("Summary statistics:")
print(pct_loss_summary, n = Inf)


#Create the Plot


bd <- ggplot(pct_loss_summary, aes(x = Compound,
                                   y = mean_loss,
                                   fill = Comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) + 
  geom_errorbar(aes(ymin = mean_loss - se_loss, 
                    ymax = mean_loss + se_loss),
                position = position_dodge(width = 0.8),
                width = 0.3,
                linewidth = 0.5,
                color = "black") +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Raw to UV"        = "#F5C98E",
      "UV to Control"    = "#D65B5A", 
      "UV to Nutri+Disp" = "#586085"
    )
  ) +
  labs(x = "PAH Compound",
       y = "Percentage Loss (%)",
       fill = "Treatment Comparison") +
  scale_y_continuous(expand = c(0, 0), limits = c(-20, 100)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y  = element_text(family = "Times New Roman", size = 12),
    axis.text.x  = element_text(family = "Times New Roman", size = 12),
    axis.title   = element_text(family = "Times New Roman", face = "bold", size = 15),
    legend.position = "right",
    legend.title = element_text(family = "Times New Roman", face = "bold", size = 15),  
    legend.text  = element_text(family = "Times New Roman", size = 12)
  )


bd



ggsave("percentage_loss_pah_with_error.png", plot = bd, 
       width = 12, height = 10, dpi = 300, bg = "white")





write.csv(pct_loss_summary, "pct_loss_summary.csv", row.names = FALSE)


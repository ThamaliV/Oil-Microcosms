library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(scales)
library(patchwork)

#loading the phase analysis df
phase_df <- read.csv( "E:/Thami-Uni/Oil Microcosm Data/phase.csv")

phase_df <- phase_df %>%
  filter(!is.na(Treatment.Group) & Treatment.Group != "") %>%
  mutate(
    Treatment.Group = str_trim(Treatment.Group),
    Phase = str_trim(Phase)
  )
# Get the original compound order from column names (excluding Treatment.Group and Phase)
original_compound_order <- setdiff(names(phase_df), c("Treatment.Group", "Phase"))


# Convert to long format
phase_long <- phase_df %>%
  pivot_longer(
    cols = -c(Treatment.Group, Phase),
    names_to = "Compound",
    values_to = "Concentration"
  )%>%
  # Replace NA with 0 instead of removing them
  mutate(Concentration = replace_na(Concentration, 0))



# DIAGNOSTIC: Check which compounds have data in both phases
diagnostic <- phase_long %>%
  group_by(Treatment.Group, Compound) %>%
  summarise(
    n_phases = n(),
    phases_present = paste(Phase[Concentration > 0], collapse = ", "),
    oil_conc = sum(Concentration[Phase == "Oil"], na.rm = TRUE),
    aqueous_conc = sum(Concentration[Phase == "Aqueous"], na.rm = TRUE),
    .groups = "drop"
  )

print("Diagnostic - Compounds by phase:")
print(diagnostic %>% filter(Compound %in% c("Acenaphthene", "Phenanthrene", "X2.4.Dimethylphenanthrene")))



# Calculate percentage for each compound within each treatment group
phase_percentage <- phase_long %>%
  group_by(Treatment.Group, Compound) %>%
  mutate(
    Total_Conc = sum(Concentration, na.rm = TRUE),
    Percentage = if_else(Total_Conc > 0, (Concentration / Total_Conc) * 100, 0)
  ) %>%
  ungroup() %>%
  filter(Total_Conc > 0)  # Only remove compounds with zero total in both phases


# Clean compound names for better display and preserve order
phase_percentage <- phase_percentage %>%
  mutate(
    Compound_Clean = str_replace_all(Compound, "\\.", " "),
    # Set factor levels to preserve original order
    Compound_Clean = factor(Compound_Clean, 
                            levels = str_replace_all(original_compound_order, "\\.", " "))
  
  )


unique(phase_percentage$Phase)
# Create the stacked bar plot


p_phase <- ggplot(phase_percentage, aes(x = Compound_Clean, y = Percentage, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  facet_wrap(~Treatment.Group, ncol = 2) +
  scale_fill_manual(
    values = c("Oil" = "#E3D26F", "Aqueous" = "#79ADDC"),  # Note: check for extra spaces in Phase names
    name = "Phase")+
  labs(
    x = "PAH Compound",
    y = "Concentration Percentage (%)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family= "Times new roman", size = 15),
    axis.text.y = element_text(family= "Times new roman", size = 15),
    axis.title = element_text(family= "Times new roman", size = 15, face = "bold"),
    strip.text = element_text(family= "Times new roman", size = 20, face = "bold"),
    legend.position = "right",
    legend.title = element_text(family= "Times new roman", size = 15, face = "bold"),
    legend.text = element_text(family= "Times new roman", size = 15)
  )

print(p_phase)

ggsave("phase.png", plot = p_phase, 
       width = 10, height = 8, dpi = 300, bg="white")





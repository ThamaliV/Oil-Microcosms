library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(RColorBrewer)

plotting_data <- read.csv("working_data_avg1.csv", stringsAsFactors = FALSE)


genus_rel <- plotting_data%>%
  group_by(Sample_base) %>%
  mutate(RelAbundance = MeanValue / sum(MeanValue)) %>%
  ungroup() %>%
 
  group_by(Genus) %>%
  mutate(mean_abundance = mean(RelAbundance)) %>%
  ungroup() %>%
  
  mutate(Genus = ifelse(mean_abundance < 0.001, "Other", Genus))

my_colors <- c(
 
 
  "Colwellia" = "#A6CEE3",
  "Pseudomonas" = "#FB9A99",
  "Sulfitobacter" = "#FDBF6F",
  "Roseivirga"   = "#00FFFF",
  "Pseudohongiella" = "#E3A1AC",
  "Pseudorhodobacter" = "#33A02C",
  "Yoonia" = "#1F78B4",
  "Maricaulis" = "#FF7F00",
  "Ulvibacter" = "#CAB2D6",
  "Pseudoalteromonas" =  "#1F7884", 
  "Sulfurimonas" = "#B2DF8A",
  "Jannaschia"  = "#B15928",
  "Croceibacter"  = "#FDB462",
  "Halomonas"    = "#FFFF99",
  "Methylotenera" = "#6A3D9A",
  "Aurantivirga" = "#CC3399",
  "Paraperlucidibaca" = "#00abff",
  "OM43 clade" = "#CC0000",
  "NS3a marine group" = "#00FF33",
  "Other" = "#000000"
  
 
)
set.seed(123)  
genus_colors <- sample(my_colors)



genera <- unique(genus_rel$Genus)
print(genera)



#This step is to summarize the data set, other wise you will get horizontal lines across the stacks randomly. 
genus_rel_summary <- genus_rel %>%
  group_by(Time, NutDisp, Genus) %>%
  summarise(RelAbundance = sum(RelAbundance), .groups= "drop")

# Factor ordering to put "Other" at the bottom
genus_rel_summary$Genus <- factor(
  genus_rel_summary$Genus,
  levels = c(
    setdiff(unique(genus_rel_summary$Genus), "Other"), "Other"
  )
)
  
genus_rel_summary$NutDisp <- factor(
  genus_rel_summary$NutDisp,
  levels = c("PRE", "CTL", "ND")
)

levels(genus_rel_summary$NutDisp)


facet_labels <- c(
  "CTL" = "Control",
  "ND"  = "Nutri + Disp",
  "PRE" = "Pre"
)


p2<- ggplot(genus_rel_summary, aes(x = Time, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", linewidth = 0.1, colour= "black", linewidth = 0.2, position = "stack" ) +
  scale_fill_manual(values= my_colors)+
  #scale_fill_manual(values = genus_color_mapping) +
  facet_grid(~ NutDisp, labeller = labeller(NutDisp = facet_labels), scales = "free_x", space = "free_x") +
  labs(x = "Time Period", y = "Relative Abundance", fill = "Genus") +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 25, face = "bold", family = "Arial"),
    #strip.text.y = element_text(size = 15, face = "bold", family = "Arial"),
    axis.title = element_text(family= "Arial", size = 25, face= "bold"),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 20, face = "bold", family = "Arial"),
    legend.text  = element_text(size = 20, family = "Arial"),
  
    legend.position = "bottom"
  )

print(p2)

ggsave("stacked_0.001.png", plot = p2, 
      width = 16, height = 8, dpi = 2000)

genus_level <- genus_rel %>%
  group_by(Genus) %>%summarise(Total_Abundance= sum(RelAbundance,na.rm=TRUE)) %>%
  #filter(!is.na(Genus)&Genus !="")
  #head(genus_level)
  
  mutate(percent_abundance = (Total_Abundance / sum(Total_Abundance)) * 100) %>%
  arrange(desc(percent_abundance)) %>%
  slice_head(n = 10)
print(genus_level)


library(vegan)
library(ggplot2)

#loading metadata

meta.in = read.csv("E:/Thami-Uni/Oil Microcosm Data/Genohub IDs Project 2599682, 3694692.csv", 
                   header = TRUE)
rownames(meta.in) = meta.in$uniq.id
meta.in = meta.in[rownames(seqtab.new),]

# Read the taxadata 
data <- read.csv("working_data_avg1.csv", stringsAsFactors = FALSE)



# Filter for CTL and ND treatments only
data_filtered <- data %>%
  filter(NutDisp %in% c("CTL", "ND")) %>%
  filter(MeanValue > 0)  # Remove zero abundance ASVs


community_matrix <- data_filtered %>%
  select(Sample_base, ASV, MeanValue) %>%
  pivot_wider(names_from = ASV, values_from = MeanValue, values_fill = 0) %>%
  column_to_rownames("Sample_base")

bray_dist <- vegdist(community_matrix, method = "bray")

set.seed(123)  # for reproducibility
nmds_result <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

nmds_points <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_points$SampleID <- rownames(nmds_points)

nmds_species <- as.data.frame(scores(nmds_result, display = "species"))
nmds_species$Taxa <- rownames(nmds_species)



rn2 <- rownames(meta.in)

# get the index positions of start and end wells
start <- which(rn2 == "plate2-E06")
end   <- which(rn2 == "plate2-F10")
meta.filtered <- meta.in[start:end, ]   # subset rows first
meta.filtered <- meta.filtered[rownames(meta.filtered) != "plate2-F02", ]

Meta_data <- meta.filtered %>%
  mutate(Sample_base = paste(Experiment,NutDisp, Time, sep = "_"))

colnames(nmds_points)
colnames(Meta_data)
head(nmds_points$SampleID)
head(Meta_data$Sample_base)




merged_points <- merge(nmds_points, Meta_data, by.x = "SampleID", by.y = "Sample_base")

nmds_result$stress

shape_map <- c(
  "1WK"   = 8,   # star
  "2.5WK" = 16,   # circle
  "4WK"   = 17,   # triangle
  "6WK"   = 15   # square
)

custom_colors <- c(
  "CTL" = "#009999",
  "ND"  = "#660066"
  
)

NMDS <- ggplot(merged_points, aes(x = NMDS1, y = NMDS2, 
                          color = NutDisp, shape = Time)) +
  geom_point(size = 3) +
  #stat_ellipse(aes(group = NutDisp), linetype = 1) +
  stat_ellipse(aes(group = NutDisp, fill = NutDisp),   # fill color
               alpha = 0.2,                              # transparency
               color = NA,                               # remove edge color if desired
               geom = "polygon") +    # make sure group matches color variable
  scale_shape_manual(values = shape_map) + 
  stat_ellipse(aes(group = NutDisp), level = 0.95, linetype = 1) +
  annotate("text", x = max(merged_points$NMDS1), 
           y = min(merged_points$NMDS2), 
           label = paste("Stress =", round(nmds_result$stress, 3)),
           hjust = -1, vjust = 4, size = 15) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors, guide = "none") +  # no duplicate legend
  
  theme_bw() +
  theme(
    
    #plot.title = element_text(family= "Arial", size = 30, face = "bold", hjust = 0.5),
    #plot.subtitle = element_text(family= "Arial", size = 30, hjust = 0.5),
    axis.title = element_text(size = 55),
    axis.text = element_text(size = 45),
    legend.position = "right",
    legend.title = element_text(size = 35, face = "bold"),  # title font size
    legend.text  = element_text(size = 35), 
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank()    # remove minor gridlines
  )+
  labs(
    x = "NMDS1", y = "NMDS2",
    shape = "Time",
    color= "Treatment Group")
print(NMDS)

nmds_result$stress

unique(merged_points$Time)

ggsave("NMDS.png", plot = NMDS, 
       width = 10, height = 7, dpi = 300, bg= "white")


permanova <- adonis2(seqtab.numeric ~ NutDisp, data = meta.filtered, method = "bray", permutations = 999)
print(permanova)

str(seqtab.filtered.new)
seqtab.numeric <- seqtab.filtered.new %>%
  select(where(is.numeric))
#seqtab.numeric <- seqtab.numeric[match(meta.filtered$SampleID, rownames(seqtab.numeric)), ]
#all(rownames(seqtab.numeric) == meta.filtered$SampleID)


permanova <- adonis2(seqtab.numeric ~ NutDisp,
                     data = meta.filtered,
                     method = "bray",
                     permutations = 999)

print(permanova)

# ANOSIM
ano <- anosim(seqtab.numeric, grouping = meta.filtered$NutDisp, distance = "bray")
print(ano)

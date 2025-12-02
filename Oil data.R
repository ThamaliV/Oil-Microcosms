
library(gplots)
#source("~/apps/FastPCA.R")
library(dendextend)
library(Rtsne)
library(vegan)
library(fastcluster)
library(ggplot2)
library(tidyverse)
library(scales)
library(viridis)
library(tibble)
library(tidyr)



seqtab <- readRDS("E:/Thami-Uni/Oil Microcosm Data/seqtab_final.rds")



colnames(seqtab)
rownames(seqtab)
colnames(
  rownames(seqtab) = gsub(pattern = "716S",replacement = "7-16S",x = rownames(seqtab)),
  rownames(seqtab) = gsub(pattern = "716S",replacement = "7-16S",x = rownames(seqtab)))
  

  
  
  sepnames.16S = separate(data=data.frame(x=rownames(seqtab)),col=x, into=c("project","primer","plate","well"),sep="[|]", fill="right")
  sepnames.16S$plate.geno = toupper(paste(sepnames.16S$platenum$p2,sepnames.16S$well,sep="_"))
  sepnames.16S$platenum = separate(data=sepnames.16S,col=plate, into=c("p1","p2","p3","p4","p5","p6"), sep="-", fill="right")
  sepnames.16S$plate.geno = toupper(paste(sepnames.16S$platenum$p2,sepnames.16S$well,sep="_"))
  seqtab.new = seqtab[grep(x=rownames(seqtab),"20250811"),]
  
  head(rownames(seqtab.new))
  seqtab.new = seqtab.new[,colSums(seqtab.new)>0]
  rownames(seqtab.new) = sepnames.16S$platenum$well[grep("20250811",rownames(seqtab))]
  
  #seqtab.new.both = cbind(seqtab.new)
  
  meta.in = read.csv("E:/Thami-Uni/Oil Microcosm Data/Genohub IDs Project 2599682, 3694692.csv", 
                     header = TRUE)
  rownames(meta.in) = meta.in$uniq.id
  meta.in = meta.in[rownames(seqtab.new),]
  
  taxout.16S = as.data.frame(readRDS("E:/Thami-Uni/Oil Microcosm Data/taxout_edit.rds")$`ref_dada2_silva_v138-2_train_set_uniq.fasta.gz`)
  taxout.16S$ASV = paste0("bASV",1:nrow(taxout.16S))
  taxout.16S$short = paste(taxout.16S$ASV,taxout.16S$Domain, taxout.16S$Order, taxout.16S$Genus)
  
  taxout.all = taxout.16S
  seqtab.new.small = seqtab.new[,colMeans(seqtab.new)>10]
  
  dim(seqtab.new.small)
  head(seqtab.new.small)
  taxnames = taxout.all[colnames(seqtab.new.small),"short"]
  
  length(seqtab.new.small)
  
  
  rn <- rownames(seqtab.new.small)
  
  # get the index positions of start and end wells
  start <- which(rn == "plate2-E06")
  end   <- which(rn == "plate2-F10")
  
  # subset rows between them
  
  seqtab.filtered <- seqtab.new.small[start:end, ]   # subset rows first
  seqtab.filtered <- seqtab.filtered[rownames(seqtab.filtered) != "plate2-F02", ]
   
  
  #setab.filtered  <- seqtab.new.small[start:end, rownames(seqtab.new.small)[start:end] != "plate2-F02"]
  
  seqtab.filtered.new <- as.data.frame(seqtab.filtered)
  
  seqtab.filtered.new <- cbind("uniq.id" = rownames(seqtab.filtered.new), seqtab.filtered.new)
  rownames(seqtab.filtered.new) <- NULL 
  

seqtab_long <- seqtab.filtered.new %>%
    pivot_longer(
      cols = -uniq.id,             # all columns except SampleID
      names_to = "Sequence",             # new column for ASV IDs
      values_to = "Abundance" )      # new column for counts
  
merged_seqtab <- seqtab_long %>%
    left_join(meta.in, by = "uniq.id")
  
taxout.all <- rownames_to_column(taxout.all, var = "Sequence")
  
working_data <- merged_seqtab %>%
    left_join(taxout.all, by = "Sequence") 


working_data <- working_data %>%
  mutate(Sample_base = paste(Experiment,NutDisp, Time, sep = "_"))


getwd()

working_data_avg <- working_data %>%
    group_by(Sample_base, Site, Time, NutDisp, ASV, Domain, Phylum, Class, Order, Family, Genus) %>%
    summarise(MeanValue = mean(Abundance, na.rm = TRUE))

write.csv(working_data_avg, "working_data_avg1.csv", row.names = FALSE)

genus_filtered <- genus_rel_summary %>%
  filter(NutDisp != "PRE")

# Create the bubble plot without facets, single Y-axis
bubble_plot <- ggplot(genus_filtered, aes(x = Time, y = Genus)) +
  geom_point(aes(size = RelAbundance, color = NutDisp), 
             alpha = 0.8) +
  facet_wrap(~ NutDisp, scales = "free_x") +
  scale_size_continuous(range = c(1, 20), 
                        name = "Relative\nAbundance") +
  scale_color_manual(values = c("CTL" = "#2E8B57",  # Sea Green
                                "ND" = "#FF6347"),  # Tomato
                     name = "Treatment") +
  labs(
       
       x = "Time",
       y = "Genus") +
  theme_classic() +
  theme(
    axis.text.x = element_text( size= 25, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size=25),
    strip.text = element_text(size = 20, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size= 25),
    legend.title = element_text(size= 25)
  )

print(bubble_plot)

ggsave("bubble_micro.png", plot = bubble_plot, 
       width = 16, height = 8, dpi = 2000)

length(unique(working_data_avg$ASV))
length(unique(working_data_avg$Genus))
length(unique(working_data_avg$Phylum))
length(unique(working_data_avg$Family))
length(unique(working_data_avg$Class))
length(unique(working_data_avg$Order))

length(seqtab.filtered.new)


col_pct <- genus_rel%>% filter(NutDisp == "PRE" ,Genus== "Colwellia")%>%
  summarise(Perc_Abundance= sum(RelAbundance,na.rm=TRUE)*100) #%>%

print(col_pct)


# T-Test Analysis for TOTAL PAH Percentage Loss

library(tidyverse)
library(flextable)
library(officer)


# Read the raw concentration data
bdloss <- read.csv("E:/Thami-Uni/Oil Microcosm Data/BDLOSS.csv")

# Clean up - remove empty rows
bdloss <- bdloss %>% 
  filter(!is.na(Treatment_type) & Treatment_type != "")

# Clean up Treatment_type (remove extra spaces)
bdloss$Treatment_type <- trimws(bdloss$Treatment_type)

# Check treatment types
print("Treatment types in data:")
print(table(bdloss$Treatment_type))



# Sum all PAH compounds for each row (replicate)
bdloss$Total_PAH <- rowSums(bdloss[, -1], na.rm = TRUE)

# View the total PAH concentrations
total_pah_data <- bdloss %>%
  select(Treatment_type, Total_PAH)

print("Total PAH concentrations by treatment:")
print(total_pah_data)


# Calculate Mean Total PAH for Baselines


mean_total_pah <- total_pah_data %>%
  group_by(Treatment_type) %>%
  summarise(mean_total = mean(Total_PAH, na.rm = TRUE),
            .groups = "drop")


print(mean_total_pah)

# Get baseline values
raw_mean <- mean_total_pah$mean_total[mean_total_pah$Treatment_type == "Raw Diesel"]
uv_mean <- mean_total_pah$mean_total[mean_total_pah$Treatment_type == "UV weathered"]



#Calculate Total % Loss for Each Replicate


total_pct_loss <- total_pah_data %>%
  mutate(
    Pct_Loss = case_when(
      Treatment_type == "UV weathered" ~ ((raw_mean - Total_PAH) / raw_mean) * 100,
      Treatment_type == "CTRL" ~ ((uv_mean - Total_PAH) / uv_mean) * 100,
      Treatment_type == "ND" ~ ((uv_mean - Total_PAH) / uv_mean) * 100,
      TRUE ~ NA_real_
    ),
    Comparison = case_when(
      Treatment_type == "UV weathered" ~ "Raw to UV",
      Treatment_type == "CTRL" ~ "UV to Control",
      Treatment_type == "ND" ~ "UV to Nutri+Disp",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Comparison))

print("Total % Loss for each replicate:")
print(total_pct_loss)

# Summary Statistics for Total % Loss


total_loss_summary <- total_pct_loss %>%
  group_by(Comparison) %>%
  summarise(
    n = n(),
    mean_loss = mean(Pct_Loss, na.rm = TRUE),
    sd_loss = sd(Pct_Loss, na.rm = TRUE),
    se_loss = sd_loss / sqrt(n),
    .groups = "drop"
  )


print(total_loss_summary)


# One-Sample T-Tests (Is % loss significantly different from 0?)

# Raw to UV

raw_to_uv_data <- total_pct_loss %>% filter(Comparison == "Raw to UV")
if(nrow(raw_to_uv_data) >= 2) {
  ttest_raw_uv <- t.test(raw_to_uv_data$Pct_Loss, mu = 0)
  print(ttest_raw_uv)
} else {
  cat("Insufficient replicates for t-test\n")
}

# UV to Control

uv_to_ctrl_data <- total_pct_loss %>% filter(Comparison == "UV to Control")
if(nrow(uv_to_ctrl_data) >= 2) {
  ttest_uv_ctrl <- t.test(uv_to_ctrl_data$Pct_Loss, mu = 0)
  print(ttest_uv_ctrl)
} else {
  cat("Insufficient replicates for t-test\n")
}

# UV to Nutri+Disp

uv_to_nd_data <- total_pct_loss %>% filter(Comparison == "UV to Nutri+Disp")
if(nrow(uv_to_nd_data) >= 2) {
  ttest_uv_nd <- t.test(uv_to_nd_data$Pct_Loss, mu = 0)
  print(ttest_uv_nd)
} else {
  cat("Insufficient replicates for t-test\n")
}


# Two-Sample T-Test: Control vs Nutri+Disp


ctrl_data <- total_pct_loss %>% filter(Comparison == "UV to Control") %>% pull(Pct_Loss)
nd_data <- total_pct_loss %>% filter(Comparison == "UV to Nutri+Disp") %>% pull(Pct_Loss)

if(length(ctrl_data) >= 2 & length(nd_data) >= 2) {
  ttest_ctrl_nd <- t.test(ctrl_data, nd_data)
  print(ttest_ctrl_nd)
} else {
  cat("Insufficient replicates for t-test\n")
}


# Create Results Summary Table


results_table <- data.frame(
  Comparison = c("Raw to UV", "UV to Control", "UV to Nutri+Disp"),
  n = c(nrow(raw_to_uv_data), nrow(uv_to_ctrl_data), nrow(uv_to_nd_data)),
  Mean_Pct_Loss = c(
    mean(raw_to_uv_data$Pct_Loss),
    mean(uv_to_ctrl_data$Pct_Loss),
    mean(uv_to_nd_data$Pct_Loss)
  ),
  SD = c(
    sd(raw_to_uv_data$Pct_Loss),
    sd(uv_to_ctrl_data$Pct_Loss),
    sd(uv_to_nd_data$Pct_Loss)
  ),
  SE = c(
    sd(raw_to_uv_data$Pct_Loss) / sqrt(nrow(raw_to_uv_data)),
    sd(uv_to_ctrl_data$Pct_Loss) / sqrt(nrow(uv_to_ctrl_data)),
    sd(uv_to_nd_data$Pct_Loss) / sqrt(nrow(uv_to_nd_data))
  ),
  t_statistic = c(
    ifelse(exists("ttest_raw_uv"), ttest_raw_uv$statistic, NA),
    ifelse(exists("ttest_uv_ctrl"), ttest_uv_ctrl$statistic, NA),
    ifelse(exists("ttest_uv_nd"), ttest_uv_nd$statistic, NA)
  ),
  p_value = c(
    ifelse(exists("ttest_raw_uv"), ttest_raw_uv$p.value, NA),
    ifelse(exists("ttest_uv_ctrl"), ttest_uv_ctrl$p.value, NA),
    ifelse(exists("ttest_uv_nd"), ttest_uv_nd$p.value, NA)
  )
)

results_table <- results_table %>%
  mutate(
    Significance = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    Mean_SE = sprintf("%.2f ± %.2f", Mean_Pct_Loss, SE)
  )


print(results_table)


write.csv(results_table, "ttest_total_pah_loss.csv", row.names = FALSE)
write.csv(total_pct_loss, "total_pah_pct_loss_data.csv", row.names = FALSE)




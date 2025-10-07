# Description: This script analyzes and visualizes the incidence rates of IBD, CD, UC, and CRC globally and in India.

library(ggplot2)
library(dplyr)
library(scales)

# Incidence data

# Global data (GIVES-21 data)

UC <- read.csv("Datasets/Incidence_data/Incidence_UC.csv")
CD <- read.csv("Datasets/Incidence_data/Incidence_CD.csv")

# avg incidence calculation
avg_IBD_global <- 4.97
avg_CRC_global <- 18.4
avg_CD_global <- mean(as.numeric(CD$incidence_rate))   # mean reported incidence
avg_UC_global <- mean(as.numeric(UC$incidence_rate))    # mean reported incidence

# India data

avg_IBD_ind  <- 2.34
avg_CRC_ind  <- 4.9
avg_CD_ind   <- CD$incidence_rate[CD$country == "India"]
avg_UC_ind   <-  UC$incidence_rate[UC$country == "India"]

global <- data.frame(
  Condition = c("IBD (ASR)", "CD (mean reported incidence)", "UC (mean reported incidence)", "CRC (ASR)"),
  Average_Incidence = c(avg_IBD_global, avg_CD_global, avg_UC_global, avg_CRC_global)
)

india <- data.frame(
  Condition = c("IBD (ASR)", "CD (mean reported incidence)", "UC (mean reported incidence)", "CRC (ASR)"),
  Average_Incidence = c(avg_IBD_ind, avg_CD_ind, avg_UC_ind, avg_CRC_ind)
)

global$Source <- "Global"
india$Source  <- "India"

plot_df <- bind_rows(global, india)

order_levels <- c("CRC (ASR)","IBD (ASR)","UC (mean reported incidence)","CD (mean reported incidence)")

plot_df$Condition <- factor(plot_df$Condition, levels = order_levels)

# Grouped bar plot

p <- ggplot(plot_df, aes(x = Condition, y = Average_Incidence, fill = Source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "black", size = 0.15) +
  geom_text(aes(label = round(Average_Incidence, 2)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = comma) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(
    title = "Global versus Indian Incidence Burden of IBD and Colorectal Cancer",
    subtitle = "Mean and age-standardized incidence rates (per 100,000 person-years, 1990–2022)",
       x = NULL,
       y = "Incidence (per 100,000 / year)",
       fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(face = "bold"))

print(p)

p + scale_y_continuous(trans = "log10", labels = scales::comma)

p <- p +
  scale_fill_manual(values = c("Global" = "#009E73", "India" = "#D55E00")) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.position = "top",
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    plot.margin = margin(10, 20, 10, 20)
  )




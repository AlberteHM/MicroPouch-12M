library(tidyverse)

# Load data 
metadata <- read_delim("data/metadata.csv") %>%  
  filter(
    project != "donor_batch",
    !(stage %in% c("treatment_5", "treatment_10", "treatment_15", "treatment_21"))) %>%
  mutate(stage = factor(stage, levels= c("inclusion", "followup_30d", "followup_1m", "followup_3m", "followup_6m", "followup_12m")))


# Function - Assign color to id 
fmt_cols <- c(
  "#5B3A2E","#854442","#9B5B59","#A06A4E","#C08E6D","#E0C6A8",
  "#B08A4F","#95825A","#5f725d","#3F5A2A","#577533","#7F9A4A",
  "#2F4B2E","#A9B98E", "#6E8B74"
)

plac_cols <- c(
  "#6FA7B0","#93C6CF","#7AA6D9","#4E79A7","#9B6B7A","#6F5068",
  "#B7A0B0","#2F5F6A","#2F4B7C","#4B8FA0","#B9DDE3","#5C86C3",
  "#1F3E5A","#A9879C","#5A3F56"
)

make_col_map <- function(id_key, fmt_cols, plac_cols,
                         group_fmt = "FMT", group_plac = "placebo",
                         id_col = "id", group_col = "group") {
  
  key <- id_key %>%
    distinct(id = .data[[id_col]], group = .data[[group_col]]) %>%
    mutate(id = as.character(id), group = as.character(group))
  
  fmt_ids  <- key %>% filter(group == group_fmt)  %>% pull(id) %>% sort()
  plac_ids <- key %>% filter(group == group_plac) %>% pull(id) %>% sort()
  
  c(
    setNames(rep(fmt_cols,  length.out = length(fmt_ids)),  fmt_ids),
    setNames(rep(plac_cols, length.out = length(plac_ids)), plac_ids)
  )
}

id_key <- metadata %>% distinct(id, group)
col_map <- make_col_map(id_key, fmt_cols, plac_cols)


# Stool frequency plot 
n_counts <- metadata %>%
  group_by(group, stage) %>%
  summarise(n = n_distinct(id), .groups = "drop")

stool_fre_plot <- ggplot(metadata, aes(x = stage, y = stool_fre, color = id)) +
  facet_grid(. ~ group, labeller = as_labeller(c("FMT" ="FMT", "placebo" = "Placebo"))) +
  geom_line(aes(group = id), alpha = 0.8, position = position_dodge(0.4), linewidth=1) +
  geom_point(aes(group = id), size = 3,  position = position_dodge(0.4)) +
  scale_x_discrete(labels = c("Baseline", "0M", "1M", "3M", "6M", "12M")) + 
  scale_color_manual(name = NULL, values = col_map) + 
  labs(x = "", y = "Stool frequency [stools/day]") + 
  scale_y_continuous(limits = c(0, NA)) + 
  theme(panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"), 
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
        axis.ticks.x = element_line(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "none", 
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10), 
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12)
  ) + 
  geom_bracket(data = subset(metadata, group == "FMT") %>%
                 distinct(group, .keep_all = TRUE),
               aes(xmin = "inclusion", xmax = "followup_1m", y.position = 26.3),
               label = "p < 0.05", label.size = 11/.pt, color = "grey20") + 
  geom_bracket(data = subset(metadata, group == "placebo") %>%
                 distinct(group, .keep_all = TRUE), 
               aes(xmin = "inclusion", xmax = "followup_3m", y.position = 22),
               label = "p < 0.05", label.size = 11/.pt, color = "grey20") + 
  expand_limits(y=28.3) +
  geom_label(
    data = n_counts,
    aes(x = stage, y = 27.5, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 10/.pt,
    label.size = 0.2,
    fill="grey97",
    alpha = 0.8
  )


# cPDAI plot 

cpdai_plot <- ggplot(metadata, aes(x = stage, y = cpdai_sum, color = id)) +
  facet_grid(. ~ group, labeller = as_labeller(c("FMT" ="FMT", "placebo" = "Placebo"))) +
  geom_line(aes(group = id), alpha = 0.8, position = position_dodge(0.4), linewidth=1) +
  geom_point(aes(group = id), size = 3,  position = position_dodge(0.4)) +
  scale_x_discrete(labels = c("Baseline", "0M", "1M", "3M", "6M", "12M")) + 
  scale_color_manual(name = NULL, values = col_map) + 
  labs(x = "", y = "Clinical PDAI") + 
  theme(panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"), 
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
        axis.ticks.x = element_line(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "none", 
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10), 
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  geom_bracket(data = subset(metadata, group == "FMT") %>%
                 distinct(group, .keep_all = TRUE), 
               aes(xmin = "inclusion", xmax = "followup_30d", y.position = 6.5),
               label = "p < 0.01", label.size = 11/.pt, color = "grey20") +
  geom_bracket(data = subset(metadata, group == "FMT") %>%
                 distinct(group, .keep_all = TRUE), 
               aes(xmin = "inclusion", xmax = "followup_1m", y.position = 6.9),
               label = "p < 0.001", label.size = 11/.pt, color = "grey20") +
  geom_bracket(data = subset(metadata, group == "FMT") %>%
                 distinct(group, .keep_all = TRUE), 
               aes(xmin = "inclusion", xmax = "followup_3m", y.position = 7.3),
               label = "p < 0.05", label.size = 11/.pt, color = "grey20") +
  geom_bracket(data = subset(metadata, group == "FMT") %>%
                 distinct(group, .keep_all = TRUE),
               aes(xmin = "inclusion", xmax = "followup_6m", y.position = 7.7),
               label = "p < 0.05", label.size = 11/.pt, color = "grey20") +
  geom_bracket(data = subset(metadata, group == "placebo") %>%
                 distinct(group, .keep_all = TRUE),
               aes(xmin = "inclusion", xmax = "followup_30d", y.position = 5.5),
               label = "p < 0.001", label.size = 11/.pt, color = "grey20") +
  geom_bracket(data = subset(metadata, group == "placebo") %>%
                 distinct(group, .keep_all = TRUE), 
               aes(xmin = "inclusion", xmax = "followup_1m", y.position = 5.9),
               label = "p < 0.05", label.size = 11/.pt, color = "grey20") + 
  geom_bracket(data = subset(metadata, group == "placebo") %>%
                 distinct(group, .keep_all = TRUE), 
               aes(xmin = "inclusion", xmax = "followup_3m", y.position = 6.3),
               label = "p < 0.001", label.size = 11/.pt, color = "grey20") + 
  scale_y_continuous(
    limits = c(0, 8.3),
    breaks = c(0, 2, 4, 6)) + 
  geom_label(
    data = n_counts,
    aes(x = stage, y = 8.0, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 10/.pt,
    label.size = 0.2,
    fill="grey97",
    alpha = 0.8
  )

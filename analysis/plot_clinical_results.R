library(tidyverse)
library(ggpubr)

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

stool_fre_stats <- read_delim("data/within_group_stoolfre_stats.tsv")

get_p_label_stool_fre <- function(xmin, xmax, group_name) {
  stool_fre_stats %>%
    filter(
      group == group_name,
      from_timepoint == xmin,
      to_timepoint == xmax
    ) %>%
    pull(p_val_bf_adj) %>%
    {
      case_when(
        length(.) == 0 ~ NA_character_,
        is.na(.) ~ NA_character_,
        . < 0.0001 ~ "p<0.0001",
        . < 0.001  ~ "p<0.001",
        . < 0.01   ~ "p<0.01",
        . < 0.05   ~ "p<0.05",
        TRUE       ~ "p>0.05"
      )
    }
}

add_sig_bracket_stool_fre <- function(xmin, xmax, group_name, y.position) {
  
  p_label <- get_p_label_stool_fre(xmin, xmax, group_name)
  
  if (is.na(p_label) || p_label == "p>0.05") {
    return(NULL)
  }
  
  geom_bracket(
    data = metadata %>%
      filter(group == group_name) %>%
      distinct(group, .keep_all = TRUE),
    aes(
      xmin = xmin,
      xmax = xmax,
      y.position = y.position
    ),
    label = p_label,
    label.size = 11 / .pt,
    family = "Helvetica",
    color = "grey20"
  )
}


# Stool frequency plot 
n_counts <- metadata %>%
  filter(!is.na(stool_fre)) %>% 
  group_by(group, stage) %>%
  summarise(n = n_distinct(id), .groups = "drop")

stool_fre_plot_front <- ggplot(metadata, aes(x = stage, y = stool_fre, color = id)) +
  facet_grid(. ~ group, labeller = as_labeller(c("FMT" ="FMT", "placebo" = "Placebo"))) +
  geom_line(aes(group = id), alpha = 0.8, position = position_dodge(0.4), linewidth=1) +
  geom_point(aes(group = id), size = 3,  position = position_dodge(0.4)) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 1.1) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "black", size = 3, shape = 18) + 
  scale_x_discrete(labels = c("Baseline", "0M", "1M", "3M", "6M", "12M")) + 
  scale_color_manual(name = NULL, values = col_map) + 
  labs(x = "", y = "Stool frequency [stools/day]") + 
  scale_y_continuous(limits = c(0, 28), 
                     breaks = c(0, 10, 20)) + 
  theme(text = element_text(family = "Helvetica"), 
        panel.background = element_rect(fill="grey97"),
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
  geom_label(
    data = n_counts,
    aes(x = stage, y = 27.5, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 10/.pt,
    label.size = 0.2,
    fill="grey97",
    alpha = 0.8
  ) + 
  add_sig_bracket_stool_fre("inclusion", "followup_30d", "FMT", 25) +
  add_sig_bracket_stool_fre("inclusion", "followup_1m", "FMT", 26.3) +
  add_sig_bracket_stool_fre("inclusion", "followup_3m", "FMT", 27) +
  add_sig_bracket_stool_fre("inclusion", "followup_6m", "FMT", 28) +
  add_sig_bracket_stool_fre("inclusion", "followup_12m", "FMT", 29) +
  add_sig_bracket_stool_fre("inclusion", "followup_30d", "placebo", 22) +
  add_sig_bracket_stool_fre("inclusion", "followup_1m", "placebo", 26) +
  add_sig_bracket_stool_fre("inclusion", "followup_3m", "placebo", 24) +
  add_sig_bracket_stool_fre("inclusion", "followup_6m", "placebo", 28) +
  add_sig_bracket_stool_fre("inclusion", "followup_12m", "placebo", 29) 

stool_fre_plot 

#ggsave("./plots/figure_2.tiff", plot=stool_fre_plot, dpi = 300, w=20.375, h=9)


# cPDAI plot 

cpdai_stats <- read_delim("data/within_group_cpdai_stats.tsv")

n_counts <- metadata %>%
  filter(!is.na(cpdai_sum)) %>% 
  group_by(group, stage) %>%
  summarise(n = n_distinct(id), .groups = "drop")

get_p_label_cpdai <- function(xmin, xmax, group_name) {
  cpdai_stats %>%
    filter(
      group == group_name,
      from_timepoint == xmin,
      to_timepoint == xmax
    ) %>%
    pull(p_val_bf_adj) %>%
    {
      case_when(
        length(.) == 0 ~ NA_character_,
        is.na(.) ~ NA_character_,
        . < 0.0001 ~ "p<0.0001",
        . < 0.001  ~ "p<0.001",
        . < 0.01   ~ "p<0.01",
        . < 0.05   ~ "p<0.05",
        TRUE       ~ "p>0.05"
      )
    }
}

add_sig_bracket_cpdai <- function(xmin, xmax, group_name, y.position) {
  
  p_label <- get_p_label_cpdai(xmin, xmax, group_name)
  
  if (is.na(p_label) || p_label == "p>0.05") {
    return(NULL)
  }
  
  geom_bracket(
    data = metadata %>%
      filter(group == group_name) %>%
      distinct(group, .keep_all = TRUE),
    aes(
      xmin = xmin,
      xmax = xmax,
      y.position = y.position
    ),
    label = p_label,
    label.size = 11 / .pt,
    family = "Helvetica",
    color = "grey20"
  )
}


cpdai_plot <- ggplot(metadata, aes(x = stage, y = cpdai_sum, color = id)) +
  facet_grid(. ~ group, labeller = as_labeller(c("FMT" ="FMT", "placebo" = "Placebo"))) +
  geom_line(aes(group = id), alpha = 0.8, position = position_dodge(0.4), linewidth=1) +
  geom_point(aes(group = id), size = 3,  position = position_dodge(0.4)) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 1.1) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "black", size = 3, shape = 18) +
  scale_x_discrete(labels = c("Baseline", "0M", "1M", "3M", "6M", "12M")) + 
  scale_color_manual(name = NULL, values = col_map) + 
  labs(x = "", y = "Clinical PDAI") + 
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill="grey97"),
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
  ) +
  add_sig_bracket_cpdai("inclusion", "followup_30d", "FMT", 6.5) +
  add_sig_bracket_cpdai("inclusion", "followup_1m", "FMT", 6.9) +
  add_sig_bracket_cpdai("inclusion", "followup_3m", "FMT", 7.3) +
  add_sig_bracket_cpdai("inclusion", "followup_6m", "FMT", 7.7) +
  add_sig_bracket_cpdai("inclusion", "followup_12m", "FMT", 11) +
  add_sig_bracket_cpdai("inclusion", "followup_30d", "placebo", 5.5) +
  add_sig_bracket_cpdai("inclusion", "followup_1m", "placebo", 5.9) +
  add_sig_bracket_cpdai("inclusion", "followup_3m", "placebo", 6.3) +
  add_sig_bracket_cpdai("inclusion", "followup_6m", "placebo", 6.7) +
  add_sig_bracket_cpdai("inclusion", "followup_12m", "placebo", 7.1) 
  
#ggsave("./plots/figure_3.tiff", plot=cpdai_plot, dpi = 300, w=20.375, h=9)

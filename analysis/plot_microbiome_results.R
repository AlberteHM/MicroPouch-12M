library(ampvis2)
library(tidyverse)
library(vegan)
library(patchwork)
# Load data 
metadata <- read_delim("data/metadata.csv") 
metaphlan <- read_delim("data/MetaPhlAn_4.1.0.txt")

  
# Functions 
## Assign color to id 
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

## Statistics
stage_levels <- c(
  "inclusion", "followup_30d", "followup_1m",
  "followup_3m", "followup_6m", "followup_12m"
)

stage_pairs <- combn(stage_levels, 2, simplify = FALSE)

run_paired_wilcox <- function(df_group, value_col) {
  value_col <- enquo(value_col)
  
  map_dfr(stage_pairs, function(pair) {
    s1 <- pair[1]
    s2 <- pair[2]
    
    dat_pair <- df_group %>%
      filter(stage %in% c(s1, s2)) %>%
      select(id, stage, value = !!value_col) %>%
      distinct()
    
    dat_wide <- dat_pair %>%
      mutate(stage = as.character(stage)) %>%
      pivot_wider(names_from = stage, values_from = value) %>%
      filter(!is.na(.data[[s1]]) & !is.na(.data[[s2]]))
    
    n_pairs <- nrow(dat_wide)
    
    warning_msg <- NA_character_
    
    if (n_pairs > 0) {

      test <- tryCatch(
        withCallingHandlers(
          wilcox.test(
            dat_wide[[s2]],
            dat_wide[[s1]],
            paired = TRUE,
            exact = FALSE,
            conf.int = TRUE,
            conf.level = 0.95
          ),
          warning = function(w) {  
            warning_msg <<- conditionMessage(w)
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) NULL
      )
      
      
      if (!is.null(test)) {
        p_val <- test$p.value
        effect_estimate <- unname(test$estimate)
        ci_lower <- unname(test$conf.int[1])
        ci_upper <- unname(test$conf.int[2])
      } else {
        p_val <- NA_real_
        effect_estimate <- NA_real_
        ci_lower <- NA_real_
        ci_upper <- NA_real_
      }
      
    } else {
      p_val <- NA_real_
      effect_estimate <- NA_real_
      ci_lower <- NA_real_
      ci_upper <- NA_real_
    }
    
    tibble(
      group = unique(df_group$group),
      comparison = paste(s1, "~", s2),
      n = n_pairs,
      effect_estimate = effect_estimate,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_val, 
      warning = warning_msg
    )
  }) %>%
    mutate(
      p_bonf_15 = p.adjust(p_value, method = "bonferroni"),
      p_bonf_5 = if_else(
        grepl("^inclusion ~ ", comparison),
        pmin(p_value * 5, 1),
        NA_real_
      ),
      sig_diff = if_else(
        p_value < 0.05 |
          p_bonf_15 < 0.05 |
          coalesce(p_bonf_5 < 0.05, FALSE),
        "yes",
        NA_character_
      )
    )
}




# Richness 
metadata_richness <- metadata %>%  
  filter(
    project != "donor_batch",
    !(stage %in% c("treatment_5", "treatment_10", "treatment_15", "treatment_21"))) %>%
  mutate(stage = factor(stage, levels= c("inclusion", "followup_30d", "followup_1m", "followup_3m", "followup_6m", "followup_12m")),
         id = factor(id))
  
richness_calc <- metaphlan %>% 
  column_to_rownames(var = "clade_name") %>% 
  t() %>%  
  specnumber() %>%  
  enframe() %>%  
  rename("sample_barcode" = "name", 
         "richness" = "value")

richness_df <- metadata_richness %>%  
  left_join(richness_calc, by="sample_barcode") %>%  
  filter(!is.na(richness))


## Paired wilcoxon test
richness_wilcox_fmt <- richness_df %>%
  filter(group == "FMT") %>%
  run_paired_wilcox(richness)

richness_wilcox_placebo <- richness_df %>%
  filter(group == "placebo") %>%
  run_paired_wilcox(richness)

## Plot
n_counts <- richness_df %>%
  filter(!is.na(richness)) %>% 
  group_by(group, stage) %>%
  summarise(n = n_distinct(id), .groups = "drop")


richness_plot <- ggplot(richness_df, aes(x = stage, y = richness, color = id)) +
  facet_grid(. ~ group, labeller = as_labeller(c("FMT" ="FMT", "placebo" = "Placebo"))) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 1.1) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "black", size = 3, shape = 18) + 
  geom_line(aes(group = id), alpha = 0.8, position = position_dodge(0.3), linewidth=1) +
  geom_point(aes(group = id), alpha = 0.9, size = 3,  position = position_dodge(0.3)) +
  scale_x_discrete(labels = c("Baseline", "0M", "1M", "3M", "6M", "12M")) + 
  scale_color_manual(name = NULL, values = col_map) + 
  labs(x = "", y = "Richness") + 
  scale_y_continuous(limits = c(0, NA)) + 
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
  expand_limits(y = 110) +
  geom_label(
    data = n_counts,
    aes(x = stage, y = 105, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 10/.pt,
    label.size = 0.2,
    fill="grey97",
    alpha = 0.8
  )


# Shannon diversity
metadata_shannon <- metadata %>%  
  filter(!(stage %in% c("treatment_5", "treatment_10", "treatment_15", "treatment_21"))) %>%
  mutate(stage = factor(stage, levels= c("inclusion", "followup_30d", "followup_1m", "followup_3m", "followup_6m", "followup_12m")),
         id = factor(id))

otutable_taxtable <- metaphlan %>%  
  mutate(OTU = paste0("OTU", row_number())) %>% 
  relocate(OTU) 

otutable <- otutable_taxtable %>% 
  select(-clade_name)

taxtable <- otutable_taxtable %>%
  select(OTU, clade_name) %>%
  mutate(clade_name = if_else(clade_name == "UNCLASSIFIED",
                              "k__UNCLASSIFIED|p__UNCLASSIFIED|c__UNCLASSIFIED|o__UNCLASSIFIED|f__UNCLASSIFIED|g__UNCLASSIFIED|s__UNCLASSIFIED", clade_name)) %>%
  separate(clade_name, sep = "\\|", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  relocate(OTU, .after = "Species") %>%
  mutate(Species = Species %>%
           sub("^s__", "", .) %>%  
           gsub("_", " ", .)) %>% 
  slice(-1)


amp_object_shannon <- amp_load(otutable = otutable, metadata = metadata_shannon, taxonomy = taxtable)
  
shannon_df <- amp_alpha_diversity(amp_object_shannon, measure = "shannon", richness = FALSE, rarefy = NULL) 

shannon_df_pt <- shannon_df %>%  
  filter(project != "donor_batch")

## Wilcoxon test
shannon_wilcox_fmt <- shannon_df_pt %>%
  filter(group == "FMT") %>%
  run_paired_wilcox(Shannon) 

shannon_wilcox_placebo <- shannon_df_pt %>%
  filter(group == "placebo") %>%
  run_paired_wilcox(Shannon) 

## Plot
n_counts <- shannon_df_pt %>%
  filter(!is.na(Shannon)) %>% 
  group_by(group, stage) %>%
  summarise(n = n_distinct(id), .groups = "drop")

shannon_plot <- ggplot(shannon_df_pt, aes(x = stage, y = Shannon,  color = id)) +
  facet_grid(. ~ group, labeller = as_labeller(c("FMT" ="FMT", "placebo" = "Placebo"))) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 1.1) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "black", size = 3, shape = 18) + 
  geom_line(aes(group = id), alpha = 0.8, position = position_dodge(0.3), linewidth=1.0) +
  geom_point(aes(group = id), alpha = 0.9, size = 3,  position = position_dodge(0.3)) +
  scale_x_discrete(labels = c("Baseline", "0M", "1M", "3M", "6M", "12M")) + 
  scale_color_manual(name = NULL, values = col_map) + 
  labs(x = "", y = "Shannon diversity") + 
  scale_y_continuous(limits = c(0, NA)) + 
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
  expand_limits(y = 4.0) +
  geom_label(
    data = n_counts,
    aes(x = stage, y = 3.8, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 10/.pt,
    label.size = 0.2,
    fill="grey97",
    alpha = 0.8
  ) 

alpha_div <- richness_plot / shannon_plot +
  plot_layout(heights = c(1, 1))

#ggsave("./plots/figure_4.tiff", plot=alpha_div, dpi = 300, width=12, height=10)


# PCA
do_metadata <- metadata %>%
  filter(grepl("do", id)) %>%
  distinct(id, .keep_all = T)

pt_metadata <- metadata %>% 
  filter(project != "donor_batch", 
         !(stage %in% c("treatment_5", "treatment_10", "treatment_15", "treatment_21"))) 
  

metadata_pca <- full_join(do_metadata, pt_metadata) %>%  
  mutate(group_ordination = case_when(
    group == "FMT" ~ "FMT",
    group == "placebo" ~ "Placebo",
    grepl("do", id) ~ "Donor",
    TRUE ~ NA_character_)) %>%  
  mutate(timepoint = if_else(!is.na(stage), stage, "Donor"), 
         timepoint = factor(timepoint, levels = c("Donor", "inclusion", "followup_30d", "followup_1m", "followup_3m", "followup_6m", "followup_12m"))) %>%  
  left_join(shannon_df %>%  
              select(sample_barcode, id, Shannon) %>%  
              filter(!is.na(Shannon)), by=c("sample_barcode", "id")) %>%  
  mutate(facet_group = if_else(group_ordination %in% c("FMT", "Placebo"), "Patients", "FMT donors")) %>%  
  mutate(group_ordination_2 = if_else(group_ordination %in% c("FMT", "Donor"), "FMT", 
                                      if_else(group_ordination == "Placebo", "Placebo", NA_character_))) %>%  
  relocate(group_ordination_2, .after="group_ordination") 


amp_object_pca <- amp_load(otutable = otutable, metadata = metadata_pca, taxonomy = taxtable)

ordination_palette <- c(
  "grey50",
  "#E8DED0",
  "#D1C3A7",
  "#B7A27C",
  "#91B7C8",
  "#5C8596",
  "#1C454F"
)


pca_plot <- amp_object_pca %>%  
  amp_ordinate(
    type = "pca",
    transform = "hellinger", 
    filter_species = 0, 
    sample_color_by = "timepoint",
    species_nlabels = 5,
    sample_shape_by = "facet_group", 
    species_label_taxonomy = "Species",
    sample_point_size = 4,
    species_plot = T) +
  scale_color_manual(name = "",
                     labels = c("Donors", "Baseline", "0M", "1M", "3M", "6M", "12M"),
                     values = ordination_palette,
                     guide = guide_legend(order = 1)) +
  scale_shape_manual(name = "Group", 
                     values = c(17, 16), 
                     guide = guide_legend(order = 2)) +
  geom_path(aes(group = id), alpha = 0.5, color = "gray30") +
  geom_point(aes(color = timepoint, shape = facet_group), size = 4) +
  theme(panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"),
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        aspect.ratio=1,
        strip.background = element_rect(fill="#d9d9d9", color="#d9d9d9"),
        strip.text = element_text(color = "black", size = 11),
        axis.title.x = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 11),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.text  = element_text(size = 11)) + 
  facet_wrap(. ~ group_ordination_2)

pca_plot

#ggsave("./plots/figure_5.svg", plot=pca_plot, dpi = 300, w=14, h=10.10)


# PERMANOVA 
### Prepare data
abund_pca <- amp_object_pca$abund %>%
  t() %>%
  as.data.frame()

### Check that metadata and abundance tables are in same order
# meta_perm <- metadata_pca %>%
#   filter(sample_barcode %in% rownames(abund_pca)) %>%
#   distinct(sample_barcode, .keep_all = TRUE)
# 
# abund_pca <- abund_pca[
#   meta_perm$sample_barcode,
#   ,
#   drop = FALSE]
# 
# stopifnot(identical(
#   rownames(abund_pca),
#   meta_perm$sample_barcode))


### Hellinger transform abundances 
abund_hell <- decostand(abund_pca, method = "hellinger")

### Euclidean distance of Hellinger-transformed abundances
dist_hell <- dist(abund_hell)

### Filter for ***FMT***
meta_perm_fmt <- meta_perm %>%  
  filter(group == "FMT") %>%  
  droplevels()

fmt_samples <- meta_perm_fmt$sample_barcode

dist_hell_fmt <- as.dist(
  as.matrix(dist_hell)[fmt_samples, fmt_samples])

### Run PERMANOVA (FMT)
set.seed(123)
permanova_fmt <- adonis2(
  dist_hell_fmt ~ timepoint,
  data = meta_perm_fmt,
  permutations = 9999,
  strata = meta_perm_fmt$id
)

### PERMANOVA assumes equal dispersions - check for that: 
disp_fmt <- betadisper(
  dist_hell_fmt,
  meta_perm_fmt$timepoint
)

permutest(disp_fmt, permutations = 9999) 

### Filter for ***Placebo***
meta_perm_placebo <- meta_perm %>%  
  filter(group == "placebo") %>%  
  droplevels()

placebo_samples <- meta_perm_placebo$sample_barcode

dist_hell_placebo <- as.dist(
  as.matrix(dist_hell)[placebo_samples, placebo_samples])

### Run PERMANOVA (placebo)
set.seed(123)
permanova_placebo <- adonis2(
  dist_hell_placebo ~ timepoint,
  data = meta_perm_placebo,
  permutations = 9999,
  strata = meta_perm_placebo$id
)

### PERMANOVA assumes equal dispersions - check for that: 
disp_placebo <- betadisper(
  dist_hell_placebo,
  meta_perm_placebo$timepoint)

permutest(disp_placebo, permutations = 9999) 



# PCA - Supplementary figure 1
pca_plot_shannon_div <- amp_object_pca %>%  
  amp_ordinate(
    type = "pca",
    transform = "hellinger", 
    filter_species = 0, 
    sample_color_by = "timepoint",
    species_nlabels = 5,
    sample_shape_by = "facet_group", 
    species_label_taxonomy = "Species",
    sample_point_size = 4,
    species_plot = T) +
  scale_color_manual(name = "",
                     labels = c("Donors", "Baseline", "0M", "1M", "3M", "6M", "12M"),
                     values = ordination_palette,
                     guide = guide_legend(order = 1)) +
  scale_shape_manual(name = "Group", 
                     values = c(17, 16), 
                     guide = guide_legend(order = 2)) +
  geom_path(aes(group = id), alpha = 0.5, color = "gray30") +
  # geom_point(aes(color = timepoint, shape = facet_group, size = Shannon)) +
  geom_point(data = function(x) x %>% filter(facet_group == "Patients"), 
             aes(color = timepoint, shape = facet_group, size = Shannon)) +
  geom_point(data = function(x) x %>% filter(facet_group == "FMT donors"),
             aes(color = timepoint, shape = facet_group), size = 4) +
  scale_size(range = c(0.5, 8),
             name = "Shannon Diversity",
             guide = guide_legend(order = 3)) +
  theme(panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"),
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        aspect.ratio=1,
        strip.background = element_rect(fill="#d9d9d9", color="#d9d9d9"),
        strip.text = element_text(color = "black", size = 11),
        axis.title.x = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 11),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.text  = element_text(size = 11)) + 
  facet_wrap(. ~ group_ordination_2)


pca_plot_cpdai <- amp_object_pca %>%  
  amp_ordinate(
    type = "pca",
    transform = "hellinger", 
    filter_species = 0, 
    sample_color_by = "timepoint",
    species_nlabels = 5,
    sample_shape_by = "facet_group", 
    species_label_taxonomy = "Species",
    sample_point_size = 4,
    species_plot = T) +
  scale_color_manual(name = "",
                     labels = c("Donors", "Baseline", "0M", "1M", "3M", "6M", "12M"), 
                     values = ordination_palette, 
                     guide = guide_legend(order = 1)) +
  scale_shape_manual(name = "Group", 
                     values = c(17, 16), 
                     guide = guide_legend(order = 2)) +
  geom_path(aes(group = id), alpha = 0.5, color = "gray30") +
  #geom_point(aes(color = timepoint, shape = facet_group, size = cpdai_sum)) +
  geom_point(data = function(x) x %>% filter(facet_group == "Patients"),
             aes(color = timepoint, shape = facet_group, size = cpdai_sum)) +
  geom_point(data = function(x) x %>% filter(facet_group == "FMT donors"),
             aes(color = timepoint, shape = facet_group),
             size = 4) + 
  scale_size(range = c(0.5, 8), 
             name = "cPDAI",
             guide = guide_legend(order = 3)) +
  theme(panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"),
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        aspect.ratio=1,
        strip.background = element_rect(fill="#d9d9d9", color="#d9d9d9"),
        strip.text = element_text(color = "black", size = 11),
        axis.title.x = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 11),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.text  = element_text(size = 11)) + 
  facet_wrap(. ~ group_ordination_2) 


suppl_fig_1 <- pca_plot_shannon_div / pca_plot_cpdai +
  plot_annotation(tag_levels = "A", tag_suffix = ")")

#ggsave("./plots/suppl_figure_1.svg", plot=suppl_fig_1, dpi = 600, w = 18.65, h = 13.41)
 

# Stacked bar plot - Supplementary figure 2
metaphlan_select <- metaphlan %>%  
  filter(str_detect(clade_name, "s__Escherichia_coli|s__Ruminococcus_gnavus|s__Faecalibacterium_prausnitzii")) %>%  
  separate(clade_name, sep = "\\|", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "clade_name")) %>%  
  select(-Kingdom:-Genus)
  
metaphlan_remaining <- metaphlan %>%  
  filter(!str_detect(clade_name, "s__Escherichia_coli|s__Ruminococcus_gnavus|s__Faecalibacterium_prausnitzii")) %>%  
  summarise(
    clade_name = "sum_remaining",
    across(-clade_name, ~ sum(.x, na.rm = TRUE)))

metaphlan_mod <- full_join(metaphlan_select, metaphlan_remaining) %>% 
  pivot_longer(
    cols = -clade_name,
    names_to = "sample_barcode",
    values_to = "relative_abundance"
  ) %>%  
  left_join(metadata %>%
              select(sample_barcode, stage, group, id), by="sample_barcode") %>%  
  mutate(group = if_else(str_detect(id, "do"), "donor", group)) %>%  
  mutate(
    stage = factor(
      stage,
      levels = c(
        "inclusion",
        "treatment_5",
        "treatment_10",
        "treatment_15",
        "treatment_21",
        "followup_30d",
        "followup_1m",
        "followup_3m",
        "followup_6m",
        "followup_12m"))) %>% 
  mutate(id = if_else(str_detect(id, "pt"), str_replace(id, "pt", "Pt"), id))

fmt_plot <- ggplot(metaphlan_mod %>% filter(group == "FMT"), aes(x=stage, relative_abundance, fill=clade_name)) +  
  geom_col(width = 0.8) +  
  facet_grid(group ~ id) +
  labs(
    x="", 
    y="Relative abundance [%]"
  ) + 
  scale_fill_manual(
    name = "Taxa",
    values = c(
      "s__Escherichia_coli" = "#854442",
      "s__Ruminococcus_gnavus" = "#428385",
      "s__Faecalibacterium_prausnitzii" = "#a3b899",
      "sum_remaining" = "#b69e8b"), 
    labels = c(
      "s__Escherichia_coli" = expression(italic("Escherichia coli")),
      "s__Ruminococcus_gnavus" = expression(italic("Ruminococcus gnavus")),
      "s__Faecalibacterium_prausnitzii" = expression(italic("Faecalibacterium prausnitzii")),
      "sum_remaining" = "Remaining taxa")) +
  scale_x_discrete(
    labels = c(
      "inclusion" = "Baseline",
      "followup_30d" = "0M",
      "followup_1m" = "1M",
      "followup_3m" = "3M",
      "followup_6m" = "6M",
      "followup_12m" = "12M")) + 
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
        axis.text.x = element_text(color = "black", size = 11, angle=45, hjust=1),
        axis.text.y = element_text(color = "black", size = 10), 
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11))

placebo_plot <- ggplot(metaphlan_mod %>% filter(group == "placebo") %>% mutate(group = if_else(str_detect(group, "placebo"), "Placebo", group)), aes(x=stage, relative_abundance, fill=clade_name)) +  
  geom_col(width = 0.8) +  
  facet_grid(group ~ id) +
  labs(
    x="", 
    y="Relative abundance [%]") + 
  scale_fill_manual(
    name = "Taxa",
    values = c(
      "s__Escherichia_coli" = "#854442",
      "s__Ruminococcus_gnavus" = "#428385",
      "s__Faecalibacterium_prausnitzii" = "#a3b899",
      "sum_remaining" = "#b69e8b"), 
    labels = c(
      "s__Escherichia_coli" = expression(italic("Escherichia coli")),
      "s__Ruminococcus_gnavus" = expression(italic("Ruminococcus gnavus")),
      "s__Faecalibacterium_prausnitzii" = expression(italic("Faecalibacterium prausnitzii")),
      "sum_remaining" = "Remaining taxa")) +
  scale_x_discrete(
    labels = c(
      "inclusion" = "Baseline",
      "followup_30d" = "0M",
      "followup_1m" = "1M",
      "followup_3m" = "3M",
      "followup_6m" = "6M",
      "followup_12m" = "12M")) + 
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"), 
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
        axis.ticks.x = element_line(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "right", 
        axis.text.x = element_text(color = "black", size = 11, angle=45, hjust=1),
        axis.text.y = element_text(color = "black", size = 10), 
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 11), 
        strip.text.y = element_text(size = 11))

suppl_2 <- (fmt_plot / placebo_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

#ggsave("./plots/suppl_figure_2.tiff", plot=suppl_2, dpi = 300, w=20, h=6)


# Sørensen similarity 
do_sim_metadata <- metadata %>%
  filter(project == "donor_batch") %>%
  select(sample_barcode, id, fecal_donation_number) %>%
  distinct(sample_barcode, .keep_all = T) %>%
  rename(do_id = id,
         do = sample_barcode,
         fecal_donation_number_do = fecal_donation_number)

pt_sim_metadata <- metadata %>%  
  filter(project == "MP") %>%  
  distinct(sample_barcode, id, stage, .keep_all = T) %>%  
  filter(!(stage %in% c("treatment_5", "treatment_10", "treatment_15", "treatment_21")))

pt_treatment_metadata <- metadata %>%  
  filter(project == "MP") %>%  
  distinct(sample_barcode, id, stage, .keep_all = T) %>% 
  filter(stage %in% c("treatment_5", "treatment_10", "treatment_15", "treatment_21")) %>%  
  filter(group == "FMT")  

pt_treatment_donor <- pt_treatment_metadata %>%
  group_by(id) %>%
  mutate(treatment_do_ids = list(unique(donor[str_detect(stage, "treatment")]))) %>%
  ungroup() %>%  
  distinct(id, treatment_do_ids)

pt_treatment_fecal_don <- pt_treatment_metadata %>%  
  group_by(id) %>%
  mutate(donor_fecal_don = paste0(donor, "-", fecal_donation_number)) %>% 
  mutate(unique_do_fecal_donation = list(unique(donor_fecal_don[str_detect(stage, "treatment")]))) %>%  
  ungroup() %>%  
  distinct(id, unique_do_fecal_donation)


sorensen <- metaphlan %>%  
  filter(clade_name != "UNCLASSIFIED") %>% 
  column_to_rownames(var = "clade_name") %>%  
  t() %>%  
  vegdist(method = "bray", binary = TRUE) %>% 
  as.matrix() %>%  
  as.data.frame()


sorensen_df_pre <- sorensen %>% 
  rownames_to_column(var = "sample_barcode") %>% 
  pivot_longer(cols = -sample_barcode, names_to = "do", values_to = "sorensen_dissimilarity") %>%  
  filter(!str_detect(sample_barcode, "Do"), 
         !str_detect(do, "Pa"), 
         !str_detect(do, "PAPR")) %>%  
  mutate(sorensen_similarity = 1-sorensen_dissimilarity) %>%
  right_join(pt_sim_metadata, by = "sample_barcode") %>%  
  left_join(do_sim_metadata, by = "do")   

sorensen_df_fmt <- sorensen_df_pre %>%  
  filter(group == "FMT") %>% 
  left_join(pt_treatment_donor, by = "id") %>%
  rowwise() %>%
  filter(do_id %in% treatment_do_ids) %>%
  ungroup() %>% 
  select(-treatment_do_ids) %>%  
  mutate(do_fecal_do_n = paste0(do_id, "-", fecal_donation_number_do)) %>%  
  left_join(pt_treatment_fecal_don, by="id") %>%  
  rowwise() %>% 
  filter(do_fecal_do_n %in% unique_do_fecal_donation) %>% 
  ungroup() %>%  
  select(-unique_do_fecal_donation, -do_fecal_do_n) %>%  
  distinct(id, do_id, stage, .keep_all = T) %>%  
  group_by(id, stage) %>%  
  mutate(median_sorensen_similarity = median(sorensen_similarity)) %>% 
  ungroup() %>%  
  distinct(id, stage, .keep_all = T) 

sorensen_df_placebo <- sorensen_df_pre %>%  
  filter(group == "placebo") %>%  
  group_by(id, stage, do_id) %>%  
  mutate(mean_sorensen_similarity_donors = mean(sorensen_similarity)) %>%  
  ungroup() %>%  
  distinct(id, stage, do_id, .keep_all = T) %>% 
  group_by(id, stage) %>%  
  mutate(median_sorensen_similarity = median(mean_sorensen_similarity_donors)) %>% 
  ungroup() %>% 
  distinct(id, stage, .keep_all = T) %>%  
  filter(!is.na(median_sorensen_similarity))


sorensen_df_combined <- full_join(sorensen_df_placebo, sorensen_df_fmt) %>%  
  mutate(stage = factor(stage, levels= c("inclusion", "followup_30d", "followup_1m", "followup_3m", "followup_6m", "followup_12m")))

## Statistics
sorensen_wilcox_fmt <- sorensen_df_combined %>%
  filter(group == "FMT") %>%
  run_paired_wilcox(median_sorensen_similarity)

sorensen_wilcox_placebo <- sorensen_df_combined %>%
  filter(group == "placebo") %>%
  run_paired_wilcox(median_sorensen_similarity) 


## Plot
n_counts <- sorensen_df_combined %>%
  filter(!is.na(median_sorensen_similarity)) %>% 
  group_by(group, stage) %>%
  summarise(n = n_distinct(id), .groups = "drop")


sorensen_plot <- ggplot(sorensen_df_combined, aes(x = stage, y = median_sorensen_similarity, color = id)) +
  facet_grid(. ~ group, labeller = as_labeller(c("FMT" ="FMT", "placebo" = "Placebo"))) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 1.1) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "black", size = 3, shape = 18) + 
  geom_line(aes(group = id), alpha = 0.8, position = position_dodge(0.3), linewidth=1.0) +
  geom_point(aes(group = id), alpha = 0.9, size = 3,  position = position_dodge(0.3)) +
  scale_x_discrete(labels = c("Baseline", "0M", "1M", "3M", "6M", "12M")) +
  scale_color_manual(name = NULL, values = col_map) + 
  labs(x = "", y = "Sørensen similarity (median)") + 
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
  scale_y_continuous(
    limits = c(0, 0.50),
    breaks = c(0.0, 0.2, 0.4)) +
  geom_label(
    data = n_counts,
    aes(x = stage, y = 0.46, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3,
    label.size = 0.2,
    fill="grey97",
    alpha = 0.8
  )


# Bray-Curtis similarity
bc <- metaphlan %>%  
  filter(clade_name != "UNCLASSIFIED") %>% 
  column_to_rownames(var = "clade_name") %>%  
  t() %>% 
  decostand(method = "hellinger") %>% 
  vegdist(method = "bray", binary = FALSE) %>% 
  as.matrix() %>%  
  as.data.frame()

bc_df_pre <- bc %>% 
  rownames_to_column(var = "sample_barcode") %>% 
  pivot_longer(cols = -sample_barcode, names_to = "do", values_to = "bray_curtis_dissimilarity") %>%  
  filter(!str_detect(sample_barcode, "Do"), 
         !str_detect(do, "Pa"), 
         !str_detect(do, "PAPR")) %>%  
  mutate(bc_similarity = 1-bray_curtis_dissimilarity) %>%
  right_join(pt_sim_metadata, by = "sample_barcode") %>%  
  left_join(do_sim_metadata, by = "do")   

bc_df_fmt <- bc_df_pre %>%  
  filter(group == "FMT") %>% 
  left_join(pt_treatment_donor, by = "id") %>%
  rowwise() %>%
  filter(do_id %in% treatment_do_ids) %>%
  ungroup() %>% 
  select(-treatment_do_ids) %>%  
  mutate(do_fecal_do_n = paste0(do_id, "-", fecal_donation_number_do)) %>%  
  left_join(pt_treatment_fecal_don, by="id") %>%  
  rowwise() %>% 
  filter(do_fecal_do_n %in% unique_do_fecal_donation) %>% 
  ungroup() %>%  
  select(-unique_do_fecal_donation, -do_fecal_do_n) %>%  
  distinct(id, do_id, stage, .keep_all = T) %>%  
  group_by(id, stage) %>%  
  mutate(median_bc_similarity = median(bc_similarity)) %>% 
  ungroup() %>%  
  distinct(id, stage, .keep_all = T) 

bc_df_placebo <- bc_df_pre %>%  
  filter(group == "placebo") %>%  
  group_by(id, stage, do_id) %>%  
  mutate(mean_bc_similarity_donors = mean(bc_similarity)) %>%  
  ungroup() %>%  
  distinct(id, stage, do_id, .keep_all = T) %>% 
  group_by(id, stage) %>%  
  mutate(median_bc_similarity = median(mean_bc_similarity_donors)) %>% 
  ungroup() %>% 
  distinct(id, stage, .keep_all = T) %>%  
  filter(!is.na(median_bc_similarity))


bc_df_combined <- full_join(bc_df_placebo, bc_df_fmt) %>%  
  mutate(stage = factor(stage, levels= c("inclusion", "followup_30d", "followup_1m", "followup_3m", "followup_6m", "followup_12m")))

## Statistics
bc_wilcox_fmt <- bc_df_combined %>%
  filter(group == "FMT") %>%
  run_paired_wilcox(median_bc_similarity) 

bc_wilcox_placebo <- bc_df_combined %>%
  filter(group == "placebo") %>%
  run_paired_wilcox(median_bc_similarity) 


## Plot
n_counts <- bc_df_combined %>%
  filter(!is.na(median_bc_similarity)) %>% 
  group_by(group, stage) %>%
  summarise(n = n_distinct(id), .groups = "drop")

bc_plot <- ggplot(bc_df_combined, aes(x = stage, y = median_bc_similarity, color = id)) +
  facet_grid(. ~ group, labeller = as_labeller(c("FMT" ="FMT", "placebo" = "Placebo"))) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 1.1) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "black", size = 3, shape = 18) + 
  geom_line(aes(group = id), alpha = 0.8, position = position_dodge(0.3), linewidth=1.0) +
  geom_point(aes(group = id), alpha = 0.9, size = 3,  position = position_dodge(0.3)) +
  scale_x_discrete(labels = c("Baseline", "0M", "1M", "3M", "6M", "12M")) +
  scale_color_manual(name = NULL, values = col_map) + 
  labs(x = "", y = "Bray-Curtis similarity (median)") + 
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
  scale_y_continuous(
    limits = c(0, 0.50),
    breaks = c(0.0, 0.2, 0.4)) +
  geom_label(
    data = n_counts,
    aes(x = stage, y = 0.46, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3,
    label.size = 0.2,
    fill="grey97",
    alpha = 0.8
  )

similarity <- sorensen_plot / bc_plot +
  plot_layout(heights = c(1, 1))


#ggsave("./plots/figure_6.tiff", plot=similarity, dpi = 300, width=12, height=10)


library(ggplot2)
library(ggalluvial)
library(dplyr)
library(ggpubr)
library(readr)
library(tidyverse)
library(ggplotify)
library(colorspace)   # for qualitative palettes
library(viridisLite)
library(patchwork) 
library(cowplot)
library(ggplotify)
library(purrr)
library(tidyr)
library(tidytext)
library(tidyr)

p <- read_delim("~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
p$ST <- as.factor(p$ST)

#Function plot
make_CS_plot <- function(combination_name, csv_path, complete_column = "combination") {  
  # --- Read & prep main dataset ---
  CS_data <- read_delim(csv_path, delim = ";", escape_double = FALSE, trim_ws = TRUE)
  CS_data <- CS_data %>%
    filter(ST != "-") %>%          # <-- remove unwanted ST
    mutate(
      ST = as.factor(ST),
      pling = trimws(pling)
    )
  # --- Prepare clean mapping of pling -> Plasmid_Lineage (avoid many-to-many)
  p_unique <- p %>%
    mutate(pling = trimws(pling)) %>%
    group_by(pling) %>%
    summarise(Plasmid_Lineage = first(Plasmid_Lineage)) %>%
    ungroup()
  
  # --- Join and preserve "Different plasmid"
  CS_data <- CS_data %>%
    left_join(p_unique, by = "pling") %>%
    mutate(
      Plasmid_Lineage = if_else(
        pling == "Different plasmid",
        "Different plasmid",
        as.character(Plasmid_Lineage)
      )
    )
  
  # --- Prepare complete/public subsets ---
  p_char <- p %>%
    filter(ST != "-") %>%          # <-- same filter here
    mutate(Plasmid_Lineage = as.character(Plasmid_Lineage))
  complete_data <- p_char %>%
    filter(.data[[complete_column]] == combination_name)  
  public_data <- CS_data %>%
    group_by(ST) %>% mutate(count_ST = n()) %>% ungroup()
  # --- Union of used lineages ---
  used_lineages <- c(complete_data$Plasmid_Lineage, public_data$Plasmid_Lineage) %>%
    na.omit() %>% unique() %>% sort()
  
  if (length(used_lineages) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste(combination_name, "- no data after filtering")))
  }
  
  has_diff <- "Different plasmid" %in% used_lineages
  normal_lineages <- setdiff(used_lineages, "Different plasmid")
  
  # --- Build color palette (no reds except for "Different plasmid")
  n_norm <- length(normal_lineages)
  base_cols <- if (n_norm > 0) {
    hues <- seq(50, 310, length.out = n_norm)     # skip red region
    hcl(h = hues, c = 100, l = 60)
  } else character(0)
  
  palette_cols <- setNames(base_cols, normal_lineages)
  if (has_diff) palette_cols <- c(palette_cols, "Different plasmid" = "#E41A1C")
  
  # --- Compute which STs have >10 isolates ---
  complete_counts <- complete_data %>%
    count(ST, name = "n") %>%
    mutate(freq = n / sum(n)) %>%
    filter(freq > 0.05) %>%
    pull(ST)
  
  public_counts <- public_data %>%
    count(ST, name = "n") %>%
    mutate(freq = n / sum(n)) %>%
    filter(freq > 0.05) %>%
    pull(ST)
  
  # --- COMPLETE panel ---
  complete_plot <- suppressWarnings({
    ggplot(complete_data, aes(axis1 = ST, axis2 = Plasmid_Lineage, y = 1)) +
      geom_alluvium(aes(fill = Plasmid_Lineage), width = 1/12, show.legend = FALSE) +
      geom_stratum(width = 1/12, fill = "grey", color = "black", show.legend = FALSE) +
      geom_text(
        stat = "stratum",
        aes(label = ifelse(stratum %in% complete_counts, as.character(after_stat(stratum)), "")),
        nudge_x = 0.5, hjust = 1,fontface="bold"
      )+ylab("Number of isolates") +
      scale_fill_manual(values = palette_cols, na.translate = FALSE) +
      scale_x_discrete(limits = c("ST", "Plasmid_Lineage"), labels = c("ST", "PL"),expand = c(0.15, 0.15)) +
      labs(title = "COMPLETE") +
      theme_classic() +
      theme(
        text = element_text(face = "bold",size=16),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold"),
        legend.position = "none"
      )
  })
  
  # --- PUBLIC panel ---
  public_plot <- suppressWarnings({
    ggplot(public_data, aes(axis1 = ST, axis2 = Plasmid_Lineage, y = 1)) +
      geom_alluvium(aes(fill = Plasmid_Lineage), width = 1/12, show.legend = FALSE) +
      geom_stratum(width = 1/12, fill = "grey", color = "black", show.legend = FALSE) +
      geom_text(
        stat = "stratum",
        aes(label = ifelse(stratum %in% public_counts, as.character(after_stat(stratum)), "")),
        nudge_x = 0.5, hjust = 1,fontface="bold"
      )+ ylab("Number of isolates") +
      scale_fill_manual(values = palette_cols, na.translate = FALSE) +
      scale_x_discrete(limits = c("ST", "Plasmid_Lineage"), labels = c("ST", "PL"),expand = c(0.15, 0.15)) +
      labs(title = "PUBLIC") +
      theme_classic() +
      theme(
        text = element_text(face = "bold",size=16),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold"),
        legend.position = "none"
      )
  })
  
  # --- Legend setup ---
  legend_labels <- used_lineages
  legend_labels[legend_labels == "Different plasmid"] <- "Different"  # rename for legend display only
  
  legend_df <- data.frame(
    Plasmid_Lineage = factor(used_lineages, levels = used_lineages),
    x = 1, y = seq_along(used_lineages)
  )
  
  legend_plot <- ggplot(legend_df, aes(x, y, fill = Plasmid_Lineage)) +
    geom_point(shape = 21, color = "black", stroke = 0.3, size = 5) +
    scale_fill_manual(
      values = palette_cols,
      name = "Plasmid lineage",
      labels = legend_labels,
      na.translate = FALSE
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 9),
      legend.text  = element_text(face = "bold", size = 8),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      legend.margin = margin(l = -10, unit = "pt")   # <-- moves legend closer to strata
    )
  
  one_legend <- cowplot::get_legend(legend_plot)
  
  # --- Arrange panels and legend ---
  # Slightly tighter spacing: smaller rel_width for the legend
  panels <- cowplot::plot_grid(complete_plot, public_plot, ncol = 2, align = "hv")
  
  final  <- cowplot::plot_grid(panels, one_legend, ncol = 2, rel_widths = c(1, 0.20))  # legend snug on right
  
  # --- Add main title ---
  final <- ggpubr::annotate_figure(
    final,
    top = ggpubr::text_grob(combination_name, face = "bold", hjust = 0.5,size=20)
  )
  
  return(final)
}

#Visualisation
CS_list <- list(
  make_CS_plot("CS4", "CS4_combine.csv"),
  make_CS_plot("CS6+CS8", "CS6_CS8_combine.csv"),
  make_CS_plot("CS7", "CS7_combine.csv"),
  make_CS_plot("CS1", "CS1_combine.csv"),
  make_CS_plot("CS5+CS6", "CS5_CS6_combine.csv"),
  make_CS_plot("CS17", "CS17_combine.csv"),
  make_CS_plot("CS19","CS19_combine.csv"),
  make_CS_plot("CS26", "CS26_combine.csv"),
  make_CS_plot("CFA/I","CFA_I_combine.csv"),
  make_CS_plot("CS27", "CS27_combine.csv"),
  make_CS_plot("CS3", "CS3_combine.csv"),
  make_CS_plot("CS20/CS28", "CS20_CS28_combine.csv"),
  make_CS_plot("CS12", "CS12_combine.csv"),
  make_CS_plot("CS13", "CS13_combine.csv"),
  make_CS_plot("CS30", "CS30_combine.csv"),
  make_CS_plot("CS14","CS14_combine.csv"),
  make_CS_plot("CS23", "CS23_combine.csv"),
  make_CS_plot("CS21","CS21_combine.csv"),
  make_CS_plot("CS8", "CS8_combine.csv"),
  make_CS_plot("CS6","CS6_combine.csv")
)
Fig4A<- ggarrange(plotlist = CS_list)
NEW_CS_list <- list(
  make_CS_plot("New_gamma-1_1+CS6", "NEW_gamma_1_CS6.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_gamma-2_1+CS6", "New_gamma2_1_CS6.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_gamma-3_1+CS6+CS8", "NEW_CF_gamma3_1_CS6_CS8.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_gamma-3_2+CS6", "NEW_gamma_3_2_CS6.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_typeIV_pili+New_kappa-1+CS6", "New_typeiv_pili_new_kappa_1_CS6.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_gamma-2_1", "New_gamma2_1.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_gamma-3_2", "NEW_CF_gamma_3_2.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_pi_1", "NEW_pi1.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_typeIV_pili+CS4", "NEW_typeiv_pili_CS4.csv", complete_column = "New_combination_easy_vis"),
  make_CS_plot("New_gamma-2_2+CS13", "New_gamma_2_2_CS13.csv", complete_column = "New_combination_easy_vis"))
  
Fig4B <- ggarrange(plotlist = NEW_CS_list, ncol = 5, nrow = 2)

Figure4 <- ggarrange(
  Fig4A,
  Fig4B,
  ncol = 1,
  heights = c(2, 1),
  labels = c("A", "B"),
  font.label = list(size = 26, face = "bold")
)

ggsave("Fig4.svg", plot = Figure4, height = 22, width = 26,bg="white")
ggsave("Fig4.png", plot = Figure4, height = 22, width = 26,bg="white")


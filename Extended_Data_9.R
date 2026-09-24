library(ggplot2)
library(ggalluvial)
library(dplyr)
library(ggpubr)
library(readr)
library(tidyverse)
library(ggplotify)
library(colorspace)
library(viridisLite)
library(patchwork) 
library(cowplot)
library(readxl)
library(ggnewscale)
library(tidytext)

VF_colors <- c(
  "Non-ETEC" = "#BDBDBD",
  "ETEC" = "#D95F02"
)

p <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv", 
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

p$ST <- as.factor(p$ST)

all_information <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/metadata_ETEC_Validation_definitive.csv", 
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

make_CS_plot <- function(combination_name, csv_path) {
  
  CS_data <- read_delim(
    csv_path,
    delim = ";",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  
  CS_data$ST <- as.factor(CS_data$ST)
  CS_data <- CS_data %>%
    mutate(pling = trimws(pling))
  
  p_unique <- p %>%
    mutate(pling = trimws(pling)) %>%
    group_by(pling) %>%
    summarise(
      Plasmid_Lineage = first(Plasmid_Lineage)
    ) %>%
    ungroup()
  
  CS_data <- CS_data %>%
    left_join(
      all_information %>% select(id, ETEC),
      by = "id"
    ) %>%
    mutate(
      ETEC = as.factor(ETEC)
    )
  
  public_data <- CS_data %>%
    left_join(
      p_unique,
      by = "pling"
    ) %>%
    mutate(
      Plasmid_Lineage = if_else(
        pling == "Different plasmid",
        "Different plasmid",
        as.character(Plasmid_Lineage)
      )
    )
  
  if (nrow(public_data) == 0) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(paste(combination_name, "- no data"))
    )
  }
  
  used_etec <- public_data$ETEC %>%
    na.omit() %>%
    unique()
  
  if (length(used_etec) == 0) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(paste(combination_name, "- no ETEC data"))
    )
  }
  
  fill_values <- c(
    VF_colors,
    "Different plasmid" = "red",
    "Other" = "grey"
  )
  
  public_counts <- public_data %>%
    count(ST, name = "n") %>%
    mutate(freq = n / sum(n)) %>%
    filter(freq > 0.05) %>%
    pull(ST)
  
  public_plot <- suppressWarnings({
    ggplot(
      public_data,
      aes(
        axis1 = ST,
        axis2 = Plasmid_Lineage,
        y = 1
      )
    ) +
      geom_alluvium(
        aes(fill = ETEC),
        width = 1/12,
        show.legend = FALSE
      ) +
      geom_stratum(
        aes(
          fill = ifelse(
            after_stat(stratum) == "Different plasmid",
            "Different plasmid",
            "Other"
          )
        ),
        width = 1/12,
        color = "black",
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = fill_values,
        na.translate = FALSE
      ) +
      geom_text(
        stat = "stratum",
        aes(
          label = ifelse(
            stratum %in% public_counts,
            as.character(after_stat(stratum)),
            ""
          )
        ),
        nudge_x = 0.5,
        hjust = 1,
        fontface = "bold"
      ) +
      ylab("Number of isolates") +
      scale_x_discrete(
        limits = c("ST", "Plasmid_Lineage"),
        expand = c(0.15, 0.15)
      ) +
      labs(title = "PUBLIC") +
      theme_classic() +
      theme(
        text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold"),
        legend.position = "none"
      )
  })
  
  final <- public_plot
  
  final <- ggpubr::annotate_figure(
    final,
    top = ggpubr::text_grob(
      combination_name,
      face = "bold",
      hjust = 0.5,
      size = 20
    )
  )
  
  return(final)
}

CS_list <- list(
  make_CS_plot("CS4", "CS4_combine.csv"),
  make_CS_plot("CS6+CS8", "CS6_CS8_combine.csv"),
  make_CS_plot("CS7", "CS7_combine.csv"),
  make_CS_plot("CS1", "CS1_combine.csv"),
  make_CS_plot("CS5+CS6", "CS5_CS6_combine.csv"),
  make_CS_plot("CS17", "CS17_combine.csv"),
  make_CS_plot("CS19", "CS19_combine.csv"),
  make_CS_plot("CS26", "CS26_combine.csv"),
  make_CS_plot("CFA/I", "CFA_I_combine.csv"),
  make_CS_plot("CS27", "CS27_combine.csv"),
  make_CS_plot("CS3", "CS3_combine.csv"),
  make_CS_plot("CS20/CS28", "CS20_CS28_combine.csv"),
  make_CS_plot("CS12", "CS12_combine.csv"),
  make_CS_plot("CS13", "CS13_combine.csv"),
  make_CS_plot("CS30", "CS30_combine.csv"),
  make_CS_plot("CS14", "CS14_combine.csv"),
  make_CS_plot("CS23", "CS23_combine.csv"),
  make_CS_plot("CS21", "CS21_combine.csv"),
  make_CS_plot("CS8", "CS8_combine.csv"),
  make_CS_plot("CS6", "CS6_combine.csv")
)

Extended_Data_9_A <- ggarrange(
  plotlist = CS_list
)


ST_ETEC <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/metadata_ETEC_Validation_definitive.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(
    id = str_trim(as.character(id)),
    ST = str_trim(as.character(ST))
  )

df_long <- ST_ETEC %>%
  pivot_longer(
    cols = `CFA/I`:CS7,
    names_to = "CF",
    values_to = "presence"
  ) %>%
  filter(presence == 1) %>%
  mutate(
    CF = str_trim(as.character(CF)),
    ST = str_trim(as.character(ST)),
    ETEC = ETEC
  )

top10_ST_per_CF <- df_long %>%
  filter(
    ST != "-",
    CF != "CS2",
    CF != "CS18"
  ) %>%
  count(
    CF,
    ST,
    name = "n_total"
  ) %>%
  group_by(CF) %>%
  slice_max(
    order_by = n_total,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  left_join(
    df_long %>%
      count(
        CF,
        ST,
        ETEC,
        name = "n"
      ),
    by = c("CF", "ST")
  )

CF_levels <- c(
  "CS4", "CS6+CS8", "CS7", "CS1", "CS5+CS6",
  "CS17", "CS19", "CS26", "CFA/I", "CS27",
  "CS3", "CS20/CS28", "CS12", "CS13", "CS30",
  "CS14", "CS23", "CS21", "CS8", "CS6"
)

top10_ST_per_CF <- top10_ST_per_CF %>%
  mutate(
    CF = factor(
      CF,
      levels = CF_levels
    )
  )

Extended_Data_9_B <- top10_ST_per_CF %>%
  ggplot(
    aes(
      x = reorder_within(ST, n_total, CF),
      y = n,
      fill = ETEC
    )
  ) +
  geom_col() +
  coord_flip() +
  facet_wrap(
    ~ CF,
    scales = "free"
  ) +
  scale_x_reordered() +
  scale_fill_manual(
    values = VF_colors,
    name = "ETEC"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold")
  ) +
  xlab("ST") +
  ylab("Number of isolates")


legend <- cowplot::get_legend(
  Extended_Data_9_B +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
)

Extended_Data_9_B_noleg <- Extended_Data_9_B +
  theme(
    legend.position = "none"
  )

plots <- ggarrange(
  Extended_Data_9_A,
  Extended_Data_9_B_noleg,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  font.label = list(
    size = 20,
    face = "bold"
  )
)

Extended_Data_9 <- cowplot::plot_grid(
  plots,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.08)
)

ggsave(
  "Extended_Data_9.png",
  plot = Extended_Data_9,
  height = 20,
  width = 18,
  dpi = 300,
  bg = "white"
)
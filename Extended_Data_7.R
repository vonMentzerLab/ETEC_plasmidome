library(tidyverse)
library(igraph)
library(ggraph)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
library(patchwork)
library(ggtree)
library(ape)
library(tidyr)
library(ggsci)
library(forcats)
library(readxl)
library(cowplot)

set.seed(123)

# Panel A pling
edges <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/pling_out_CF_plasmids/all_plasmids_distances_CF.tsv",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

nodes <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  sub("\\.fasta$", "", x, ignore.case = TRUE)
}

edges <- edges %>%
  mutate(
    plasmid_1 = clean_id(plasmid_1),
    plasmid_2 = clean_id(plasmid_2)
  )

nodes <- nodes %>%
  mutate(name = clean_id(plasmid_id))

all_nodes <- tibble(
  name = unique(c(edges$plasmid_1, edges$plasmid_2))
)

vertices <- all_nodes %>%
  left_join(nodes, by = "name")

g <- graph_from_data_frame(
  d = edges,
  vertices = vertices,
  directed = FALSE
)

cf_palette <- c(
  "NEW_Gamma-2_1"                = "#E69F00",
  "NEW_Pi_1"                     = "#56B4E9",
  "NEW_Gamma-2_2"                = "#009E73",
  "NEW_Type_IV_pili+NEW_Kappa_1" = "#F0E442",
  "NEW_Gamma-3_2"                = "#0072B2",
  "NEW_Gamma-3_1"                = "#D55E00",
  "NEW_Gamma-1_1"                = "#CC79A7",
  "NEW_Type_IV_pili"             = "#000000"
)

layout_coords <- create_layout(g, layout = "fr")

FigureS14 <- ggraph(layout_coords) +
  geom_edge_link(
    alpha = 0.05,
    color = "grey60"
  ) +
  geom_node_point(
    aes(color = NEW_CF),
    data = layout_coords[
      !is.na(layout_coords$NEW_CF) &
        layout_coords$NEW_CF != "Known CF",
    ]
  ) +
  scale_color_manual(
    values = cf_palette,
    name = "New putative CFs",
    na.translate = FALSE
  ) +
  theme_void() +
  theme(
    text = element_text(face = "bold")
  )


# Panel B subcomunity
p <- nodes

sc <- p %>%
  filter(
    pling != "Hub plasmids",
    !is.na(NEW_CF),
    NEW_CF != "Known CF"
  ) %>%
  group_by(
    pling,
    NEW_CF
  ) %>%
  summarise(
    n = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    total_n = sum(n)
  ) %>%
  ungroup()

FigS15A <- sc %>%
  ggplot(
    aes(
      x = reorder(pling, total_n),
      y = n,
      fill = NEW_CF
    )
  ) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = cf_palette,
    name = "New putative CFs"
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    text = element_text(face = "bold")
  ) +
  xlab("Plasmid subcommunity") +
  ylab("Number of plasmids")

#Panel C Lineage
PL <- p %>%
  filter(
    Plasmid_distance_containment_singleton != "Singleton",
    !is.na(NEW_CF),
    NEW_CF != "Known CF"
  ) %>%
  group_by(
    Plasmid_Lineage,
    NEW_CF
  ) %>%
  summarise(
    n = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    total_n = sum(n)
  ) %>%
  ungroup()

FigS15B <- PL %>%
  ggplot(
    aes(
      x = reorder(Plasmid_Lineage, total_n),
      y = n,
      fill = NEW_CF
    )
  ) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = cf_palette,
    name = "New putative CFs"
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    text = element_text(face = "bold")
  ) +
  xlab("Plasmid lineage") +
  ylab("Number of plasmids")


# Panel D
links_df <- read_excel(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/contigs/links_plasmid_genome.xlsx"
)

genome_tree <- read.tree(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/ETEC_maximum_likelihood.nwk"
)

plasmid_tree <- read.tree(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/plasmid_tree_NJ_containment_distance.nwk"
)

clean <- function(x) {
  x %>%
    sub(
      "\\.f(ast)?a$|\\.fna$|\\.fa$|\\.fasta$",
      "",
      .,
      ignore.case = TRUE
    ) %>%
    sub("_fasta$", "", .) %>%
    sub("\\s+$", "", .)
}

genome_tree$tip.label <- clean(genome_tree$tip.label)
plasmid_tree$tip.label <- clean(plasmid_tree$tip.label)
links_df$plasmid <- clean(links_df$plasmid)
links_df$genome <- clean(links_df$genome)

p_genome <- ggtree(genome_tree) +
  theme_tree()

g_coords <- p_genome$data %>%
  filter(isTip) %>%
  select(label, x, y)

colnames(g_coords) <- c(
  "genome",
  "x_genome",
  "y_genome"
)

p_plasmid <- ggtree(plasmid_tree) +
  theme_tree()

offset <- max(p_genome$data$x) + 2

p_plasmid$data$x <- offset - p_plasmid$data$x

p_coords <- p_plasmid$data %>%
  filter(isTip) %>%
  select(label, x, y)

colnames(p_coords) <- c(
  "plasmid",
  "x_plasmid",
  "y_plasmid"
)

link_coords <- links_df %>%
  left_join(
    g_coords,
    by = "genome"
  ) %>%
  left_join(
    p_coords,
    by = "plasmid"
  ) %>%
  filter(
    !is.na(x_genome),
    !is.na(y_genome),
    !is.na(x_plasmid),
    !is.na(y_plasmid),
    !(x_genome == x_plasmid & y_genome == y_plasmid)
  )

facet_order <- c(
  "New_gamma-1_1+CS6",
  "New_gamma-2_1",
  "New_gamma-2_1+CS6",
  "New_gamma-2_2+CS13",
  "New_gamma-3_1+CS6+CS8",
  "New_gamma-3_2",
  "New_gamma-3_2+CS6",
  "New_pi_1",
  "New_typeIV_pili+CS4",
  "New_typeIV_pili+New_kappa-1+CS6"
)

link_coords_filtered <- link_coords %>%
  count(NEW_CF) %>%
  filter(n >= 4) %>%
  filter(
    NEW_CF != "SAME",
    !is.na(NEW_CF)
  ) %>%
  inner_join(
    link_coords,
    by = "NEW_CF"
  ) %>%
  filter(
    NEW_CF %in% facet_order
  ) %>%
  mutate(
    NEW_CF = factor(
      NEW_CF,
      levels = facet_order
    )
  )

FigureS16 <- ggplot() +
  geom_tree(
    data = p_genome$data,
    aes(
      x = x,
      y = y
    )
  ) +
  geom_tree(
    data = p_plasmid$data,
    aes(
      x = x,
      y = y
    )
  ) +
  geom_curve(
    data = link_coords_filtered,
    aes(
      x = x_genome,
      y = y_genome,
      xend = x_plasmid,
      yend = y_plasmid,
      color = as.factor(Mobility)
    ),
    curvature = 0,
    alpha = 0.5,
    linewidth = 0.25
  ) +
  scale_color_manual(
    values = c(
      "#E83F6F",
      "#FFBF00",
      "#2274A5"
    )
  ) +
  theme_void() +
  labs(
    colour = "Plasmid mobility"
  ) +
  theme(
    text = element_text(face = "bold")
  ) +
  facet_wrap(
    ~ NEW_CF,
    nrow = 2,
    ncol = 5,
    drop = FALSE
  )+theme(legend.position = "bottom")


# Shared CF legend
cf_legend <- cowplot::get_legend(
  FigS15A +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
)

FigureS14_noleg <- FigureS14 +
  theme(
    legend.position = "none"
  )

FigS15A_noleg <- FigS15A +
  theme(
    legend.position = "none"
  )

FigS15B_noleg <- FigS15B +
  theme(
    legend.position = "none"
  )


# A, B, C
# A
A_row <- cowplot::plot_grid(
  FigureS14_noleg,
  cf_legend,
  ncol = 2,
  rel_widths = c(1, 0.16)
)

A_row <- ggarrange(
  A_row,
  labels = "A",
  font.label = list(
    size = 20,
    face = "bold"
  )
)

# B and C
BC_row <- ggarrange(
  FigS15A_noleg,
  FigS15B_noleg,
  labels = c("B", "C"),
  ncol = 2,
  nrow = 1,
  font.label = list(
    size = 20,
    face = "bold"
  )
)

# Extended Data 7
Extended_Data_7 <- ggarrange(
  A_row,
  BC_row,
  FigureS16,
  labels = c("", "", "D"),
  ncol = 1,
  nrow = 3,
  heights = c(1, 1, 1.2),
  font.label = list(
    size = 20,
    face = "bold"
  )
)

ggsave(
  "Extended_Data_7.png",
  plot = Extended_Data_7,
  width = 12,
  height = 16,
  dpi = 300,
  bg = "white"
)
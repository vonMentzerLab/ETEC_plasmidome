## -------------------------
## 0) Libraries
## -------------------------
library(tidyverse)
library(igraph)
library(ggraph)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(pheatmap)

set.seed(123)


edges_file <- "~/Desktop/Project_plasmidome/Assemblies_good_quality/pling_out_CF_plasmids/all_plasmids_distances_CF.tsv"
nodes_file <- "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv"
distances_file <- "~/Downloads/containment_distance_CF.csv"

edges <- read_delim(
  edges_file,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE
)

nodes_raw <- read_delim(
  nodes_file,
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE
)

distances <- read_delim(
  distances_file,
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE
)

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  sub("\\.fasta$", "", x, ignore.case = TRUE)
}

first_existing_col <- function(df, candidates) {
  out <- candidates[candidates %in% names(df)][1]
  if (is.na(out)) NA_character_ else out
}

node_id_col <- first_existing_col(nodes_raw, c("plasmid_id", "name", "plasmid", "id"))
if (is.na(node_id_col)) {
  stop("Could not find plasmid ID column in nodes file.")
}

comb_col <- first_existing_col(nodes_raw, c("combination", "CF_plasmid", "CF", "Colonisation_factor"))
if (is.na(comb_col)) {
  stop("Could not find CF/combination column in nodes file.")
}

lineage_col <- first_existing_col(nodes_raw, c("Plasmid_Lineage", "plasmid_lineage", "Lineage"))
if (is.na(lineage_col)) {
  stop("Could not find Plasmid_Lineage column in nodes file.")
}

dist_col_edges <- first_existing_col(edges, c("distance", "dist", "containment_distance", "mean_distance", "weight"))
if (is.na(dist_col_edges)) {
  stop("Could not find pairwise distance column in edges file.")
}

dist_col_distances <- first_existing_col(distances, c("distance", "dist", "containment_distance", "mean_distance", "weight"))
if (is.na(dist_col_distances)) {
  stop("Could not find pairwise distance column in distances file.")
}

if (!all(c("plasmid_1", "plasmid_2") %in% names(edges))) {
  stop("Edges file must contain columns named plasmid_1 and plasmid_2.")
}

if (!all(c("plasmid_1", "plasmid_2") %in% names(distances))) {
  stop("Distances file must contain columns named plasmid_1 and plasmid_2.")
}

edges <- edges %>%
  mutate(
    plasmid_1 = clean_id(plasmid_1),
    plasmid_2 = clean_id(plasmid_2),
    distance  = as.numeric(.data[[dist_col_edges]])
  ) %>%
  select(plasmid_1, plasmid_2, distance, everything())

nodes <- nodes_raw %>%
  mutate(
    plasmid_id      = clean_id(.data[[node_id_col]]),
    combination     = as.character(.data[[comb_col]]),
    Plasmid_Lineage = as.character(.data[[lineage_col]])
  ) %>%
  distinct(plasmid_id, .keep_all = TRUE)

distances <- distances %>%
  mutate(
    plasmid_1 = clean_id(plasmid_1),
    plasmid_2 = clean_id(plasmid_2),
    distance  = as.numeric(.data[[dist_col_distances]])
  ) %>%
  select(plasmid_1, plasmid_2, distance, everything())

p <- nodes

#Fig1A
all_nodes <- tibble(name = unique(c(edges$plasmid_1, edges$plasmid_2)))

vertices <- all_nodes %>%
  left_join(
    nodes %>% rename(name = plasmid_id),
    by = "name"
  )

g <- graph_from_data_frame(
  d = edges %>% select(from = plasmid_1, to = plasmid_2, distance),
  vertices = vertices,
  directed = FALSE
)

no_beige_palette <- c(
  "#440154", "#482878", "#3E4A89", "#31688E", "#26828E", "#1F9E89", "#35B779", "#6DCD59",
  "#000004", "#1B0B42", "#4B0C6B", "#781C6D", "#A52C60", "#CF4446", "#ED6925", "#FB9A06",
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
  "black", "#BCBD22"
)

layout_coords <- create_layout(g, layout = "fr")

Figure1A <- ggraph(layout_coords) +
  geom_edge_link(alpha = 0.05, color = "grey60") +
  geom_node_point(
    aes(color = combination),
    data = layout_coords[!is.na(layout_coords$combination), ]
  ) +
  scale_color_manual(values = no_beige_palette, na.translate = FALSE) +
  theme_void() +
  labs(color = "Colonisation factors") +
  theme(
    text = element_text(face = "bold"),
    plot.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 9, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.spacing.y = unit(0.12, "cm")
  )

#Figure1B
bad_ids <- c(
  "50497_1#116_5_length=10785_depth=0.52x",
  "GCA_001264175.1_ASM126417v1_genomic_LFZH01000001.1"
)

distances_clean <- distances %>%
  filter(!(plasmid_1 %in% bad_ids | plasmid_2 %in% bad_ids))

dist_meta <- distances_clean %>%
  left_join(
    nodes %>% select(plasmid_id, combination, Plasmid_Lineage),
    by = c("plasmid_1" = "plasmid_id")
  ) %>%
  rename(
    combination1 = combination,
    Plasmid_Lineage1 = Plasmid_Lineage
  ) %>%
  left_join(
    nodes %>% select(plasmid_id, combination, Plasmid_Lineage),
    by = c("plasmid_2" = "plasmid_id")
  ) %>%
  rename(
    combination2 = combination,
    Plasmid_Lineage2 = Plasmid_Lineage
  )

combi_mean <- dist_meta %>%
  group_by(combination1, combination2) %>%
  summarise(mean_distance = mean(distance, na.rm = TRUE), .groups = "drop")

combi_mean <- combi_mean %>%
  filter(!is.na(mean_distance))

combi_wide <- combi_mean %>%
  filter(!is.na(combination2)) %>%
  pivot_wider(names_from = combination2, values_from = mean_distance)

m1 <- as.matrix(combi_wide[, -1])
rownames(m1) <- combi_wide$combination1

Fig1B <- pheatmap(
  m1,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  border_color = "grey70",
  color = colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdBu"))(50),
  number_color = "black"
)

Fig1B_gg <- as.ggplot(Fig1B)

#Fig1C
cf_counts <- p %>%
  filter(!is.na(combination), !is.na(Plasmid_Lineage)) %>%
  group_by(combination) %>%
  summarise(
    n_samples  = n(),
    n_lineages = n_distinct(Plasmid_Lineage),
    max_lin_n  = max(as.integer(table(Plasmid_Lineage))),
    .groups = "drop"
  ) %>%
  mutate(max_lin_perc = 100 * max_lin_n / n_samples)

cf_mean_lineage_dist <- dist_meta %>%
  filter(
    !is.na(combination1), !is.na(combination2),
    !is.na(Plasmid_Lineage1), !is.na(Plasmid_Lineage2),
    !is.na(distance),
    combination1 == combination2,
    Plasmid_Lineage1 != Plasmid_Lineage2
  ) %>%
  group_by(combination = combination1) %>%
  summarise(
    mean_lineage_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- cf_counts %>%
  left_join(cf_mean_lineage_dist, by = "combination") %>%
  mutate(
    mean_lineage_distance = ifelse(is.na(mean_lineage_distance), 0, mean_lineage_distance)
  )

Fig1C <- plot_df %>%
  filter(
    n_samples >= 5,
    combination != "Only Toxin",
    combination != "Only Non-Canonical"
  ) %>%
  ggplot(aes(
    x = n_lineages,
    y = mean_lineage_distance,
    size = n_samples
  )) +
  geom_point(
    aes(fill = max_lin_perc),
    alpha = 0.6,
    shape = 21
  ) +
  geom_text_repel(
    aes(label = combination),
    size = 3,
    color = "black",
    fontface = "bold",
    box.padding = 0.7,
    point.padding = 0.4,
    force = 2.5,
    max.overlaps = Inf,
    max.iter = 20000,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.size = 0.3
  ) +
  scale_fill_gradientn(
    colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
    limits = c(0, 100),
    name = "% in most common lineage"
  ) +
  scale_size_area(
    max_size = 20,
    name = "Number of plasmids"
  ) +
  scale_x_continuous(
    breaks = 1:max(plot_df$n_lineages, na.rm = TRUE)
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  xlab("Number of plasmid lineages") +
  ylab("Mean distance between plasmid lineages")
## -------------------------
## 12) Assemble final figure
## -------------------------
blank <- ggplot() + theme_void()

TopA_centered <- ggarrange(
  blank,
  Figure1A,
  blank,
  ncol = 3,
  widths = c(0.25, 3.5, 0.25),
  labels = c("", "A", ""),
  font.label = list(size = 20, face = "bold"),
  label.x = 0, label.y = 1
)

Bottom_BC <- ggarrange(
  Fig1B_gg,
  Fig1C,
  ncol = 2,
  widths = c(1.25, 1),
  align = "hv",
  labels = c("B", "C"),
  font.label = list(size = 20, face = "bold"),
  label.x = 0, label.y = 1
)

Figure1_final <- ggarrange(
  TopA_centered,
  Bottom_BC,
  ncol = 1,
  heights = c(1.05, 1.00)
)

ggsave(
  filename = "Figure1_final_centeredA.png",
  plot = Figure1_final,
  width = 16,
  height = 16,
  units = "in",
  dpi = 300,bg="white"
)
#Modided later in Inkscape

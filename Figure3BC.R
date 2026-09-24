#Figure 3BC
library(dplyr)
library(tidyr)
library(pheatmap)
library(readr)
library(igraph)
library(ggraph)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggraph)
library(readr)
set.seed(123)
#Load datasets
library(dplyr)
library(readr)
library(igraph)
library(ggpubr)
library(ggplotify)  
library(pheatmap)
library(viridisLite)

p <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv",
  delim = ";", escape_double = FALSE, trim_ws = TRUE
)
p$Other_CF[p$Other_CF=="None"]<-"Alone"
m1 <- p %>%
  filter(!is.na(Other_CF), Other_CF != "",
         !is.na(NEW_CF), NEW_CF != "NA") %>%
  count(NEW_CF, Other_CF) %>%
  group_by(NEW_CF) %>%
  mutate(ratio = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = Other_CF, values_from = ratio, values_fill = 0) %>%
  as.data.frame()

rownames(m1) <- m1$NEW_CF
m1$NEW_CF <- NULL

m1 <- as.matrix(m1)
storage.mode(m1) <- "numeric"
m1 <- m1[, c("Alone", setdiff(colnames(m1), "Alone"))]

Fig3B <- pheatmap(
  m1,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = round(m1, 2),
  border_color = "grey70",
  color = colorRampPalette(c("white", "red"))(50),
  breaks = seq(0, 1, length.out = 51),
  number_color = "black",
  fontsize_row = 15,
  fontsize_col = 15,
  fontsize_number = 12,
  angle_col = 90,
  silent = TRUE
)

make_bold <- function(g) {
  if (!is.null(g$gp)) {
    g$gp$fontface <- "bold"
  }
  if (!is.null(g$children)) {
    g$children <- lapply(g$children, make_bold)
  }
  g
}

Fig3B$gtable$grobs <- lapply(Fig3B$gtable$grobs, make_bold)

grid.newpage()
grid.draw(Fig3B$gtable)


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

#Fig 3C
edges <-read_delim("~/Desktop/Project_plasmidome/Assemblies_good_quality/pling_out_CF_plasmids/all_plasmids_distances_CF.tsv", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)
nodes <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv",
  delim = ";", escape_double = FALSE, trim_ws = TRUE
)

# Consistent ID cleaning for BOTH tables
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

# Build node list from edges (all unique IDs that appear in edges)
all_nodes <- tibble(name = unique(c(edges$plasmid_1, edges$plasmid_2)))

# LEFT join so we keep every vertex referenced by edges (even if metadata missing)
vertices <- all_nodes %>%
  left_join(nodes, by = "name")

set.seed(123)
g <- graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)


# Calculate layout once and save it

layout_coords <- create_layout(g, layout = "fr")
cf_palette <- c(
  "CS6 Negative"                 = "#7F7F7F",  # grey
  "CS6 only"                     = "#000000",  # black
  "NEW_Type_IV_pili+NEW_Kappa_1" = "#E69F00",  # orange
  "CS5"                          = "#56B4E9",  # sky blue
  "NEW_Gamma-3_2"                = "#009E73",  # bluish green
  "NEW_Gamma-1_1"                = "#F0E442",  # yellow
  "NEW_Gamma-2_1"                = "#0072B2",  # blue
  "CS8+NEW_Gamma-3_1"            = "#D55E00"   # vermillion
)
#Plot
Figure_NEW_CF <- ggraph(layout_coords) +
  geom_edge_link(alpha=0.05, color="grey60") +
  geom_node_point(aes(color = Other_CF), 
                  data = layout_coords[!is.na(layout_coords$Other_CF), ]) +  theme_void() +
  labs(color = "Other CF") +
  scale_color_manual(values = cf_palette)+
  theme(text = element_text(face="bold"))
ggsave("Figure3C.png",plot=Figure_NEW_CF,dpi=300,bg="white",width = 7.5,height = 6)
#Modified in Inkscape

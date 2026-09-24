library(data.table)
library(ape)
library(ggplot2)

tree <- read.tree(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/ETEC_maximum_likelihood.nwk"
)

meta <- fread(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv",
  sep = ";",
  header = TRUE,
  showProgress = FALSE
)

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\.fasta$", "", x, ignore.case = TRUE)
  x <- gsub("^`|`$", "", x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("#", "_", x)
  x <- gsub("\\.", "_", x)
  x
}

median_layer <- function() {
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 2.5
  )
}

pairwise_distances_from_idx <- function(idx, D, max_tree_dist) {
  k <- length(idx)
  if (k < 2L) return(numeric(0))
  
  subm <- D[idx, idx, drop = FALSE]
  vals <- subm[upper.tri(subm, diag = FALSE)]
  vals / max_tree_dist
}

build_box_dt <- function(group_dt, group_col, D, max_tree_dist) {
  group_col <- as.character(group_col)
  
  res <- vector("list", nrow(group_dt))
  
  for (i in seq_len(nrow(group_dt))) {
    idx <- group_dt$idx[[i]]
    k <- length(idx)
    
    if (k < 2L) {
      res[[i]] <- data.table(
        combination = group_dt$combination[i],
        group_value = group_dt[[group_col]][i],
        distance = 0,
        n = group_dt$n[i],
        n_genomes = k,
        Distribution = "Singleton"
      )
    } else {
      vals <- pairwise_distances_from_idx(idx, D, max_tree_dist)
      max_breadth <- max(vals)
      dist_lab <- if (max_breadth < 0.1) "Restricted" else "Distributed"
      
      res[[i]] <- data.table(
        combination = group_dt$combination[i],
        group_value = group_dt[[group_col]][i],
        distance = vals,
        n = group_dt$n[i],
        n_genomes = k,
        Distribution = dist_lab
      )
    }
  }
  
  out <- rbindlist(res, use.names = TRUE)
  setnames(out, "group_value", group_col)
  out[]
}


tree$tip.label <- clean_id(tree$tip.label)
tree$node.label <- NULL
tree <- keep.tip(tree, unique(tree$tip.label))

D <- cophenetic(tree)
storage.mode(D) <- "double"

max_tree_dist <- max(D, na.rm = TRUE)

tree_ids <- rownames(D)
idx_map <- setNames(seq_along(tree_ids), tree_ids)

rm(tree)
gc()

cat("Tips in tree:", length(tree_ids), "\n")
cat("Distance matrix size:", format(object.size(D), units = "auto"), "\n")
cat("Max tree distance:", max_tree_dist, "\n")

meta[, Genome := clean_id(Genome)]

meta_pairs <- unique(
  meta[
    !is.na(Genome) & Genome %in% tree_ids,
    .(Genome, combination, pling, Plasmid_Lineage)
  ]
)

facet_order <- c(
  "CS4","CS6+CS8","CS7","CS1","CS5+CS6",
  "CS17","CS19","CS26","CFA/I","CS27",
  "CS3","CS20/CS28","CS12","CS13","CS30",
  "CS14","CS23","CS21","CS8","CS6"
)

meta_pairs[, combination := factor(combination, levels = facet_order)]
meta_pairs[, idx := unname(idx_map[Genome])]


group_pling <- meta_pairs[
  !is.na(pling) &
    pling != "" &
    pling != "Hub plasmids" &
    combination != "Only toxin" &
    combination != "Only New" &
    combination != "PCFO71",
  .(
    idx = list(sort(unique(idx))),
    n = .N
  ),
  by = .(combination, pling)
]

group_lineage <- meta_pairs[
  !is.na(Plasmid_Lineage) &
    Plasmid_Lineage != "" &
    Plasmid_Lineage != "Hub plasmids" &
    combination != "Only toxin" &
    combination != "Only New" &
    combination != "PCFO71",
  .(
    idx = list(sort(unique(idx))),
    n = .N
  ),
  by = .(combination, Plasmid_Lineage)
]


df_box_pling <- build_box_dt(
  group_dt = group_pling,
  group_col = "pling",
  D = D,
  max_tree_dist = max_tree_dist
)

df_box_lineage <- build_box_dt(
  group_dt = group_lineage,
  group_col = "Plasmid_Lineage",
  D = D,
  max_tree_dist = max_tree_dist
)

df_box_pling[, combination := factor(as.character(combination), levels = facet_order)]
df_box_lineage[, combination := factor(as.character(combination), levels = facet_order)]

#Visualize
# Colorblind-friendly palette
distribution_colors <- c(
  "Singleton" = "#999999",
  "Distributed" = "#E69F00",
  "Restricted" = "#0072B2"
)

# Visualize
A_Lineage <- ggplot(
  df_box_lineage,
  aes(y = Plasmid_Lineage, x = distance, fill = Distribution)
) +
  geom_boxplot() +
  facet_wrap(~ combination, scales = "free_y") +
  median_layer() +
  scale_fill_manual(values = distribution_colors) +
  theme_classic() +
  theme(
    text = element_text(face = "bold")
  ) +
  labs(
    x = "Normalized phylogenetic breadth",
    y = "Plasmid lineage",
    fill = "Distribution"
  )

B_pling <- ggplot(
  df_box_pling,
  aes(y = pling, x = distance, fill = Distribution)
) +
  geom_boxplot() +
  facet_wrap(~ combination, scales = "free_y") +
  median_layer() +
  scale_fill_manual(values = distribution_colors) +
  theme_classic() +
  theme(
    text = element_text(face = "bold")
  ) +
  labs(
    x = "Normalized phylogenetic breadth",
    y = "Pling",
    fill = "Distribution"
  )

# Get shared legend
legend <- cowplot::get_legend(
 A_Lineage + theme(legend.position = "right")
)

# Remove legends
A <- A_Lineage + theme(legend.position = "none")
B <-B_pling+ theme(legend.position = "none")

# Align plotting areas
aligned_plots <- cowplot::align_plots(
  A,
  B,
  align = "v",
  axis = "lr"
)

# Combine A and B
plots <- cowplot::plot_grid(
  aligned_plots[[1]],
  aligned_plots[[2]],
  ncol = 1,
  labels = c("A", "B"),
  rel_heights = c(1, 1)
)

# Add shared legend
Extended_Data_6 <- cowplot::plot_grid(
  plots,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.10)
)

ggsave(
  "Extended_Data_6.png",
  plot = Extended_Data_6,
  width = 14,
  height = 17,
  dpi=300,
  
  bg = "white"
)
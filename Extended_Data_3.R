library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

# Load data
p <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv",
  delim = ";", escape_double = FALSE, trim_ws = TRUE
)

# Colour palette
no_beige_palette <- c(
  "#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
  "#1F9E89", "#35B779", "#6DCD59",
  "#000004", "#1B0B42", "#4B0C6B", "#781C6D", "#A52C60",
  "#CF4446", "#ED6925", "#FB9A06",
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "black", "#BCBD22"
)

no_beige_palette <- c(
  "#440154", "#482878", "#3E4A89", "#31688E", "#26828E", "#1F9E89", "#35B779", "#6DCD59",
  "#000004", "#1B0B42", "#4B0C6B", "#781C6D", "#A52C60", "#CF4446", "#ED6925", "#FB9A06",
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
  "black", "#BCBD22"
)

combination_levels <- sort(unique(na.omit(p$combination)))

CF_colours <- setNames(
  no_beige_palette[seq_along(combination_levels)],
  combination_levels
)
# A - Pling subcommunities
sc <- p %>%
  filter(pling != "Hub plasmids") %>%
  group_by(pling, combination) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pling) %>%
  mutate(total_n = sum(n)) %>%
  filter(total_n >= 5) %>%
  ungroup()

Extended_data_3A <- sc %>%
  ggplot(aes(x = reorder(pling, total_n), y = n, fill = combination)) +
  geom_col() +
  scale_fill_manual(
    values = CF_colours,
    name = "Colonisation factor"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    x = "Plasmid subcommunity",
    y = "Number of plasmids"
  ) +
  theme(text = element_text(face = "bold"))


# B - Plasmid lineages
PL <- p %>%
  filter(Plasmid_distance_containment_singleton != "Singleton") %>%
  group_by(Plasmid_Lineage, combination) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Plasmid_Lineage) %>%
  mutate(total_n = sum(n)) %>%
  filter(total_n >= 5) %>%
  ungroup()

Extended_data_3B <- PL %>%
  ggplot(aes(
    x = reorder(Plasmid_Lineage, total_n),
    y = n,
    fill = combination
  )) +
  geom_col() +
  scale_fill_manual(
    values = CF_colours,
    name = "Colonisation factor"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    x = "Plasmid lineage",
    y = "Number of plasmids"
  ) +
  theme(text = element_text(face = "bold"))


# Get legend
legend <- get_legend(
  Extended_data_3A +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
)

# Remove legends from plots
A <- Extended_data_3A + theme(legend.position = "none")
B <- Extended_data_3B + theme(legend.position = "none")

# Combine A and B
plots <- plot_grid(
  A, B,
  labels = c("A", "B"),
  ncol = 2,
  align = "h",
  axis = "tb"
)

# Add shared legend
Extended_data_3 <- plot_grid(
  plots,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.22)
)

# Save
ggsave(
  "Extended_data_3.png",
  Extended_data_3,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)

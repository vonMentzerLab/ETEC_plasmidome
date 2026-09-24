# Extended Data Figure 2

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(cowplot)
library(rsvg)
library(grid)
library(igraph)
library(ggraph)

options(scipen = 999)

#Load data

d <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/all_plasmid_information_DEF.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

edges <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/pling_out_all_plasmids/all_plasmids_distances.tsv",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
d <- d %>%
  mutate(
    VF = case_when(
      ETEC == "ETEC" ~ "Virulence factor",
      ETEC == "Non-ETEC" ~ "Non Virulence factor",
      TRUE ~ NA_character_
    )
  )

d$VF <- factor(
  d$VF,
  levels = c("Non Virulence factor", "Virulence factor")
)

VF_colors <- c(
  "Non Virulence factor" = "#BDBDBD",
  "Virulence factor" = "#D95F02"
)

# Panel A: Plasmid Size

Extended_data_2A <- ggplot(
  d,
  aes(x = size, fill = VF)
) +
  geom_histogram(
    color = "black",
    linewidth = 0.3,
    bins = 30,
    position = "stack"
  ) +
  scale_fill_manual(values = VF_colors) +
  scale_x_continuous(
    breaks = seq(0, 300000, by = 25000)
  ) +
  labs(
    x = "Plasmid size (bp)",
    y = "Number of plasmids",
    fill = "Virulence factor"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# Panel B: MObility
d <- d %>%
  mutate(
    predicted_mobility = case_when(
      predicted_mobility == "conjugative" ~ "Conjugative",
      predicted_mobility == "mobilizable" ~ "Mobilizable",
      predicted_mobility == "non-mobilizable" ~ "Non-mobilizable",
      TRUE ~ predicted_mobility
    )
  )

mobility <- d %>%
  group_by(predicted_mobility, VF) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

Extended_data_2B <- ggplot(
  mobility,
  aes(
    x = predicted_mobility,
    y = n,
    fill = VF
  )
) +
  geom_col(
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(values = VF_colors) +
  labs(
    x = "Mobility",
    y = "Number of plasmids",
    fill = "Virulence factor"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold")
  )


#Panel C Replicon Type
rep <- d %>%
  group_by(VF) %>%
  summarise(
    `No Replicon` = sum(`No Replicon`, na.rm = TRUE),
    `Col-like` = sum(`Col-like`, na.rm = TRUE),
    IncF = sum(IncF, na.rm = TRUE),
    IncH = sum(IncH, na.rm = TRUE),
    IncI = sum(IncI, na.rm = TRUE),
    `IncB/O/K/Z` = sum(`IncB/O/K/Z_2`, na.rm = TRUE),
    IncN = sum(IncN, na.rm = TRUE),
    IncP = sum(IncP, na.rm = TRUE),
    IncQ = sum(IncQ, na.rm = TRUE),
    IncR = sum(IncR, na.rm = TRUE),
    IncX = sum(IncX, na.rm = TRUE),
    IncY = sum(IncY, na.rm = TRUE),
    p0111 = sum(p0111, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -VF,
    names_to = "Replicon",
    values_to = "value"
  )

# Order by total abundance rather than separately within VF groups
rep_order <- rep %>%
  group_by(Replicon) %>%
  summarise(total = sum(value), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Replicon)

rep$Replicon <- factor(
  rep$Replicon,
  levels = rep_order
)

Extended_data_2C <- ggplot(
  rep,
  aes(
    x = Replicon,
    y = value,
    fill = VF
  )
) +
  geom_col(
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(values = VF_colors) +
  labs(
    x = "Plasmid replicon",
    y = "Number of plasmids",
    fill = "Virulence factor"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


#Pling Network
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

nodes <- d %>%
  mutate(
    name = clean_id(plasmid_id)
  )

all_nodes <- tibble(
  name = unique(
    c(
      edges$plasmid_1,
      edges$plasmid_2
    )
  )
)

vertices <- all_nodes %>%
  left_join(
    nodes,
    by = "name"
  )

set.seed(123)

g <- graph_from_data_frame(
  d = edges,
  vertices = vertices,
  directed = FALSE
)

layout_coords <- create_layout(
  g,
  layout = "fr"
)


#Pling communities
pling_community <- nodes %>%
  filter(
    !is.na(community),
    community != "hub plasmid"
  ) %>%
  group_by(
    community,
    VF
  ) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

community_keep <- pling_community %>%
  group_by(community) %>%
  summarise(
    total = sum(n),
    .groups = "drop"
  ) %>%
  filter(total >= 5)

pling_community_plot <- pling_community %>%
  semi_join(
    community_keep,
    by = "community"
  )

Extended_data_2D <- ggplot(
  pling_community_plot,
  aes(
    x = reorder(community, n, FUN = sum),
    y = n,
    fill = VF
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = VF_colors) +
  labs(
    x = "Pling community",
    y = "Number of plasmids",
    fill = "Virulence factor"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold")
  )


#Panel E: pling subcommunities
pling_subcommunity_0 <- nodes %>%
  filter(
    as.character(community) == "0",
    !is.na(subcommunity)
  ) %>%
  group_by(
    subcommunity,
    VF
  ) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

# Keep subcommunities with at least 5 plasmids in total
subcommunity_keep <- pling_subcommunity_0 %>%
  group_by(subcommunity) %>%
  summarise(
    total = sum(n),
    .groups = "drop"
  ) %>%
  filter(total >= 5)

pling_subcommunity_plot <- pling_subcommunity_0 %>%
  semi_join(
    subcommunity_keep,
    by = "subcommunity"
  )

Extended_data_2E <- ggplot(
  pling_subcommunity_plot,
  aes(
    x = reorder(subcommunity, n, FUN = sum),
    y = n,
    fill = VF
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = VF_colors) +
  labs(
    x = "Pling subcommunity",
    y = "Number of plasmids",
    fill = "Virulence factor"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold")
  )


#Panel F pling
Extended_data_2F <- ggraph(layout_coords) +
  geom_edge_link(
    alpha = 0.05,
    color = "grey60"
  ) +
  geom_node_point(
    aes(color = VF),
    data = layout_coords[
      !is.na(layout_coords$VF),
    ],
    size = 1.5
  ) +
  scale_color_manual(values = VF_colors) +
  theme_void() +
  labs(
    color = "Virulence factor"
  ) +
  theme(
    text = element_text(face = "bold"),
    legend.position = "right"
  )



# Extract ONE legend from a fill-based plot
shared_legend <- cowplot::get_legend(
  Extended_data_2A +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
)

# Remove legends from ALL panels
A <- Extended_data_2A + theme(legend.position = "none")
B <- Extended_data_2B + theme(legend.position = "none")
C <- Extended_data_2C + theme(legend.position = "none")
D <- Extended_data_2D + theme(legend.position = "none")
E <- Extended_data_2E + theme(legend.position = "none")
F <- Extended_data_2F + theme(legend.position = "none")



top_row <- plot_grid(
  A,
  B,
  C,
  labels = c("A", "B", "C"),
  ncol = 3,
  align = "h",
  axis = "tb",
  rel_widths = c(1, 1, 1)
)



middle_row <- plot_grid(
  D,
  E,
  labels = c("D", "E"),
  ncol = 2,
  align = "h",
  axis = "tb",
  rel_widths = c(1, 1)
)


main_figure <- plot_grid(
  top_row,
  middle_row,
  F,
  labels = c("", "", "F"),
  ncol = 1,
  rel_heights = c(1, 1.4, 1)
)

Extended_data_2 <- plot_grid(
  main_figure,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.12)
)


# SAVE figure

ggsave(
  "Extended_data_2.png",
  plot = Extended_data_2,
  width = 16,
  height = 18,
  bg = "white"
)


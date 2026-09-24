library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# Import dataset
p <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Check if plasmids containing CFs also contain toxin
toxin <- p %>%
  filter(
    combination != "Only Toxin",
    combination != "Only Non-Canonical",
    pling != "Hub plasmids"
  )

# Make a column with the combination of toxins in plasmids
toxin <- toxin %>%
  rowwise() %>%
  mutate(
    Toxin_Combination = {
      toxins <- c(
        if (LT == 1) "LT",
        if (STah == 1) "STah",
        if (STap == 1) "STap"
      )
      
      if (length(toxins) == 0) "None" else paste(toxins, collapse = "+")
    }
  ) %>%
  ungroup()

# Total genomes per combination x pling
totals <- toxin %>%
  group_by(combination, pling) %>%
  summarise(total = n(), .groups = "drop")

pling_label <- totals %>%
  mutate(pling_label = paste0(pling, " (n=", total, ")"))

# Counts per toxin combination, excluding None from the stacks
Toxin_presence <- toxin %>%
  filter(Toxin_Combination != "None") %>%
  group_by(combination, pling, Toxin_Combination) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(totals, by = c("combination", "pling")) %>%
  mutate(freq = n / total * 100)

# Ensure plings that are all-None still show up as 0-height bars
all_pling <- distinct(toxin, combination, pling)

Toxin_presence_plot <- all_pling %>%
  left_join(Toxin_presence, by = c("combination", "pling")) %>%
  mutate(
    freq = replace_na(freq, 0),
    Toxin_Combination = replace_na(Toxin_Combination, "No toxin (0%)")
  ) %>%
  left_join(
    pling_label %>% select(combination, pling, pling_label),
    by = c("combination", "pling")
  )

# Color-blind friendly Okabe-Ito palette
okabe_ito <- c(
  "LT" = "#0072B2",
  "STah" = "#E69F00",
  "STap" = "#009E73",
  "LT+STah" = "#CC79A7",
  "LT+STap" = "#56B4E9",
  "STah+STap" = "#D55E00",
  "LT+STah+STap" = "#F0E442"
)

present_toxins <- intersect(
  names(okabe_ito),
  unique(Toxin_presence_plot$Toxin_Combination)
)

# Plot
Extended_data_4 <- Toxin_presence_plot %>%
  ggplot(aes(x = pling_label, y = freq, fill = Toxin_Combination)) +
  geom_col() +
  facet_wrap(vars(combination), scales = "free_y") +
  coord_flip() +
  scale_fill_manual(
    values = okabe_ito[present_toxins],
    breaks = present_toxins,
    drop = FALSE
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Plasmid subcommunity",
    y = "Percentage (%)",
    fill = "Toxin profile"
  )

# Save
ggsave(
  "Extended_Data_4.png",
  plot = Extended_data_4,
  width = 12,
  height = 11,
  dpi = 300,
  bg = "white"
)


library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

library(tidyverse)
library(readr)
library(ggplot2)
library(tidytext)
library(ggpubr)

# -----------------------------
# Read main CF metadata
# -----------------------------

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

# -----------------------------
# Long format for known CFs
# -----------------------------

df_long <- ST_ETEC %>%
  pivot_longer(
    cols = `CFA/I`:CS7,
    names_to = "CF",
    values_to = "presence"
  ) %>%
  filter(presence == 1) %>%
  mutate(
    CF = str_trim(as.character(CF)),
    ST = str_trim(as.character(ST))
  )

# -----------------------------
# Read global MLST table
# ETEC column: 1 = ETEC, 0 = Non-ETEC
# -----------------------------

ETEC_all <- read_delim("~/Desktop/Project_plasmidome/Assemblies_good_quality/E_coli_ATB_mlst_all.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(
    ST = str_trim(as.character(ST)),
    ETEC = as.numeric(ETEC)
  ) %>%
  filter(!is.na(ST), ST != "-", ST != "")

# -----------------------------
# Ratio ETEC per ST
# mean(ETEC) = proportion of ETEC isolates per ST
# -----------------------------

ratio_ETEC_per_ST <- ETEC_all %>%
  group_by(ST) %>%
  summarise(
    ratio_ETEC = mean(ETEC, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Top 10 STs per known CF
# -----------------------------

top10_ST_per_CF <- df_long %>%
  filter(!is.na(ST), ST != "-", ST != "") %>%
  count(CF, ST, name = "n_total") %>%
  group_by(CF) %>%
  slice_max(order_by = n_total, n = 10, with_ties = FALSE) %>%
  ungroup()

ratio_df <- top10_ST_per_CF %>%
  left_join(ratio_ETEC_per_ST, by = "ST")

CF_levels <- c(
  "CS4", "CS6+CS8", "CS7", "CS1", "CS5+CS6", "CS17", "CS19", "CS26", "CFA/I", "CS27",
  "CS3", "CS20/CS28", "CS12", "CS13", "CS30", "CS14", "CS23", "CS21", "CS8", "CS6"
)

ratio_df2 <- ratio_df %>%
  mutate(
    CF = factor(CF, levels = CF_levels),
    ST = str_trim(as.character(ST))
  )

# -----------------------------
# Plot Figure 5A
# -----------------------------

Fig5A <- ratio_df2 %>%
  filter(CF != "CS2", CF != "CS18") %>%
  ggplot(
    aes(
      x = reorder_within(ST, n_total, CF),
      y = n_total,
      fill = ratio_ETEC
    )
  ) +
  geom_col(color = "black", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(~ CF, scales = "free") +
  scale_x_reordered() +
  scale_fill_gradient2(
    low = "white",
    high = "red",
    limits = c(0, 1),
    oob = scales::squish,
    na.value = "grey80"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.spacing.x = unit(0.6, "cm")
  ) +
  labs(
    x = "ST",
    y = "Number of isolates",
    fill = "Ratio ETEC"
  )

# -----------------------------
# NEW CF files
# -----------------------------

NEW_CF_files <- tibble::tibble(
  CF = c(
    "New_gamma-1_1+CS6",
    "New_gamma-2_1",
    "New_gamma-2_1+CS6",
    "New_gamma-2_2+CS13",
    "New_gamma-3_1+CS6+CS8",
    "New_gamma-3_2",
    "New_gamma-3_2+CS6",
    "New_typeIV_pili+CS4",
    "New_typeIV_pili+New_kappa-1+CS6",
    "New_pi_1"
  ),
  file = c(
    "NEW_gamma_1_CS6.csv",
    "New_gamma2_1.csv",
    "New_gamma2_1_CS6.csv",
    "New_gamma_2_2_CS13.csv",
    "NEW_CF_gamma3_1_CS6_CS8.csv",
    "NEW_CF_gamma_3_2.csv",
    "NEW_gamma_3_2_CS6.csv",
    "NEW_typeiv_pili_CS4.csv",
    "New_typeiv_pili_new_kappa_1_CS6.csv",
    "NEW_pi1.csv"
  )
)

NEW_CF_levels <- NEW_CF_files$CF

# -----------------------------
# Function to read NEW files
# -----------------------------

read_new_cf_file <- function(file) {
  read_delim(
    file,
    delim = ";",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE,
    col_types = cols(
      .default = col_guess(),
      ST = col_character(),
      id = col_character(),
      ETEC = col_character()
    )
  )
}

# -----------------------------
# Read NEW files
# -----------------------------

NEW_df_long <- NEW_CF_files %>%
  mutate(data = map(file, read_new_cf_file)) %>%
  unnest(data) %>%
  mutate(
    CF = str_trim(as.character(CF)),
    id = str_trim(as.character(id)),
    ST = str_trim(as.character(ST)),
    ETEC = as.character(ETEC)
  ) %>%
  filter(
    !is.na(id), id != "",
    !is.na(ST), ST != "-", ST != "",
    !is.na(CF), CF != ""
  )

# -----------------------------
# Top 10 STs per NEW CF
# -----------------------------

top10_ST_per_NEW_CF <- NEW_df_long %>%
  count(CF, ST, name = "n_total") %>%
  group_by(CF) %>%
  slice_max(order_by = n_total, n = 10, with_ties = FALSE) %>%
  ungroup()

ratio_NEW_df <- top10_ST_per_NEW_CF %>%
  left_join(ratio_ETEC_per_ST, by = "ST")

ratio_NEW_df2 <- ratio_NEW_df %>%
  mutate(
    CF = factor(CF, levels = NEW_CF_levels),
    ST = str_trim(as.character(ST))
  )

# -----------------------------
# Plot Figure 5B
# -----------------------------

Fig5B <- ratio_NEW_df2 %>%
  ggplot(
    aes(
      x = reorder_within(ST, n_total, CF),
      y = n_total,
      fill = ratio_ETEC
    )
  ) +
  geom_col(color = "black", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(
    ~ CF,
    scales = "free",
    ncol = 5,
    labeller = labeller(CF = label_wrap_gen(width = 16))
  ) +
  scale_x_reordered() +
  scale_fill_gradient2(
    low = "white",
    high = "red",
    limits = c(0, 1),
    oob = scales::squish,
    na.value = "grey80"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold", size = 20),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.spacing.x = unit(0.6, "cm")
  ) +
  labs(
    x = "ST",
    y = "Number of isolates",
    fill = "Ratio ETEC"
  )

# -----------------------------
# Combine Figure 5
# -----------------------------

Figure5 <- ggarrange(
  Fig5A,
  Fig5B,
  ncol = 1,
  heights = c(2, 1),
  labels = c("A", "B"),
  font.label = list(size = 26, face = "bold"),
  common.legend = TRUE,
  legend = "right"
)

ggsave(
  "Figure5_trial.png",
  plot = Figure5,
  width = 26,
  height = 18.5,
  bg = "white"
)
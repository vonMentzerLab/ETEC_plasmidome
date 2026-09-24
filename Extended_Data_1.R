library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(cowplot)
library(rsvg)
library(grid)

# Import dataset
d <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/ETEC_VF_final_summary.tsv", 
  delim = "\t",
  escape_double = FALSE, 
  trim_ws = TRUE
)

# A - Temporal distribution
year <- d %>%
  group_by(Year, Origin) %>%
  summarise(n = n())

Yearg <- year %>%
  ggplot(aes(x = Year, y = n, fill = Origin)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(face = "bold")) +
  ylab("Number of isolates")

# B - Geographical distribution
continent <- d %>%
  group_by(Continent, Origin) %>%
  summarise(n = n())

continentg <- continent %>%
  filter(!is.na(Continent)) %>%
  ggplot(aes(x = reorder(Continent, n), y = n, fill = Origin)) +
  geom_bar(stat = "identity") +
  xlab("Continent") +
  ylab("Number of isolates") +
  theme_classic() +
  theme(text = element_text(face = "bold")) +
  coord_flip()

# C - Virulence factor frequency
CS <- d %>%
  summarise(across(eatA:stx2, mean), n = n())

CS <- CS %>%
  pivot_longer((1:30), names_to = "VF")

CS$Percentage <- CS$value * 100
CS$n <- (CS$Percentage * 535) / 100

write_csv(CS, "CS_freq.csv")

CS_freq <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/CS_freq_mod.csv",
  delim = ";", 
  escape_double = FALSE,
  trim_ws = TRUE
)

colours_Plasmid <- c("red", "black")

CS_freq2 <- CS_freq %>%
  arrange(Class, Percentage) %>%
  mutate(VF = factor(VF, levels = unique(VF)))

CS_freq2$Location[CS_freq2$Location == "Enterotoxin"] <- "Plasmid"
CS_freq2$Location[CS_freq2$Location == "Shiga-toxin"] <- "Chromosome"

class_lines <- CS_freq2 %>%
  group_by(Class) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    end = cumsum(n),
    start = lag(end, default = 0) + 1,
    mid = (start + end) / 2
  )

sep_positions <- class_lines$end[-nrow(class_lines)] + 0.5

FigS2 <- ggplot(
  CS_freq2,
  aes(
    x = VF,
    y = Percentage,
    color = Location
  )
) +
  geom_point(alpha = 0.8, size = 3) +
  geom_vline(
    xintercept = sep_positions,
    colour = "grey50"
  ) +
  scale_color_manual(values = colours_Plasmid) +
  theme_classic() +
  xlab("Virulence factor") +
  ylab("Percentage (%)") +
  theme(
    text = element_text(face = "bold"),
    legend.position = c(0.8, 0.25),
    legend.background = element_rect(
      fill = alpha("white", 0.6)
    ),
    plot.margin = margin(3, 3, 3, 30)
  ) +
  coord_flip()

# Extended Data 1
AB <- ggarrange(
  Yearg,
  continentg,
  labels = c("A", "B"),
  ncol = 2,
  common.legend = TRUE,
  legend = "right"
)

Extended_Data_1 <- ggarrange(
  AB,
  FigS2,
  labels = c("", "C"),
  ncol = 1,
  nrow = 2,
  heights = c(0.8, 1)
)

ggsave(
  "Extended_Data_1.png",
  plot = Extended_Data_1,
  width = 8,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Extended_Data_1.svg",
  plot = Extended_Data_1,
  width = 8,
  height = 8,
  bg = "white"
)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(countrycode)
library(cowplot)
library(ggpubr)
library(grid)

p <- read_delim(
  "~/Desktop/Project_plasmidome/Assemblies_good_quality/Plasmid_CF_summary_DEF.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

CS_combinations <- c(
  "CS1","CS4","CS7","CS6+CS8","CS5+CS6","CFA/I","CS26","CS12","CS3","CS17",
  "CS19","CS27","CS14","CS20/CS28","CS23","CS6","CS21","CS8","CS13","CS30"
)

year_limits <- c(
  floor(min(p$Year, na.rm = TRUE) / 5) * 5,
  ceiling(max(p$Year, na.rm = TRUE) / 5) * 5
)

p_country <- p %>%
  filter(combination %in% CS_combinations) %>%
  mutate(
    pling = str_trim(as.character(pling)),
    Country = str_trim(as.character(Country)),
    Country = ifelse(Country == "unknown" | is.na(Country), NA, Country),
    Country = countrycode(Country, "country.name", "iso2c"),
    Plasmid_Lineage = as.character(Plasmid_Lineage)
  ) %>%
  filter(
    !is.na(Country),
    !is.na(Plasmid_Lineage),
    !is.na(pling)
  )

all_countries <- sort(unique(p_country$Country))


add_panel_title <- function(p, title) {
  p +
    ggtitle(title) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5,
        margin = margin(b = 6)
      )
    )
}


get_palette_for_combination <- function(combination_name, p) {
  
  used_lineages <- p %>%
    filter(
      combination == combination_name,
      pling != "Hub plasmids"
    ) %>%
    pull(Plasmid_Lineage) %>%
    as.character() %>%
    na.omit() %>%
    unique() %>%
    sort()
  
  if (length(used_lineages) == 0) {
    return(
      list(
        cols = NULL,
        lineages = character(0)
      )
    )
  }
  
  cols <- setNames(
    hcl(
      h = seq(50, 310, length.out = length(used_lineages)),
      c = 100,
      l = 60
    ),
    used_lineages
  )
  
  list(
    cols = cols,
    lineages = used_lineages
  )
}


make_year_panel <- function(combination_name, p, year_limits, palette_cols) {
  
  df_sum <- p %>%
    filter(pling != "Hub plasmids") %>%
    mutate(
      pling = str_trim(as.character(pling))
    ) %>%
    filter(combination == combination_name) %>%
    group_by(
      pling,
      Plasmid_Lineage
    ) %>%
    summarise(
      first_year = suppressWarnings(min(Year, na.rm = TRUE)),
      last_year = suppressWarnings(max(Year, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    filter(
      is.finite(first_year),
      is.finite(last_year)
    )
  
  if (nrow(df_sum) == 0) {
    return(
      ggplot() +
        theme_void() +
        ggtitle("No YEAR data")
    )
  }
  
  ggplot(df_sum, aes(y = pling)) +
    geom_segment(
      aes(
        x = first_year,
        xend = last_year,
        yend = pling
      ),
      linewidth = 1,
      colour = "grey40"
    ) +
    geom_point(
      aes(
        x = first_year,
        colour = Plasmid_Lineage
      ),
      size = 2.6
    ) +
    geom_point(
      aes(
        x = last_year,
        colour = Plasmid_Lineage
      ),
      size = 2.6
    ) +
    scale_colour_manual(
      values = palette_cols,
      na.translate = FALSE
    ) +
    scale_x_continuous(
      limits = year_limits,
      breaks = seq(
        year_limits[1],
        year_limits[2],
        by = 5
      ),
      expand = c(0, 0)
    ) +
    labs(
      x = "",
      y = NULL
    ) +
    theme_classic() +
    theme(
      text = element_text(face = "bold"),
      legend.position = "none",
      axis.line.x = element_line(linewidth = 0.8),
      axis.ticks.x = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text.x = element_text(
        size = 9,
        angle = 90
      ),
      axis.title.x = element_text(size = 10)
    )
}


make_country_panel <- function(combination_name, p_country, all_countries, palette_cols) {
  
  df_tile <- p_country %>%
    filter(pling != "Hub plasmids") %>%
    filter(combination == combination_name) %>%
    distinct(
      pling,
      Country,
      Plasmid_Lineage
    ) %>%
    group_by(
      pling,
      Country
    ) %>%
    summarise(
      Plasmid_Lineage = first(Plasmid_Lineage),
      .groups = "drop"
    ) %>%
    mutate(
      Country = factor(
        Country,
        levels = all_countries
      ),
      pling = factor(pling)
    ) %>%
    tidyr::complete(
      pling,
      Country = all_countries,
      fill = list(
        Plasmid_Lineage = NA_character_
      )
    )
  
  if (nrow(df_tile) == 0) {
    return(
      ggplot() +
        theme_void() +
        ggtitle("No COUNTRY data")
    )
  }
  
  ggplot(
    df_tile,
    aes(
      x = Country,
      y = pling,
      fill = Plasmid_Lineage
    )
  ) +
    geom_tile(
      color = "grey90",
      linewidth = 0.2
    ) +
    scale_fill_manual(
      values = palette_cols,
      na.translate = FALSE
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic() +
    theme(
      text = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}


year_panel_list <- lapply(
  CS_combinations,
  function(cf) {
    
    pal <- get_palette_for_combination(
      cf,
      p
    )
    
    if (is.null(pal$cols) || length(pal$cols) == 0) {
      return(
        ggplot() +
          theme_void() +
          ggtitle(cf)
      )
    }
    
    make_year_panel(
      combination_name = cf,
      p = p,
      year_limits = year_limits,
      palette_cols = pal$cols
    ) +
      ggtitle(cf) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 12,
          hjust = 0.5
        )
      )
  }
)


country_panel_list <- lapply(
  CS_combinations,
  function(cf) {
    
    pal <- get_palette_for_combination(
      cf,
      p
    )
    
    if (is.null(pal$cols) || length(pal$cols) == 0) {
      return(
        ggplot() +
          theme_void() +
          ggtitle(cf)
      )
    }
    
    make_country_panel(
      combination_name = cf,
      p_country = p_country,
      all_countries = all_countries,
      palette_cols = pal$cols
    ) +
      ggtitle(cf) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 12,
          hjust = 0.5
        )
      )
  }
)


Extended_Data_5A <- ggarrange(
  plotlist = year_panel_list,
  ncol = 5,
  nrow = 4
)

Extended_Data_5B <- ggarrange(
  plotlist = country_panel_list,
  ncol = 5,
  nrow = 4
)

Extended_Data_5 <- ggarrange(
  Extended_Data_5A,
  Extended_Data_5B,
  ncol = 1,
  nrow = 2,
  labels = c("A", "B"),
  font.label = list(
    size = 20,
    face = "bold"
  ),
  heights = c(1, 1)
)

ggsave(
  "Extended_Data_5.png",
  plot = Extended_Data_5,
  width = 20,
  height =20,
  bg = "white",
  dpi=300
)

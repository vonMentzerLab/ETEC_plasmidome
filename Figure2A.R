#Figure2A ETEC (cophylo)
library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(forcats)
library(readxl)
#Import dataset
#Import datasets and trees
links_df <- read_excel("~/Desktop/Project_plasmidome/Assemblies_good_quality/contigs/links_plasmid_genome.xlsx")
#We upload a core genome tree (based on core gene alignment) and a mash plasmid tree
genome_tree <- read.tree("~/Desktop/Project_plasmidome/Assemblies_good_quality/ETEC_maximum_likelihood.nwk")
plasmid_tree <- read.tree("~/Desktop/Project_plasmidome/Assemblies_good_quality/plasmid_tree_NJ_containment_distance.nwk")
#Clean names (remove fasta)
clean <- function(x) {
  x %>%
    sub("\\.f(ast)?a$|\\.fna$|\\.fa$|\\.fasta$", "", ., ignore.case = TRUE) %>%  # remove extensions
    sub("_fasta$", "", .) %>%                                                   # remove _fasta
    sub("\\s+$", "", .)                                                         # trim trailing spaces
}

genome_tree$tip.label <- clean(genome_tree$tip.label)
plasmid_tree$tip.label <- clean(plasmid_tree$tip.label)
links_df$plasmid<-clean(links_df$plasmid)
links_df$genome<-clean(links_df$genome)

p_genome <- ggtree(genome_tree,size=0.3) + theme_tree()
#Genome tip coords
g_coords <- p_genome$data %>% filter(isTip) %>% select(label, x, y)
colnames(g_coords) <- c("genome", "x_genome", "y_genome")
#Plot plasmid tree
p_plasmid <- ggtree(plasmid_tree,size=0.3) + theme_tree()
#Flip horizontally and shift right
offset <- max(p_genome$data$x) + 2
p_plasmid$data$x <- offset - p_plasmid$data$x
#Plasmid coords
p_coords <- p_plasmid$data %>% filter(isTip) %>% select(label, x, y)
colnames(p_coords) <- c("plasmid", "x_plasmid", "y_plasmid")

#Merge coordinates
link_coords <- links_df %>%
  left_join(g_coords, by = "genome") %>%
  left_join(p_coords, by = "plasmid") %>%
  filter(
    !is.na(x_genome), !is.na(y_genome),
    !is.na(x_plasmid), !is.na(y_plasmid),
    !(x_genome == x_plasmid & y_genome == y_plasmid)
  )
#Tanglegram
# Keep only combinations with 5 or more rows
facet_order <- c(
  "CS4","CS6+CS8","CS7","CS1","CS5+CS6",
  "CS17","CS19","CS26","CFA/I","CS27",
  "CS3","CS20/CS28","CS12","CS13","CS30",
  "CS14","CS23","CS21","CS8","CS6"
)
link_coords_filtered <- link_coords %>%
  count(combination) %>%
  filter(n >= 5) %>%
  filter(combination!="Only New"&combination!="Only toxin")%>%
  inner_join(link_coords, by = "combination") %>%
  mutate(combination = factor(combination, levels = facet_order))

cophyloCF <- ggplot() +
  geom_tree(data = p_genome$data, aes(x = x, y = y), linewidth = 0.3) +
  geom_tree(data = p_plasmid$data, aes(x = x, y = y), linewidth = 0.3) +
  geom_curve(
    data = link_coords_filtered,
    aes(
      x = x_genome, y = y_genome,
      xend = x_plasmid, yend = y_plasmid,
      color = as.factor(Mobility)
    ),
    curvature = 0,
    alpha = 0.5,
    linewidth = 0.25
  ) +
  scale_color_manual(
    values = c("#E83F6F", "#FFBF00", "#2274A5")
  ) +
  theme_void() +
  labs(colour = "Plasmid mobility") +
  theme(
    text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  facet_wrap(
    ~ combination,
    nrow = 4,
    ncol = 5,
    drop = FALSE
  )
ggsave("Fig2_modify.png",plot=cophyloCF,dpi=300,width = 7.5,height = 6,bg="white")
ggsave("Fig2_modify.svg",plot=cophyloCF,dpi=300,width = 7.5,height = 6,bg="white")
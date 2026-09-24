#Figure S17
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
#Import dataset with ETEC frequency
d <- read_delim("~/Desktop/Project_plasmidome/Assemblies_good_quality/NEW_CF_stats.csv", 
                           delim = ";", escape_double = FALSE, trim_ws = TRUE)
#Barplots ETEC
VF_colors <- c(
  "Non-ETEC" = "#BDBDBD",
  "ETEC" = "#D95F02"
)

Extended_Data_8<-d%>%ggplot(aes(x=NEW_CF,y=n,fill=ETEC))+
  geom_bar(stat="identity")+coord_flip()+
  theme_classic()+theme(text = element_text(face="bold"))+scale_fill_manual(values =VF_colors )+
  xlab("New putative CF operon")+ylab("Number of isolates")

ggsave(
  "Extended_data_8.png",
  plot = Extended_Data_8,
  width = 6,
  height = 5,
  bg = "white"
)



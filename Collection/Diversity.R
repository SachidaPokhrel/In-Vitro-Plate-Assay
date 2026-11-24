#load necessary library
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)
library(multcomp)
library(multcompView)
library(dplyr)
library(car)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb18 <- c(
  "#000000", # black
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  
  "#882255", # wine
  "#44AA99", # teal
  "#117733", # green
  "#332288", # indigo
  "#88CCEE", # light blue
  "#DDCC77", # sand
  "#CC6677", # rose
  "#AA4499", # magenta
  "#999933", # olive
  "#661100"  # brown
)

collection <- read.csv("Collection/YeastCollection.csv", na.strings = "na") 

collection2 <- collection[collection$Date.Collected == "2021",]

diversity <- ggplot(collection2,
       aes(x = Method,
           fill = Genus)) +
  geom_bar(width = 0.8) +
  labs(x = "",
       y = "Count") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.spacing = unit(0.05, "lines")
  ) +
  scale_fill_manual(values = cb18)

library(dplyr)

richness <- collection2 %>%
  group_by(Method, Tree, Leaf) %>%
  summarise(
    Richness = n_distinct(ManuscriptName_phylotype),
    .groups = "drop"
  )

richness.plot <- ggplot(richness,
       aes(Method, Richness)) +
  
  stat_summary(
    fun = mean,
    geom = "col",
    fill = "grey40",
    color = "black",
    width = 0.75
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.25
  ) +
  
  geom_point(
    position = position_jitter(width = 0.08),
    size = 2
  ) +
  
  theme_classic() +
  labs(
    x = NULL,
    y = "Richness per leaf"
  ) +
  ggpubr::stat_compare_means(
    method = "t.test",
    label = "p.format"
  ) 


ggarrange(diversity, richness.plot, labels = "auto")



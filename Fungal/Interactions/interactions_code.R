filename1 = "Fungal/Interactions/Interactions Metadata.csv"
Interactions = read.csv(filename1, na.strings = "na", header = TRUE)

#unlink("C:/Users/beatr/AppData/Local/R/win-library/4.3/00LOCK", recursive = TRUE)
library(ggplot2)
#install.packages("ggtext", dependencies = TRUE, INSTALL_opts = '--no-lock')
library(ggtext)

# color blind pallet used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000") 
tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

Interactions2 <- Interactions %>%
  subset(Pathogen != "F98")

Interactions2$Average <- mean(Interactions2$Rating)

ggplot(Interactions2, aes(x = as.factor(Yeast), y = Rating, fill = Pathogen)) + geom_bar(stat = 'identity')

xlabels <- c("*M. elongata* NVP5", "*G. irregulare* OEO_O28", "*F. graminearum* PH-1", "*R. solani* R_12")

ggplot(Interactions2, aes(x = as.factor(Yeast), y = Rating, fill = Pathogen)) + 
  stat_summary(fun.y=mean,geom="bar", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) + 
  ylab("Rating") + 
  xlab("Yeast")+ 
  scale_x_discrete(label = xlabels) +
  geom_jitter(width = 0.1, alpha = 0.8) + 
  #scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +
  #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
  theme_classic() +
  theme(legend.position="none")

ggplot(Interactions2, aes(x=Pathogen, y = Label, fill = as.factor(Rating))) + 
  geom_tile() +
  scale_fill_manual(values = ibm.cbb) +
  guides(fill = guide_legend(title = "Rating")) +
  labs(title = "Interactions Data", y = "Taxonomy")

ggplot(Interactions, aes(x=Pathogen, y = Label, fill = as.factor(Average))) + 
  geom_tile() +
  scale_fill_manual(values = tol.cbb) +
  guides(fill = guide_legend(title = "Average Rating")) +
  labs(title = "Interactions Data", y = "Taxonomy")

ggplot(Interactions2, aes(x=Pathogen, y = Label, fill = as.factor(Average))) + 
  geom_tile() +
  scale_fill_manual(values = tol.cbb) +
  guides(fill = guide_legend(title = "Average Rating")) +
  labs(title = "Interactions Data", y = "Taxonomy")

yeast_fungal <- ggplot(Interactions2, aes(x=Pathogen, y = Isolate, fill = CombinedAvg)) + 
  geom_tile() +
  scale_x_discrete(labels = c(expression(italic("M. elongata") ~ "NVP5"), expression(italic("G. irregulare") ~ "OEO_O28"),
                              expression(italic("F. graminearum") ~ "PH-1"), expression(italic("R. solani") ~ "R_12"))) +
#  scale_fill_manual(values = tol.cbb) +
  guides(fill = guide_legend(title = "Average Rating")) +
  labs(title = "Yeast-Fungi Interactions", y = "Taxonomy") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_gradientn(
    colours = c(
      "#56B4E9",  # light blue
      "#009E73",  # green
      "#F0E442",  # yellow
      "#E69F00"   # orange
    )
  )


filename2 = "Fungal/Interactions/Fusarium_Weight.csv"
Fusarium = read.csv(filename2, na.strings = "na", header = TRUE)

fusarim_wetweight <- ggplot(Fusarium, aes(x = Treatment, y = WetWeight)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.05) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_discrete(labels = c(
    "F3/Fus"   = "EMM_F3 + PH-1",
    "Fus alone" = "PH-1"
  )) +
  #ggpubr::stat_compare_means(
  #  method = "t.test",
  #  label = "p.format"
  #) +
  xlab("") +
  ylab("Biomass (mg)")

t.test(
  WetWeight ~ Treatment,
  data = Fusarium,
  alternative = "less"
)

ggpubr::ggarrange(yeast_fungal, fusarim_wetweight, labels = "auto", nrow = 2)





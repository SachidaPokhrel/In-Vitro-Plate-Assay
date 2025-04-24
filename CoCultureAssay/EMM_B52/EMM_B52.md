## Things to know before Analysis

***This is applied for all the dataset for Co-culture Assay***

- One .csv file is loaded for one bacteria interacting with multiple
  yeast isolates
- There are three R chunks for each bacterial data.
- There are 16 different yeast interacting with one bacteria. So, to
  make it more clear in the plot and to make it less cumbersome, these
  yeast isolates are further categorised into 4 Yeast classes.
- There are 8 colonies of yeast and 8 colonies of bacteria at different
  distances in one plate. Data is taken only for bacterial colony size
  which reflects the impact of yeast on bacteria.
- Plots have multiple line graphs, thus the significant letters are not
  assigned through code and will be added later after exporting the
  plot.

``` r
#load necessary library
library(ggplot2)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.4     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.4     ✔ tidyr     1.3.1
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggpubr)
library(rstatix)
```

    ## 
    ## Attaching package: 'rstatix'
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(emmeans)
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
library(multcomp)
```

    ## Loading required package: mvtnorm
    ## Loading required package: survival
    ## Loading required package: TH.data
    ## Loading required package: MASS
    ## 
    ## Attaching package: 'MASS'
    ## 
    ## The following object is masked from 'package:rstatix':
    ## 
    ##     select
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select
    ## 
    ## 
    ## Attaching package: 'TH.data'
    ## 
    ## The following object is masked from 'package:MASS':
    ## 
    ##     geyser

``` r
library(multcompView)
library(dplyr)
library(car)
```

    ## Loading required package: carData
    ## 
    ## Attaching package: 'car'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     some

``` r
library(nlme)
```

    ## 
    ## Attaching package: 'nlme'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

## ***Methylobacterium platani* EMM_B52**

Plot is generated using loop around the 4 different classes of yeast,
coming up with 4 plots as an output which will be combined in one plot.

``` r
#read data
B52 <- read.csv("CoCultureAssay/CoCultureAssayData/MergedB52.csv", na.strings = "na") 
#load cbb color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
str(B52)
```

    ## 'data.frame':    1728 obs. of  9 variables:
    ##  $ Bacteria               : chr  "B52" "B52" "B52" "B52" ...
    ##  $ Yeast                  : chr  "EMM_F3" "EMM_F3" "EMM_F3" "EMM_F3" ...
    ##  $ Class                  : chr  "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" ...
    ##  $ Replication            : int  1 1 1 1 1 1 1 1 2 2 ...
    ##  $ DAI                    : int  2 2 2 2 2 2 2 2 2 2 ...
    ##  $ distance_to_yeast      : num  0 11.4 17.4 24.9 31.8 ...
    ##  $ colony_diameter        : num  6.21 5.75 5.75 6.15 5.39 6.14 6.39 6.39 6.33 6.33 ...
    ##  $ colony_diameter_control: num  6.34 5.63 5.67 8.45 7.26 6.88 7.6 6.69 7.01 7.02 ...
    ##  $ increase               : num  0 0 0 0 0 0 0 0 0 0 ...

``` r
#round off of the distance to yeast so that it is visibly a categorical variable
B52$distance_to_yeast <- round(B52$distance_to_yeast, 0)
#change into categorical dataset
B52$Replication=as.factor(B52$Replication)
B52$Yeast=as.factor(B52$Yeast)
B52$DAI=as.factor(B52$DAI)
B52$distance_to_yeast=as.factor(B52$distance_to_yeast)

#remove the data for the distance zero
B52.no.contact <- B52[B52$distance_to_yeast != 0, ]
head(B52.no.contact)
```

    ##   Bacteria  Yeast           Class Replication DAI distance_to_yeast
    ## 2      B52 EMM_F3 Dothideomycetes           1   2                11
    ## 3      B52 EMM_F3 Dothideomycetes           1   2                17
    ## 4      B52 EMM_F3 Dothideomycetes           1   2                25
    ## 5      B52 EMM_F3 Dothideomycetes           1   2                32
    ## 6      B52 EMM_F3 Dothideomycetes           1   2                41
    ## 7      B52 EMM_F3 Dothideomycetes           1   2                48
    ##   colony_diameter colony_diameter_control increase
    ## 2            5.75                    5.63        0
    ## 3            5.75                    5.67        0
    ## 4            6.15                    8.45        0
    ## 5            5.39                    7.26        0
    ## 6            6.14                    6.88        0
    ## 7            6.39                    7.60        0

``` r
# remove data first data since it is a reference for increase in colony size and donot follow the assumption of linear model i.e. normal distribution
B52.no.1st.data <- B52.no.contact[B52.no.contact$DAI != "2",]

#for plot purposes we used Day 6 data
B52.Day6 <- B52.no.contact[B52.no.contact$DAI == "6",]


plot_listB52 <- list()

# Get unique classes
unique_class <- unique(B52.Day6$Class[B52.Day6$Class != "Control"])

# Loop through each class
for (i in seq_along(unique_class)) {
  
  # Filter for Control and the current class
  data <- B52.Day6 %>%
    filter(Class %in% c("Control", unique_class[i]))
  
  # Generate the plot
  p <- ggplot(data, aes(x = distance_to_yeast, y = increase, group = Yeast, color = Yeast)) +
    stat_summary(fun = mean, geom = "line") +   
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +   
    ylab("Increase in colony diameter (mm)") +   
    xlab("Distance to Yeast (mm)") +   
    theme_classic() +
    scale_color_manual(values = cbbPalette) +
    ggtitle(unique_class[i]) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Store the plot
  plot_listB52[[i]] <- p
}

combined_plotB52 <- ggarrange(plotlist = plot_listB52, ncol = 2, nrow = 2, common.legend = FALSE)
```

    ## Warning: Removed 7 rows containing non-finite outside the scale range
    ## (`stat_summary()`).
    ## Removed 7 rows containing non-finite outside the scale range
    ## (`stat_summary()`).

``` r
# Annotation for the title
final_plotB52 <- annotate_figure(combined_plotB52,
                                  top = text_grob(
    expression("Impact on growth of"~italic("Methylobacterium platanii")~"EMM_B52 by Yeast"), face = "bold", color = "Blue2", size = 14, hjust = 0.5))

print(final_plotB52)
```

![](EMM_B52_files/figure-gfm/Plot%20for%20EMM_B52-1.png)<!-- -->

### Stats *Methylobacterium platanii* EMM_B52

We are using linear mixed model. Our dependent variable or y is increase
(increase in colony diameter from 1st data) and independent variables
are different Yeast isolates, days after inoculation (DAI), and distance
to yeast which is the distance between the yeast and bacterial colony in
the plate. Each plate is replicated 3 times.

``` r
#filter data to remove 1st day data since the first data is taken as a base to measure the increase in colony size to rule out the variability that is caused by the drop inoculation. So, initially the increase in the colony diameter for 1st data for all colony is "0" that violates the assumption of normality, thus we remove that from analysis. This would be similar for all the bacterial isolates.
B52.no.1st.data <- B52.no.contact[B52.no.contact$DAI != "2",]
B52try <- lme(increase~distance_to_yeast*DAI*Yeast, data = B52.no.1st.data, random = ~1|Replication, na.action = na.omit)
anova(B52try)
```

    ##                             numDF denDF   F-value p-value
    ## (Intercept)                     1   760 1635.4062  <.0001
    ## distance_to_yeast               6   760    3.1752  0.0044
    ## DAI                             2   760  317.0915  <.0001
    ## Yeast                          16   760   29.5620  <.0001
    ## distance_to_yeast:DAI          12   760    0.7432  0.7093
    ## distance_to_yeast:Yeast        96   760    1.4432  0.0054
    ## DAI:Yeast                      32   760    1.6741  0.0118
    ## distance_to_yeast:DAI:Yeast   192   760    0.2265  1.0000

``` r
resultsB52=lme(increase~DAI+distance_to_yeast*Yeast+DAI*Yeast, data = B52.no.1st.data, random = ~1|Replication, na.action = na.omit)
summary(resultsB52)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: B52.no.1st.data 
    ##        AIC     BIC    logLik
    ##   2597.509 3352.85 -1143.755
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept)  Residual
    ## StdDev:  0.07065429 0.6584648
    ## 
    ## Fixed effects:  increase ~ DAI + distance_to_yeast * Yeast + DAI * Yeast 
    ##                                       Value Std.Error  DF   t-value p-value
    ## (Intercept)                       2.0408730 0.1806480 964 11.297510  0.0000
    ## DAI6                              0.8897619 0.1436888 964  6.192285  0.0000
    ## DAI8                              1.2659524 0.1436888 964  8.810377  0.0000
    ## distance_to_yeast17              -0.1511111 0.2194883 964 -0.688470  0.4913
    ## distance_to_yeast25              -0.4616667 0.2194883 964 -2.103378  0.0357
    ## distance_to_yeast32              -0.1516667 0.2194883 964 -0.691001  0.4897
    ## distance_to_yeast41              -0.0794444 0.2194883 964 -0.361953  0.7175
    ## distance_to_yeast48              -0.2594444 0.2194883 964 -1.182042  0.2375
    ## distance_to_yeast55              -0.1127778 0.2194883 964 -0.513821  0.6075
    ## YeastEMM_F3                      -1.1065873 0.3048100 964 -3.630417  0.0003
    ## YeastEMM_F34                     -0.5810317 0.3048100 964 -1.906210  0.0569
    ## YeastEMM_F47                     -1.7605556 0.3048100 964 -5.775912  0.0000
    ## YeastEMM_F48                     -1.2870635 0.3048100 964 -4.222511  0.0000
    ## YeastEMM_F49                     -1.6959524 0.3048100 964 -5.563966  0.0000
    ## YeastEMM_F5                      -0.8683306 0.3050193 964 -2.846806  0.0045
    ## YeastEMM_F63                     -1.0616667 0.3048100 964 -3.483045  0.0005
    ## YeastEMM_F64                     -0.1964286 0.3048100 964 -0.644430  0.5195
    ## YeastEMM_F65                     -0.6578571 0.3048100 964 -2.158253  0.0312
    ## YeastEMM_F66                     -0.6946145 0.3235848 964 -2.146623  0.0321
    ## YeastEMM_F70                     -0.1294444 0.3048100 964 -0.424673  0.6712
    ## YeastEMM_F89                     -1.1051587 0.3048100 964 -3.625730  0.0003
    ## YeastF7                          -0.1864286 0.3048100 964 -0.611622  0.5409
    ## YeastSP_F14                      -1.3942063 0.3048100 964 -4.574018  0.0000
    ## YeastZAN_F3                      -0.9919841 0.3048100 964 -3.254435  0.0012
    ## YeastZAN_F4                      -0.9311905 0.3048100 964 -3.054987  0.0023
    ## distance_to_yeast17:YeastEMM_F3  -0.0466667 0.3801648 964 -0.122754  0.9023
    ## distance_to_yeast25:YeastEMM_F3   0.2872222 0.3801648 964  0.755520  0.4501
    ## distance_to_yeast32:YeastEMM_F3   0.3205556 0.3801648 964  0.843202  0.3993
    ## distance_to_yeast41:YeastEMM_F3   0.3638889 0.3801648 964  0.957187  0.3387
    ## distance_to_yeast48:YeastEMM_F3   0.1527778 0.3801648 964  0.401872  0.6879
    ## distance_to_yeast55:YeastEMM_F3   0.3650000 0.3801648 964  0.960110  0.3372
    ## distance_to_yeast17:YeastEMM_F34  0.0444444 0.3801648 964  0.116908  0.9070
    ## distance_to_yeast25:YeastEMM_F34  0.0494444 0.3801648 964  0.130061  0.8965
    ## distance_to_yeast32:YeastEMM_F34  0.0505556 0.3801648 964  0.132983  0.8942
    ## distance_to_yeast41:YeastEMM_F34  0.1650000 0.3801648 964  0.434022  0.6644
    ## distance_to_yeast48:YeastEMM_F34  0.2861111 0.3801648 964  0.752598  0.4519
    ## distance_to_yeast55:YeastEMM_F34 -0.3683333 0.3801648 964 -0.968878  0.3328
    ## distance_to_yeast17:YeastEMM_F47  0.4677778 0.3801648 964  1.230460  0.2188
    ## distance_to_yeast25:YeastEMM_F47  1.1350000 0.3801648 964  2.985547  0.0029
    ## distance_to_yeast32:YeastEMM_F47  0.3772222 0.3801648 964  0.992260  0.3213
    ## distance_to_yeast41:YeastEMM_F47  0.4583333 0.3801648 964  1.205617  0.2283
    ## distance_to_yeast48:YeastEMM_F47  0.4627778 0.3801648 964  1.217308  0.2238
    ## distance_to_yeast55:YeastEMM_F47  0.1694444 0.3801648 964  0.445713  0.6559
    ## distance_to_yeast17:YeastEMM_F48 -0.0722222 0.3801648 964 -0.189976  0.8494
    ## distance_to_yeast25:YeastEMM_F48  0.1838889 0.3801648 964  0.483708  0.6287
    ## distance_to_yeast32:YeastEMM_F48  0.2894444 0.3801648 964  0.761366  0.4466
    ## distance_to_yeast41:YeastEMM_F48 -0.0627778 0.3801648 964 -0.165133  0.8689
    ## distance_to_yeast48:YeastEMM_F48  0.4138889 0.3801648 964  1.088709  0.2766
    ## distance_to_yeast55:YeastEMM_F48  0.1172222 0.3801648 964  0.308346  0.7579
    ## distance_to_yeast17:YeastEMM_F49  0.4377778 0.3801648 964  1.151547  0.2498
    ## distance_to_yeast25:YeastEMM_F49  0.8316667 0.3801648 964  2.187648  0.0289
    ## distance_to_yeast32:YeastEMM_F49  0.4838889 0.3801648 964  1.272840  0.2034
    ## distance_to_yeast41:YeastEMM_F49  0.7116667 0.3801648 964  1.871995  0.0615
    ## distance_to_yeast48:YeastEMM_F49  0.9350000 0.3801648 964  2.459460  0.0141
    ## distance_to_yeast55:YeastEMM_F49  0.4916667 0.3801648 964  1.293299  0.1962
    ## distance_to_yeast17:YeastEMM_F5   0.4544444 0.3801648 964  1.195388  0.2322
    ## distance_to_yeast25:YeastEMM_F5   0.5594444 0.3801648 964  1.471584  0.1415
    ## distance_to_yeast32:YeastEMM_F5   0.4950000 0.3801648 964  1.302067  0.1932
    ## distance_to_yeast41:YeastEMM_F5   0.7172222 0.3801648 964  1.886609  0.0595
    ## distance_to_yeast48:YeastEMM_F5   0.7905556 0.3801648 964  2.079507  0.0378
    ## distance_to_yeast55:YeastEMM_F5   0.1116478 0.3883053 964  0.287526  0.7738
    ## distance_to_yeast17:YeastEMM_F63 -0.0144444 0.3801648 964 -0.037995  0.9697
    ## distance_to_yeast25:YeastEMM_F63  0.4905556 0.3801648 964  1.290376  0.1972
    ## distance_to_yeast32:YeastEMM_F63  0.1016667 0.3801648 964  0.267428  0.7892
    ## distance_to_yeast41:YeastEMM_F63 -0.3350000 0.3801648 964 -0.881197  0.3784
    ## distance_to_yeast48:YeastEMM_F63  0.6305556 0.3801648 964  1.658637  0.0975
    ## distance_to_yeast55:YeastEMM_F63  0.2116667 0.3801648 964  0.556776  0.5778
    ## distance_to_yeast17:YeastEMM_F64 -0.5144444 0.3801648 964 -1.353214  0.1763
    ## distance_to_yeast25:YeastEMM_F64 -0.8783333 0.3801648 964 -2.310401  0.0211
    ## distance_to_yeast32:YeastEMM_F64 -1.2361111 0.3801648 964 -3.251514  0.0012
    ## distance_to_yeast41:YeastEMM_F64 -0.9838889 0.3801648 964 -2.588059  0.0098
    ## distance_to_yeast48:YeastEMM_F64 -0.4472222 0.3801648 964 -1.176390  0.2397
    ## distance_to_yeast55:YeastEMM_F64 -0.9283333 0.3801648 964 -2.441923  0.0148
    ## distance_to_yeast17:YeastEMM_F65 -0.5288889 0.3801648 964 -1.391209  0.1645
    ## distance_to_yeast25:YeastEMM_F65 -0.3916667 0.3801648 964 -1.030255  0.3031
    ## distance_to_yeast32:YeastEMM_F65  0.2916667 0.3801648 964  0.767211  0.4431
    ## distance_to_yeast41:YeastEMM_F65 -0.9927778 0.3801648 964 -2.611441  0.0092
    ## distance_to_yeast48:YeastEMM_F65 -0.1205556 0.3801648 964 -0.317114  0.7512
    ## distance_to_yeast55:YeastEMM_F65 -0.1027778 0.3801648 964 -0.270351  0.7869
    ## distance_to_yeast17:YeastEMM_F66  0.5696825 0.4147938 964  1.373411  0.1699
    ## distance_to_yeast25:YeastEMM_F66 -0.2969048 0.4147938 964 -0.715789  0.4743
    ## distance_to_yeast32:YeastEMM_F66 -0.0640476 0.4147938 964 -0.154408  0.8773
    ## distance_to_yeast41:YeastEMM_F66 -0.0719841 0.4147938 964 -0.173542  0.8623
    ## distance_to_yeast48:YeastEMM_F66  0.4751587 0.4147938 964  1.145530  0.2523
    ## distance_to_yeast55:YeastEMM_F66  0.3670635 0.4147938 964  0.884930  0.3764
    ## distance_to_yeast17:YeastEMM_F70 -0.7944444 0.3801648 964 -2.089737  0.0369
    ## distance_to_yeast25:YeastEMM_F70 -0.2727778 0.3801648 964 -0.717525  0.4732
    ## distance_to_yeast32:YeastEMM_F70 -0.9416667 0.3801648 964 -2.476996  0.0134
    ## distance_to_yeast41:YeastEMM_F70 -0.0638889 0.3801648 964 -0.168056  0.8666
    ## distance_to_yeast48:YeastEMM_F70 -0.2861111 0.3801648 964 -0.752598  0.4519
    ## distance_to_yeast55:YeastEMM_F70 -0.6450000 0.3801648 964 -1.696633  0.0901
    ## distance_to_yeast17:YeastEMM_F89  0.0255556 0.3801648 964  0.067222  0.9464
    ## distance_to_yeast25:YeastEMM_F89  0.3872222 0.3801648 964  1.018564  0.3087
    ## distance_to_yeast32:YeastEMM_F89  0.5261111 0.3801648 964  1.383903  0.1667
    ## distance_to_yeast41:YeastEMM_F89  0.0350000 0.3801648 964  0.092065  0.9267
    ## distance_to_yeast48:YeastEMM_F89  0.2738889 0.3801648 964  0.720448  0.4714
    ## distance_to_yeast55:YeastEMM_F89 -0.2816667 0.3801648 964 -0.740907  0.4589
    ## distance_to_yeast17:YeastF7       0.1022222 0.3801648 964  0.268889  0.7881
    ## distance_to_yeast25:YeastF7       0.0450000 0.3801648 964  0.118370  0.9058
    ## distance_to_yeast32:YeastF7      -0.0805556 0.3801648 964 -0.211896  0.8322
    ## distance_to_yeast41:YeastF7      -0.5416667 0.3801648 964 -1.424821  0.1545
    ## distance_to_yeast48:YeastF7      -0.3338889 0.3801648 964 -0.878274  0.3800
    ## distance_to_yeast55:YeastF7      -0.7461111 0.3801648 964 -1.962599  0.0500
    ## distance_to_yeast17:YeastSP_F14   0.6000000 0.3801648 964  1.578263  0.1148
    ## distance_to_yeast25:YeastSP_F14   0.3027778 0.3801648 964  0.796438  0.4260
    ## distance_to_yeast32:YeastSP_F14   0.4494444 0.3801648 964  1.182236  0.2374
    ## distance_to_yeast41:YeastSP_F14   0.8316667 0.3801648 964  2.187648  0.0289
    ## distance_to_yeast48:YeastSP_F14   0.4683333 0.3801648 964  1.231922  0.2183
    ## distance_to_yeast55:YeastSP_F14   0.0472222 0.3801648 964  0.124215  0.9012
    ## distance_to_yeast17:YeastZAN_F3   0.2044444 0.3801648 964  0.537778  0.5909
    ## distance_to_yeast25:YeastZAN_F3   0.4883333 0.3801648 964  1.284531  0.1993
    ## distance_to_yeast32:YeastZAN_F3   0.4527778 0.3801648 964  1.191004  0.2339
    ## distance_to_yeast41:YeastZAN_F3  -0.1927778 0.3801648 964 -0.507090  0.6122
    ## distance_to_yeast48:YeastZAN_F3   0.1772222 0.3801648 964  0.466172  0.6412
    ## distance_to_yeast55:YeastZAN_F3  -0.3927778 0.3801648 964 -1.033178  0.3018
    ## distance_to_yeast17:YeastZAN_F4   0.0222222 0.3801648 964  0.058454  0.9534
    ## distance_to_yeast25:YeastZAN_F4   0.4494444 0.3801648 964  1.182236  0.2374
    ## distance_to_yeast32:YeastZAN_F4   0.0416667 0.3801648 964  0.109602  0.9127
    ## distance_to_yeast41:YeastZAN_F4  -0.1150000 0.3801648 964 -0.302500  0.7623
    ## distance_to_yeast48:YeastZAN_F4  -0.1738889 0.3801648 964 -0.457404  0.6475
    ## distance_to_yeast55:YeastZAN_F4  -0.6094444 0.3801648 964 -1.603106  0.1092
    ## DAI6:YeastEMM_F3                  0.2669048 0.2488763 964  1.072440  0.2838
    ## DAI8:YeastEMM_F3                  0.4011905 0.2488763 964  1.612008  0.1073
    ## DAI6:YeastEMM_F34                 0.2154762 0.2488763 964  0.865796  0.3868
    ## DAI8:YeastEMM_F34                 0.2926190 0.2488763 964  1.175761  0.2400
    ## DAI6:YeastEMM_F47                -0.2816667 0.2488763 964 -1.131754  0.2580
    ## DAI8:YeastEMM_F47                -0.3250000 0.2488763 964 -1.305870  0.1919
    ## DAI6:YeastEMM_F48                -0.2288095 0.2488763 964 -0.919371  0.3581
    ## DAI8:YeastEMM_F48                -0.2550000 0.2488763 964 -1.024605  0.3058
    ## DAI6:YeastEMM_F49                -0.3311905 0.2488763 964 -1.330743  0.1836
    ## DAI8:YeastEMM_F49                -0.3892857 0.2488763 964 -1.564174  0.1181
    ## DAI6:YeastEMM_F5                  0.2830952 0.2488763 964  1.137494  0.2556
    ## DAI8:YeastEMM_F5                  0.2668967 0.2511741 964  1.062596  0.2882
    ## DAI6:YeastEMM_F63                 0.5030952 0.2488763 964  2.021467  0.0435
    ## DAI8:YeastEMM_F63                 0.3902381 0.2488763 964  1.568000  0.1172
    ## DAI6:YeastEMM_F64                 0.4283333 0.2488763 964  1.721069  0.0856
    ## DAI8:YeastEMM_F64                 0.3459524 0.2488763 964  1.390058  0.1648
    ## DAI6:YeastEMM_F65                -0.2135714 0.2488763 964 -0.858143  0.3910
    ## DAI8:YeastEMM_F65                -0.0945238 0.2488763 964 -0.379802  0.7042
    ## DAI6:YeastEMM_F66                -0.6872156 0.2691156 964 -2.553608  0.0108
    ## DAI8:YeastEMM_F66                -0.1798347 0.2691156 964 -0.668243  0.5041
    ## DAI6:YeastEMM_F70                -0.1459524 0.2488763 964 -0.586446  0.5577
    ## DAI8:YeastEMM_F70                 0.1392857 0.2488763 964  0.559658  0.5758
    ## DAI6:YeastEMM_F89                 0.1988095 0.2488763 964  0.798829  0.4246
    ## DAI8:YeastEMM_F89                 0.2183333 0.2488763 964  0.877277  0.3806
    ## DAI6:YeastF7                      0.1673810 0.2488763 964  0.672547  0.5014
    ## DAI8:YeastF7                      0.2602381 0.2488763 964  1.045652  0.2960
    ## DAI6:YeastSP_F14                 -0.6507143 0.2488763 964 -2.614609  0.0091
    ## DAI8:YeastSP_F14                 -0.5616667 0.2488763 964 -2.256811  0.0242
    ## DAI6:YeastZAN_F3                  0.3607143 0.2488763 964  1.449372  0.1476
    ## DAI8:YeastZAN_F3                  0.3102381 0.2488763 964  1.246555  0.2129
    ## DAI6:YeastZAN_F4                 -0.1969048 0.2488763 964 -0.791175  0.4290
    ## DAI8:YeastZAN_F4                 -0.3578571 0.2488763 964 -1.437892  0.1508
    ##  Correlation: 
    ##                                  (Intr) DAI6   DAI8   ds__17 ds__25 ds__32
    ## DAI6                             -0.398                                   
    ## DAI8                             -0.398  0.500                            
    ## distance_to_yeast17              -0.608  0.000  0.000                     
    ## distance_to_yeast25              -0.608  0.000  0.000  0.500              
    ## distance_to_yeast32              -0.608  0.000  0.000  0.500  0.500       
    ## distance_to_yeast41              -0.608  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast48              -0.608  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast55              -0.608  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F3                      -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F34                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F47                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F48                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F49                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F5                      -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F63                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F64                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F65                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F66                     -0.530  0.222  0.222  0.339  0.339  0.339
    ## YeastEMM_F70                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastEMM_F89                     -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastF7                          -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastSP_F14                      -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastZAN_F3                      -0.562  0.236  0.236  0.360  0.360  0.360
    ## YeastZAN_F4                      -0.562  0.236  0.236  0.360  0.360  0.360
    ## distance_to_yeast17:YeastEMM_F3   0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F3   0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F3   0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F3   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F3   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F3   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F34  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F34  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F34  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F34  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F34  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F34  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F47  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F47  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F47  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F47  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F47  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F47  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F48  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F48  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F48  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F48  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F48  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F48  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F49  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F49  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F49  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F49  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F49  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F49  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F5   0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F5   0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F5   0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F5   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F5   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F5   0.343  0.000  0.000 -0.283 -0.283 -0.283
    ## distance_to_yeast17:YeastEMM_F63  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F63  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F63  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F63  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F63  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F63  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F64  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F64  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F64  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F64  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F64  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F64  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F65  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F65  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F65  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F65  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F65  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F65  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F66  0.321  0.000  0.000 -0.529 -0.265 -0.265
    ## distance_to_yeast25:YeastEMM_F66  0.321  0.000  0.000 -0.265 -0.529 -0.265
    ## distance_to_yeast32:YeastEMM_F66  0.321  0.000  0.000 -0.265 -0.265 -0.529
    ## distance_to_yeast41:YeastEMM_F66  0.321  0.000  0.000 -0.265 -0.265 -0.265
    ## distance_to_yeast48:YeastEMM_F66  0.321  0.000  0.000 -0.265 -0.265 -0.265
    ## distance_to_yeast55:YeastEMM_F66  0.321  0.000  0.000 -0.265 -0.265 -0.265
    ## distance_to_yeast17:YeastEMM_F70  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F70  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F70  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F70  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F70  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F70  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastEMM_F89  0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastEMM_F89  0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastEMM_F89  0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastEMM_F89  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastEMM_F89  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastEMM_F89  0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastF7       0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastF7       0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastF7       0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastF7       0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastF7       0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastF7       0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastSP_F14   0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastSP_F14   0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastSP_F14   0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastSP_F14   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastSP_F14   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastSP_F14   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastZAN_F3   0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastZAN_F3   0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastZAN_F3   0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastZAN_F3   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastZAN_F3   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastZAN_F3   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast17:YeastZAN_F4   0.351  0.000  0.000 -0.577 -0.289 -0.289
    ## distance_to_yeast25:YeastZAN_F4   0.351  0.000  0.000 -0.289 -0.577 -0.289
    ## distance_to_yeast32:YeastZAN_F4   0.351  0.000  0.000 -0.289 -0.289 -0.577
    ## distance_to_yeast41:YeastZAN_F4   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast48:YeastZAN_F4   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## distance_to_yeast55:YeastZAN_F4   0.351  0.000  0.000 -0.289 -0.289 -0.289
    ## DAI6:YeastEMM_F3                  0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F3                  0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F34                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F34                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F47                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F47                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F48                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F48                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F49                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F49                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F5                  0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F5                  0.228 -0.286 -0.572  0.000  0.000  0.000
    ## DAI6:YeastEMM_F63                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F63                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F64                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F64                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F65                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F65                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F66                 0.212 -0.534 -0.267  0.000  0.000  0.000
    ## DAI8:YeastEMM_F66                 0.212 -0.267 -0.534  0.000  0.000  0.000
    ## DAI6:YeastEMM_F70                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F70                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastEMM_F89                 0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastEMM_F89                 0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastF7                      0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastF7                      0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastSP_F14                  0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastSP_F14                  0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastZAN_F3                  0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastZAN_F3                  0.230 -0.289 -0.577  0.000  0.000  0.000
    ## DAI6:YeastZAN_F4                  0.230 -0.577 -0.289  0.000  0.000  0.000
    ## DAI8:YeastZAN_F4                  0.230 -0.289 -0.577  0.000  0.000  0.000
    ##                                  ds__41 ds__48 ds__55 YsEMM_F3 YEMM_F34
    ## DAI6                                                                   
    ## DAI8                                                                   
    ## distance_to_yeast17                                                    
    ## distance_to_yeast25                                                    
    ## distance_to_yeast32                                                    
    ## distance_to_yeast41                                                    
    ## distance_to_yeast48               0.500                                
    ## distance_to_yeast55               0.500  0.500                         
    ## YeastEMM_F3                       0.360  0.360  0.360                  
    ## YeastEMM_F34                      0.360  0.360  0.360  0.333           
    ## YeastEMM_F47                      0.360  0.360  0.360  0.333    0.333  
    ## YeastEMM_F48                      0.360  0.360  0.360  0.333    0.333  
    ## YeastEMM_F49                      0.360  0.360  0.360  0.333    0.333  
    ## YeastEMM_F5                       0.360  0.360  0.360  0.333    0.333  
    ## YeastEMM_F63                      0.360  0.360  0.360  0.333    0.333  
    ## YeastEMM_F64                      0.360  0.360  0.360  0.333    0.333  
    ## YeastEMM_F65                      0.360  0.360  0.360  0.333    0.333  
    ## YeastEMM_F66                      0.339  0.339  0.339  0.314    0.314  
    ## YeastEMM_F70                      0.360  0.360  0.360  0.333    0.333  
    ## YeastEMM_F89                      0.360  0.360  0.360  0.333    0.333  
    ## YeastF7                           0.360  0.360  0.360  0.333    0.333  
    ## YeastSP_F14                       0.360  0.360  0.360  0.333    0.333  
    ## YeastZAN_F3                       0.360  0.360  0.360  0.333    0.333  
    ## YeastZAN_F4                       0.360  0.360  0.360  0.333    0.333  
    ## distance_to_yeast17:YeastEMM_F3  -0.289 -0.289 -0.289 -0.624   -0.208  
    ## distance_to_yeast25:YeastEMM_F3  -0.289 -0.289 -0.289 -0.624   -0.208  
    ## distance_to_yeast32:YeastEMM_F3  -0.289 -0.289 -0.289 -0.624   -0.208  
    ## distance_to_yeast41:YeastEMM_F3  -0.577 -0.289 -0.289 -0.624   -0.208  
    ## distance_to_yeast48:YeastEMM_F3  -0.289 -0.577 -0.289 -0.624   -0.208  
    ## distance_to_yeast55:YeastEMM_F3  -0.289 -0.289 -0.577 -0.624   -0.208  
    ## distance_to_yeast17:YeastEMM_F34 -0.289 -0.289 -0.289 -0.208   -0.624  
    ## distance_to_yeast25:YeastEMM_F34 -0.289 -0.289 -0.289 -0.208   -0.624  
    ## distance_to_yeast32:YeastEMM_F34 -0.289 -0.289 -0.289 -0.208   -0.624  
    ## distance_to_yeast41:YeastEMM_F34 -0.577 -0.289 -0.289 -0.208   -0.624  
    ## distance_to_yeast48:YeastEMM_F34 -0.289 -0.577 -0.289 -0.208   -0.624  
    ## distance_to_yeast55:YeastEMM_F34 -0.289 -0.289 -0.577 -0.208   -0.624  
    ## distance_to_yeast17:YeastEMM_F47 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F47 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F47 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F47 -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F47 -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F47 -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastEMM_F48 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F48 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F48 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F48 -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F48 -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F48 -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastEMM_F49 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F49 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F49 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F49 -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F49 -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F49 -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastEMM_F5  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F5  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F5  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F5  -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F5  -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F5  -0.283 -0.283 -0.565 -0.204   -0.204  
    ## distance_to_yeast17:YeastEMM_F63 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F63 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F63 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F63 -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F63 -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F63 -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastEMM_F64 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F64 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F64 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F64 -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F64 -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F64 -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastEMM_F65 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F65 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F65 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F65 -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F65 -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F65 -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastEMM_F66 -0.265 -0.265 -0.265 -0.191   -0.191  
    ## distance_to_yeast25:YeastEMM_F66 -0.265 -0.265 -0.265 -0.191   -0.191  
    ## distance_to_yeast32:YeastEMM_F66 -0.265 -0.265 -0.265 -0.191   -0.191  
    ## distance_to_yeast41:YeastEMM_F66 -0.529 -0.265 -0.265 -0.191   -0.191  
    ## distance_to_yeast48:YeastEMM_F66 -0.265 -0.529 -0.265 -0.191   -0.191  
    ## distance_to_yeast55:YeastEMM_F66 -0.265 -0.265 -0.529 -0.191   -0.191  
    ## distance_to_yeast17:YeastEMM_F70 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F70 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F70 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F70 -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F70 -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F70 -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastEMM_F89 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastEMM_F89 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastEMM_F89 -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastEMM_F89 -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastEMM_F89 -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastEMM_F89 -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastF7      -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastF7      -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastF7      -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastF7      -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastF7      -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastF7      -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastSP_F14  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastSP_F14  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastSP_F14  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastSP_F14  -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastSP_F14  -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastSP_F14  -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastZAN_F3  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastZAN_F3  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastZAN_F3  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastZAN_F3  -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastZAN_F3  -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastZAN_F3  -0.289 -0.289 -0.577 -0.208   -0.208  
    ## distance_to_yeast17:YeastZAN_F4  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast25:YeastZAN_F4  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast32:YeastZAN_F4  -0.289 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast41:YeastZAN_F4  -0.577 -0.289 -0.289 -0.208   -0.208  
    ## distance_to_yeast48:YeastZAN_F4  -0.289 -0.577 -0.289 -0.208   -0.208  
    ## distance_to_yeast55:YeastZAN_F4  -0.289 -0.289 -0.577 -0.208   -0.208  
    ## DAI6:YeastEMM_F3                  0.000  0.000  0.000 -0.408   -0.136  
    ## DAI8:YeastEMM_F3                  0.000  0.000  0.000 -0.408   -0.136  
    ## DAI6:YeastEMM_F34                 0.000  0.000  0.000 -0.136   -0.408  
    ## DAI8:YeastEMM_F34                 0.000  0.000  0.000 -0.136   -0.408  
    ## DAI6:YeastEMM_F47                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F47                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastEMM_F48                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F48                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastEMM_F49                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F49                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastEMM_F5                  0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F5                  0.000  0.000  0.000 -0.135   -0.135  
    ## DAI6:YeastEMM_F63                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F63                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastEMM_F64                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F64                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastEMM_F65                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F65                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastEMM_F66                 0.000  0.000  0.000 -0.126   -0.126  
    ## DAI8:YeastEMM_F66                 0.000  0.000  0.000 -0.126   -0.126  
    ## DAI6:YeastEMM_F70                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F70                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastEMM_F89                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastEMM_F89                 0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastF7                      0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastF7                      0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastSP_F14                  0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastSP_F14                  0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastZAN_F3                  0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastZAN_F3                  0.000  0.000  0.000 -0.136   -0.136  
    ## DAI6:YeastZAN_F4                  0.000  0.000  0.000 -0.136   -0.136  
    ## DAI8:YeastZAN_F4                  0.000  0.000  0.000 -0.136   -0.136  
    ##                                  YEMM_F47 YEMM_F48 YEMM_F49 YEMM_F5 YEMM_F63
    ## DAI6                                                                        
    ## DAI8                                                                        
    ## distance_to_yeast17                                                         
    ## distance_to_yeast25                                                         
    ## distance_to_yeast32                                                         
    ## distance_to_yeast41                                                         
    ## distance_to_yeast48                                                         
    ## distance_to_yeast55                                                         
    ## YeastEMM_F3                                                                 
    ## YeastEMM_F34                                                                
    ## YeastEMM_F47                                                                
    ## YeastEMM_F48                      0.333                                     
    ## YeastEMM_F49                      0.333    0.333                            
    ## YeastEMM_F5                       0.333    0.333    0.333                   
    ## YeastEMM_F63                      0.333    0.333    0.333    0.333          
    ## YeastEMM_F64                      0.333    0.333    0.333    0.333   0.333  
    ## YeastEMM_F65                      0.333    0.333    0.333    0.333   0.333  
    ## YeastEMM_F66                      0.314    0.314    0.314    0.314   0.314  
    ## YeastEMM_F70                      0.333    0.333    0.333    0.333   0.333  
    ## YeastEMM_F89                      0.333    0.333    0.333    0.333   0.333  
    ## YeastF7                           0.333    0.333    0.333    0.333   0.333  
    ## YeastSP_F14                       0.333    0.333    0.333    0.333   0.333  
    ## YeastZAN_F3                       0.333    0.333    0.333    0.333   0.333  
    ## YeastZAN_F4                       0.333    0.333    0.333    0.333   0.333  
    ## distance_to_yeast17:YeastEMM_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastEMM_F34 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F34 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F34 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F34 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F34 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F34 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastEMM_F47 -0.624   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F47 -0.624   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F47 -0.624   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F47 -0.624   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F47 -0.624   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F47 -0.624   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastEMM_F48 -0.208   -0.624   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F48 -0.208   -0.624   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F48 -0.208   -0.624   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F48 -0.208   -0.624   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F48 -0.208   -0.624   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F48 -0.208   -0.624   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastEMM_F49 -0.208   -0.208   -0.624   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F49 -0.208   -0.208   -0.624   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F49 -0.208   -0.208   -0.624   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F49 -0.208   -0.208   -0.624   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F49 -0.208   -0.208   -0.624   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F49 -0.208   -0.208   -0.624   -0.208  -0.208  
    ## distance_to_yeast17:YeastEMM_F5  -0.208   -0.208   -0.208   -0.623  -0.208  
    ## distance_to_yeast25:YeastEMM_F5  -0.208   -0.208   -0.208   -0.623  -0.208  
    ## distance_to_yeast32:YeastEMM_F5  -0.208   -0.208   -0.208   -0.623  -0.208  
    ## distance_to_yeast41:YeastEMM_F5  -0.208   -0.208   -0.208   -0.623  -0.208  
    ## distance_to_yeast48:YeastEMM_F5  -0.208   -0.208   -0.208   -0.623  -0.208  
    ## distance_to_yeast55:YeastEMM_F5  -0.204   -0.204   -0.204   -0.618  -0.204  
    ## distance_to_yeast17:YeastEMM_F63 -0.208   -0.208   -0.208   -0.208  -0.624  
    ## distance_to_yeast25:YeastEMM_F63 -0.208   -0.208   -0.208   -0.208  -0.624  
    ## distance_to_yeast32:YeastEMM_F63 -0.208   -0.208   -0.208   -0.208  -0.624  
    ## distance_to_yeast41:YeastEMM_F63 -0.208   -0.208   -0.208   -0.208  -0.624  
    ## distance_to_yeast48:YeastEMM_F63 -0.208   -0.208   -0.208   -0.208  -0.624  
    ## distance_to_yeast55:YeastEMM_F63 -0.208   -0.208   -0.208   -0.208  -0.624  
    ## distance_to_yeast17:YeastEMM_F64 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F64 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F64 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F64 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F64 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F64 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastEMM_F65 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F65 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F65 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F65 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F65 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F65 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastEMM_F66 -0.191   -0.191   -0.191   -0.190  -0.191  
    ## distance_to_yeast25:YeastEMM_F66 -0.191   -0.191   -0.191   -0.190  -0.191  
    ## distance_to_yeast32:YeastEMM_F66 -0.191   -0.191   -0.191   -0.190  -0.191  
    ## distance_to_yeast41:YeastEMM_F66 -0.191   -0.191   -0.191   -0.190  -0.191  
    ## distance_to_yeast48:YeastEMM_F66 -0.191   -0.191   -0.191   -0.190  -0.191  
    ## distance_to_yeast55:YeastEMM_F66 -0.191   -0.191   -0.191   -0.190  -0.191  
    ## distance_to_yeast17:YeastEMM_F70 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F70 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F70 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F70 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F70 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F70 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastEMM_F89 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastEMM_F89 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastEMM_F89 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastEMM_F89 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastEMM_F89 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastEMM_F89 -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastF7      -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastF7      -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastF7      -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastF7      -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastF7      -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastF7      -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastSP_F14  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastSP_F14  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastSP_F14  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastSP_F14  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastSP_F14  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastSP_F14  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastZAN_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastZAN_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastZAN_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastZAN_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastZAN_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastZAN_F3  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast17:YeastZAN_F4  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast25:YeastZAN_F4  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast32:YeastZAN_F4  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast41:YeastZAN_F4  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast48:YeastZAN_F4  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## distance_to_yeast55:YeastZAN_F4  -0.208   -0.208   -0.208   -0.208  -0.208  
    ## DAI6:YeastEMM_F3                 -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastEMM_F3                 -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastEMM_F34                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastEMM_F34                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastEMM_F47                -0.408   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastEMM_F47                -0.408   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastEMM_F48                -0.136   -0.408   -0.136   -0.136  -0.136  
    ## DAI8:YeastEMM_F48                -0.136   -0.408   -0.136   -0.136  -0.136  
    ## DAI6:YeastEMM_F49                -0.136   -0.136   -0.408   -0.136  -0.136  
    ## DAI8:YeastEMM_F49                -0.136   -0.136   -0.408   -0.136  -0.136  
    ## DAI6:YeastEMM_F5                 -0.136   -0.136   -0.136   -0.408  -0.136  
    ## DAI8:YeastEMM_F5                 -0.135   -0.135   -0.135   -0.409  -0.135  
    ## DAI6:YeastEMM_F63                -0.136   -0.136   -0.136   -0.136  -0.408  
    ## DAI8:YeastEMM_F63                -0.136   -0.136   -0.136   -0.136  -0.408  
    ## DAI6:YeastEMM_F64                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastEMM_F64                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastEMM_F65                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastEMM_F65                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastEMM_F66                -0.126   -0.126   -0.126   -0.126  -0.126  
    ## DAI8:YeastEMM_F66                -0.126   -0.126   -0.126   -0.126  -0.126  
    ## DAI6:YeastEMM_F70                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastEMM_F70                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastEMM_F89                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastEMM_F89                -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastF7                     -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastF7                     -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastSP_F14                 -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastSP_F14                 -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastZAN_F3                 -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastZAN_F3                 -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI6:YeastZAN_F4                 -0.136   -0.136   -0.136   -0.136  -0.136  
    ## DAI8:YeastZAN_F4                 -0.136   -0.136   -0.136   -0.136  -0.136  
    ##                                  YEMM_F64 YEMM_F65 YEMM_F66 YEMM_F7 YEMM_F8
    ## DAI6                                                                       
    ## DAI8                                                                       
    ## distance_to_yeast17                                                        
    ## distance_to_yeast25                                                        
    ## distance_to_yeast32                                                        
    ## distance_to_yeast41                                                        
    ## distance_to_yeast48                                                        
    ## distance_to_yeast55                                                        
    ## YeastEMM_F3                                                                
    ## YeastEMM_F34                                                               
    ## YeastEMM_F47                                                               
    ## YeastEMM_F48                                                               
    ## YeastEMM_F49                                                               
    ## YeastEMM_F5                                                                
    ## YeastEMM_F63                                                               
    ## YeastEMM_F64                                                               
    ## YeastEMM_F65                      0.333                                    
    ## YeastEMM_F66                      0.314    0.314                           
    ## YeastEMM_F70                      0.333    0.333    0.314                  
    ## YeastEMM_F89                      0.333    0.333    0.314    0.333         
    ## YeastF7                           0.333    0.333    0.314    0.333   0.333 
    ## YeastSP_F14                       0.333    0.333    0.314    0.333   0.333 
    ## YeastZAN_F3                       0.333    0.333    0.314    0.333   0.333 
    ## YeastZAN_F4                       0.333    0.333    0.314    0.333   0.333 
    ## distance_to_yeast17:YeastEMM_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastEMM_F34 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F34 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F34 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F34 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F34 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F34 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastEMM_F47 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F47 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F47 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F47 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F47 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F47 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastEMM_F48 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F48 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F48 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F48 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F48 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F48 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastEMM_F49 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F49 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F49 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F49 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F49 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F49 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastEMM_F5  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F5  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F5  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F5  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F5  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F5  -0.204   -0.204   -0.192   -0.204  -0.204 
    ## distance_to_yeast17:YeastEMM_F63 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F63 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F63 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F63 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F63 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F63 -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastEMM_F64 -0.624   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F64 -0.624   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F64 -0.624   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F64 -0.624   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F64 -0.624   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F64 -0.624   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastEMM_F65 -0.208   -0.624   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastEMM_F65 -0.208   -0.624   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastEMM_F65 -0.208   -0.624   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastEMM_F65 -0.208   -0.624   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastEMM_F65 -0.208   -0.624   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastEMM_F65 -0.208   -0.624   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastEMM_F66 -0.191   -0.191   -0.641   -0.191  -0.191 
    ## distance_to_yeast25:YeastEMM_F66 -0.191   -0.191   -0.641   -0.191  -0.191 
    ## distance_to_yeast32:YeastEMM_F66 -0.191   -0.191   -0.641   -0.191  -0.191 
    ## distance_to_yeast41:YeastEMM_F66 -0.191   -0.191   -0.641   -0.191  -0.191 
    ## distance_to_yeast48:YeastEMM_F66 -0.191   -0.191   -0.641   -0.191  -0.191 
    ## distance_to_yeast55:YeastEMM_F66 -0.191   -0.191   -0.641   -0.191  -0.191 
    ## distance_to_yeast17:YeastEMM_F70 -0.208   -0.208   -0.196   -0.624  -0.208 
    ## distance_to_yeast25:YeastEMM_F70 -0.208   -0.208   -0.196   -0.624  -0.208 
    ## distance_to_yeast32:YeastEMM_F70 -0.208   -0.208   -0.196   -0.624  -0.208 
    ## distance_to_yeast41:YeastEMM_F70 -0.208   -0.208   -0.196   -0.624  -0.208 
    ## distance_to_yeast48:YeastEMM_F70 -0.208   -0.208   -0.196   -0.624  -0.208 
    ## distance_to_yeast55:YeastEMM_F70 -0.208   -0.208   -0.196   -0.624  -0.208 
    ## distance_to_yeast17:YeastEMM_F89 -0.208   -0.208   -0.196   -0.208  -0.624 
    ## distance_to_yeast25:YeastEMM_F89 -0.208   -0.208   -0.196   -0.208  -0.624 
    ## distance_to_yeast32:YeastEMM_F89 -0.208   -0.208   -0.196   -0.208  -0.624 
    ## distance_to_yeast41:YeastEMM_F89 -0.208   -0.208   -0.196   -0.208  -0.624 
    ## distance_to_yeast48:YeastEMM_F89 -0.208   -0.208   -0.196   -0.208  -0.624 
    ## distance_to_yeast55:YeastEMM_F89 -0.208   -0.208   -0.196   -0.208  -0.624 
    ## distance_to_yeast17:YeastF7      -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastF7      -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastF7      -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastF7      -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastF7      -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastF7      -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastSP_F14  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastSP_F14  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastSP_F14  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastSP_F14  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastSP_F14  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastSP_F14  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastZAN_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastZAN_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastZAN_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastZAN_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastZAN_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastZAN_F3  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast17:YeastZAN_F4  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast25:YeastZAN_F4  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast32:YeastZAN_F4  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast41:YeastZAN_F4  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast48:YeastZAN_F4  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## distance_to_yeast55:YeastZAN_F4  -0.208   -0.208   -0.196   -0.208  -0.208 
    ## DAI6:YeastEMM_F3                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F3                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastEMM_F34                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F34                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastEMM_F47                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F47                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastEMM_F48                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F48                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastEMM_F49                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F49                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastEMM_F5                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F5                 -0.135   -0.135   -0.127   -0.135  -0.135 
    ## DAI6:YeastEMM_F63                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F63                -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastEMM_F64                -0.408   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F64                -0.408   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastEMM_F65                -0.136   -0.408   -0.128   -0.136  -0.136 
    ## DAI8:YeastEMM_F65                -0.136   -0.408   -0.128   -0.136  -0.136 
    ## DAI6:YeastEMM_F66                -0.126   -0.126   -0.356   -0.126  -0.126 
    ## DAI8:YeastEMM_F66                -0.126   -0.126   -0.356   -0.126  -0.126 
    ## DAI6:YeastEMM_F70                -0.136   -0.136   -0.128   -0.408  -0.136 
    ## DAI8:YeastEMM_F70                -0.136   -0.136   -0.128   -0.408  -0.136 
    ## DAI6:YeastEMM_F89                -0.136   -0.136   -0.128   -0.136  -0.408 
    ## DAI8:YeastEMM_F89                -0.136   -0.136   -0.128   -0.136  -0.408 
    ## DAI6:YeastF7                     -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastF7                     -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastSP_F14                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastSP_F14                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastZAN_F3                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastZAN_F3                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI6:YeastZAN_F4                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ## DAI8:YeastZAN_F4                 -0.136   -0.136   -0.128   -0.136  -0.136 
    ##                                  YestF7 YSP_F1 YZAN_F3 YZAN_F4 ds__17:YEMM_F3
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                       0.333                                      
    ## YeastZAN_F3                       0.333  0.333                               
    ## YeastZAN_F4                       0.333  0.333  0.333                        
    ## distance_to_yeast17:YeastEMM_F3  -0.208 -0.208 -0.208  -0.208                
    ## distance_to_yeast25:YeastEMM_F3  -0.208 -0.208 -0.208  -0.208   0.500        
    ## distance_to_yeast32:YeastEMM_F3  -0.208 -0.208 -0.208  -0.208   0.500        
    ## distance_to_yeast41:YeastEMM_F3  -0.208 -0.208 -0.208  -0.208   0.500        
    ## distance_to_yeast48:YeastEMM_F3  -0.208 -0.208 -0.208  -0.208   0.500        
    ## distance_to_yeast55:YeastEMM_F3  -0.208 -0.208 -0.208  -0.208   0.500        
    ## distance_to_yeast17:YeastEMM_F34 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F34 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F34 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F34 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F34 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F34 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastEMM_F47 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F47 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F47 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F47 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F47 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F47 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastEMM_F48 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F48 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F48 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F48 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F48 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F48 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastEMM_F49 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F49 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F49 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F49 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F49 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F49 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastEMM_F5  -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F5  -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F5  -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F5  -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F5  -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F5  -0.204 -0.204 -0.204  -0.204   0.163        
    ## distance_to_yeast17:YeastEMM_F63 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F63 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F63 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F63 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F63 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F63 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastEMM_F64 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F64 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F64 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F64 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F64 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F64 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastEMM_F65 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F65 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F65 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F65 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F65 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F65 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastEMM_F66 -0.191 -0.191 -0.191  -0.191   0.306        
    ## distance_to_yeast25:YeastEMM_F66 -0.191 -0.191 -0.191  -0.191   0.153        
    ## distance_to_yeast32:YeastEMM_F66 -0.191 -0.191 -0.191  -0.191   0.153        
    ## distance_to_yeast41:YeastEMM_F66 -0.191 -0.191 -0.191  -0.191   0.153        
    ## distance_to_yeast48:YeastEMM_F66 -0.191 -0.191 -0.191  -0.191   0.153        
    ## distance_to_yeast55:YeastEMM_F66 -0.191 -0.191 -0.191  -0.191   0.153        
    ## distance_to_yeast17:YeastEMM_F70 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F70 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F70 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F70 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F70 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F70 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastEMM_F89 -0.208 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastEMM_F89 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastEMM_F89 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastEMM_F89 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastEMM_F89 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastEMM_F89 -0.208 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastF7      -0.624 -0.208 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastF7      -0.624 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastF7      -0.624 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastF7      -0.624 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastF7      -0.624 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastF7      -0.624 -0.208 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastSP_F14  -0.208 -0.624 -0.208  -0.208   0.333        
    ## distance_to_yeast25:YeastSP_F14  -0.208 -0.624 -0.208  -0.208   0.167        
    ## distance_to_yeast32:YeastSP_F14  -0.208 -0.624 -0.208  -0.208   0.167        
    ## distance_to_yeast41:YeastSP_F14  -0.208 -0.624 -0.208  -0.208   0.167        
    ## distance_to_yeast48:YeastSP_F14  -0.208 -0.624 -0.208  -0.208   0.167        
    ## distance_to_yeast55:YeastSP_F14  -0.208 -0.624 -0.208  -0.208   0.167        
    ## distance_to_yeast17:YeastZAN_F3  -0.208 -0.208 -0.624  -0.208   0.333        
    ## distance_to_yeast25:YeastZAN_F3  -0.208 -0.208 -0.624  -0.208   0.167        
    ## distance_to_yeast32:YeastZAN_F3  -0.208 -0.208 -0.624  -0.208   0.167        
    ## distance_to_yeast41:YeastZAN_F3  -0.208 -0.208 -0.624  -0.208   0.167        
    ## distance_to_yeast48:YeastZAN_F3  -0.208 -0.208 -0.624  -0.208   0.167        
    ## distance_to_yeast55:YeastZAN_F3  -0.208 -0.208 -0.624  -0.208   0.167        
    ## distance_to_yeast17:YeastZAN_F4  -0.208 -0.208 -0.208  -0.624   0.333        
    ## distance_to_yeast25:YeastZAN_F4  -0.208 -0.208 -0.208  -0.624   0.167        
    ## distance_to_yeast32:YeastZAN_F4  -0.208 -0.208 -0.208  -0.624   0.167        
    ## distance_to_yeast41:YeastZAN_F4  -0.208 -0.208 -0.208  -0.624   0.167        
    ## distance_to_yeast48:YeastZAN_F4  -0.208 -0.208 -0.208  -0.624   0.167        
    ## distance_to_yeast55:YeastZAN_F4  -0.208 -0.208 -0.208  -0.624   0.167        
    ## DAI6:YeastEMM_F3                 -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F3                 -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F34                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F34                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F47                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F47                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F48                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F48                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F49                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F49                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F5                 -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F5                 -0.135 -0.135 -0.135  -0.135   0.000        
    ## DAI6:YeastEMM_F63                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F63                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F64                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F64                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F65                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F65                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F66                -0.126 -0.126 -0.126  -0.126   0.000        
    ## DAI8:YeastEMM_F66                -0.126 -0.126 -0.126  -0.126   0.000        
    ## DAI6:YeastEMM_F70                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F70                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastEMM_F89                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastEMM_F89                -0.136 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastF7                     -0.408 -0.136 -0.136  -0.136   0.000        
    ## DAI8:YeastF7                     -0.408 -0.136 -0.136  -0.136   0.000        
    ## DAI6:YeastSP_F14                 -0.136 -0.408 -0.136  -0.136   0.000        
    ## DAI8:YeastSP_F14                 -0.136 -0.408 -0.136  -0.136   0.000        
    ## DAI6:YeastZAN_F3                 -0.136 -0.136 -0.408  -0.136   0.000        
    ## DAI8:YeastZAN_F3                 -0.136 -0.136 -0.408  -0.136   0.000        
    ## DAI6:YeastZAN_F4                 -0.136 -0.136 -0.136  -0.408   0.000        
    ## DAI8:YeastZAN_F4                 -0.136 -0.136 -0.136  -0.408   0.000        
    ##                                  ds__25:YEMM_F3 ds__32:YEMM_F3 ds__41:YEMM_F3
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3   0.500                                      
    ## distance_to_yeast41:YeastEMM_F3   0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F3   0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F3   0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F34  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F34  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F34  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F34  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F34  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F34  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F47  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F47  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F47  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F48  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F48  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F48  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F49  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F49  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F49  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.163          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast25:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast48:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.167          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  ds__48:YEMM_F3 ds__55:YEMM_F3 d__17:YEMM_F34
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3   0.500                                      
    ## distance_to_yeast17:YeastEMM_F34  0.167          0.167                       
    ## distance_to_yeast25:YeastEMM_F34  0.167          0.167          0.500        
    ## distance_to_yeast32:YeastEMM_F34  0.167          0.167          0.500        
    ## distance_to_yeast41:YeastEMM_F34  0.167          0.167          0.500        
    ## distance_to_yeast48:YeastEMM_F34  0.333          0.167          0.500        
    ## distance_to_yeast55:YeastEMM_F34  0.167          0.333          0.500        
    ## distance_to_yeast17:YeastEMM_F47  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F47  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F47  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F48  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F48  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F48  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F49  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F49  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F49  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.326          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast25:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast48:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast25:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast48:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.333          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__25:YEMM_F34 d__32:YEMM_F34 d__41:YEMM_F34
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34  0.500                                      
    ## distance_to_yeast41:YeastEMM_F34  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F34  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F34  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F47  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F47  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F47  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F47  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F48  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F48  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F48  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F49  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F49  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F49  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.163          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast25:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast48:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.167          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__48:YEMM_F34 d__55:YEMM_F34 d__17:YEMM_F47
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34  0.500                                      
    ## distance_to_yeast17:YeastEMM_F47  0.167          0.167                       
    ## distance_to_yeast25:YeastEMM_F47  0.167          0.167          0.500        
    ## distance_to_yeast32:YeastEMM_F47  0.167          0.167          0.500        
    ## distance_to_yeast41:YeastEMM_F47  0.167          0.167          0.500        
    ## distance_to_yeast48:YeastEMM_F47  0.333          0.167          0.500        
    ## distance_to_yeast55:YeastEMM_F47  0.167          0.333          0.500        
    ## distance_to_yeast17:YeastEMM_F48  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F48  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F48  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F49  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F49  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F49  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.326          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast25:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast48:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast25:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast48:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.333          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__25:YEMM_F47 d__32:YEMM_F47 d__41:YEMM_F47
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47  0.500                                      
    ## distance_to_yeast41:YeastEMM_F47  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F47  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F47  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F48  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F48  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F48  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F48  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F49  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F49  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F49  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.163          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast25:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast48:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.167          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__48:YEMM_F47 d__55:YEMM_F47 d__17:YEMM_F48
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47  0.500                                      
    ## distance_to_yeast17:YeastEMM_F48  0.167          0.167                       
    ## distance_to_yeast25:YeastEMM_F48  0.167          0.167          0.500        
    ## distance_to_yeast32:YeastEMM_F48  0.167          0.167          0.500        
    ## distance_to_yeast41:YeastEMM_F48  0.167          0.167          0.500        
    ## distance_to_yeast48:YeastEMM_F48  0.333          0.167          0.500        
    ## distance_to_yeast55:YeastEMM_F48  0.167          0.333          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F49  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F49  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.326          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast25:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast48:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast25:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast48:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.333          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__25:YEMM_F48 d__32:YEMM_F48 d__41:YEMM_F48
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48  0.500                                      
    ## distance_to_yeast41:YeastEMM_F48  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F48  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F48  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F49  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F49  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F49  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F49  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.163          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast25:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast48:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.167          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__48:YEMM_F48 d__55:YEMM_F48 d__17:YEMM_F49
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48  0.500                                      
    ## distance_to_yeast17:YeastEMM_F49  0.167          0.167                       
    ## distance_to_yeast25:YeastEMM_F49  0.167          0.167          0.500        
    ## distance_to_yeast32:YeastEMM_F49  0.167          0.167          0.500        
    ## distance_to_yeast41:YeastEMM_F49  0.167          0.167          0.500        
    ## distance_to_yeast48:YeastEMM_F49  0.333          0.167          0.500        
    ## distance_to_yeast55:YeastEMM_F49  0.167          0.333          0.500        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.326          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast25:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast48:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast25:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast48:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.333          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__25:YEMM_F49 d__32:YEMM_F49 d__41:YEMM_F49
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49  0.500                                      
    ## distance_to_yeast41:YeastEMM_F49  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F49  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F49  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F5   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F5   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.163          0.163        
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F63  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast25:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast48:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.167          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__48:YEMM_F49 d__55:YEMM_F49 d__17:YEMM_F5
    ## DAI6                                                                        
    ## DAI8                                                                        
    ## distance_to_yeast17                                                         
    ## distance_to_yeast25                                                         
    ## distance_to_yeast32                                                         
    ## distance_to_yeast41                                                         
    ## distance_to_yeast48                                                         
    ## distance_to_yeast55                                                         
    ## YeastEMM_F3                                                                 
    ## YeastEMM_F34                                                                
    ## YeastEMM_F47                                                                
    ## YeastEMM_F48                                                                
    ## YeastEMM_F49                                                                
    ## YeastEMM_F5                                                                 
    ## YeastEMM_F63                                                                
    ## YeastEMM_F64                                                                
    ## YeastEMM_F65                                                                
    ## YeastEMM_F66                                                                
    ## YeastEMM_F70                                                                
    ## YeastEMM_F89                                                                
    ## YeastF7                                                                     
    ## YeastSP_F14                                                                 
    ## YeastZAN_F3                                                                 
    ## YeastZAN_F4                                                                 
    ## distance_to_yeast17:YeastEMM_F3                                             
    ## distance_to_yeast25:YeastEMM_F3                                             
    ## distance_to_yeast32:YeastEMM_F3                                             
    ## distance_to_yeast41:YeastEMM_F3                                             
    ## distance_to_yeast48:YeastEMM_F3                                             
    ## distance_to_yeast55:YeastEMM_F3                                             
    ## distance_to_yeast17:YeastEMM_F34                                            
    ## distance_to_yeast25:YeastEMM_F34                                            
    ## distance_to_yeast32:YeastEMM_F34                                            
    ## distance_to_yeast41:YeastEMM_F34                                            
    ## distance_to_yeast48:YeastEMM_F34                                            
    ## distance_to_yeast55:YeastEMM_F34                                            
    ## distance_to_yeast17:YeastEMM_F47                                            
    ## distance_to_yeast25:YeastEMM_F47                                            
    ## distance_to_yeast32:YeastEMM_F47                                            
    ## distance_to_yeast41:YeastEMM_F47                                            
    ## distance_to_yeast48:YeastEMM_F47                                            
    ## distance_to_yeast55:YeastEMM_F47                                            
    ## distance_to_yeast17:YeastEMM_F48                                            
    ## distance_to_yeast25:YeastEMM_F48                                            
    ## distance_to_yeast32:YeastEMM_F48                                            
    ## distance_to_yeast41:YeastEMM_F48                                            
    ## distance_to_yeast48:YeastEMM_F48                                            
    ## distance_to_yeast55:YeastEMM_F48                                            
    ## distance_to_yeast17:YeastEMM_F49                                            
    ## distance_to_yeast25:YeastEMM_F49                                            
    ## distance_to_yeast32:YeastEMM_F49                                            
    ## distance_to_yeast41:YeastEMM_F49                                            
    ## distance_to_yeast48:YeastEMM_F49                                            
    ## distance_to_yeast55:YeastEMM_F49  0.500                                     
    ## distance_to_yeast17:YeastEMM_F5   0.167          0.167                      
    ## distance_to_yeast25:YeastEMM_F5   0.167          0.167          0.500       
    ## distance_to_yeast32:YeastEMM_F5   0.167          0.167          0.500       
    ## distance_to_yeast41:YeastEMM_F5   0.167          0.167          0.500       
    ## distance_to_yeast48:YeastEMM_F5   0.333          0.167          0.500       
    ## distance_to_yeast55:YeastEMM_F5   0.163          0.326          0.490       
    ## distance_to_yeast17:YeastEMM_F63  0.167          0.167          0.333       
    ## distance_to_yeast25:YeastEMM_F63  0.167          0.167          0.167       
    ## distance_to_yeast32:YeastEMM_F63  0.167          0.167          0.167       
    ## distance_to_yeast41:YeastEMM_F63  0.167          0.167          0.167       
    ## distance_to_yeast48:YeastEMM_F63  0.333          0.167          0.167       
    ## distance_to_yeast55:YeastEMM_F63  0.167          0.333          0.167       
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.333       
    ## distance_to_yeast25:YeastEMM_F64  0.167          0.167          0.167       
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.167          0.167       
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.167       
    ## distance_to_yeast48:YeastEMM_F64  0.333          0.167          0.167       
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.333          0.167       
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.333       
    ## distance_to_yeast25:YeastEMM_F65  0.167          0.167          0.167       
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.167          0.167       
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.167       
    ## distance_to_yeast48:YeastEMM_F65  0.333          0.167          0.167       
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.333          0.167       
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.306       
    ## distance_to_yeast25:YeastEMM_F66  0.153          0.153          0.153       
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.153          0.153       
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.153       
    ## distance_to_yeast48:YeastEMM_F66  0.306          0.153          0.153       
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.306          0.153       
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.333       
    ## distance_to_yeast25:YeastEMM_F70  0.167          0.167          0.167       
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.167          0.167       
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.167       
    ## distance_to_yeast48:YeastEMM_F70  0.333          0.167          0.167       
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.333          0.167       
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.333       
    ## distance_to_yeast25:YeastEMM_F89  0.167          0.167          0.167       
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.167          0.167       
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.167       
    ## distance_to_yeast48:YeastEMM_F89  0.333          0.167          0.167       
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.333          0.167       
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.333       
    ## distance_to_yeast25:YeastF7       0.167          0.167          0.167       
    ## distance_to_yeast32:YeastF7       0.167          0.167          0.167       
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.167       
    ## distance_to_yeast48:YeastF7       0.333          0.167          0.167       
    ## distance_to_yeast55:YeastF7       0.167          0.333          0.167       
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.333       
    ## distance_to_yeast25:YeastSP_F14   0.167          0.167          0.167       
    ## distance_to_yeast32:YeastSP_F14   0.167          0.167          0.167       
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.167       
    ## distance_to_yeast48:YeastSP_F14   0.333          0.167          0.167       
    ## distance_to_yeast55:YeastSP_F14   0.167          0.333          0.167       
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.333       
    ## distance_to_yeast25:YeastZAN_F3   0.167          0.167          0.167       
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.167          0.167       
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.167       
    ## distance_to_yeast48:YeastZAN_F3   0.333          0.167          0.167       
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.333          0.167       
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.333       
    ## distance_to_yeast25:YeastZAN_F4   0.167          0.167          0.167       
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.167          0.167       
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.167       
    ## distance_to_yeast48:YeastZAN_F4   0.333          0.167          0.167       
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.333          0.167       
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000       
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000       
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000       
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000       
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000       
    ## DAI6:YeastF7                      0.000          0.000          0.000       
    ## DAI8:YeastF7                      0.000          0.000          0.000       
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000       
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000       
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000       
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000       
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000       
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000       
    ##                                  d__25:YEMM_F5 d__32:YEMM_F5 d__41:YEMM_F5
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5   0.500                                   
    ## distance_to_yeast41:YeastEMM_F5   0.500         0.500                     
    ## distance_to_yeast48:YeastEMM_F5   0.500         0.500         0.500       
    ## distance_to_yeast55:YeastEMM_F5   0.490         0.490         0.490       
    ## distance_to_yeast17:YeastEMM_F63  0.167         0.167         0.167       
    ## distance_to_yeast25:YeastEMM_F63  0.333         0.167         0.167       
    ## distance_to_yeast32:YeastEMM_F63  0.167         0.333         0.167       
    ## distance_to_yeast41:YeastEMM_F63  0.167         0.167         0.333       
    ## distance_to_yeast48:YeastEMM_F63  0.167         0.167         0.167       
    ## distance_to_yeast55:YeastEMM_F63  0.167         0.167         0.167       
    ## distance_to_yeast17:YeastEMM_F64  0.167         0.167         0.167       
    ## distance_to_yeast25:YeastEMM_F64  0.333         0.167         0.167       
    ## distance_to_yeast32:YeastEMM_F64  0.167         0.333         0.167       
    ## distance_to_yeast41:YeastEMM_F64  0.167         0.167         0.333       
    ## distance_to_yeast48:YeastEMM_F64  0.167         0.167         0.167       
    ## distance_to_yeast55:YeastEMM_F64  0.167         0.167         0.167       
    ## distance_to_yeast17:YeastEMM_F65  0.167         0.167         0.167       
    ## distance_to_yeast25:YeastEMM_F65  0.333         0.167         0.167       
    ## distance_to_yeast32:YeastEMM_F65  0.167         0.333         0.167       
    ## distance_to_yeast41:YeastEMM_F65  0.167         0.167         0.333       
    ## distance_to_yeast48:YeastEMM_F65  0.167         0.167         0.167       
    ## distance_to_yeast55:YeastEMM_F65  0.167         0.167         0.167       
    ## distance_to_yeast17:YeastEMM_F66  0.153         0.153         0.153       
    ## distance_to_yeast25:YeastEMM_F66  0.306         0.153         0.153       
    ## distance_to_yeast32:YeastEMM_F66  0.153         0.306         0.153       
    ## distance_to_yeast41:YeastEMM_F66  0.153         0.153         0.306       
    ## distance_to_yeast48:YeastEMM_F66  0.153         0.153         0.153       
    ## distance_to_yeast55:YeastEMM_F66  0.153         0.153         0.153       
    ## distance_to_yeast17:YeastEMM_F70  0.167         0.167         0.167       
    ## distance_to_yeast25:YeastEMM_F70  0.333         0.167         0.167       
    ## distance_to_yeast32:YeastEMM_F70  0.167         0.333         0.167       
    ## distance_to_yeast41:YeastEMM_F70  0.167         0.167         0.333       
    ## distance_to_yeast48:YeastEMM_F70  0.167         0.167         0.167       
    ## distance_to_yeast55:YeastEMM_F70  0.167         0.167         0.167       
    ## distance_to_yeast17:YeastEMM_F89  0.167         0.167         0.167       
    ## distance_to_yeast25:YeastEMM_F89  0.333         0.167         0.167       
    ## distance_to_yeast32:YeastEMM_F89  0.167         0.333         0.167       
    ## distance_to_yeast41:YeastEMM_F89  0.167         0.167         0.333       
    ## distance_to_yeast48:YeastEMM_F89  0.167         0.167         0.167       
    ## distance_to_yeast55:YeastEMM_F89  0.167         0.167         0.167       
    ## distance_to_yeast17:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast25:YeastF7       0.333         0.167         0.167       
    ## distance_to_yeast32:YeastF7       0.167         0.333         0.167       
    ## distance_to_yeast41:YeastF7       0.167         0.167         0.333       
    ## distance_to_yeast48:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast55:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast17:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastSP_F14   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastSP_F14   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastSP_F14   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast17:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastZAN_F3   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F3   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastZAN_F3   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast17:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastZAN_F4   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F4   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastZAN_F4   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F4   0.167         0.167         0.167       
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI6:YeastF7                      0.000         0.000         0.000       
    ## DAI8:YeastF7                      0.000         0.000         0.000       
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000       
    ##                                  d__48:YEMM_F5 d__55:YEMM_F5 d__17:YEMM_F63
    ## DAI6                                                                       
    ## DAI8                                                                       
    ## distance_to_yeast17                                                        
    ## distance_to_yeast25                                                        
    ## distance_to_yeast32                                                        
    ## distance_to_yeast41                                                        
    ## distance_to_yeast48                                                        
    ## distance_to_yeast55                                                        
    ## YeastEMM_F3                                                                
    ## YeastEMM_F34                                                               
    ## YeastEMM_F47                                                               
    ## YeastEMM_F48                                                               
    ## YeastEMM_F49                                                               
    ## YeastEMM_F5                                                                
    ## YeastEMM_F63                                                               
    ## YeastEMM_F64                                                               
    ## YeastEMM_F65                                                               
    ## YeastEMM_F66                                                               
    ## YeastEMM_F70                                                               
    ## YeastEMM_F89                                                               
    ## YeastF7                                                                    
    ## YeastSP_F14                                                                
    ## YeastZAN_F3                                                                
    ## YeastZAN_F4                                                                
    ## distance_to_yeast17:YeastEMM_F3                                            
    ## distance_to_yeast25:YeastEMM_F3                                            
    ## distance_to_yeast32:YeastEMM_F3                                            
    ## distance_to_yeast41:YeastEMM_F3                                            
    ## distance_to_yeast48:YeastEMM_F3                                            
    ## distance_to_yeast55:YeastEMM_F3                                            
    ## distance_to_yeast17:YeastEMM_F34                                           
    ## distance_to_yeast25:YeastEMM_F34                                           
    ## distance_to_yeast32:YeastEMM_F34                                           
    ## distance_to_yeast41:YeastEMM_F34                                           
    ## distance_to_yeast48:YeastEMM_F34                                           
    ## distance_to_yeast55:YeastEMM_F34                                           
    ## distance_to_yeast17:YeastEMM_F47                                           
    ## distance_to_yeast25:YeastEMM_F47                                           
    ## distance_to_yeast32:YeastEMM_F47                                           
    ## distance_to_yeast41:YeastEMM_F47                                           
    ## distance_to_yeast48:YeastEMM_F47                                           
    ## distance_to_yeast55:YeastEMM_F47                                           
    ## distance_to_yeast17:YeastEMM_F48                                           
    ## distance_to_yeast25:YeastEMM_F48                                           
    ## distance_to_yeast32:YeastEMM_F48                                           
    ## distance_to_yeast41:YeastEMM_F48                                           
    ## distance_to_yeast48:YeastEMM_F48                                           
    ## distance_to_yeast55:YeastEMM_F48                                           
    ## distance_to_yeast17:YeastEMM_F49                                           
    ## distance_to_yeast25:YeastEMM_F49                                           
    ## distance_to_yeast32:YeastEMM_F49                                           
    ## distance_to_yeast41:YeastEMM_F49                                           
    ## distance_to_yeast48:YeastEMM_F49                                           
    ## distance_to_yeast55:YeastEMM_F49                                           
    ## distance_to_yeast17:YeastEMM_F5                                            
    ## distance_to_yeast25:YeastEMM_F5                                            
    ## distance_to_yeast32:YeastEMM_F5                                            
    ## distance_to_yeast41:YeastEMM_F5                                            
    ## distance_to_yeast48:YeastEMM_F5                                            
    ## distance_to_yeast55:YeastEMM_F5   0.490                                    
    ## distance_to_yeast17:YeastEMM_F63  0.167         0.163                      
    ## distance_to_yeast25:YeastEMM_F63  0.167         0.163         0.500        
    ## distance_to_yeast32:YeastEMM_F63  0.167         0.163         0.500        
    ## distance_to_yeast41:YeastEMM_F63  0.167         0.163         0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.333         0.163         0.500        
    ## distance_to_yeast55:YeastEMM_F63  0.167         0.326         0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.167         0.163         0.333        
    ## distance_to_yeast25:YeastEMM_F64  0.167         0.163         0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167         0.163         0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167         0.163         0.167        
    ## distance_to_yeast48:YeastEMM_F64  0.333         0.163         0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167         0.326         0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167         0.163         0.333        
    ## distance_to_yeast25:YeastEMM_F65  0.167         0.163         0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167         0.163         0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167         0.163         0.167        
    ## distance_to_yeast48:YeastEMM_F65  0.333         0.163         0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167         0.326         0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153         0.150         0.306        
    ## distance_to_yeast25:YeastEMM_F66  0.153         0.150         0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153         0.150         0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153         0.150         0.153        
    ## distance_to_yeast48:YeastEMM_F66  0.306         0.150         0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153         0.299         0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167         0.163         0.333        
    ## distance_to_yeast25:YeastEMM_F70  0.167         0.163         0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167         0.163         0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167         0.163         0.167        
    ## distance_to_yeast48:YeastEMM_F70  0.333         0.163         0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167         0.326         0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167         0.163         0.333        
    ## distance_to_yeast25:YeastEMM_F89  0.167         0.163         0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167         0.163         0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167         0.163         0.167        
    ## distance_to_yeast48:YeastEMM_F89  0.333         0.163         0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167         0.326         0.167        
    ## distance_to_yeast17:YeastF7       0.167         0.163         0.333        
    ## distance_to_yeast25:YeastF7       0.167         0.163         0.167        
    ## distance_to_yeast32:YeastF7       0.167         0.163         0.167        
    ## distance_to_yeast41:YeastF7       0.167         0.163         0.167        
    ## distance_to_yeast48:YeastF7       0.333         0.163         0.167        
    ## distance_to_yeast55:YeastF7       0.167         0.326         0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167         0.163         0.333        
    ## distance_to_yeast25:YeastSP_F14   0.167         0.163         0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167         0.163         0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167         0.163         0.167        
    ## distance_to_yeast48:YeastSP_F14   0.333         0.163         0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167         0.326         0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167         0.163         0.333        
    ## distance_to_yeast25:YeastZAN_F3   0.167         0.163         0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167         0.163         0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167         0.163         0.167        
    ## distance_to_yeast48:YeastZAN_F3   0.333         0.163         0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167         0.326         0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167         0.163         0.333        
    ## distance_to_yeast25:YeastZAN_F4   0.167         0.163         0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167         0.163         0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167         0.163         0.167        
    ## distance_to_yeast48:YeastZAN_F4   0.333         0.163         0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167         0.326         0.167        
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000        
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000        
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000        
    ## DAI8:YeastEMM_F5                  0.000         0.027         0.000        
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000        
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000        
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000        
    ## DAI6:YeastF7                      0.000         0.000         0.000        
    ## DAI8:YeastF7                      0.000         0.000         0.000        
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000        
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000        
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000        
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000        
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000        
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000        
    ##                                  d__25:YEMM_F63 d__32:YEMM_F63 d__41:YEMM_F63
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63  0.500                                      
    ## distance_to_yeast41:YeastEMM_F63  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F63  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F64  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast25:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast48:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.167          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__48:YEMM_F63 d__55:YEMM_F63 d__17:YEMM_F64
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63  0.500                                      
    ## distance_to_yeast17:YeastEMM_F64  0.167          0.167                       
    ## distance_to_yeast25:YeastEMM_F64  0.167          0.167          0.500        
    ## distance_to_yeast32:YeastEMM_F64  0.167          0.167          0.500        
    ## distance_to_yeast41:YeastEMM_F64  0.167          0.167          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.333          0.167          0.500        
    ## distance_to_yeast55:YeastEMM_F64  0.167          0.333          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast25:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast48:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast25:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast48:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.333          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__25:YEMM_F64 d__32:YEMM_F64 d__41:YEMM_F64
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63                                             
    ## distance_to_yeast17:YeastEMM_F64                                             
    ## distance_to_yeast25:YeastEMM_F64                                             
    ## distance_to_yeast32:YeastEMM_F64  0.500                                      
    ## distance_to_yeast41:YeastEMM_F64  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F64  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F65  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast25:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast48:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.167          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__48:YEMM_F64 d__55:YEMM_F64 d__17:YEMM_F65
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63                                             
    ## distance_to_yeast17:YeastEMM_F64                                             
    ## distance_to_yeast25:YeastEMM_F64                                             
    ## distance_to_yeast32:YeastEMM_F64                                             
    ## distance_to_yeast41:YeastEMM_F64                                             
    ## distance_to_yeast48:YeastEMM_F64                                             
    ## distance_to_yeast55:YeastEMM_F64  0.500                                      
    ## distance_to_yeast17:YeastEMM_F65  0.167          0.167                       
    ## distance_to_yeast25:YeastEMM_F65  0.167          0.167          0.500        
    ## distance_to_yeast32:YeastEMM_F65  0.167          0.167          0.500        
    ## distance_to_yeast41:YeastEMM_F65  0.167          0.167          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.333          0.167          0.500        
    ## distance_to_yeast55:YeastEMM_F65  0.167          0.333          0.500        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast25:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast48:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast25:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast48:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast25:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast48:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast25:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast48:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.333          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__25:YEMM_F65 d__32:YEMM_F65 d__41:YEMM_F65
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63                                             
    ## distance_to_yeast17:YeastEMM_F64                                             
    ## distance_to_yeast25:YeastEMM_F64                                             
    ## distance_to_yeast32:YeastEMM_F64                                             
    ## distance_to_yeast41:YeastEMM_F64                                             
    ## distance_to_yeast48:YeastEMM_F64                                             
    ## distance_to_yeast55:YeastEMM_F64                                             
    ## distance_to_yeast17:YeastEMM_F65                                             
    ## distance_to_yeast25:YeastEMM_F65                                             
    ## distance_to_yeast32:YeastEMM_F65  0.500                                      
    ## distance_to_yeast41:YeastEMM_F65  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F65  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F66  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F70  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast25:YeastEMM_F89  0.333          0.167          0.167        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.333          0.167        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.333        
    ## distance_to_yeast48:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.167          0.167        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast25:YeastF7       0.333          0.167          0.167        
    ## distance_to_yeast32:YeastF7       0.167          0.333          0.167        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.333        
    ## distance_to_yeast48:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast55:YeastF7       0.167          0.167          0.167        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastSP_F14   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F3   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.167          0.167        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast25:YeastZAN_F4   0.333          0.167          0.167        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.333          0.167        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.333        
    ## distance_to_yeast48:YeastZAN_F4   0.167          0.167          0.167        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.167          0.167        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__48:YEMM_F65 d__55:YEMM_F65 d__17:YEMM_F66
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63                                             
    ## distance_to_yeast17:YeastEMM_F64                                             
    ## distance_to_yeast25:YeastEMM_F64                                             
    ## distance_to_yeast32:YeastEMM_F64                                             
    ## distance_to_yeast41:YeastEMM_F64                                             
    ## distance_to_yeast48:YeastEMM_F64                                             
    ## distance_to_yeast55:YeastEMM_F64                                             
    ## distance_to_yeast17:YeastEMM_F65                                             
    ## distance_to_yeast25:YeastEMM_F65                                             
    ## distance_to_yeast32:YeastEMM_F65                                             
    ## distance_to_yeast41:YeastEMM_F65                                             
    ## distance_to_yeast48:YeastEMM_F65                                             
    ## distance_to_yeast55:YeastEMM_F65  0.500                                      
    ## distance_to_yeast17:YeastEMM_F66  0.153          0.153                       
    ## distance_to_yeast25:YeastEMM_F66  0.153          0.153          0.500        
    ## distance_to_yeast32:YeastEMM_F66  0.153          0.153          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.153          0.153          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.306          0.153          0.500        
    ## distance_to_yeast55:YeastEMM_F66  0.153          0.306          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.167          0.167          0.306        
    ## distance_to_yeast25:YeastEMM_F70  0.167          0.167          0.153        
    ## distance_to_yeast32:YeastEMM_F70  0.167          0.167          0.153        
    ## distance_to_yeast41:YeastEMM_F70  0.167          0.167          0.153        
    ## distance_to_yeast48:YeastEMM_F70  0.333          0.167          0.153        
    ## distance_to_yeast55:YeastEMM_F70  0.167          0.333          0.153        
    ## distance_to_yeast17:YeastEMM_F89  0.167          0.167          0.306        
    ## distance_to_yeast25:YeastEMM_F89  0.167          0.167          0.153        
    ## distance_to_yeast32:YeastEMM_F89  0.167          0.167          0.153        
    ## distance_to_yeast41:YeastEMM_F89  0.167          0.167          0.153        
    ## distance_to_yeast48:YeastEMM_F89  0.333          0.167          0.153        
    ## distance_to_yeast55:YeastEMM_F89  0.167          0.333          0.153        
    ## distance_to_yeast17:YeastF7       0.167          0.167          0.306        
    ## distance_to_yeast25:YeastF7       0.167          0.167          0.153        
    ## distance_to_yeast32:YeastF7       0.167          0.167          0.153        
    ## distance_to_yeast41:YeastF7       0.167          0.167          0.153        
    ## distance_to_yeast48:YeastF7       0.333          0.167          0.153        
    ## distance_to_yeast55:YeastF7       0.167          0.333          0.153        
    ## distance_to_yeast17:YeastSP_F14   0.167          0.167          0.306        
    ## distance_to_yeast25:YeastSP_F14   0.167          0.167          0.153        
    ## distance_to_yeast32:YeastSP_F14   0.167          0.167          0.153        
    ## distance_to_yeast41:YeastSP_F14   0.167          0.167          0.153        
    ## distance_to_yeast48:YeastSP_F14   0.333          0.167          0.153        
    ## distance_to_yeast55:YeastSP_F14   0.167          0.333          0.153        
    ## distance_to_yeast17:YeastZAN_F3   0.167          0.167          0.306        
    ## distance_to_yeast25:YeastZAN_F3   0.167          0.167          0.153        
    ## distance_to_yeast32:YeastZAN_F3   0.167          0.167          0.153        
    ## distance_to_yeast41:YeastZAN_F3   0.167          0.167          0.153        
    ## distance_to_yeast48:YeastZAN_F3   0.333          0.167          0.153        
    ## distance_to_yeast55:YeastZAN_F3   0.167          0.333          0.153        
    ## distance_to_yeast17:YeastZAN_F4   0.167          0.167          0.306        
    ## distance_to_yeast25:YeastZAN_F4   0.167          0.167          0.153        
    ## distance_to_yeast32:YeastZAN_F4   0.167          0.167          0.153        
    ## distance_to_yeast41:YeastZAN_F4   0.167          0.167          0.153        
    ## distance_to_yeast48:YeastZAN_F4   0.333          0.167          0.153        
    ## distance_to_yeast55:YeastZAN_F4   0.167          0.333          0.153        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__25:YEMM_F66 d__32:YEMM_F66 d__41:YEMM_F66
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63                                             
    ## distance_to_yeast17:YeastEMM_F64                                             
    ## distance_to_yeast25:YeastEMM_F64                                             
    ## distance_to_yeast32:YeastEMM_F64                                             
    ## distance_to_yeast41:YeastEMM_F64                                             
    ## distance_to_yeast48:YeastEMM_F64                                             
    ## distance_to_yeast55:YeastEMM_F64                                             
    ## distance_to_yeast17:YeastEMM_F65                                             
    ## distance_to_yeast25:YeastEMM_F65                                             
    ## distance_to_yeast32:YeastEMM_F65                                             
    ## distance_to_yeast41:YeastEMM_F65                                             
    ## distance_to_yeast48:YeastEMM_F65                                             
    ## distance_to_yeast55:YeastEMM_F65                                             
    ## distance_to_yeast17:YeastEMM_F66                                             
    ## distance_to_yeast25:YeastEMM_F66                                             
    ## distance_to_yeast32:YeastEMM_F66  0.500                                      
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F66  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F70  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F70  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F70  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F70  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F70  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastEMM_F89  0.153          0.153          0.153        
    ## distance_to_yeast25:YeastEMM_F89  0.306          0.153          0.153        
    ## distance_to_yeast32:YeastEMM_F89  0.153          0.306          0.153        
    ## distance_to_yeast41:YeastEMM_F89  0.153          0.153          0.306        
    ## distance_to_yeast48:YeastEMM_F89  0.153          0.153          0.153        
    ## distance_to_yeast55:YeastEMM_F89  0.153          0.153          0.153        
    ## distance_to_yeast17:YeastF7       0.153          0.153          0.153        
    ## distance_to_yeast25:YeastF7       0.306          0.153          0.153        
    ## distance_to_yeast32:YeastF7       0.153          0.306          0.153        
    ## distance_to_yeast41:YeastF7       0.153          0.153          0.306        
    ## distance_to_yeast48:YeastF7       0.153          0.153          0.153        
    ## distance_to_yeast55:YeastF7       0.153          0.153          0.153        
    ## distance_to_yeast17:YeastSP_F14   0.153          0.153          0.153        
    ## distance_to_yeast25:YeastSP_F14   0.306          0.153          0.153        
    ## distance_to_yeast32:YeastSP_F14   0.153          0.306          0.153        
    ## distance_to_yeast41:YeastSP_F14   0.153          0.153          0.306        
    ## distance_to_yeast48:YeastSP_F14   0.153          0.153          0.153        
    ## distance_to_yeast55:YeastSP_F14   0.153          0.153          0.153        
    ## distance_to_yeast17:YeastZAN_F3   0.153          0.153          0.153        
    ## distance_to_yeast25:YeastZAN_F3   0.306          0.153          0.153        
    ## distance_to_yeast32:YeastZAN_F3   0.153          0.306          0.153        
    ## distance_to_yeast41:YeastZAN_F3   0.153          0.153          0.306        
    ## distance_to_yeast48:YeastZAN_F3   0.153          0.153          0.153        
    ## distance_to_yeast55:YeastZAN_F3   0.153          0.153          0.153        
    ## distance_to_yeast17:YeastZAN_F4   0.153          0.153          0.153        
    ## distance_to_yeast25:YeastZAN_F4   0.306          0.153          0.153        
    ## distance_to_yeast32:YeastZAN_F4   0.153          0.306          0.153        
    ## distance_to_yeast41:YeastZAN_F4   0.153          0.153          0.306        
    ## distance_to_yeast48:YeastZAN_F4   0.153          0.153          0.153        
    ## distance_to_yeast55:YeastZAN_F4   0.153          0.153          0.153        
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000        
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000        
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000        
    ## DAI6:YeastF7                      0.000          0.000          0.000        
    ## DAI8:YeastF7                      0.000          0.000          0.000        
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000        
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000        
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000        
    ##                                  d__48:YEMM_F66 d__55:YEMM_F66 d__17:YEMM_F7
    ## DAI6                                                                        
    ## DAI8                                                                        
    ## distance_to_yeast17                                                         
    ## distance_to_yeast25                                                         
    ## distance_to_yeast32                                                         
    ## distance_to_yeast41                                                         
    ## distance_to_yeast48                                                         
    ## distance_to_yeast55                                                         
    ## YeastEMM_F3                                                                 
    ## YeastEMM_F34                                                                
    ## YeastEMM_F47                                                                
    ## YeastEMM_F48                                                                
    ## YeastEMM_F49                                                                
    ## YeastEMM_F5                                                                 
    ## YeastEMM_F63                                                                
    ## YeastEMM_F64                                                                
    ## YeastEMM_F65                                                                
    ## YeastEMM_F66                                                                
    ## YeastEMM_F70                                                                
    ## YeastEMM_F89                                                                
    ## YeastF7                                                                     
    ## YeastSP_F14                                                                 
    ## YeastZAN_F3                                                                 
    ## YeastZAN_F4                                                                 
    ## distance_to_yeast17:YeastEMM_F3                                             
    ## distance_to_yeast25:YeastEMM_F3                                             
    ## distance_to_yeast32:YeastEMM_F3                                             
    ## distance_to_yeast41:YeastEMM_F3                                             
    ## distance_to_yeast48:YeastEMM_F3                                             
    ## distance_to_yeast55:YeastEMM_F3                                             
    ## distance_to_yeast17:YeastEMM_F34                                            
    ## distance_to_yeast25:YeastEMM_F34                                            
    ## distance_to_yeast32:YeastEMM_F34                                            
    ## distance_to_yeast41:YeastEMM_F34                                            
    ## distance_to_yeast48:YeastEMM_F34                                            
    ## distance_to_yeast55:YeastEMM_F34                                            
    ## distance_to_yeast17:YeastEMM_F47                                            
    ## distance_to_yeast25:YeastEMM_F47                                            
    ## distance_to_yeast32:YeastEMM_F47                                            
    ## distance_to_yeast41:YeastEMM_F47                                            
    ## distance_to_yeast48:YeastEMM_F47                                            
    ## distance_to_yeast55:YeastEMM_F47                                            
    ## distance_to_yeast17:YeastEMM_F48                                            
    ## distance_to_yeast25:YeastEMM_F48                                            
    ## distance_to_yeast32:YeastEMM_F48                                            
    ## distance_to_yeast41:YeastEMM_F48                                            
    ## distance_to_yeast48:YeastEMM_F48                                            
    ## distance_to_yeast55:YeastEMM_F48                                            
    ## distance_to_yeast17:YeastEMM_F49                                            
    ## distance_to_yeast25:YeastEMM_F49                                            
    ## distance_to_yeast32:YeastEMM_F49                                            
    ## distance_to_yeast41:YeastEMM_F49                                            
    ## distance_to_yeast48:YeastEMM_F49                                            
    ## distance_to_yeast55:YeastEMM_F49                                            
    ## distance_to_yeast17:YeastEMM_F5                                             
    ## distance_to_yeast25:YeastEMM_F5                                             
    ## distance_to_yeast32:YeastEMM_F5                                             
    ## distance_to_yeast41:YeastEMM_F5                                             
    ## distance_to_yeast48:YeastEMM_F5                                             
    ## distance_to_yeast55:YeastEMM_F5                                             
    ## distance_to_yeast17:YeastEMM_F63                                            
    ## distance_to_yeast25:YeastEMM_F63                                            
    ## distance_to_yeast32:YeastEMM_F63                                            
    ## distance_to_yeast41:YeastEMM_F63                                            
    ## distance_to_yeast48:YeastEMM_F63                                            
    ## distance_to_yeast55:YeastEMM_F63                                            
    ## distance_to_yeast17:YeastEMM_F64                                            
    ## distance_to_yeast25:YeastEMM_F64                                            
    ## distance_to_yeast32:YeastEMM_F64                                            
    ## distance_to_yeast41:YeastEMM_F64                                            
    ## distance_to_yeast48:YeastEMM_F64                                            
    ## distance_to_yeast55:YeastEMM_F64                                            
    ## distance_to_yeast17:YeastEMM_F65                                            
    ## distance_to_yeast25:YeastEMM_F65                                            
    ## distance_to_yeast32:YeastEMM_F65                                            
    ## distance_to_yeast41:YeastEMM_F65                                            
    ## distance_to_yeast48:YeastEMM_F65                                            
    ## distance_to_yeast55:YeastEMM_F65                                            
    ## distance_to_yeast17:YeastEMM_F66                                            
    ## distance_to_yeast25:YeastEMM_F66                                            
    ## distance_to_yeast32:YeastEMM_F66                                            
    ## distance_to_yeast41:YeastEMM_F66                                            
    ## distance_to_yeast48:YeastEMM_F66                                            
    ## distance_to_yeast55:YeastEMM_F66  0.500                                     
    ## distance_to_yeast17:YeastEMM_F70  0.153          0.153                      
    ## distance_to_yeast25:YeastEMM_F70  0.153          0.153          0.500       
    ## distance_to_yeast32:YeastEMM_F70  0.153          0.153          0.500       
    ## distance_to_yeast41:YeastEMM_F70  0.153          0.153          0.500       
    ## distance_to_yeast48:YeastEMM_F70  0.306          0.153          0.500       
    ## distance_to_yeast55:YeastEMM_F70  0.153          0.306          0.500       
    ## distance_to_yeast17:YeastEMM_F89  0.153          0.153          0.333       
    ## distance_to_yeast25:YeastEMM_F89  0.153          0.153          0.167       
    ## distance_to_yeast32:YeastEMM_F89  0.153          0.153          0.167       
    ## distance_to_yeast41:YeastEMM_F89  0.153          0.153          0.167       
    ## distance_to_yeast48:YeastEMM_F89  0.306          0.153          0.167       
    ## distance_to_yeast55:YeastEMM_F89  0.153          0.306          0.167       
    ## distance_to_yeast17:YeastF7       0.153          0.153          0.333       
    ## distance_to_yeast25:YeastF7       0.153          0.153          0.167       
    ## distance_to_yeast32:YeastF7       0.153          0.153          0.167       
    ## distance_to_yeast41:YeastF7       0.153          0.153          0.167       
    ## distance_to_yeast48:YeastF7       0.306          0.153          0.167       
    ## distance_to_yeast55:YeastF7       0.153          0.306          0.167       
    ## distance_to_yeast17:YeastSP_F14   0.153          0.153          0.333       
    ## distance_to_yeast25:YeastSP_F14   0.153          0.153          0.167       
    ## distance_to_yeast32:YeastSP_F14   0.153          0.153          0.167       
    ## distance_to_yeast41:YeastSP_F14   0.153          0.153          0.167       
    ## distance_to_yeast48:YeastSP_F14   0.306          0.153          0.167       
    ## distance_to_yeast55:YeastSP_F14   0.153          0.306          0.167       
    ## distance_to_yeast17:YeastZAN_F3   0.153          0.153          0.333       
    ## distance_to_yeast25:YeastZAN_F3   0.153          0.153          0.167       
    ## distance_to_yeast32:YeastZAN_F3   0.153          0.153          0.167       
    ## distance_to_yeast41:YeastZAN_F3   0.153          0.153          0.167       
    ## distance_to_yeast48:YeastZAN_F3   0.306          0.153          0.167       
    ## distance_to_yeast55:YeastZAN_F3   0.153          0.306          0.167       
    ## distance_to_yeast17:YeastZAN_F4   0.153          0.153          0.333       
    ## distance_to_yeast25:YeastZAN_F4   0.153          0.153          0.167       
    ## distance_to_yeast32:YeastZAN_F4   0.153          0.153          0.167       
    ## distance_to_yeast41:YeastZAN_F4   0.153          0.153          0.167       
    ## distance_to_yeast48:YeastZAN_F4   0.306          0.153          0.167       
    ## distance_to_yeast55:YeastZAN_F4   0.153          0.306          0.167       
    ## DAI6:YeastEMM_F3                  0.000          0.000          0.000       
    ## DAI8:YeastEMM_F3                  0.000          0.000          0.000       
    ## DAI6:YeastEMM_F34                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F34                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F47                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F47                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F48                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F48                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F49                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F49                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F5                  0.000          0.000          0.000       
    ## DAI8:YeastEMM_F5                  0.000          0.000          0.000       
    ## DAI6:YeastEMM_F63                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F63                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F64                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F64                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F65                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F65                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F66                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F66                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F70                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F70                 0.000          0.000          0.000       
    ## DAI6:YeastEMM_F89                 0.000          0.000          0.000       
    ## DAI8:YeastEMM_F89                 0.000          0.000          0.000       
    ## DAI6:YeastF7                      0.000          0.000          0.000       
    ## DAI8:YeastF7                      0.000          0.000          0.000       
    ## DAI6:YeastSP_F14                  0.000          0.000          0.000       
    ## DAI8:YeastSP_F14                  0.000          0.000          0.000       
    ## DAI6:YeastZAN_F3                  0.000          0.000          0.000       
    ## DAI8:YeastZAN_F3                  0.000          0.000          0.000       
    ## DAI6:YeastZAN_F4                  0.000          0.000          0.000       
    ## DAI8:YeastZAN_F4                  0.000          0.000          0.000       
    ##                                  d__25:YEMM_F7 d__32:YEMM_F7 d__41:YEMM_F7
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70  0.500                                   
    ## distance_to_yeast41:YeastEMM_F70  0.500         0.500                     
    ## distance_to_yeast48:YeastEMM_F70  0.500         0.500         0.500       
    ## distance_to_yeast55:YeastEMM_F70  0.500         0.500         0.500       
    ## distance_to_yeast17:YeastEMM_F89  0.167         0.167         0.167       
    ## distance_to_yeast25:YeastEMM_F89  0.333         0.167         0.167       
    ## distance_to_yeast32:YeastEMM_F89  0.167         0.333         0.167       
    ## distance_to_yeast41:YeastEMM_F89  0.167         0.167         0.333       
    ## distance_to_yeast48:YeastEMM_F89  0.167         0.167         0.167       
    ## distance_to_yeast55:YeastEMM_F89  0.167         0.167         0.167       
    ## distance_to_yeast17:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast25:YeastF7       0.333         0.167         0.167       
    ## distance_to_yeast32:YeastF7       0.167         0.333         0.167       
    ## distance_to_yeast41:YeastF7       0.167         0.167         0.333       
    ## distance_to_yeast48:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast55:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast17:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastSP_F14   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastSP_F14   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastSP_F14   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast17:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastZAN_F3   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F3   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastZAN_F3   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast17:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastZAN_F4   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F4   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastZAN_F4   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F4   0.167         0.167         0.167       
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI6:YeastF7                      0.000         0.000         0.000       
    ## DAI8:YeastF7                      0.000         0.000         0.000       
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000       
    ##                                  d__48:YEMM_F7 d__55:YEMM_F7 d__17:YEMM_F8
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70  0.500                                   
    ## distance_to_yeast17:YeastEMM_F89  0.167         0.167                     
    ## distance_to_yeast25:YeastEMM_F89  0.167         0.167         0.500       
    ## distance_to_yeast32:YeastEMM_F89  0.167         0.167         0.500       
    ## distance_to_yeast41:YeastEMM_F89  0.167         0.167         0.500       
    ## distance_to_yeast48:YeastEMM_F89  0.333         0.167         0.500       
    ## distance_to_yeast55:YeastEMM_F89  0.167         0.333         0.500       
    ## distance_to_yeast17:YeastF7       0.167         0.167         0.333       
    ## distance_to_yeast25:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast32:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast41:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast48:YeastF7       0.333         0.167         0.167       
    ## distance_to_yeast55:YeastF7       0.167         0.333         0.167       
    ## distance_to_yeast17:YeastSP_F14   0.167         0.167         0.333       
    ## distance_to_yeast25:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast32:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast41:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast48:YeastSP_F14   0.333         0.167         0.167       
    ## distance_to_yeast55:YeastSP_F14   0.167         0.333         0.167       
    ## distance_to_yeast17:YeastZAN_F3   0.167         0.167         0.333       
    ## distance_to_yeast25:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast41:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast48:YeastZAN_F3   0.333         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F3   0.167         0.333         0.167       
    ## distance_to_yeast17:YeastZAN_F4   0.167         0.167         0.333       
    ## distance_to_yeast25:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast41:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast48:YeastZAN_F4   0.333         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F4   0.167         0.333         0.167       
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI6:YeastF7                      0.000         0.000         0.000       
    ## DAI8:YeastF7                      0.000         0.000         0.000       
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000       
    ##                                  d__25:YEMM_F8 d__32:YEMM_F8 d__41:YEMM_F8
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89  0.500                                   
    ## distance_to_yeast41:YeastEMM_F89  0.500         0.500                     
    ## distance_to_yeast48:YeastEMM_F89  0.500         0.500         0.500       
    ## distance_to_yeast55:YeastEMM_F89  0.500         0.500         0.500       
    ## distance_to_yeast17:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast25:YeastF7       0.333         0.167         0.167       
    ## distance_to_yeast32:YeastF7       0.167         0.333         0.167       
    ## distance_to_yeast41:YeastF7       0.167         0.167         0.333       
    ## distance_to_yeast48:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast55:YeastF7       0.167         0.167         0.167       
    ## distance_to_yeast17:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastSP_F14   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastSP_F14   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastSP_F14   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastSP_F14   0.167         0.167         0.167       
    ## distance_to_yeast17:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastZAN_F3   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F3   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastZAN_F3   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F3   0.167         0.167         0.167       
    ## distance_to_yeast17:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastZAN_F4   0.333         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F4   0.167         0.333         0.167       
    ## distance_to_yeast41:YeastZAN_F4   0.167         0.167         0.333       
    ## distance_to_yeast48:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F4   0.167         0.167         0.167       
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI6:YeastF7                      0.000         0.000         0.000       
    ## DAI8:YeastF7                      0.000         0.000         0.000       
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000       
    ##                                  d__48:YEMM_F8 d__55:YEMM_F8 d__17:YF d__25:YF
    ## DAI6                                                                          
    ## DAI8                                                                          
    ## distance_to_yeast17                                                           
    ## distance_to_yeast25                                                           
    ## distance_to_yeast32                                                           
    ## distance_to_yeast41                                                           
    ## distance_to_yeast48                                                           
    ## distance_to_yeast55                                                           
    ## YeastEMM_F3                                                                   
    ## YeastEMM_F34                                                                  
    ## YeastEMM_F47                                                                  
    ## YeastEMM_F48                                                                  
    ## YeastEMM_F49                                                                  
    ## YeastEMM_F5                                                                   
    ## YeastEMM_F63                                                                  
    ## YeastEMM_F64                                                                  
    ## YeastEMM_F65                                                                  
    ## YeastEMM_F66                                                                  
    ## YeastEMM_F70                                                                  
    ## YeastEMM_F89                                                                  
    ## YeastF7                                                                       
    ## YeastSP_F14                                                                   
    ## YeastZAN_F3                                                                   
    ## YeastZAN_F4                                                                   
    ## distance_to_yeast17:YeastEMM_F3                                               
    ## distance_to_yeast25:YeastEMM_F3                                               
    ## distance_to_yeast32:YeastEMM_F3                                               
    ## distance_to_yeast41:YeastEMM_F3                                               
    ## distance_to_yeast48:YeastEMM_F3                                               
    ## distance_to_yeast55:YeastEMM_F3                                               
    ## distance_to_yeast17:YeastEMM_F34                                              
    ## distance_to_yeast25:YeastEMM_F34                                              
    ## distance_to_yeast32:YeastEMM_F34                                              
    ## distance_to_yeast41:YeastEMM_F34                                              
    ## distance_to_yeast48:YeastEMM_F34                                              
    ## distance_to_yeast55:YeastEMM_F34                                              
    ## distance_to_yeast17:YeastEMM_F47                                              
    ## distance_to_yeast25:YeastEMM_F47                                              
    ## distance_to_yeast32:YeastEMM_F47                                              
    ## distance_to_yeast41:YeastEMM_F47                                              
    ## distance_to_yeast48:YeastEMM_F47                                              
    ## distance_to_yeast55:YeastEMM_F47                                              
    ## distance_to_yeast17:YeastEMM_F48                                              
    ## distance_to_yeast25:YeastEMM_F48                                              
    ## distance_to_yeast32:YeastEMM_F48                                              
    ## distance_to_yeast41:YeastEMM_F48                                              
    ## distance_to_yeast48:YeastEMM_F48                                              
    ## distance_to_yeast55:YeastEMM_F48                                              
    ## distance_to_yeast17:YeastEMM_F49                                              
    ## distance_to_yeast25:YeastEMM_F49                                              
    ## distance_to_yeast32:YeastEMM_F49                                              
    ## distance_to_yeast41:YeastEMM_F49                                              
    ## distance_to_yeast48:YeastEMM_F49                                              
    ## distance_to_yeast55:YeastEMM_F49                                              
    ## distance_to_yeast17:YeastEMM_F5                                               
    ## distance_to_yeast25:YeastEMM_F5                                               
    ## distance_to_yeast32:YeastEMM_F5                                               
    ## distance_to_yeast41:YeastEMM_F5                                               
    ## distance_to_yeast48:YeastEMM_F5                                               
    ## distance_to_yeast55:YeastEMM_F5                                               
    ## distance_to_yeast17:YeastEMM_F63                                              
    ## distance_to_yeast25:YeastEMM_F63                                              
    ## distance_to_yeast32:YeastEMM_F63                                              
    ## distance_to_yeast41:YeastEMM_F63                                              
    ## distance_to_yeast48:YeastEMM_F63                                              
    ## distance_to_yeast55:YeastEMM_F63                                              
    ## distance_to_yeast17:YeastEMM_F64                                              
    ## distance_to_yeast25:YeastEMM_F64                                              
    ## distance_to_yeast32:YeastEMM_F64                                              
    ## distance_to_yeast41:YeastEMM_F64                                              
    ## distance_to_yeast48:YeastEMM_F64                                              
    ## distance_to_yeast55:YeastEMM_F64                                              
    ## distance_to_yeast17:YeastEMM_F65                                              
    ## distance_to_yeast25:YeastEMM_F65                                              
    ## distance_to_yeast32:YeastEMM_F65                                              
    ## distance_to_yeast41:YeastEMM_F65                                              
    ## distance_to_yeast48:YeastEMM_F65                                              
    ## distance_to_yeast55:YeastEMM_F65                                              
    ## distance_to_yeast17:YeastEMM_F66                                              
    ## distance_to_yeast25:YeastEMM_F66                                              
    ## distance_to_yeast32:YeastEMM_F66                                              
    ## distance_to_yeast41:YeastEMM_F66                                              
    ## distance_to_yeast48:YeastEMM_F66                                              
    ## distance_to_yeast55:YeastEMM_F66                                              
    ## distance_to_yeast17:YeastEMM_F70                                              
    ## distance_to_yeast25:YeastEMM_F70                                              
    ## distance_to_yeast32:YeastEMM_F70                                              
    ## distance_to_yeast41:YeastEMM_F70                                              
    ## distance_to_yeast48:YeastEMM_F70                                              
    ## distance_to_yeast55:YeastEMM_F70                                              
    ## distance_to_yeast17:YeastEMM_F89                                              
    ## distance_to_yeast25:YeastEMM_F89                                              
    ## distance_to_yeast32:YeastEMM_F89                                              
    ## distance_to_yeast41:YeastEMM_F89                                              
    ## distance_to_yeast48:YeastEMM_F89                                              
    ## distance_to_yeast55:YeastEMM_F89  0.500                                       
    ## distance_to_yeast17:YeastF7       0.167         0.167                         
    ## distance_to_yeast25:YeastF7       0.167         0.167         0.500           
    ## distance_to_yeast32:YeastF7       0.167         0.167         0.500    0.500  
    ## distance_to_yeast41:YeastF7       0.167         0.167         0.500    0.500  
    ## distance_to_yeast48:YeastF7       0.333         0.167         0.500    0.500  
    ## distance_to_yeast55:YeastF7       0.167         0.333         0.500    0.500  
    ## distance_to_yeast17:YeastSP_F14   0.167         0.167         0.333    0.167  
    ## distance_to_yeast25:YeastSP_F14   0.167         0.167         0.167    0.333  
    ## distance_to_yeast32:YeastSP_F14   0.167         0.167         0.167    0.167  
    ## distance_to_yeast41:YeastSP_F14   0.167         0.167         0.167    0.167  
    ## distance_to_yeast48:YeastSP_F14   0.333         0.167         0.167    0.167  
    ## distance_to_yeast55:YeastSP_F14   0.167         0.333         0.167    0.167  
    ## distance_to_yeast17:YeastZAN_F3   0.167         0.167         0.333    0.167  
    ## distance_to_yeast25:YeastZAN_F3   0.167         0.167         0.167    0.333  
    ## distance_to_yeast32:YeastZAN_F3   0.167         0.167         0.167    0.167  
    ## distance_to_yeast41:YeastZAN_F3   0.167         0.167         0.167    0.167  
    ## distance_to_yeast48:YeastZAN_F3   0.333         0.167         0.167    0.167  
    ## distance_to_yeast55:YeastZAN_F3   0.167         0.333         0.167    0.167  
    ## distance_to_yeast17:YeastZAN_F4   0.167         0.167         0.333    0.167  
    ## distance_to_yeast25:YeastZAN_F4   0.167         0.167         0.167    0.333  
    ## distance_to_yeast32:YeastZAN_F4   0.167         0.167         0.167    0.167  
    ## distance_to_yeast41:YeastZAN_F4   0.167         0.167         0.167    0.167  
    ## distance_to_yeast48:YeastZAN_F4   0.333         0.167         0.167    0.167  
    ## distance_to_yeast55:YeastZAN_F4   0.167         0.333         0.167    0.167  
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000    0.000  
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000    0.000  
    ## DAI6:YeastF7                      0.000         0.000         0.000    0.000  
    ## DAI8:YeastF7                      0.000         0.000         0.000    0.000  
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000    0.000  
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000    0.000  
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000    0.000  
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000    0.000  
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000    0.000  
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000    0.000  
    ##                                  d__32:YF d__41:YF d__48:YF d__55:YF d__17:YS
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63                                             
    ## distance_to_yeast17:YeastEMM_F64                                             
    ## distance_to_yeast25:YeastEMM_F64                                             
    ## distance_to_yeast32:YeastEMM_F64                                             
    ## distance_to_yeast41:YeastEMM_F64                                             
    ## distance_to_yeast48:YeastEMM_F64                                             
    ## distance_to_yeast55:YeastEMM_F64                                             
    ## distance_to_yeast17:YeastEMM_F65                                             
    ## distance_to_yeast25:YeastEMM_F65                                             
    ## distance_to_yeast32:YeastEMM_F65                                             
    ## distance_to_yeast41:YeastEMM_F65                                             
    ## distance_to_yeast48:YeastEMM_F65                                             
    ## distance_to_yeast55:YeastEMM_F65                                             
    ## distance_to_yeast17:YeastEMM_F66                                             
    ## distance_to_yeast25:YeastEMM_F66                                             
    ## distance_to_yeast32:YeastEMM_F66                                             
    ## distance_to_yeast41:YeastEMM_F66                                             
    ## distance_to_yeast48:YeastEMM_F66                                             
    ## distance_to_yeast55:YeastEMM_F66                                             
    ## distance_to_yeast17:YeastEMM_F70                                             
    ## distance_to_yeast25:YeastEMM_F70                                             
    ## distance_to_yeast32:YeastEMM_F70                                             
    ## distance_to_yeast41:YeastEMM_F70                                             
    ## distance_to_yeast48:YeastEMM_F70                                             
    ## distance_to_yeast55:YeastEMM_F70                                             
    ## distance_to_yeast17:YeastEMM_F89                                             
    ## distance_to_yeast25:YeastEMM_F89                                             
    ## distance_to_yeast32:YeastEMM_F89                                             
    ## distance_to_yeast41:YeastEMM_F89                                             
    ## distance_to_yeast48:YeastEMM_F89                                             
    ## distance_to_yeast55:YeastEMM_F89                                             
    ## distance_to_yeast17:YeastF7                                                  
    ## distance_to_yeast25:YeastF7                                                  
    ## distance_to_yeast32:YeastF7                                                  
    ## distance_to_yeast41:YeastF7       0.500                                      
    ## distance_to_yeast48:YeastF7       0.500    0.500                             
    ## distance_to_yeast55:YeastF7       0.500    0.500    0.500                    
    ## distance_to_yeast17:YeastSP_F14   0.167    0.167    0.167    0.167           
    ## distance_to_yeast25:YeastSP_F14   0.167    0.167    0.167    0.167    0.500  
    ## distance_to_yeast32:YeastSP_F14   0.333    0.167    0.167    0.167    0.500  
    ## distance_to_yeast41:YeastSP_F14   0.167    0.333    0.167    0.167    0.500  
    ## distance_to_yeast48:YeastSP_F14   0.167    0.167    0.333    0.167    0.500  
    ## distance_to_yeast55:YeastSP_F14   0.167    0.167    0.167    0.333    0.500  
    ## distance_to_yeast17:YeastZAN_F3   0.167    0.167    0.167    0.167    0.333  
    ## distance_to_yeast25:YeastZAN_F3   0.167    0.167    0.167    0.167    0.167  
    ## distance_to_yeast32:YeastZAN_F3   0.333    0.167    0.167    0.167    0.167  
    ## distance_to_yeast41:YeastZAN_F3   0.167    0.333    0.167    0.167    0.167  
    ## distance_to_yeast48:YeastZAN_F3   0.167    0.167    0.333    0.167    0.167  
    ## distance_to_yeast55:YeastZAN_F3   0.167    0.167    0.167    0.333    0.167  
    ## distance_to_yeast17:YeastZAN_F4   0.167    0.167    0.167    0.167    0.333  
    ## distance_to_yeast25:YeastZAN_F4   0.167    0.167    0.167    0.167    0.167  
    ## distance_to_yeast32:YeastZAN_F4   0.333    0.167    0.167    0.167    0.167  
    ## distance_to_yeast41:YeastZAN_F4   0.167    0.333    0.167    0.167    0.167  
    ## distance_to_yeast48:YeastZAN_F4   0.167    0.167    0.333    0.167    0.167  
    ## distance_to_yeast55:YeastZAN_F4   0.167    0.167    0.167    0.333    0.167  
    ## DAI6:YeastEMM_F3                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F3                  0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F34                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F34                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F47                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F47                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F48                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F48                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F49                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F49                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F5                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F5                  0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F63                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F63                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F64                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F64                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F65                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F65                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F66                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F66                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F70                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F70                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F89                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F89                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastF7                      0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastF7                      0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastSP_F14                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastSP_F14                  0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastZAN_F3                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastZAN_F3                  0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastZAN_F4                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastZAN_F4                  0.000    0.000    0.000    0.000    0.000  
    ##                                  d__25:YS d__32:YS d__41:YS d__48:YS d__55:YS
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17                                                          
    ## distance_to_yeast25                                                          
    ## distance_to_yeast32                                                          
    ## distance_to_yeast41                                                          
    ## distance_to_yeast48                                                          
    ## distance_to_yeast55                                                          
    ## YeastEMM_F3                                                                  
    ## YeastEMM_F34                                                                 
    ## YeastEMM_F47                                                                 
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastF7                                                                      
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34                                             
    ## distance_to_yeast32:YeastEMM_F34                                             
    ## distance_to_yeast41:YeastEMM_F34                                             
    ## distance_to_yeast48:YeastEMM_F34                                             
    ## distance_to_yeast55:YeastEMM_F34                                             
    ## distance_to_yeast17:YeastEMM_F47                                             
    ## distance_to_yeast25:YeastEMM_F47                                             
    ## distance_to_yeast32:YeastEMM_F47                                             
    ## distance_to_yeast41:YeastEMM_F47                                             
    ## distance_to_yeast48:YeastEMM_F47                                             
    ## distance_to_yeast55:YeastEMM_F47                                             
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48                                             
    ## distance_to_yeast17:YeastEMM_F49                                             
    ## distance_to_yeast25:YeastEMM_F49                                             
    ## distance_to_yeast32:YeastEMM_F49                                             
    ## distance_to_yeast41:YeastEMM_F49                                             
    ## distance_to_yeast48:YeastEMM_F49                                             
    ## distance_to_yeast55:YeastEMM_F49                                             
    ## distance_to_yeast17:YeastEMM_F5                                              
    ## distance_to_yeast25:YeastEMM_F5                                              
    ## distance_to_yeast32:YeastEMM_F5                                              
    ## distance_to_yeast41:YeastEMM_F5                                              
    ## distance_to_yeast48:YeastEMM_F5                                              
    ## distance_to_yeast55:YeastEMM_F5                                              
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63                                             
    ## distance_to_yeast17:YeastEMM_F64                                             
    ## distance_to_yeast25:YeastEMM_F64                                             
    ## distance_to_yeast32:YeastEMM_F64                                             
    ## distance_to_yeast41:YeastEMM_F64                                             
    ## distance_to_yeast48:YeastEMM_F64                                             
    ## distance_to_yeast55:YeastEMM_F64                                             
    ## distance_to_yeast17:YeastEMM_F65                                             
    ## distance_to_yeast25:YeastEMM_F65                                             
    ## distance_to_yeast32:YeastEMM_F65                                             
    ## distance_to_yeast41:YeastEMM_F65                                             
    ## distance_to_yeast48:YeastEMM_F65                                             
    ## distance_to_yeast55:YeastEMM_F65                                             
    ## distance_to_yeast17:YeastEMM_F66                                             
    ## distance_to_yeast25:YeastEMM_F66                                             
    ## distance_to_yeast32:YeastEMM_F66                                             
    ## distance_to_yeast41:YeastEMM_F66                                             
    ## distance_to_yeast48:YeastEMM_F66                                             
    ## distance_to_yeast55:YeastEMM_F66                                             
    ## distance_to_yeast17:YeastEMM_F70                                             
    ## distance_to_yeast25:YeastEMM_F70                                             
    ## distance_to_yeast32:YeastEMM_F70                                             
    ## distance_to_yeast41:YeastEMM_F70                                             
    ## distance_to_yeast48:YeastEMM_F70                                             
    ## distance_to_yeast55:YeastEMM_F70                                             
    ## distance_to_yeast17:YeastEMM_F89                                             
    ## distance_to_yeast25:YeastEMM_F89                                             
    ## distance_to_yeast32:YeastEMM_F89                                             
    ## distance_to_yeast41:YeastEMM_F89                                             
    ## distance_to_yeast48:YeastEMM_F89                                             
    ## distance_to_yeast55:YeastEMM_F89                                             
    ## distance_to_yeast17:YeastF7                                                  
    ## distance_to_yeast25:YeastF7                                                  
    ## distance_to_yeast32:YeastF7                                                  
    ## distance_to_yeast41:YeastF7                                                  
    ## distance_to_yeast48:YeastF7                                                  
    ## distance_to_yeast55:YeastF7                                                  
    ## distance_to_yeast17:YeastSP_F14                                              
    ## distance_to_yeast25:YeastSP_F14                                              
    ## distance_to_yeast32:YeastSP_F14   0.500                                      
    ## distance_to_yeast41:YeastSP_F14   0.500    0.500                             
    ## distance_to_yeast48:YeastSP_F14   0.500    0.500    0.500                    
    ## distance_to_yeast55:YeastSP_F14   0.500    0.500    0.500    0.500           
    ## distance_to_yeast17:YeastZAN_F3   0.167    0.167    0.167    0.167    0.167  
    ## distance_to_yeast25:YeastZAN_F3   0.333    0.167    0.167    0.167    0.167  
    ## distance_to_yeast32:YeastZAN_F3   0.167    0.333    0.167    0.167    0.167  
    ## distance_to_yeast41:YeastZAN_F3   0.167    0.167    0.333    0.167    0.167  
    ## distance_to_yeast48:YeastZAN_F3   0.167    0.167    0.167    0.333    0.167  
    ## distance_to_yeast55:YeastZAN_F3   0.167    0.167    0.167    0.167    0.333  
    ## distance_to_yeast17:YeastZAN_F4   0.167    0.167    0.167    0.167    0.167  
    ## distance_to_yeast25:YeastZAN_F4   0.333    0.167    0.167    0.167    0.167  
    ## distance_to_yeast32:YeastZAN_F4   0.167    0.333    0.167    0.167    0.167  
    ## distance_to_yeast41:YeastZAN_F4   0.167    0.167    0.333    0.167    0.167  
    ## distance_to_yeast48:YeastZAN_F4   0.167    0.167    0.167    0.333    0.167  
    ## distance_to_yeast55:YeastZAN_F4   0.167    0.167    0.167    0.167    0.333  
    ## DAI6:YeastEMM_F3                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F3                  0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F34                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F34                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F47                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F47                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F48                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F48                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F49                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F49                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F5                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F5                  0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F63                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F63                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F64                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F64                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F65                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F65                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F66                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F66                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F70                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F70                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastEMM_F89                 0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastEMM_F89                 0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastF7                      0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastF7                      0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastSP_F14                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastSP_F14                  0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastZAN_F3                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastZAN_F3                  0.000    0.000    0.000    0.000    0.000  
    ## DAI6:YeastZAN_F4                  0.000    0.000    0.000    0.000    0.000  
    ## DAI8:YeastZAN_F4                  0.000    0.000    0.000    0.000    0.000  
    ##                                  d__17:YZAN_F3 d__25:YZAN_F3 d__32:YZAN_F3
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3   0.500                                   
    ## distance_to_yeast32:YeastZAN_F3   0.500         0.500                     
    ## distance_to_yeast41:YeastZAN_F3   0.500         0.500         0.500       
    ## distance_to_yeast48:YeastZAN_F3   0.500         0.500         0.500       
    ## distance_to_yeast55:YeastZAN_F3   0.500         0.500         0.500       
    ## distance_to_yeast17:YeastZAN_F4   0.333         0.167         0.167       
    ## distance_to_yeast25:YeastZAN_F4   0.167         0.333         0.167       
    ## distance_to_yeast32:YeastZAN_F4   0.167         0.167         0.333       
    ## distance_to_yeast41:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast48:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast55:YeastZAN_F4   0.167         0.167         0.167       
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI6:YeastF7                      0.000         0.000         0.000       
    ## DAI8:YeastF7                      0.000         0.000         0.000       
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000       
    ##                                  d__41:YZAN_F3 d__48:YZAN_F3 d__55:YZAN_F3
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3                                           
    ## distance_to_yeast32:YeastZAN_F3                                           
    ## distance_to_yeast41:YeastZAN_F3                                           
    ## distance_to_yeast48:YeastZAN_F3   0.500                                   
    ## distance_to_yeast55:YeastZAN_F3   0.500         0.500                     
    ## distance_to_yeast17:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast25:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast32:YeastZAN_F4   0.167         0.167         0.167       
    ## distance_to_yeast41:YeastZAN_F4   0.333         0.167         0.167       
    ## distance_to_yeast48:YeastZAN_F4   0.167         0.333         0.167       
    ## distance_to_yeast55:YeastZAN_F4   0.167         0.167         0.333       
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI6:YeastF7                      0.000         0.000         0.000       
    ## DAI8:YeastF7                      0.000         0.000         0.000       
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000       
    ##                                  d__17:YZAN_F4 d__25:YZAN_F4 d__32:YZAN_F4
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3                                           
    ## distance_to_yeast32:YeastZAN_F3                                           
    ## distance_to_yeast41:YeastZAN_F3                                           
    ## distance_to_yeast48:YeastZAN_F3                                           
    ## distance_to_yeast55:YeastZAN_F3                                           
    ## distance_to_yeast17:YeastZAN_F4                                           
    ## distance_to_yeast25:YeastZAN_F4   0.500                                   
    ## distance_to_yeast32:YeastZAN_F4   0.500         0.500                     
    ## distance_to_yeast41:YeastZAN_F4   0.500         0.500         0.500       
    ## distance_to_yeast48:YeastZAN_F4   0.500         0.500         0.500       
    ## distance_to_yeast55:YeastZAN_F4   0.500         0.500         0.500       
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI6:YeastF7                      0.000         0.000         0.000       
    ## DAI8:YeastF7                      0.000         0.000         0.000       
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000       
    ##                                  d__41:YZAN_F4 d__48:YZAN_F4 d__55:YZAN_F4
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3                                           
    ## distance_to_yeast32:YeastZAN_F3                                           
    ## distance_to_yeast41:YeastZAN_F3                                           
    ## distance_to_yeast48:YeastZAN_F3                                           
    ## distance_to_yeast55:YeastZAN_F3                                           
    ## distance_to_yeast17:YeastZAN_F4                                           
    ## distance_to_yeast25:YeastZAN_F4                                           
    ## distance_to_yeast32:YeastZAN_F4                                           
    ## distance_to_yeast41:YeastZAN_F4                                           
    ## distance_to_yeast48:YeastZAN_F4   0.500                                   
    ## distance_to_yeast55:YeastZAN_F4   0.500         0.500                     
    ## DAI6:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F3                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F34                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F47                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F48                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F49                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI8:YeastEMM_F5                  0.000         0.000         0.000       
    ## DAI6:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F63                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F64                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F65                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F66                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F70                 0.000         0.000         0.000       
    ## DAI6:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI8:YeastEMM_F89                 0.000         0.000         0.000       
    ## DAI6:YeastF7                      0.000         0.000         0.000       
    ## DAI8:YeastF7                      0.000         0.000         0.000       
    ## DAI6:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI8:YeastSP_F14                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F3                  0.000         0.000         0.000       
    ## DAI6:YeastZAN_F4                  0.000         0.000         0.000       
    ## DAI8:YeastZAN_F4                  0.000         0.000         0.000       
    ##                                  DAI6:YsEMM_F3 DAI8:YsEMM_F3 DAI6:YEMM_F34
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3                                           
    ## distance_to_yeast32:YeastZAN_F3                                           
    ## distance_to_yeast41:YeastZAN_F3                                           
    ## distance_to_yeast48:YeastZAN_F3                                           
    ## distance_to_yeast55:YeastZAN_F3                                           
    ## distance_to_yeast17:YeastZAN_F4                                           
    ## distance_to_yeast25:YeastZAN_F4                                           
    ## distance_to_yeast32:YeastZAN_F4                                           
    ## distance_to_yeast41:YeastZAN_F4                                           
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## DAI6:YeastEMM_F3                                                          
    ## DAI8:YeastEMM_F3                  0.500                                   
    ## DAI6:YeastEMM_F34                 0.333         0.167                     
    ## DAI8:YeastEMM_F34                 0.167         0.333         0.500       
    ## DAI6:YeastEMM_F47                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F47                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F48                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F48                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F49                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F49                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F5                  0.333         0.167         0.333       
    ## DAI8:YeastEMM_F5                  0.165         0.330         0.165       
    ## DAI6:YeastEMM_F63                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F63                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F64                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F64                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F65                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F65                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F66                 0.308         0.154         0.308       
    ## DAI8:YeastEMM_F66                 0.154         0.308         0.154       
    ## DAI6:YeastEMM_F70                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F70                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F89                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F89                 0.167         0.333         0.167       
    ## DAI6:YeastF7                      0.333         0.167         0.333       
    ## DAI8:YeastF7                      0.167         0.333         0.167       
    ## DAI6:YeastSP_F14                  0.333         0.167         0.333       
    ## DAI8:YeastSP_F14                  0.167         0.333         0.167       
    ## DAI6:YeastZAN_F3                  0.333         0.167         0.333       
    ## DAI8:YeastZAN_F3                  0.167         0.333         0.167       
    ## DAI6:YeastZAN_F4                  0.333         0.167         0.333       
    ## DAI8:YeastZAN_F4                  0.167         0.333         0.167       
    ##                                  DAI8:YEMM_F34 DAI6:YEMM_F47 DAI8:YEMM_F47
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3                                           
    ## distance_to_yeast32:YeastZAN_F3                                           
    ## distance_to_yeast41:YeastZAN_F3                                           
    ## distance_to_yeast48:YeastZAN_F3                                           
    ## distance_to_yeast55:YeastZAN_F3                                           
    ## distance_to_yeast17:YeastZAN_F4                                           
    ## distance_to_yeast25:YeastZAN_F4                                           
    ## distance_to_yeast32:YeastZAN_F4                                           
    ## distance_to_yeast41:YeastZAN_F4                                           
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## DAI6:YeastEMM_F3                                                          
    ## DAI8:YeastEMM_F3                                                          
    ## DAI6:YeastEMM_F34                                                         
    ## DAI8:YeastEMM_F34                                                         
    ## DAI6:YeastEMM_F47                 0.167                                   
    ## DAI8:YeastEMM_F47                 0.333         0.500                     
    ## DAI6:YeastEMM_F48                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F48                 0.333         0.167         0.333       
    ## DAI6:YeastEMM_F49                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F49                 0.333         0.167         0.333       
    ## DAI6:YeastEMM_F5                  0.167         0.333         0.167       
    ## DAI8:YeastEMM_F5                  0.330         0.165         0.330       
    ## DAI6:YeastEMM_F63                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F63                 0.333         0.167         0.333       
    ## DAI6:YeastEMM_F64                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F64                 0.333         0.167         0.333       
    ## DAI6:YeastEMM_F65                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F65                 0.333         0.167         0.333       
    ## DAI6:YeastEMM_F66                 0.154         0.308         0.154       
    ## DAI8:YeastEMM_F66                 0.308         0.154         0.308       
    ## DAI6:YeastEMM_F70                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F70                 0.333         0.167         0.333       
    ## DAI6:YeastEMM_F89                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F89                 0.333         0.167         0.333       
    ## DAI6:YeastF7                      0.167         0.333         0.167       
    ## DAI8:YeastF7                      0.333         0.167         0.333       
    ## DAI6:YeastSP_F14                  0.167         0.333         0.167       
    ## DAI8:YeastSP_F14                  0.333         0.167         0.333       
    ## DAI6:YeastZAN_F3                  0.167         0.333         0.167       
    ## DAI8:YeastZAN_F3                  0.333         0.167         0.333       
    ## DAI6:YeastZAN_F4                  0.167         0.333         0.167       
    ## DAI8:YeastZAN_F4                  0.333         0.167         0.333       
    ##                                  DAI6:YEMM_F48 DAI8:YEMM_F48 DAI6:YEMM_F49
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3                                           
    ## distance_to_yeast32:YeastZAN_F3                                           
    ## distance_to_yeast41:YeastZAN_F3                                           
    ## distance_to_yeast48:YeastZAN_F3                                           
    ## distance_to_yeast55:YeastZAN_F3                                           
    ## distance_to_yeast17:YeastZAN_F4                                           
    ## distance_to_yeast25:YeastZAN_F4                                           
    ## distance_to_yeast32:YeastZAN_F4                                           
    ## distance_to_yeast41:YeastZAN_F4                                           
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## DAI6:YeastEMM_F3                                                          
    ## DAI8:YeastEMM_F3                                                          
    ## DAI6:YeastEMM_F34                                                         
    ## DAI8:YeastEMM_F34                                                         
    ## DAI6:YeastEMM_F47                                                         
    ## DAI8:YeastEMM_F47                                                         
    ## DAI6:YeastEMM_F48                                                         
    ## DAI8:YeastEMM_F48                 0.500                                   
    ## DAI6:YeastEMM_F49                 0.333         0.167                     
    ## DAI8:YeastEMM_F49                 0.167         0.333         0.500       
    ## DAI6:YeastEMM_F5                  0.333         0.167         0.333       
    ## DAI8:YeastEMM_F5                  0.165         0.330         0.165       
    ## DAI6:YeastEMM_F63                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F63                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F64                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F64                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F65                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F65                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F66                 0.308         0.154         0.308       
    ## DAI8:YeastEMM_F66                 0.154         0.308         0.154       
    ## DAI6:YeastEMM_F70                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F70                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F89                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F89                 0.167         0.333         0.167       
    ## DAI6:YeastF7                      0.333         0.167         0.333       
    ## DAI8:YeastF7                      0.167         0.333         0.167       
    ## DAI6:YeastSP_F14                  0.333         0.167         0.333       
    ## DAI8:YeastSP_F14                  0.167         0.333         0.167       
    ## DAI6:YeastZAN_F3                  0.333         0.167         0.333       
    ## DAI8:YeastZAN_F3                  0.167         0.333         0.167       
    ## DAI6:YeastZAN_F4                  0.333         0.167         0.333       
    ## DAI8:YeastZAN_F4                  0.167         0.333         0.167       
    ##                                  DAI8:YEMM_F49 DAI6:YEMM_F5 DAI8:YEMM_F5
    ## DAI6                                                                    
    ## DAI8                                                                    
    ## distance_to_yeast17                                                     
    ## distance_to_yeast25                                                     
    ## distance_to_yeast32                                                     
    ## distance_to_yeast41                                                     
    ## distance_to_yeast48                                                     
    ## distance_to_yeast55                                                     
    ## YeastEMM_F3                                                             
    ## YeastEMM_F34                                                            
    ## YeastEMM_F47                                                            
    ## YeastEMM_F48                                                            
    ## YeastEMM_F49                                                            
    ## YeastEMM_F5                                                             
    ## YeastEMM_F63                                                            
    ## YeastEMM_F64                                                            
    ## YeastEMM_F65                                                            
    ## YeastEMM_F66                                                            
    ## YeastEMM_F70                                                            
    ## YeastEMM_F89                                                            
    ## YeastF7                                                                 
    ## YeastSP_F14                                                             
    ## YeastZAN_F3                                                             
    ## YeastZAN_F4                                                             
    ## distance_to_yeast17:YeastEMM_F3                                         
    ## distance_to_yeast25:YeastEMM_F3                                         
    ## distance_to_yeast32:YeastEMM_F3                                         
    ## distance_to_yeast41:YeastEMM_F3                                         
    ## distance_to_yeast48:YeastEMM_F3                                         
    ## distance_to_yeast55:YeastEMM_F3                                         
    ## distance_to_yeast17:YeastEMM_F34                                        
    ## distance_to_yeast25:YeastEMM_F34                                        
    ## distance_to_yeast32:YeastEMM_F34                                        
    ## distance_to_yeast41:YeastEMM_F34                                        
    ## distance_to_yeast48:YeastEMM_F34                                        
    ## distance_to_yeast55:YeastEMM_F34                                        
    ## distance_to_yeast17:YeastEMM_F47                                        
    ## distance_to_yeast25:YeastEMM_F47                                        
    ## distance_to_yeast32:YeastEMM_F47                                        
    ## distance_to_yeast41:YeastEMM_F47                                        
    ## distance_to_yeast48:YeastEMM_F47                                        
    ## distance_to_yeast55:YeastEMM_F47                                        
    ## distance_to_yeast17:YeastEMM_F48                                        
    ## distance_to_yeast25:YeastEMM_F48                                        
    ## distance_to_yeast32:YeastEMM_F48                                        
    ## distance_to_yeast41:YeastEMM_F48                                        
    ## distance_to_yeast48:YeastEMM_F48                                        
    ## distance_to_yeast55:YeastEMM_F48                                        
    ## distance_to_yeast17:YeastEMM_F49                                        
    ## distance_to_yeast25:YeastEMM_F49                                        
    ## distance_to_yeast32:YeastEMM_F49                                        
    ## distance_to_yeast41:YeastEMM_F49                                        
    ## distance_to_yeast48:YeastEMM_F49                                        
    ## distance_to_yeast55:YeastEMM_F49                                        
    ## distance_to_yeast17:YeastEMM_F5                                         
    ## distance_to_yeast25:YeastEMM_F5                                         
    ## distance_to_yeast32:YeastEMM_F5                                         
    ## distance_to_yeast41:YeastEMM_F5                                         
    ## distance_to_yeast48:YeastEMM_F5                                         
    ## distance_to_yeast55:YeastEMM_F5                                         
    ## distance_to_yeast17:YeastEMM_F63                                        
    ## distance_to_yeast25:YeastEMM_F63                                        
    ## distance_to_yeast32:YeastEMM_F63                                        
    ## distance_to_yeast41:YeastEMM_F63                                        
    ## distance_to_yeast48:YeastEMM_F63                                        
    ## distance_to_yeast55:YeastEMM_F63                                        
    ## distance_to_yeast17:YeastEMM_F64                                        
    ## distance_to_yeast25:YeastEMM_F64                                        
    ## distance_to_yeast32:YeastEMM_F64                                        
    ## distance_to_yeast41:YeastEMM_F64                                        
    ## distance_to_yeast48:YeastEMM_F64                                        
    ## distance_to_yeast55:YeastEMM_F64                                        
    ## distance_to_yeast17:YeastEMM_F65                                        
    ## distance_to_yeast25:YeastEMM_F65                                        
    ## distance_to_yeast32:YeastEMM_F65                                        
    ## distance_to_yeast41:YeastEMM_F65                                        
    ## distance_to_yeast48:YeastEMM_F65                                        
    ## distance_to_yeast55:YeastEMM_F65                                        
    ## distance_to_yeast17:YeastEMM_F66                                        
    ## distance_to_yeast25:YeastEMM_F66                                        
    ## distance_to_yeast32:YeastEMM_F66                                        
    ## distance_to_yeast41:YeastEMM_F66                                        
    ## distance_to_yeast48:YeastEMM_F66                                        
    ## distance_to_yeast55:YeastEMM_F66                                        
    ## distance_to_yeast17:YeastEMM_F70                                        
    ## distance_to_yeast25:YeastEMM_F70                                        
    ## distance_to_yeast32:YeastEMM_F70                                        
    ## distance_to_yeast41:YeastEMM_F70                                        
    ## distance_to_yeast48:YeastEMM_F70                                        
    ## distance_to_yeast55:YeastEMM_F70                                        
    ## distance_to_yeast17:YeastEMM_F89                                        
    ## distance_to_yeast25:YeastEMM_F89                                        
    ## distance_to_yeast32:YeastEMM_F89                                        
    ## distance_to_yeast41:YeastEMM_F89                                        
    ## distance_to_yeast48:YeastEMM_F89                                        
    ## distance_to_yeast55:YeastEMM_F89                                        
    ## distance_to_yeast17:YeastF7                                             
    ## distance_to_yeast25:YeastF7                                             
    ## distance_to_yeast32:YeastF7                                             
    ## distance_to_yeast41:YeastF7                                             
    ## distance_to_yeast48:YeastF7                                             
    ## distance_to_yeast55:YeastF7                                             
    ## distance_to_yeast17:YeastSP_F14                                         
    ## distance_to_yeast25:YeastSP_F14                                         
    ## distance_to_yeast32:YeastSP_F14                                         
    ## distance_to_yeast41:YeastSP_F14                                         
    ## distance_to_yeast48:YeastSP_F14                                         
    ## distance_to_yeast55:YeastSP_F14                                         
    ## distance_to_yeast17:YeastZAN_F3                                         
    ## distance_to_yeast25:YeastZAN_F3                                         
    ## distance_to_yeast32:YeastZAN_F3                                         
    ## distance_to_yeast41:YeastZAN_F3                                         
    ## distance_to_yeast48:YeastZAN_F3                                         
    ## distance_to_yeast55:YeastZAN_F3                                         
    ## distance_to_yeast17:YeastZAN_F4                                         
    ## distance_to_yeast25:YeastZAN_F4                                         
    ## distance_to_yeast32:YeastZAN_F4                                         
    ## distance_to_yeast41:YeastZAN_F4                                         
    ## distance_to_yeast48:YeastZAN_F4                                         
    ## distance_to_yeast55:YeastZAN_F4                                         
    ## DAI6:YeastEMM_F3                                                        
    ## DAI8:YeastEMM_F3                                                        
    ## DAI6:YeastEMM_F34                                                       
    ## DAI8:YeastEMM_F34                                                       
    ## DAI6:YeastEMM_F47                                                       
    ## DAI8:YeastEMM_F47                                                       
    ## DAI6:YeastEMM_F48                                                       
    ## DAI8:YeastEMM_F48                                                       
    ## DAI6:YeastEMM_F49                                                       
    ## DAI8:YeastEMM_F49                                                       
    ## DAI6:YeastEMM_F5                  0.167                                 
    ## DAI8:YeastEMM_F5                  0.330         0.495                   
    ## DAI6:YeastEMM_F63                 0.167         0.333        0.165      
    ## DAI8:YeastEMM_F63                 0.333         0.167        0.330      
    ## DAI6:YeastEMM_F64                 0.167         0.333        0.165      
    ## DAI8:YeastEMM_F64                 0.333         0.167        0.330      
    ## DAI6:YeastEMM_F65                 0.167         0.333        0.165      
    ## DAI8:YeastEMM_F65                 0.333         0.167        0.330      
    ## DAI6:YeastEMM_F66                 0.154         0.308        0.153      
    ## DAI8:YeastEMM_F66                 0.308         0.154        0.306      
    ## DAI6:YeastEMM_F70                 0.167         0.333        0.165      
    ## DAI8:YeastEMM_F70                 0.333         0.167        0.330      
    ## DAI6:YeastEMM_F89                 0.167         0.333        0.165      
    ## DAI8:YeastEMM_F89                 0.333         0.167        0.330      
    ## DAI6:YeastF7                      0.167         0.333        0.165      
    ## DAI8:YeastF7                      0.333         0.167        0.330      
    ## DAI6:YeastSP_F14                  0.167         0.333        0.165      
    ## DAI8:YeastSP_F14                  0.333         0.167        0.330      
    ## DAI6:YeastZAN_F3                  0.167         0.333        0.165      
    ## DAI8:YeastZAN_F3                  0.333         0.167        0.330      
    ## DAI6:YeastZAN_F4                  0.167         0.333        0.165      
    ## DAI8:YeastZAN_F4                  0.333         0.167        0.330      
    ##                                  DAI6:YEMM_F63 DAI8:YEMM_F63 DAI6:YEMM_F64
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3                                           
    ## distance_to_yeast32:YeastZAN_F3                                           
    ## distance_to_yeast41:YeastZAN_F3                                           
    ## distance_to_yeast48:YeastZAN_F3                                           
    ## distance_to_yeast55:YeastZAN_F3                                           
    ## distance_to_yeast17:YeastZAN_F4                                           
    ## distance_to_yeast25:YeastZAN_F4                                           
    ## distance_to_yeast32:YeastZAN_F4                                           
    ## distance_to_yeast41:YeastZAN_F4                                           
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## DAI6:YeastEMM_F3                                                          
    ## DAI8:YeastEMM_F3                                                          
    ## DAI6:YeastEMM_F34                                                         
    ## DAI8:YeastEMM_F34                                                         
    ## DAI6:YeastEMM_F47                                                         
    ## DAI8:YeastEMM_F47                                                         
    ## DAI6:YeastEMM_F48                                                         
    ## DAI8:YeastEMM_F48                                                         
    ## DAI6:YeastEMM_F49                                                         
    ## DAI8:YeastEMM_F49                                                         
    ## DAI6:YeastEMM_F5                                                          
    ## DAI8:YeastEMM_F5                                                          
    ## DAI6:YeastEMM_F63                                                         
    ## DAI8:YeastEMM_F63                 0.500                                   
    ## DAI6:YeastEMM_F64                 0.333         0.167                     
    ## DAI8:YeastEMM_F64                 0.167         0.333         0.500       
    ## DAI6:YeastEMM_F65                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F65                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F66                 0.308         0.154         0.308       
    ## DAI8:YeastEMM_F66                 0.154         0.308         0.154       
    ## DAI6:YeastEMM_F70                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F70                 0.167         0.333         0.167       
    ## DAI6:YeastEMM_F89                 0.333         0.167         0.333       
    ## DAI8:YeastEMM_F89                 0.167         0.333         0.167       
    ## DAI6:YeastF7                      0.333         0.167         0.333       
    ## DAI8:YeastF7                      0.167         0.333         0.167       
    ## DAI6:YeastSP_F14                  0.333         0.167         0.333       
    ## DAI8:YeastSP_F14                  0.167         0.333         0.167       
    ## DAI6:YeastZAN_F3                  0.333         0.167         0.333       
    ## DAI8:YeastZAN_F3                  0.167         0.333         0.167       
    ## DAI6:YeastZAN_F4                  0.333         0.167         0.333       
    ## DAI8:YeastZAN_F4                  0.167         0.333         0.167       
    ##                                  DAI8:YEMM_F64 DAI6:YEMM_F65 DAI8:YEMM_F65
    ## DAI6                                                                      
    ## DAI8                                                                      
    ## distance_to_yeast17                                                       
    ## distance_to_yeast25                                                       
    ## distance_to_yeast32                                                       
    ## distance_to_yeast41                                                       
    ## distance_to_yeast48                                                       
    ## distance_to_yeast55                                                       
    ## YeastEMM_F3                                                               
    ## YeastEMM_F34                                                              
    ## YeastEMM_F47                                                              
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
    ## YeastF7                                                                   
    ## YeastSP_F14                                                               
    ## YeastZAN_F3                                                               
    ## YeastZAN_F4                                                               
    ## distance_to_yeast17:YeastEMM_F3                                           
    ## distance_to_yeast25:YeastEMM_F3                                           
    ## distance_to_yeast32:YeastEMM_F3                                           
    ## distance_to_yeast41:YeastEMM_F3                                           
    ## distance_to_yeast48:YeastEMM_F3                                           
    ## distance_to_yeast55:YeastEMM_F3                                           
    ## distance_to_yeast17:YeastEMM_F34                                          
    ## distance_to_yeast25:YeastEMM_F34                                          
    ## distance_to_yeast32:YeastEMM_F34                                          
    ## distance_to_yeast41:YeastEMM_F34                                          
    ## distance_to_yeast48:YeastEMM_F34                                          
    ## distance_to_yeast55:YeastEMM_F34                                          
    ## distance_to_yeast17:YeastEMM_F47                                          
    ## distance_to_yeast25:YeastEMM_F47                                          
    ## distance_to_yeast32:YeastEMM_F47                                          
    ## distance_to_yeast41:YeastEMM_F47                                          
    ## distance_to_yeast48:YeastEMM_F47                                          
    ## distance_to_yeast55:YeastEMM_F47                                          
    ## distance_to_yeast17:YeastEMM_F48                                          
    ## distance_to_yeast25:YeastEMM_F48                                          
    ## distance_to_yeast32:YeastEMM_F48                                          
    ## distance_to_yeast41:YeastEMM_F48                                          
    ## distance_to_yeast48:YeastEMM_F48                                          
    ## distance_to_yeast55:YeastEMM_F48                                          
    ## distance_to_yeast17:YeastEMM_F49                                          
    ## distance_to_yeast25:YeastEMM_F49                                          
    ## distance_to_yeast32:YeastEMM_F49                                          
    ## distance_to_yeast41:YeastEMM_F49                                          
    ## distance_to_yeast48:YeastEMM_F49                                          
    ## distance_to_yeast55:YeastEMM_F49                                          
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5                                           
    ## distance_to_yeast41:YeastEMM_F5                                           
    ## distance_to_yeast48:YeastEMM_F5                                           
    ## distance_to_yeast55:YeastEMM_F5                                           
    ## distance_to_yeast17:YeastEMM_F63                                          
    ## distance_to_yeast25:YeastEMM_F63                                          
    ## distance_to_yeast32:YeastEMM_F63                                          
    ## distance_to_yeast41:YeastEMM_F63                                          
    ## distance_to_yeast48:YeastEMM_F63                                          
    ## distance_to_yeast55:YeastEMM_F63                                          
    ## distance_to_yeast17:YeastEMM_F64                                          
    ## distance_to_yeast25:YeastEMM_F64                                          
    ## distance_to_yeast32:YeastEMM_F64                                          
    ## distance_to_yeast41:YeastEMM_F64                                          
    ## distance_to_yeast48:YeastEMM_F64                                          
    ## distance_to_yeast55:YeastEMM_F64                                          
    ## distance_to_yeast17:YeastEMM_F65                                          
    ## distance_to_yeast25:YeastEMM_F65                                          
    ## distance_to_yeast32:YeastEMM_F65                                          
    ## distance_to_yeast41:YeastEMM_F65                                          
    ## distance_to_yeast48:YeastEMM_F65                                          
    ## distance_to_yeast55:YeastEMM_F65                                          
    ## distance_to_yeast17:YeastEMM_F66                                          
    ## distance_to_yeast25:YeastEMM_F66                                          
    ## distance_to_yeast32:YeastEMM_F66                                          
    ## distance_to_yeast41:YeastEMM_F66                                          
    ## distance_to_yeast48:YeastEMM_F66                                          
    ## distance_to_yeast55:YeastEMM_F66                                          
    ## distance_to_yeast17:YeastEMM_F70                                          
    ## distance_to_yeast25:YeastEMM_F70                                          
    ## distance_to_yeast32:YeastEMM_F70                                          
    ## distance_to_yeast41:YeastEMM_F70                                          
    ## distance_to_yeast48:YeastEMM_F70                                          
    ## distance_to_yeast55:YeastEMM_F70                                          
    ## distance_to_yeast17:YeastEMM_F89                                          
    ## distance_to_yeast25:YeastEMM_F89                                          
    ## distance_to_yeast32:YeastEMM_F89                                          
    ## distance_to_yeast41:YeastEMM_F89                                          
    ## distance_to_yeast48:YeastEMM_F89                                          
    ## distance_to_yeast55:YeastEMM_F89                                          
    ## distance_to_yeast17:YeastF7                                               
    ## distance_to_yeast25:YeastF7                                               
    ## distance_to_yeast32:YeastF7                                               
    ## distance_to_yeast41:YeastF7                                               
    ## distance_to_yeast48:YeastF7                                               
    ## distance_to_yeast55:YeastF7                                               
    ## distance_to_yeast17:YeastSP_F14                                           
    ## distance_to_yeast25:YeastSP_F14                                           
    ## distance_to_yeast32:YeastSP_F14                                           
    ## distance_to_yeast41:YeastSP_F14                                           
    ## distance_to_yeast48:YeastSP_F14                                           
    ## distance_to_yeast55:YeastSP_F14                                           
    ## distance_to_yeast17:YeastZAN_F3                                           
    ## distance_to_yeast25:YeastZAN_F3                                           
    ## distance_to_yeast32:YeastZAN_F3                                           
    ## distance_to_yeast41:YeastZAN_F3                                           
    ## distance_to_yeast48:YeastZAN_F3                                           
    ## distance_to_yeast55:YeastZAN_F3                                           
    ## distance_to_yeast17:YeastZAN_F4                                           
    ## distance_to_yeast25:YeastZAN_F4                                           
    ## distance_to_yeast32:YeastZAN_F4                                           
    ## distance_to_yeast41:YeastZAN_F4                                           
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## DAI6:YeastEMM_F3                                                          
    ## DAI8:YeastEMM_F3                                                          
    ## DAI6:YeastEMM_F34                                                         
    ## DAI8:YeastEMM_F34                                                         
    ## DAI6:YeastEMM_F47                                                         
    ## DAI8:YeastEMM_F47                                                         
    ## DAI6:YeastEMM_F48                                                         
    ## DAI8:YeastEMM_F48                                                         
    ## DAI6:YeastEMM_F49                                                         
    ## DAI8:YeastEMM_F49                                                         
    ## DAI6:YeastEMM_F5                                                          
    ## DAI8:YeastEMM_F5                                                          
    ## DAI6:YeastEMM_F63                                                         
    ## DAI8:YeastEMM_F63                                                         
    ## DAI6:YeastEMM_F64                                                         
    ## DAI8:YeastEMM_F64                                                         
    ## DAI6:YeastEMM_F65                 0.167                                   
    ## DAI8:YeastEMM_F65                 0.333         0.500                     
    ## DAI6:YeastEMM_F66                 0.154         0.308         0.154       
    ## DAI8:YeastEMM_F66                 0.308         0.154         0.308       
    ## DAI6:YeastEMM_F70                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F70                 0.333         0.167         0.333       
    ## DAI6:YeastEMM_F89                 0.167         0.333         0.167       
    ## DAI8:YeastEMM_F89                 0.333         0.167         0.333       
    ## DAI6:YeastF7                      0.167         0.333         0.167       
    ## DAI8:YeastF7                      0.333         0.167         0.333       
    ## DAI6:YeastSP_F14                  0.167         0.333         0.167       
    ## DAI8:YeastSP_F14                  0.333         0.167         0.333       
    ## DAI6:YeastZAN_F3                  0.167         0.333         0.167       
    ## DAI8:YeastZAN_F3                  0.333         0.167         0.333       
    ## DAI6:YeastZAN_F4                  0.167         0.333         0.167       
    ## DAI8:YeastZAN_F4                  0.333         0.167         0.333       
    ##                                  DAI6:YEMM_F66 DAI8:YEMM_F66 DAI6:YEMM_F7
    ## DAI6                                                                     
    ## DAI8                                                                     
    ## distance_to_yeast17                                                      
    ## distance_to_yeast25                                                      
    ## distance_to_yeast32                                                      
    ## distance_to_yeast41                                                      
    ## distance_to_yeast48                                                      
    ## distance_to_yeast55                                                      
    ## YeastEMM_F3                                                              
    ## YeastEMM_F34                                                             
    ## YeastEMM_F47                                                             
    ## YeastEMM_F48                                                             
    ## YeastEMM_F49                                                             
    ## YeastEMM_F5                                                              
    ## YeastEMM_F63                                                             
    ## YeastEMM_F64                                                             
    ## YeastEMM_F65                                                             
    ## YeastEMM_F66                                                             
    ## YeastEMM_F70                                                             
    ## YeastEMM_F89                                                             
    ## YeastF7                                                                  
    ## YeastSP_F14                                                              
    ## YeastZAN_F3                                                              
    ## YeastZAN_F4                                                              
    ## distance_to_yeast17:YeastEMM_F3                                          
    ## distance_to_yeast25:YeastEMM_F3                                          
    ## distance_to_yeast32:YeastEMM_F3                                          
    ## distance_to_yeast41:YeastEMM_F3                                          
    ## distance_to_yeast48:YeastEMM_F3                                          
    ## distance_to_yeast55:YeastEMM_F3                                          
    ## distance_to_yeast17:YeastEMM_F34                                         
    ## distance_to_yeast25:YeastEMM_F34                                         
    ## distance_to_yeast32:YeastEMM_F34                                         
    ## distance_to_yeast41:YeastEMM_F34                                         
    ## distance_to_yeast48:YeastEMM_F34                                         
    ## distance_to_yeast55:YeastEMM_F34                                         
    ## distance_to_yeast17:YeastEMM_F47                                         
    ## distance_to_yeast25:YeastEMM_F47                                         
    ## distance_to_yeast32:YeastEMM_F47                                         
    ## distance_to_yeast41:YeastEMM_F47                                         
    ## distance_to_yeast48:YeastEMM_F47                                         
    ## distance_to_yeast55:YeastEMM_F47                                         
    ## distance_to_yeast17:YeastEMM_F48                                         
    ## distance_to_yeast25:YeastEMM_F48                                         
    ## distance_to_yeast32:YeastEMM_F48                                         
    ## distance_to_yeast41:YeastEMM_F48                                         
    ## distance_to_yeast48:YeastEMM_F48                                         
    ## distance_to_yeast55:YeastEMM_F48                                         
    ## distance_to_yeast17:YeastEMM_F49                                         
    ## distance_to_yeast25:YeastEMM_F49                                         
    ## distance_to_yeast32:YeastEMM_F49                                         
    ## distance_to_yeast41:YeastEMM_F49                                         
    ## distance_to_yeast48:YeastEMM_F49                                         
    ## distance_to_yeast55:YeastEMM_F49                                         
    ## distance_to_yeast17:YeastEMM_F5                                          
    ## distance_to_yeast25:YeastEMM_F5                                          
    ## distance_to_yeast32:YeastEMM_F5                                          
    ## distance_to_yeast41:YeastEMM_F5                                          
    ## distance_to_yeast48:YeastEMM_F5                                          
    ## distance_to_yeast55:YeastEMM_F5                                          
    ## distance_to_yeast17:YeastEMM_F63                                         
    ## distance_to_yeast25:YeastEMM_F63                                         
    ## distance_to_yeast32:YeastEMM_F63                                         
    ## distance_to_yeast41:YeastEMM_F63                                         
    ## distance_to_yeast48:YeastEMM_F63                                         
    ## distance_to_yeast55:YeastEMM_F63                                         
    ## distance_to_yeast17:YeastEMM_F64                                         
    ## distance_to_yeast25:YeastEMM_F64                                         
    ## distance_to_yeast32:YeastEMM_F64                                         
    ## distance_to_yeast41:YeastEMM_F64                                         
    ## distance_to_yeast48:YeastEMM_F64                                         
    ## distance_to_yeast55:YeastEMM_F64                                         
    ## distance_to_yeast17:YeastEMM_F65                                         
    ## distance_to_yeast25:YeastEMM_F65                                         
    ## distance_to_yeast32:YeastEMM_F65                                         
    ## distance_to_yeast41:YeastEMM_F65                                         
    ## distance_to_yeast48:YeastEMM_F65                                         
    ## distance_to_yeast55:YeastEMM_F65                                         
    ## distance_to_yeast17:YeastEMM_F66                                         
    ## distance_to_yeast25:YeastEMM_F66                                         
    ## distance_to_yeast32:YeastEMM_F66                                         
    ## distance_to_yeast41:YeastEMM_F66                                         
    ## distance_to_yeast48:YeastEMM_F66                                         
    ## distance_to_yeast55:YeastEMM_F66                                         
    ## distance_to_yeast17:YeastEMM_F70                                         
    ## distance_to_yeast25:YeastEMM_F70                                         
    ## distance_to_yeast32:YeastEMM_F70                                         
    ## distance_to_yeast41:YeastEMM_F70                                         
    ## distance_to_yeast48:YeastEMM_F70                                         
    ## distance_to_yeast55:YeastEMM_F70                                         
    ## distance_to_yeast17:YeastEMM_F89                                         
    ## distance_to_yeast25:YeastEMM_F89                                         
    ## distance_to_yeast32:YeastEMM_F89                                         
    ## distance_to_yeast41:YeastEMM_F89                                         
    ## distance_to_yeast48:YeastEMM_F89                                         
    ## distance_to_yeast55:YeastEMM_F89                                         
    ## distance_to_yeast17:YeastF7                                              
    ## distance_to_yeast25:YeastF7                                              
    ## distance_to_yeast32:YeastF7                                              
    ## distance_to_yeast41:YeastF7                                              
    ## distance_to_yeast48:YeastF7                                              
    ## distance_to_yeast55:YeastF7                                              
    ## distance_to_yeast17:YeastSP_F14                                          
    ## distance_to_yeast25:YeastSP_F14                                          
    ## distance_to_yeast32:YeastSP_F14                                          
    ## distance_to_yeast41:YeastSP_F14                                          
    ## distance_to_yeast48:YeastSP_F14                                          
    ## distance_to_yeast55:YeastSP_F14                                          
    ## distance_to_yeast17:YeastZAN_F3                                          
    ## distance_to_yeast25:YeastZAN_F3                                          
    ## distance_to_yeast32:YeastZAN_F3                                          
    ## distance_to_yeast41:YeastZAN_F3                                          
    ## distance_to_yeast48:YeastZAN_F3                                          
    ## distance_to_yeast55:YeastZAN_F3                                          
    ## distance_to_yeast17:YeastZAN_F4                                          
    ## distance_to_yeast25:YeastZAN_F4                                          
    ## distance_to_yeast32:YeastZAN_F4                                          
    ## distance_to_yeast41:YeastZAN_F4                                          
    ## distance_to_yeast48:YeastZAN_F4                                          
    ## distance_to_yeast55:YeastZAN_F4                                          
    ## DAI6:YeastEMM_F3                                                         
    ## DAI8:YeastEMM_F3                                                         
    ## DAI6:YeastEMM_F34                                                        
    ## DAI8:YeastEMM_F34                                                        
    ## DAI6:YeastEMM_F47                                                        
    ## DAI8:YeastEMM_F47                                                        
    ## DAI6:YeastEMM_F48                                                        
    ## DAI8:YeastEMM_F48                                                        
    ## DAI6:YeastEMM_F49                                                        
    ## DAI8:YeastEMM_F49                                                        
    ## DAI6:YeastEMM_F5                                                         
    ## DAI8:YeastEMM_F5                                                         
    ## DAI6:YeastEMM_F63                                                        
    ## DAI8:YeastEMM_F63                                                        
    ## DAI6:YeastEMM_F64                                                        
    ## DAI8:YeastEMM_F64                                                        
    ## DAI6:YeastEMM_F65                                                        
    ## DAI8:YeastEMM_F65                                                        
    ## DAI6:YeastEMM_F66                                                        
    ## DAI8:YeastEMM_F66                 0.430                                  
    ## DAI6:YeastEMM_F70                 0.308         0.154                    
    ## DAI8:YeastEMM_F70                 0.154         0.308         0.500      
    ## DAI6:YeastEMM_F89                 0.308         0.154         0.333      
    ## DAI8:YeastEMM_F89                 0.154         0.308         0.167      
    ## DAI6:YeastF7                      0.308         0.154         0.333      
    ## DAI8:YeastF7                      0.154         0.308         0.167      
    ## DAI6:YeastSP_F14                  0.308         0.154         0.333      
    ## DAI8:YeastSP_F14                  0.154         0.308         0.167      
    ## DAI6:YeastZAN_F3                  0.308         0.154         0.333      
    ## DAI8:YeastZAN_F3                  0.154         0.308         0.167      
    ## DAI6:YeastZAN_F4                  0.308         0.154         0.333      
    ## DAI8:YeastZAN_F4                  0.154         0.308         0.167      
    ##                                  DAI8:YEMM_F7 DAI6:YEMM_F8 DAI8:YEMM_F8 DAI6:YF
    ## DAI6                                                                           
    ## DAI8                                                                           
    ## distance_to_yeast17                                                            
    ## distance_to_yeast25                                                            
    ## distance_to_yeast32                                                            
    ## distance_to_yeast41                                                            
    ## distance_to_yeast48                                                            
    ## distance_to_yeast55                                                            
    ## YeastEMM_F3                                                                    
    ## YeastEMM_F34                                                                   
    ## YeastEMM_F47                                                                   
    ## YeastEMM_F48                                                                   
    ## YeastEMM_F49                                                                   
    ## YeastEMM_F5                                                                    
    ## YeastEMM_F63                                                                   
    ## YeastEMM_F64                                                                   
    ## YeastEMM_F65                                                                   
    ## YeastEMM_F66                                                                   
    ## YeastEMM_F70                                                                   
    ## YeastEMM_F89                                                                   
    ## YeastF7                                                                        
    ## YeastSP_F14                                                                    
    ## YeastZAN_F3                                                                    
    ## YeastZAN_F4                                                                    
    ## distance_to_yeast17:YeastEMM_F3                                                
    ## distance_to_yeast25:YeastEMM_F3                                                
    ## distance_to_yeast32:YeastEMM_F3                                                
    ## distance_to_yeast41:YeastEMM_F3                                                
    ## distance_to_yeast48:YeastEMM_F3                                                
    ## distance_to_yeast55:YeastEMM_F3                                                
    ## distance_to_yeast17:YeastEMM_F34                                               
    ## distance_to_yeast25:YeastEMM_F34                                               
    ## distance_to_yeast32:YeastEMM_F34                                               
    ## distance_to_yeast41:YeastEMM_F34                                               
    ## distance_to_yeast48:YeastEMM_F34                                               
    ## distance_to_yeast55:YeastEMM_F34                                               
    ## distance_to_yeast17:YeastEMM_F47                                               
    ## distance_to_yeast25:YeastEMM_F47                                               
    ## distance_to_yeast32:YeastEMM_F47                                               
    ## distance_to_yeast41:YeastEMM_F47                                               
    ## distance_to_yeast48:YeastEMM_F47                                               
    ## distance_to_yeast55:YeastEMM_F47                                               
    ## distance_to_yeast17:YeastEMM_F48                                               
    ## distance_to_yeast25:YeastEMM_F48                                               
    ## distance_to_yeast32:YeastEMM_F48                                               
    ## distance_to_yeast41:YeastEMM_F48                                               
    ## distance_to_yeast48:YeastEMM_F48                                               
    ## distance_to_yeast55:YeastEMM_F48                                               
    ## distance_to_yeast17:YeastEMM_F49                                               
    ## distance_to_yeast25:YeastEMM_F49                                               
    ## distance_to_yeast32:YeastEMM_F49                                               
    ## distance_to_yeast41:YeastEMM_F49                                               
    ## distance_to_yeast48:YeastEMM_F49                                               
    ## distance_to_yeast55:YeastEMM_F49                                               
    ## distance_to_yeast17:YeastEMM_F5                                                
    ## distance_to_yeast25:YeastEMM_F5                                                
    ## distance_to_yeast32:YeastEMM_F5                                                
    ## distance_to_yeast41:YeastEMM_F5                                                
    ## distance_to_yeast48:YeastEMM_F5                                                
    ## distance_to_yeast55:YeastEMM_F5                                                
    ## distance_to_yeast17:YeastEMM_F63                                               
    ## distance_to_yeast25:YeastEMM_F63                                               
    ## distance_to_yeast32:YeastEMM_F63                                               
    ## distance_to_yeast41:YeastEMM_F63                                               
    ## distance_to_yeast48:YeastEMM_F63                                               
    ## distance_to_yeast55:YeastEMM_F63                                               
    ## distance_to_yeast17:YeastEMM_F64                                               
    ## distance_to_yeast25:YeastEMM_F64                                               
    ## distance_to_yeast32:YeastEMM_F64                                               
    ## distance_to_yeast41:YeastEMM_F64                                               
    ## distance_to_yeast48:YeastEMM_F64                                               
    ## distance_to_yeast55:YeastEMM_F64                                               
    ## distance_to_yeast17:YeastEMM_F65                                               
    ## distance_to_yeast25:YeastEMM_F65                                               
    ## distance_to_yeast32:YeastEMM_F65                                               
    ## distance_to_yeast41:YeastEMM_F65                                               
    ## distance_to_yeast48:YeastEMM_F65                                               
    ## distance_to_yeast55:YeastEMM_F65                                               
    ## distance_to_yeast17:YeastEMM_F66                                               
    ## distance_to_yeast25:YeastEMM_F66                                               
    ## distance_to_yeast32:YeastEMM_F66                                               
    ## distance_to_yeast41:YeastEMM_F66                                               
    ## distance_to_yeast48:YeastEMM_F66                                               
    ## distance_to_yeast55:YeastEMM_F66                                               
    ## distance_to_yeast17:YeastEMM_F70                                               
    ## distance_to_yeast25:YeastEMM_F70                                               
    ## distance_to_yeast32:YeastEMM_F70                                               
    ## distance_to_yeast41:YeastEMM_F70                                               
    ## distance_to_yeast48:YeastEMM_F70                                               
    ## distance_to_yeast55:YeastEMM_F70                                               
    ## distance_to_yeast17:YeastEMM_F89                                               
    ## distance_to_yeast25:YeastEMM_F89                                               
    ## distance_to_yeast32:YeastEMM_F89                                               
    ## distance_to_yeast41:YeastEMM_F89                                               
    ## distance_to_yeast48:YeastEMM_F89                                               
    ## distance_to_yeast55:YeastEMM_F89                                               
    ## distance_to_yeast17:YeastF7                                                    
    ## distance_to_yeast25:YeastF7                                                    
    ## distance_to_yeast32:YeastF7                                                    
    ## distance_to_yeast41:YeastF7                                                    
    ## distance_to_yeast48:YeastF7                                                    
    ## distance_to_yeast55:YeastF7                                                    
    ## distance_to_yeast17:YeastSP_F14                                                
    ## distance_to_yeast25:YeastSP_F14                                                
    ## distance_to_yeast32:YeastSP_F14                                                
    ## distance_to_yeast41:YeastSP_F14                                                
    ## distance_to_yeast48:YeastSP_F14                                                
    ## distance_to_yeast55:YeastSP_F14                                                
    ## distance_to_yeast17:YeastZAN_F3                                                
    ## distance_to_yeast25:YeastZAN_F3                                                
    ## distance_to_yeast32:YeastZAN_F3                                                
    ## distance_to_yeast41:YeastZAN_F3                                                
    ## distance_to_yeast48:YeastZAN_F3                                                
    ## distance_to_yeast55:YeastZAN_F3                                                
    ## distance_to_yeast17:YeastZAN_F4                                                
    ## distance_to_yeast25:YeastZAN_F4                                                
    ## distance_to_yeast32:YeastZAN_F4                                                
    ## distance_to_yeast41:YeastZAN_F4                                                
    ## distance_to_yeast48:YeastZAN_F4                                                
    ## distance_to_yeast55:YeastZAN_F4                                                
    ## DAI6:YeastEMM_F3                                                               
    ## DAI8:YeastEMM_F3                                                               
    ## DAI6:YeastEMM_F34                                                              
    ## DAI8:YeastEMM_F34                                                              
    ## DAI6:YeastEMM_F47                                                              
    ## DAI8:YeastEMM_F47                                                              
    ## DAI6:YeastEMM_F48                                                              
    ## DAI8:YeastEMM_F48                                                              
    ## DAI6:YeastEMM_F49                                                              
    ## DAI8:YeastEMM_F49                                                              
    ## DAI6:YeastEMM_F5                                                               
    ## DAI8:YeastEMM_F5                                                               
    ## DAI6:YeastEMM_F63                                                              
    ## DAI8:YeastEMM_F63                                                              
    ## DAI6:YeastEMM_F64                                                              
    ## DAI8:YeastEMM_F64                                                              
    ## DAI6:YeastEMM_F65                                                              
    ## DAI8:YeastEMM_F65                                                              
    ## DAI6:YeastEMM_F66                                                              
    ## DAI8:YeastEMM_F66                                                              
    ## DAI6:YeastEMM_F70                                                              
    ## DAI8:YeastEMM_F70                                                              
    ## DAI6:YeastEMM_F89                 0.167                                        
    ## DAI8:YeastEMM_F89                 0.333        0.500                           
    ## DAI6:YeastF7                      0.167        0.333        0.167              
    ## DAI8:YeastF7                      0.333        0.167        0.333        0.500 
    ## DAI6:YeastSP_F14                  0.167        0.333        0.167        0.333 
    ## DAI8:YeastSP_F14                  0.333        0.167        0.333        0.167 
    ## DAI6:YeastZAN_F3                  0.167        0.333        0.167        0.333 
    ## DAI8:YeastZAN_F3                  0.333        0.167        0.333        0.167 
    ## DAI6:YeastZAN_F4                  0.167        0.333        0.167        0.333 
    ## DAI8:YeastZAN_F4                  0.333        0.167        0.333        0.167 
    ##                                  DAI8:YF DAI6:YS DAI8:YS DAI6:YZAN_F3
    ## DAI6                                                                 
    ## DAI8                                                                 
    ## distance_to_yeast17                                                  
    ## distance_to_yeast25                                                  
    ## distance_to_yeast32                                                  
    ## distance_to_yeast41                                                  
    ## distance_to_yeast48                                                  
    ## distance_to_yeast55                                                  
    ## YeastEMM_F3                                                          
    ## YeastEMM_F34                                                         
    ## YeastEMM_F47                                                         
    ## YeastEMM_F48                                                         
    ## YeastEMM_F49                                                         
    ## YeastEMM_F5                                                          
    ## YeastEMM_F63                                                         
    ## YeastEMM_F64                                                         
    ## YeastEMM_F65                                                         
    ## YeastEMM_F66                                                         
    ## YeastEMM_F70                                                         
    ## YeastEMM_F89                                                         
    ## YeastF7                                                              
    ## YeastSP_F14                                                          
    ## YeastZAN_F3                                                          
    ## YeastZAN_F4                                                          
    ## distance_to_yeast17:YeastEMM_F3                                      
    ## distance_to_yeast25:YeastEMM_F3                                      
    ## distance_to_yeast32:YeastEMM_F3                                      
    ## distance_to_yeast41:YeastEMM_F3                                      
    ## distance_to_yeast48:YeastEMM_F3                                      
    ## distance_to_yeast55:YeastEMM_F3                                      
    ## distance_to_yeast17:YeastEMM_F34                                     
    ## distance_to_yeast25:YeastEMM_F34                                     
    ## distance_to_yeast32:YeastEMM_F34                                     
    ## distance_to_yeast41:YeastEMM_F34                                     
    ## distance_to_yeast48:YeastEMM_F34                                     
    ## distance_to_yeast55:YeastEMM_F34                                     
    ## distance_to_yeast17:YeastEMM_F47                                     
    ## distance_to_yeast25:YeastEMM_F47                                     
    ## distance_to_yeast32:YeastEMM_F47                                     
    ## distance_to_yeast41:YeastEMM_F47                                     
    ## distance_to_yeast48:YeastEMM_F47                                     
    ## distance_to_yeast55:YeastEMM_F47                                     
    ## distance_to_yeast17:YeastEMM_F48                                     
    ## distance_to_yeast25:YeastEMM_F48                                     
    ## distance_to_yeast32:YeastEMM_F48                                     
    ## distance_to_yeast41:YeastEMM_F48                                     
    ## distance_to_yeast48:YeastEMM_F48                                     
    ## distance_to_yeast55:YeastEMM_F48                                     
    ## distance_to_yeast17:YeastEMM_F49                                     
    ## distance_to_yeast25:YeastEMM_F49                                     
    ## distance_to_yeast32:YeastEMM_F49                                     
    ## distance_to_yeast41:YeastEMM_F49                                     
    ## distance_to_yeast48:YeastEMM_F49                                     
    ## distance_to_yeast55:YeastEMM_F49                                     
    ## distance_to_yeast17:YeastEMM_F5                                      
    ## distance_to_yeast25:YeastEMM_F5                                      
    ## distance_to_yeast32:YeastEMM_F5                                      
    ## distance_to_yeast41:YeastEMM_F5                                      
    ## distance_to_yeast48:YeastEMM_F5                                      
    ## distance_to_yeast55:YeastEMM_F5                                      
    ## distance_to_yeast17:YeastEMM_F63                                     
    ## distance_to_yeast25:YeastEMM_F63                                     
    ## distance_to_yeast32:YeastEMM_F63                                     
    ## distance_to_yeast41:YeastEMM_F63                                     
    ## distance_to_yeast48:YeastEMM_F63                                     
    ## distance_to_yeast55:YeastEMM_F63                                     
    ## distance_to_yeast17:YeastEMM_F64                                     
    ## distance_to_yeast25:YeastEMM_F64                                     
    ## distance_to_yeast32:YeastEMM_F64                                     
    ## distance_to_yeast41:YeastEMM_F64                                     
    ## distance_to_yeast48:YeastEMM_F64                                     
    ## distance_to_yeast55:YeastEMM_F64                                     
    ## distance_to_yeast17:YeastEMM_F65                                     
    ## distance_to_yeast25:YeastEMM_F65                                     
    ## distance_to_yeast32:YeastEMM_F65                                     
    ## distance_to_yeast41:YeastEMM_F65                                     
    ## distance_to_yeast48:YeastEMM_F65                                     
    ## distance_to_yeast55:YeastEMM_F65                                     
    ## distance_to_yeast17:YeastEMM_F66                                     
    ## distance_to_yeast25:YeastEMM_F66                                     
    ## distance_to_yeast32:YeastEMM_F66                                     
    ## distance_to_yeast41:YeastEMM_F66                                     
    ## distance_to_yeast48:YeastEMM_F66                                     
    ## distance_to_yeast55:YeastEMM_F66                                     
    ## distance_to_yeast17:YeastEMM_F70                                     
    ## distance_to_yeast25:YeastEMM_F70                                     
    ## distance_to_yeast32:YeastEMM_F70                                     
    ## distance_to_yeast41:YeastEMM_F70                                     
    ## distance_to_yeast48:YeastEMM_F70                                     
    ## distance_to_yeast55:YeastEMM_F70                                     
    ## distance_to_yeast17:YeastEMM_F89                                     
    ## distance_to_yeast25:YeastEMM_F89                                     
    ## distance_to_yeast32:YeastEMM_F89                                     
    ## distance_to_yeast41:YeastEMM_F89                                     
    ## distance_to_yeast48:YeastEMM_F89                                     
    ## distance_to_yeast55:YeastEMM_F89                                     
    ## distance_to_yeast17:YeastF7                                          
    ## distance_to_yeast25:YeastF7                                          
    ## distance_to_yeast32:YeastF7                                          
    ## distance_to_yeast41:YeastF7                                          
    ## distance_to_yeast48:YeastF7                                          
    ## distance_to_yeast55:YeastF7                                          
    ## distance_to_yeast17:YeastSP_F14                                      
    ## distance_to_yeast25:YeastSP_F14                                      
    ## distance_to_yeast32:YeastSP_F14                                      
    ## distance_to_yeast41:YeastSP_F14                                      
    ## distance_to_yeast48:YeastSP_F14                                      
    ## distance_to_yeast55:YeastSP_F14                                      
    ## distance_to_yeast17:YeastZAN_F3                                      
    ## distance_to_yeast25:YeastZAN_F3                                      
    ## distance_to_yeast32:YeastZAN_F3                                      
    ## distance_to_yeast41:YeastZAN_F3                                      
    ## distance_to_yeast48:YeastZAN_F3                                      
    ## distance_to_yeast55:YeastZAN_F3                                      
    ## distance_to_yeast17:YeastZAN_F4                                      
    ## distance_to_yeast25:YeastZAN_F4                                      
    ## distance_to_yeast32:YeastZAN_F4                                      
    ## distance_to_yeast41:YeastZAN_F4                                      
    ## distance_to_yeast48:YeastZAN_F4                                      
    ## distance_to_yeast55:YeastZAN_F4                                      
    ## DAI6:YeastEMM_F3                                                     
    ## DAI8:YeastEMM_F3                                                     
    ## DAI6:YeastEMM_F34                                                    
    ## DAI8:YeastEMM_F34                                                    
    ## DAI6:YeastEMM_F47                                                    
    ## DAI8:YeastEMM_F47                                                    
    ## DAI6:YeastEMM_F48                                                    
    ## DAI8:YeastEMM_F48                                                    
    ## DAI6:YeastEMM_F49                                                    
    ## DAI8:YeastEMM_F49                                                    
    ## DAI6:YeastEMM_F5                                                     
    ## DAI8:YeastEMM_F5                                                     
    ## DAI6:YeastEMM_F63                                                    
    ## DAI8:YeastEMM_F63                                                    
    ## DAI6:YeastEMM_F64                                                    
    ## DAI8:YeastEMM_F64                                                    
    ## DAI6:YeastEMM_F65                                                    
    ## DAI8:YeastEMM_F65                                                    
    ## DAI6:YeastEMM_F66                                                    
    ## DAI8:YeastEMM_F66                                                    
    ## DAI6:YeastEMM_F70                                                    
    ## DAI8:YeastEMM_F70                                                    
    ## DAI6:YeastEMM_F89                                                    
    ## DAI8:YeastEMM_F89                                                    
    ## DAI6:YeastF7                                                         
    ## DAI8:YeastF7                                                         
    ## DAI6:YeastSP_F14                  0.167                              
    ## DAI8:YeastSP_F14                  0.333   0.500                      
    ## DAI6:YeastZAN_F3                  0.167   0.333   0.167              
    ## DAI8:YeastZAN_F3                  0.333   0.167   0.333   0.500      
    ## DAI6:YeastZAN_F4                  0.167   0.333   0.167   0.333      
    ## DAI8:YeastZAN_F4                  0.333   0.167   0.333   0.167      
    ##                                  DAI8:YZAN_F3 DAI6:YZAN_F4
    ## DAI6                                                      
    ## DAI8                                                      
    ## distance_to_yeast17                                       
    ## distance_to_yeast25                                       
    ## distance_to_yeast32                                       
    ## distance_to_yeast41                                       
    ## distance_to_yeast48                                       
    ## distance_to_yeast55                                       
    ## YeastEMM_F3                                               
    ## YeastEMM_F34                                              
    ## YeastEMM_F47                                              
    ## YeastEMM_F48                                              
    ## YeastEMM_F49                                              
    ## YeastEMM_F5                                               
    ## YeastEMM_F63                                              
    ## YeastEMM_F64                                              
    ## YeastEMM_F65                                              
    ## YeastEMM_F66                                              
    ## YeastEMM_F70                                              
    ## YeastEMM_F89                                              
    ## YeastF7                                                   
    ## YeastSP_F14                                               
    ## YeastZAN_F3                                               
    ## YeastZAN_F4                                               
    ## distance_to_yeast17:YeastEMM_F3                           
    ## distance_to_yeast25:YeastEMM_F3                           
    ## distance_to_yeast32:YeastEMM_F3                           
    ## distance_to_yeast41:YeastEMM_F3                           
    ## distance_to_yeast48:YeastEMM_F3                           
    ## distance_to_yeast55:YeastEMM_F3                           
    ## distance_to_yeast17:YeastEMM_F34                          
    ## distance_to_yeast25:YeastEMM_F34                          
    ## distance_to_yeast32:YeastEMM_F34                          
    ## distance_to_yeast41:YeastEMM_F34                          
    ## distance_to_yeast48:YeastEMM_F34                          
    ## distance_to_yeast55:YeastEMM_F34                          
    ## distance_to_yeast17:YeastEMM_F47                          
    ## distance_to_yeast25:YeastEMM_F47                          
    ## distance_to_yeast32:YeastEMM_F47                          
    ## distance_to_yeast41:YeastEMM_F47                          
    ## distance_to_yeast48:YeastEMM_F47                          
    ## distance_to_yeast55:YeastEMM_F47                          
    ## distance_to_yeast17:YeastEMM_F48                          
    ## distance_to_yeast25:YeastEMM_F48                          
    ## distance_to_yeast32:YeastEMM_F48                          
    ## distance_to_yeast41:YeastEMM_F48                          
    ## distance_to_yeast48:YeastEMM_F48                          
    ## distance_to_yeast55:YeastEMM_F48                          
    ## distance_to_yeast17:YeastEMM_F49                          
    ## distance_to_yeast25:YeastEMM_F49                          
    ## distance_to_yeast32:YeastEMM_F49                          
    ## distance_to_yeast41:YeastEMM_F49                          
    ## distance_to_yeast48:YeastEMM_F49                          
    ## distance_to_yeast55:YeastEMM_F49                          
    ## distance_to_yeast17:YeastEMM_F5                           
    ## distance_to_yeast25:YeastEMM_F5                           
    ## distance_to_yeast32:YeastEMM_F5                           
    ## distance_to_yeast41:YeastEMM_F5                           
    ## distance_to_yeast48:YeastEMM_F5                           
    ## distance_to_yeast55:YeastEMM_F5                           
    ## distance_to_yeast17:YeastEMM_F63                          
    ## distance_to_yeast25:YeastEMM_F63                          
    ## distance_to_yeast32:YeastEMM_F63                          
    ## distance_to_yeast41:YeastEMM_F63                          
    ## distance_to_yeast48:YeastEMM_F63                          
    ## distance_to_yeast55:YeastEMM_F63                          
    ## distance_to_yeast17:YeastEMM_F64                          
    ## distance_to_yeast25:YeastEMM_F64                          
    ## distance_to_yeast32:YeastEMM_F64                          
    ## distance_to_yeast41:YeastEMM_F64                          
    ## distance_to_yeast48:YeastEMM_F64                          
    ## distance_to_yeast55:YeastEMM_F64                          
    ## distance_to_yeast17:YeastEMM_F65                          
    ## distance_to_yeast25:YeastEMM_F65                          
    ## distance_to_yeast32:YeastEMM_F65                          
    ## distance_to_yeast41:YeastEMM_F65                          
    ## distance_to_yeast48:YeastEMM_F65                          
    ## distance_to_yeast55:YeastEMM_F65                          
    ## distance_to_yeast17:YeastEMM_F66                          
    ## distance_to_yeast25:YeastEMM_F66                          
    ## distance_to_yeast32:YeastEMM_F66                          
    ## distance_to_yeast41:YeastEMM_F66                          
    ## distance_to_yeast48:YeastEMM_F66                          
    ## distance_to_yeast55:YeastEMM_F66                          
    ## distance_to_yeast17:YeastEMM_F70                          
    ## distance_to_yeast25:YeastEMM_F70                          
    ## distance_to_yeast32:YeastEMM_F70                          
    ## distance_to_yeast41:YeastEMM_F70                          
    ## distance_to_yeast48:YeastEMM_F70                          
    ## distance_to_yeast55:YeastEMM_F70                          
    ## distance_to_yeast17:YeastEMM_F89                          
    ## distance_to_yeast25:YeastEMM_F89                          
    ## distance_to_yeast32:YeastEMM_F89                          
    ## distance_to_yeast41:YeastEMM_F89                          
    ## distance_to_yeast48:YeastEMM_F89                          
    ## distance_to_yeast55:YeastEMM_F89                          
    ## distance_to_yeast17:YeastF7                               
    ## distance_to_yeast25:YeastF7                               
    ## distance_to_yeast32:YeastF7                               
    ## distance_to_yeast41:YeastF7                               
    ## distance_to_yeast48:YeastF7                               
    ## distance_to_yeast55:YeastF7                               
    ## distance_to_yeast17:YeastSP_F14                           
    ## distance_to_yeast25:YeastSP_F14                           
    ## distance_to_yeast32:YeastSP_F14                           
    ## distance_to_yeast41:YeastSP_F14                           
    ## distance_to_yeast48:YeastSP_F14                           
    ## distance_to_yeast55:YeastSP_F14                           
    ## distance_to_yeast17:YeastZAN_F3                           
    ## distance_to_yeast25:YeastZAN_F3                           
    ## distance_to_yeast32:YeastZAN_F3                           
    ## distance_to_yeast41:YeastZAN_F3                           
    ## distance_to_yeast48:YeastZAN_F3                           
    ## distance_to_yeast55:YeastZAN_F3                           
    ## distance_to_yeast17:YeastZAN_F4                           
    ## distance_to_yeast25:YeastZAN_F4                           
    ## distance_to_yeast32:YeastZAN_F4                           
    ## distance_to_yeast41:YeastZAN_F4                           
    ## distance_to_yeast48:YeastZAN_F4                           
    ## distance_to_yeast55:YeastZAN_F4                           
    ## DAI6:YeastEMM_F3                                          
    ## DAI8:YeastEMM_F3                                          
    ## DAI6:YeastEMM_F34                                         
    ## DAI8:YeastEMM_F34                                         
    ## DAI6:YeastEMM_F47                                         
    ## DAI8:YeastEMM_F47                                         
    ## DAI6:YeastEMM_F48                                         
    ## DAI8:YeastEMM_F48                                         
    ## DAI6:YeastEMM_F49                                         
    ## DAI8:YeastEMM_F49                                         
    ## DAI6:YeastEMM_F5                                          
    ## DAI8:YeastEMM_F5                                          
    ## DAI6:YeastEMM_F63                                         
    ## DAI8:YeastEMM_F63                                         
    ## DAI6:YeastEMM_F64                                         
    ## DAI8:YeastEMM_F64                                         
    ## DAI6:YeastEMM_F65                                         
    ## DAI8:YeastEMM_F65                                         
    ## DAI6:YeastEMM_F66                                         
    ## DAI8:YeastEMM_F66                                         
    ## DAI6:YeastEMM_F70                                         
    ## DAI8:YeastEMM_F70                                         
    ## DAI6:YeastEMM_F89                                         
    ## DAI8:YeastEMM_F89                                         
    ## DAI6:YeastF7                                              
    ## DAI8:YeastF7                                              
    ## DAI6:YeastSP_F14                                          
    ## DAI8:YeastSP_F14                                          
    ## DAI6:YeastZAN_F3                                          
    ## DAI8:YeastZAN_F3                                          
    ## DAI6:YeastZAN_F4                  0.167                   
    ## DAI8:YeastZAN_F4                  0.333        0.500      
    ## 
    ## Standardized Within-Group Residuals:
    ##          Min           Q1          Med           Q3          Max 
    ## -4.221965716 -0.573871473  0.004985585  0.616073766  3.265101420 
    ## 
    ## Number of Observations: 1119
    ## Number of Groups: 3

``` r
#Anova(resultsB5)
anova(resultsB52)
```

    ##                         numDF denDF   F-value p-value
    ## (Intercept)                 1   964 1633.0371  <.0001
    ## DAI                         2   964  376.1636  <.0001
    ## distance_to_yeast           6   964    3.7083  0.0012
    ## Yeast                      16   964   35.0552  <.0001
    ## distance_to_yeast:Yeast    96   964    1.7051  0.0001
    ## DAI:Yeast                  32   964    1.9843  0.0010

### Loop for running analysis for each day separately for EMM_B52

Days after inoculation (DAI) as a factor is always significantly
impacting the growth. Also, our plot will represent data for Day 6 thus
we want relevant stats and comparison on Day 6 to present in the plot.
So, loop was made for each Day data and removing DAI from the model and
keeping rest of it present.

``` r
# Unique days
day <- unique(B52.no.1st.data$DAI)
# Initialize empty lists to store results
B52_Day <- list()
meansB52_Day <- list()
resultsB52_Day <- list()
model_summary_B52 <- list()
anova_table_B52 <- list()
Results_lsmeansB52 <- list()
filtered_comparisonsB52 <- list()

for (i in seq_along(day)) {
  # Filter data for the current day
  B52_Day[[i]] <- filter(B52.no.1st.data, DAI == day[i])
  # Fit the mixed effects model
  resultsB52_Day[[i]] <- lme(increase ~ distance_to_yeast * Yeast,
                             data = B52_Day[[i]],
                             random = ~1 | Replication,
                             na.action = na.omit)
  model_summary_B52[[i]] <- summary(resultsB52_Day[[i]])
  anova_table_B52[[i]] <- anova(resultsB52_Day[[i]])
  # Estimated marginal means
  meansB52_Day[[i]] <- emmeans(resultsB52_Day[[i]], ~ Yeast | distance_to_yeast)
  # Compact letter display with Tukey adjustment
  Results_lsmeansB52[[i]] <- multcomp::cld(meansB52_Day[[i]], alpha = 0.05,
                                           Letters = letters, reversed = TRUE,
                                           details = TRUE)
  # Print lsmeans results
  #print(Results_lsmeansB52[[i]])
  
  # Filter comparisons involving "Control"
  filtered_comparisonsB52[[i]] <- Results_lsmeansB52[[i]]$comparisons %>%
    filter(grepl("Control", contrast))
  # Print filtered contrasts
  print(filtered_comparisonsB52[[i]])
}
```

    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49    1.423 0.424 257   3.358  0.0751
    ##  Control - EMM_F47    1.323 0.424 257   3.122  0.1427
    ##  Control - EMM_F63    1.057 0.424 257   2.493  0.5100
    ##  Control - EMM_F48    1.017 0.424 257   2.398  0.5812
    ##  Control - EMM_F89    0.963 0.424 257   2.272  0.6745
    ##  Control - SP_F14     0.923 0.424 257   2.178  0.7402
    ##  Control - EMM_F66    0.863 0.424 257   2.037  0.8269
    ##  Control - EMM_F3     0.823 0.424 257   1.942  0.8746
    ##  Control - ZAN_F3     0.780 0.424 257   1.840  0.9161
    ##  Control - ZAN_F4     0.470 0.424 257   1.109  0.9996
    ##  Control - EMM_F70    0.450 0.424 257   1.062  0.9997
    ##  Control - EMM_F5     0.437 0.424 257   1.030  0.9998
    ##  Control - EMM_F65    0.380 0.424 257   0.896  1.0000
    ##  Control - EMM_F34    0.250 0.424 257   0.590  1.0000
    ##  Control - EMM_F64    0.060 0.424 257   0.142  1.0000
    ##  Control - F7         0.020 0.424 257   0.047  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65    1.380 0.424 257   3.255  0.1002
    ##  Control - EMM_F47    1.350 0.424 257   3.185  0.1212
    ##  Control - EMM_F48    1.253 0.424 257   2.957  0.2127
    ##  Control - EMM_F49    1.143 0.424 257   2.697  0.3635
    ##  Control - EMM_F3     1.120 0.424 257   2.642  0.4011
    ##  Control - EMM_F63    1.120 0.424 257   2.642  0.4011
    ##  Control - ZAN_F4     0.960 0.424 257   2.265  0.6801
    ##  Control - ZAN_F3     0.957 0.424 257   2.257  0.6858
    ##  Control - EMM_F89    0.950 0.424 257   2.241  0.6969
    ##  Control - SP_F14     0.843 0.424 257   1.989  0.8518
    ##  Control - EMM_F34    0.833 0.424 257   1.966  0.8635
    ##  Control - EMM_F64    0.827 0.424 257   1.950  0.8709
    ##  Control - EMM_F70    0.497 0.424 257   1.172  0.9991
    ##  Control - EMM_F5     0.420 0.424 257   0.991  0.9999
    ##  Control - EMM_F66    0.380 0.424 257   0.896  1.0000
    ##  Control - F7         0.250 0.424 257   0.590  1.0000
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14     1.655 0.424 257   3.904  0.0127
    ##  Control - EMM_F65    1.415 0.424 257   3.338  0.0795
    ##  Control - EMM_F48    1.365 0.424 257   3.220  0.1103
    ##  Control - EMM_F66    1.168 0.424 257   2.756  0.3251
    ##  Control - EMM_F89    1.142 0.424 257   2.693  0.3661
    ##  Control - EMM_F3     1.105 0.424 257   2.607  0.4262
    ##  Control - EMM_F49    1.088 0.424 257   2.567  0.4546
    ##  Control - EMM_F64    1.082 0.424 257   2.552  0.4661
    ##  Control - EMM_F47    1.075 0.424 257   2.536  0.4777
    ##  Control - ZAN_F3     1.025 0.424 257   2.418  0.5663
    ##  Control - EMM_F63    0.995 0.424 257   2.347  0.6196
    ##  Control - EMM_F34    0.968 0.424 257   2.284  0.6659
    ##  Control - ZAN_F4     0.962 0.424 257   2.268  0.6773
    ##  Control - EMM_F5     0.788 0.424 257   1.860  0.9089
    ##  Control - EMM_F70    0.652 0.424 257   1.537  0.9828
    ##  Control - F7         0.588 0.424 257   1.388  0.9940
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49    1.562 0.424 257   3.684  0.0272
    ##  Control - EMM_F64    1.478 0.424 257   3.487  0.0510
    ##  Control - EMM_F47    1.438 0.424 257   3.393  0.0677
    ##  Control - EMM_F48    1.072 0.424 257   2.528  0.4836
    ##  Control - EMM_F63    1.015 0.424 257   2.394  0.5841
    ##  Control - ZAN_F4     0.965 0.424 257   2.276  0.6716
    ##  Control - EMM_F89    0.898 0.424 257   2.119  0.7783
    ##  Control - EMM_F3     0.872 0.424 257   2.056  0.8159
    ##  Control - SP_F14     0.825 0.424 257   1.946  0.8728
    ##  Control - EMM_F70    0.712 0.424 257   1.679  0.9608
    ##  Control - ZAN_F3     0.615 0.424 257   1.451  0.9904
    ##  Control - EMM_F65    0.572 0.424 257   1.349  0.9956
    ##  Control - EMM_F34    0.522 0.424 257   1.231  0.9984
    ##  Control - EMM_F66    0.412 0.424 257   0.971  0.9999
    ##  Control - EMM_F5     0.318 0.424 257   0.751  1.0000
    ##  Control - F7         0.218 0.424 257   0.515  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65    1.752 0.424 257   4.132  0.0055
    ##  Control - EMM_F48    1.478 0.424 257   3.487  0.0510
    ##  Control - EMM_F47    1.455 0.424 257   3.432  0.0603
    ##  Control - EMM_F63    1.435 0.424 257   3.385  0.0693
    ##  Control - EMM_F64    1.355 0.424 257   3.196  0.1175
    ##  Control - ZAN_F3     1.212 0.424 257   2.858  0.2641
    ##  Control - ZAN_F4     1.145 0.424 257   2.701  0.3608
    ##  Control - EMM_F89    1.142 0.424 257   2.693  0.3661
    ##  Control - EMM_F49    1.142 0.424 257   2.693  0.3661
    ##  Control - EMM_F3     0.792 0.424 257   1.867  0.9060
    ##  Control - F7         0.748 0.424 257   1.765  0.9398
    ##  Control - SP_F14     0.678 0.424 257   1.600  0.9747
    ##  Control - EMM_F66    0.532 0.424 257   1.254  0.9981
    ##  Control - EMM_F70    0.488 0.424 257   1.152  0.9993
    ##  Control - EMM_F5     0.422 0.424 257   0.995  0.9999
    ##  Control - EMM_F34    0.342 0.424 257   0.806  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47    1.405 0.424 257   3.314  0.0850
    ##  Control - ZAN_F4     1.075 0.424 257   2.536  0.4777
    ##  Control - SP_F14     0.995 0.424 257   2.347  0.6196
    ##  Control - EMM_F3     0.962 0.424 257   2.268  0.6773
    ##  Control - EMM_F64    0.845 0.424 257   1.993  0.8498
    ##  Control - EMM_F48    0.825 0.424 257   1.946  0.8728
    ##  Control - EMM_F65    0.808 0.424 257   1.907  0.8902
    ##  Control - EMM_F89    0.762 0.424 257   1.797  0.9305
    ##  Control - EMM_F49    0.632 0.424 257   1.490  0.9874
    ##  Control - ZAN_F3     0.625 0.424 257   1.474  0.9886
    ##  Control - F7         0.395 0.424 257   0.932  1.0000
    ##  Control - EMM_F70    0.368 0.424 257   0.869  1.0000
    ##  Control - EMM_F63    0.188 0.424 257   0.444  1.0000
    ##  Control - EMM_F34    0.188 0.424 257   0.444  1.0000
    ##  Control - EMM_F5     0.055 0.424 257   0.130  1.0000
    ##  Control - EMM_F66    0.035 0.424 257   0.083  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F4     1.327 0.424 257   3.129  0.1399
    ##  Control - EMM_F47    1.207 0.424 257   2.846  0.2707
    ##  Control - SP_F14     1.140 0.424 257   2.689  0.3687
    ##  Control - EMM_F48    1.130 0.424 257   2.666  0.3848
    ##  Control - ZAN_F3     0.993 0.424 257   2.343  0.6225
    ##  Control - EMM_F49    0.990 0.424 257   2.335  0.6284
    ##  Control - EMM_F89    0.913 0.424 257   2.154  0.7558
    ##  Control - EMM_F70    0.743 0.424 257   1.753  0.9431
    ##  Control - EMM_F34    0.737 0.424 257   1.738  0.9472
    ##  Control - EMM_F64    0.717 0.424 257   1.691  0.9583
    ##  Control - F7         0.640 0.424 257   1.510  0.9856
    ##  Control - EMM_F3     0.630 0.424 257   1.486  0.9877
    ##  Control - EMM_F63    0.537 0.424 257   1.266  0.9978
    ##  Control - EMM_F5     0.510 0.424 257   1.203  0.9988
    ##  Control - EMM_F66    0.493 0.424 257   1.164  0.9992
    ##  Control - EMM_F65    0.143 0.424 257   0.338  1.0000
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 17 estimates 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  2.31500 0.550 250   4.212  0.0040
    ##  Control - EMM_F49  2.20500 0.550 250   4.012  0.0086
    ##  Control - SP_F14   2.07833 0.550 250   3.782  0.0196
    ##  Control - EMM_F48  1.64500 0.550 250   2.993  0.1956
    ##  Control - ZAN_F4   1.47833 0.550 250   2.690  0.3684
    ##  Control - EMM_F65  1.13167 0.550 250   2.059  0.8141
    ##  Control - EMM_F66  1.12439 0.635 250   1.771  0.9382
    ##  Control - EMM_F89  1.06167 0.550 250   1.932  0.8792
    ##  Control - EMM_F3   1.02500 0.550 250   1.865  0.9068
    ##  Control - ZAN_F3   0.82167 0.550 250   1.495  0.9869
    ##  Control - EMM_F5   0.78167 0.550 250   1.422  0.9922
    ##  Control - EMM_F63  0.65500 0.550 250   1.192  0.9989
    ##  Control - EMM_F34  0.41500 0.550 250   0.755  1.0000
    ##  Control - EMM_F70  0.16833 0.550 250   0.306  1.0000
    ##  Control - F7       0.11833 0.550 250   0.215  1.0000
    ##  EMM_F64 - Control  0.14500 0.550 250   0.264  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F48  1.65167 0.550 250   3.005  0.1901
    ##  Control - EMM_F47  1.59167 0.550 250   2.896  0.2436
    ##  Control - EMM_F49  1.54833 0.550 250   2.817  0.2878
    ##  Control - EMM_F65  1.32500 0.550 250   2.411  0.5717
    ##  Control - SP_F14   1.31833 0.550 250   2.399  0.5808
    ##  Control - ZAN_F4   1.15167 0.550 250   2.096  0.7927
    ##  Control - EMM_F70  1.06833 0.550 250   1.944  0.8737
    ##  Control - EMM_F89  0.93833 0.550 250   1.707  0.9545
    ##  Control - EMM_F3   0.87833 0.550 250   1.598  0.9749
    ##  Control - EMM_F63  0.65167 0.550 250   1.186  0.9990
    ##  Control - EMM_F66  0.62605 0.635 250   0.986  0.9999
    ##  Control - ZAN_F3   0.37167 0.550 250   0.676  1.0000
    ##  Control - EMM_F34  0.30833 0.550 250   0.561  1.0000
    ##  Control - EMM_F5   0.19167 0.550 250   0.349  1.0000
    ##  Control - EMM_F64  0.12167 0.550 250   0.221  1.0000
    ##  F7 - Control       0.11833 0.550 250   0.215  1.0000
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   1.48000 0.550 250   2.693  0.3664
    ##  Control - EMM_F66  1.35105 0.635 250   2.128  0.7729
    ##  Control - EMM_F48  1.08333 0.550 250   1.971  0.8608
    ##  Control - EMM_F49  1.05667 0.550 250   1.923  0.8832
    ##  Control - EMM_F65  1.03667 0.550 250   1.886  0.8985
    ##  Control - EMM_F64  0.65667 0.550 250   1.195  0.9989
    ##  Control - EMM_F47  0.64667 0.550 250   1.177  0.9991
    ##  Control - EMM_F70  0.45333 0.550 250   0.825  1.0000
    ##  Control - ZAN_F4   0.43333 0.550 250   0.788  1.0000
    ##  Control - EMM_F3   0.30667 0.550 250   0.558  1.0000
    ##  Control - EMM_F89  0.29333 0.550 250   0.534  1.0000
    ##  EMM_F34 - Control  0.00333 0.550 250   0.006  1.0000
    ##  F7 - Control       0.22000 0.550 250   0.400  1.0000
    ##  EMM_F63 - Control  0.22333 0.550 250   0.406  1.0000
    ##  ZAN_F3 - Control   0.30000 0.550 250   0.546  1.0000
    ##  EMM_F5 - Control   0.36667 0.550 250   0.667  1.0000
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F66  1.77439 0.635 250   2.794  0.3016
    ##  Control - SP_F14   1.70833 0.550 250   3.108  0.1478
    ##  Control - EMM_F47  1.62167 0.550 250   2.951  0.2157
    ##  Control - EMM_F49  1.50833 0.550 250   2.744  0.3326
    ##  Control - EMM_F70  1.42833 0.550 250   2.599  0.4318
    ##  Control - EMM_F48  1.13833 0.550 250   2.071  0.8071
    ##  Control - EMM_F64  1.11833 0.550 250   2.035  0.8278
    ##  Control - ZAN_F4   1.08167 0.550 250   1.968  0.8623
    ##  Control - EMM_F63  0.54500 0.550 250   0.992  0.9999
    ##  Control - EMM_F65  0.48833 0.550 250   0.889  1.0000
    ##  Control - EMM_F3   0.44167 0.550 250   0.804  1.0000
    ##  Control - EMM_F34  0.25833 0.550 250   0.470  1.0000
    ##  Control - EMM_F89  0.12833 0.550 250   0.234  1.0000
    ##  Control - ZAN_F3   0.10833 0.550 250   0.197  1.0000
    ##  Control - F7       0.07167 0.550 250   0.130  1.0000
    ##  Control - EMM_F5   0.06833 0.550 250   0.124  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65  1.93500 0.550 250   3.521  0.0461
    ##  Control - EMM_F47  1.65167 0.550 250   3.005  0.1901
    ##  Control - EMM_F48  1.60500 0.550 250   2.920  0.2310
    ##  Control - EMM_F66  1.57439 0.635 250   2.479  0.5200
    ##  Control - EMM_F49  1.34833 0.550 250   2.453  0.5396
    ##  Control - ZAN_F4   1.31500 0.550 250   2.393  0.5854
    ##  Control - SP_F14   1.29500 0.550 250   2.356  0.6127
    ##  Control - ZAN_F3   0.95500 0.550 250   1.738  0.9472
    ##  Control - EMM_F63  0.92500 0.550 250   1.683  0.9599
    ##  Control - EMM_F89  0.91500 0.550 250   1.665  0.9636
    ##  Control - EMM_F64  0.73500 0.550 250   1.337  0.9960
    ##  Control - F7       0.60500 0.550 250   1.101  0.9996
    ##  Control - EMM_F3   0.51500 0.550 250   0.937  1.0000
    ##  Control - EMM_F34  0.46500 0.550 250   0.846  1.0000
    ##  Control - EMM_F70  0.24833 0.550 250   0.452  1.0000
    ##  EMM_F5 - Control   0.13167 0.550 250   0.240  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   1.64833 0.550 250   2.999  0.1929
    ##  Control - EMM_F47  1.48500 0.550 250   2.702  0.3603
    ##  Control - EMM_F48  1.26167 0.550 250   2.296  0.6576
    ##  Control - ZAN_F4   1.24167 0.550 250   2.259  0.6839
    ##  Control - EMM_F49  1.07833 0.550 250   1.962  0.8652
    ##  Control - EMM_F66  0.93772 0.635 250   1.477  0.9884
    ##  Control - EMM_F65  0.90833 0.550 250   1.653  0.9659
    ##  Control - EMM_F3   0.75500 0.550 250   1.374  0.9946
    ##  Control - ZAN_F3   0.65500 0.550 250   1.192  0.9989
    ##  Control - EMM_F89  0.61500 0.550 250   1.119  0.9995
    ##  Control - EMM_F70  0.47833 0.550 250   0.870  1.0000
    ##  Control - F7       0.46833 0.550 250   0.852  1.0000
    ##  Control - EMM_F34  0.09167 0.550 250   0.167  1.0000
    ##  Control - EMM_F64  0.04500 0.550 250   0.082  1.0000
    ##  EMM_F63 - Control  0.05500 0.550 250   0.100  1.0000
    ##  EMM_F5 - Control   0.06833 0.550 250   0.124  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   2.08667 0.550 250   3.797  0.0186
    ##  Control - EMM_F47  1.91333 0.550 250   3.481  0.0521
    ##  Control - ZAN_F4   1.58000 0.550 250   2.875  0.2551
    ##  Control - EMM_F49  1.55333 0.550 250   2.826  0.2825
    ##  Control - EMM_F89  1.42667 0.550 250   2.596  0.4340
    ##  Control - EMM_F48  1.35667 0.550 250   2.469  0.5282
    ##  Control - EMM_F66  1.29772 0.635 250   2.044  0.8229
    ##  Control - EMM_F65  1.12000 0.550 250   2.038  0.8261
    ##  Control - EMM_F70  1.08667 0.550 250   1.977  0.8578
    ##  Control - ZAN_F3   1.07000 0.550 250   1.947  0.8723
    ##  Control - EMM_F64  0.83333 0.550 250   1.516  0.9849
    ##  Control - EMM_F34  0.79667 0.550 250   1.450  0.9904
    ##  Control - F7       0.76333 0.550 250   1.389  0.9939
    ##  Control - EMM_F3   0.51333 0.550 250   0.934  1.0000
    ##  Control - EMM_F5   0.49333 0.550 250   0.898  1.0000
    ##  Control - EMM_F63  0.32667 0.550 250   0.594  1.0000
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 17 estimates 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   2.39333 0.542 249   4.415  0.0018
    ##  Control - EMM_F47  2.25000 0.542 249   4.151  0.0051
    ##  Control - EMM_F49  2.18000 0.542 249   4.022  0.0083
    ##  Control - EMM_F48  1.68333 0.542 249   3.106  0.1489
    ##  Control - ZAN_F4   1.40000 0.542 249   2.583  0.4434
    ##  Control - EMM_F89  0.87333 0.542 249   1.611  0.9730
    ##  Control - EMM_F5   0.83667 0.542 249   1.544  0.9820
    ##  Control - EMM_F3   0.80333 0.542 249   1.482  0.9880
    ##  Control - EMM_F66  0.77772 0.626 249   1.242  0.9983
    ##  Control - EMM_F65  0.77000 0.542 249   1.421  0.9923
    ##  Control - ZAN_F3   0.70333 0.542 249   1.298  0.9971
    ##  Control - EMM_F63  0.58000 0.542 249   1.070  0.9997
    ##  Control - EMM_F34  0.57000 0.542 249   1.052  0.9998
    ##  F7 - Control       0.00667 0.542 249   0.012  1.0000
    ##  EMM_F64 - Control  0.10000 0.542 249   0.184  1.0000
    ##  EMM_F70 - Control  0.22333 0.542 249   0.412  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49  1.80333 0.542 249   3.327  0.0822
    ##  Control - EMM_F48  1.65667 0.542 249   3.056  0.1682
    ##  Control - EMM_F47  1.54333 0.542 249   2.847  0.2705
    ##  Control - SP_F14   1.43333 0.542 249   2.644  0.3996
    ##  Control - EMM_F70  1.21333 0.542 249   2.238  0.6987
    ##  Control - ZAN_F4   1.17000 0.542 249   2.159  0.7531
    ##  Control - EMM_F65  1.16333 0.542 249   2.146  0.7610
    ##  Control - EMM_F89  0.93333 0.542 249   1.722  0.9511
    ##  Control - EMM_F3   0.79333 0.542 249   1.464  0.9894
    ##  Control - EMM_F63  0.56333 0.542 249   1.039  0.9998
    ##  Control - EMM_F64  0.41000 0.542 249   0.756  1.0000
    ##  Control - ZAN_F3   0.36333 0.542 249   0.670  1.0000
    ##  Control - EMM_F66  0.12272 0.626 249   0.196  1.0000
    ##  Control - EMM_F5   0.08000 0.542 249   0.148  1.0000
    ##  EMM_F34 - Control  0.04000 0.542 249   0.074  1.0000
    ##  F7 - Control       0.30667 0.542 249   0.566  1.0000
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F66  1.45772 0.626 249   2.328  0.6337
    ##  Control - SP_F14   1.35167 0.542 249   2.494  0.5093
    ##  Control - EMM_F48  1.34500 0.542 249   2.481  0.5185
    ##  Control - EMM_F49  1.16833 0.542 249   2.155  0.7551
    ##  Control - EMM_F65  1.00500 0.542 249   1.854  0.9109
    ##  Control - EMM_F47  0.76167 0.542 249   1.405  0.9931
    ##  Control - EMM_F64  0.71167 0.542 249   1.313  0.9967
    ##  Control - ZAN_F4   0.60500 0.542 249   1.116  0.9995
    ##  Control - EMM_F3   0.37833 0.542 249   0.698  1.0000
    ##  Control - EMM_F89  0.30167 0.542 249   0.557  1.0000
    ##  Control - EMM_F34  0.12167 0.542 249   0.224  1.0000
    ##  Control - ZAN_F3   0.11500 0.542 249   0.212  1.0000
    ##  Control - EMM_F70  0.10833 0.542 249   0.200  1.0000
    ##  Control - EMM_F63  0.04833 0.542 249   0.089  1.0000
    ##  EMM_F5 - Control   0.04500 0.542 249   0.083  1.0000
    ##  F7 - Control       0.37167 0.542 249   0.686  1.0000
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  1.69667 0.542 249   3.130  0.1399
    ##  Control - SP_F14   1.51333 0.542 249   2.792  0.3030
    ##  Control - EMM_F49  1.28667 0.542 249   2.374  0.5996
    ##  Control - EMM_F48  1.26667 0.542 249   2.337  0.6272
    ##  Control - ZAN_F4   1.17667 0.542 249   2.171  0.7450
    ##  Control - EMM_F66  1.13772 0.626 249   1.817  0.9238
    ##  Control - EMM_F70  1.08000 0.542 249   1.992  0.8502
    ##  Control - EMM_F64  0.92667 0.542 249   1.710  0.9540
    ##  Control - EMM_F63  0.42667 0.542 249   0.787  1.0000
    ##  Control - EMM_F3   0.37667 0.542 249   0.695  1.0000
    ##  Control - EMM_F65  0.34667 0.542 249   0.640  1.0000
    ##  Control - EMM_F34  0.30333 0.542 249   0.560  1.0000
    ##  Control - EMM_F89  0.29333 0.542 249   0.541  1.0000
    ##  Control - ZAN_F3   0.22333 0.542 249   0.412  1.0000
    ##  Control - EMM_F5   0.18333 0.542 249   0.338  1.0000
    ##  Control - F7       0.08333 0.542 249   0.154  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65  1.57333 0.542 249   2.903  0.2402
    ##  Control - EMM_F48  1.45000 0.542 249   2.675  0.3784
    ##  Control - EMM_F47  1.40667 0.542 249   2.595  0.4345
    ##  Control - ZAN_F4   1.23333 0.542 249   2.275  0.6723
    ##  Control - EMM_F66  1.21105 0.626 249   1.934  0.8781
    ##  Control - EMM_F49  1.18333 0.542 249   2.183  0.7368
    ##  Control - EMM_F63  0.93667 0.542 249   1.728  0.9496
    ##  Control - SP_F14   0.92667 0.542 249   1.710  0.9540
    ##  Control - EMM_F89  0.73667 0.542 249   1.359  0.9952
    ##  Control - ZAN_F3   0.71667 0.542 249   1.322  0.9964
    ##  Control - EMM_F64  0.67667 0.542 249   1.248  0.9982
    ##  Control - F7       0.40333 0.542 249   0.744  1.0000
    ##  Control - EMM_F3   0.25333 0.542 249   0.467  1.0000
    ##  EMM_F34 - Control  0.06667 0.542 249   0.123  1.0000
    ##  EMM_F70 - Control  0.15000 0.542 249   0.277  1.0000
    ##  EMM_F5 - Control   0.38667 0.542 249   0.713  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  1.61000 0.542 249   2.970  0.2063
    ##  Control - ZAN_F4   1.55333 0.542 249   2.866  0.2601
    ##  Control - SP_F14   1.34667 0.542 249   2.484  0.5162
    ##  Control - EMM_F49  1.29333 0.542 249   2.386  0.5904
    ##  Control - EMM_F48  1.01667 0.542 249   1.876  0.9027
    ##  Control - EMM_F65  0.92667 0.542 249   1.710  0.9540
    ##  Control - EMM_F89  0.70000 0.542 249   1.291  0.9973
    ##  Control - EMM_F66  0.58772 0.626 249   0.939  0.9999
    ##  Control - ZAN_F3   0.49333 0.542 249   0.910  1.0000
    ##  Control - EMM_F3   0.47667 0.542 249   0.879  1.0000
    ##  Control - EMM_F70  0.40667 0.542 249   0.750  1.0000
    ##  Control - F7       0.27000 0.542 249   0.498  1.0000
    ##  Control - EMM_F64  0.26667 0.542 249   0.492  1.0000
    ##  Control - EMM_F63  0.26667 0.542 249   0.492  1.0000
    ##  Control - EMM_F34  0.09667 0.542 249   0.178  1.0000
    ##  EMM_F5 - Control   0.30333 0.542 249   0.560  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F4   2.27000 0.542 249   4.188  0.0044
    ##  Control - EMM_F47  2.26000 0.542 249   4.169  0.0048
    ##  Control - SP_F14   2.02667 0.542 249   3.739  0.0227
    ##  Control - EMM_F49  1.79000 0.542 249   3.302  0.0881
    ##  Control - EMM_F48  1.50667 0.542 249   2.780  0.3105
    ##  Control - ZAN_F3   1.42000 0.542 249   2.620  0.4169
    ##  Control - EMM_F89  1.40333 0.542 249   2.589  0.4389
    ##  Control - EMM_F65  1.32667 0.542 249   2.448  0.5440
    ##  Control - EMM_F64  1.05000 0.542 249   1.937  0.8768
    ##  Control - F7       0.96667 0.542 249   1.783  0.9345
    ##  Control - EMM_F34  0.80667 0.542 249   1.488  0.9875
    ##  Control - EMM_F63  0.79333 0.542 249   1.464  0.9894
    ##  Control - EMM_F5   0.67105 0.626 249   1.072  0.9997
    ##  Control - EMM_F70  0.50000 0.542 249   0.922  1.0000
    ##  Control - EMM_F3   0.41333 0.542 249   0.763  1.0000
    ##  EMM_F66 - Control  0.21895 0.626 249   0.350  1.0000
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 17 estimates

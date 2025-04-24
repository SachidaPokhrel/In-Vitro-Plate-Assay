### Things to know before Analysis

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

## ***Curtobacterium flaccumfaciens* EMM_B44**

Plot is generated using loop around the 4 different classes of yeast,
coming up with 4 plots as an output which will be combined in one plot.

``` r
#read data
B44 <- read.csv("CoCultureAssay/CoCultureAssayData/2024-08-09_PeaceAssay_B44.csv", na.strings = "na") 


#load cbb color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
str(B44)
```

    ## 'data.frame':    1536 obs. of  8 variables:
    ##  $ Bacteria         : chr  "B44" "B44" "B44" "B44" ...
    ##  $ Yeast            : chr  "EMM_F3" "EMM_F3" "EMM_F3" "EMM_F3" ...
    ##  $ Class            : chr  "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" ...
    ##  $ Replication      : int  1 1 1 1 1 1 1 1 2 2 ...
    ##  $ DAI              : int  2 2 2 2 2 2 2 2 2 2 ...
    ##  $ distance_to_yeast: num  0 11.4 17.4 24.9 31.8 ...
    ##  $ colony_diameter  : num  7.21 6.62 7.08 7.4 7.59 7.73 7.56 7.65 6.89 6.98 ...
    ##  $ increase         : num  0 0 0 0 0 0 0 0 0 0 ...

``` r
#round off of the distance to yeast so that it is visibly a categorical variable
B44$distance_to_yeast <- round(B44$distance_to_yeast, 0)
#change into categorical dataset
B44$Replication=as.factor(B44$Replication)
B44$Yeast=as.factor(B44$Yeast)
B44$DAI=as.factor(B44$DAI)
B44$distance_to_yeast=as.factor(B44$distance_to_yeast)

#remove the data for the distance between yeast and bacteria at zero to remove misleading details
B44.no.contact <- B44[B44$distance_to_yeast != 0, ]
head(B44.no.contact)
```

    ##   Bacteria  Yeast           Class Replication DAI distance_to_yeast
    ## 2      B44 EMM_F3 Dothideomycetes           1   2                11
    ## 3      B44 EMM_F3 Dothideomycetes           1   2                17
    ## 4      B44 EMM_F3 Dothideomycetes           1   2                25
    ## 5      B44 EMM_F3 Dothideomycetes           1   2                32
    ## 6      B44 EMM_F3 Dothideomycetes           1   2                41
    ## 7      B44 EMM_F3 Dothideomycetes           1   2                48
    ##   colony_diameter increase
    ## 2            6.62        0
    ## 3            7.08        0
    ## 4            7.40        0
    ## 5            7.59        0
    ## 6            7.73        0
    ## 7            7.56        0

``` r
#remove first data since it is a reference for increase in colony size and donot follow the assumption of linear model i.e. normal distributuion
B44.no.1st.data <-B44.no.contact[B44.no.contact$DAI != "2",]

# For plots purpose, I am only taking Day 6 since it has visible difference compared to other 
B44.Day6 <- B44.no.contact[B44.no.contact$DAI == "6",]

#making a list object to store all the plots
plot_listB44 <- list()

# Get unique classes from the data
unique_class <- unique(B44.Day6$Class[B44.Day6$Class != "Control"])

# Loop through each class
for (i in seq_along(unique_class)) {
  
  # Filter for Control and the current class for using it in one iteration
  data <- B44.Day6 %>%
    filter(Class %in% c("Control", unique_class[i]))
  
  # Generate the plot
  p <- ggplot(data, aes(x = distance_to_yeast, y = increase, group = Yeast, color = Yeast)) +
    stat_summary(fun = mean, geom = "line") +   #making a line graph
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +   #standard error in plot
    ylab("Increase in colony diameter (mm)") +   
    xlab("Distance to Yeast (mm)") +   
    theme_classic() +
    scale_color_manual(values = cbbPalette) + #providing color blind palette
    ggtitle(unique_class[i]) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Store the plot in the list
  plot_listB44[[i]] <- p
}

#combining each iteration into one plot for one bacteria or dataset
combined_plotB44 <- ggarrange(plotlist = plot_listB44, ncol = 2, nrow = 2, common.legend = FALSE) #using the list generated from the loop, use it to make a combined plot without common legend.
# Annotation for the title
final_plotB44 <- annotate_figure(combined_plotB44,
                                 top = text_grob(
    expression("Impact on growth of"~italic("Curtobacterium flaccumfaciens")~"EMM_B44 by Yeast"), color = "Blue2", face = "bold", size = 14, hjust = 0.5))

#display final plot
print(final_plotB44)
```

![](EMM_B44_files/figure-gfm/Plot%20for%20EMM_B44-1.png)<!-- -->

### Stats for *Curtobacterium flaccumfaciens* EMM_B44

We are using linear mixed model. Our dependent variable or y is increase
(increase in colony diameter from 1st data) and independent variables
are different Yeast isolates, days after inoculation (DAI), and distance
to yeast which is the distance between the yeast and bacterial colony in
the plate. Each plate is replicated 3 times.

***This applies for all the dataset for Co-Culture Plate Assay***

``` r
#filter data to remove 1st day data since the first data is taken as a base to measure the increase in colony size to rule out the variability that is caused by the drop inoculation. So, initially the increase in the colony diameter for 1st data for all colony is "0" that violates the assumption of normality, thus we remove that from analysis. This would be similar for all the bacterial isolates.

B44.no.1st.data <- B44.no.contact[B44.no.contact$DAI != "2",]
#We check interaction for all before cutting it down to the significant one.
B44try=lme(increase~DAI*distance_to_yeast*Yeast, data = B44.no.1st.data, random = ~1|Replication, na.action = na.omit)
anova(B44try)
```

    ##                             numDF denDF  F-value p-value
    ## (Intercept)                     1   670 90.89249  <.0001
    ## DAI                             2   670 30.61220  <.0001
    ## distance_to_yeast               6   670  8.56900  <.0001
    ## Yeast                          15   670  6.44533  <.0001
    ## DAI:distance_to_yeast          12   670  0.39830  0.9644
    ## DAI:Yeast                      30   670  0.54957  0.9768
    ## distance_to_yeast:Yeast        90   670  4.61931  <.0001
    ## DAI:distance_to_yeast:Yeast   180   670  0.42207  1.0000

``` r
#removed interaction terms that were not significant.
resultsB44=lme(increase~DAI+distance_to_yeast*Yeast, data = B44.no.1st.data, random = ~1|Replication, na.action = na.omit)
anova(B44try, resultsB44) #check if complex model was better or not
```

    ## Warning in anova.lme(B44try, resultsB44): fitted objects with different fixed
    ## effects. REML comparisons are not meaningful.

    ##            Model  df       AIC      BIC    logLik   Test  L.Ratio p-value
    ## B44try         1 338 1304.5791 2829.046 -314.2896                        
    ## resultsB44     2 116  697.3477 1253.650 -232.6738 1 vs 2 163.2314  0.9989

``` r
summary(resultsB44) #gives summary of the linear model
```

    ## Linear mixed-effects model fit by REML
    ##   Data: B44.no.1st.data 
    ##        AIC     BIC    logLik
    ##   697.3477 1253.65 -232.6739
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept)  Residual
    ## StdDev:  0.07293941 0.2709505
    ## 
    ## Fixed effects:  increase ~ DAI + distance_to_yeast * Yeast 
    ##                                       Value  Std.Error  DF    t-value p-value
    ## (Intercept)                       0.3254861 0.10038017 892   3.242534  0.0012
    ## DAI6                              0.0807143 0.02090428 892   3.861136  0.0001
    ## DAI8                              0.1761607 0.02090428 892   8.427015  0.0000
    ## distance_to_yeast17               0.1077778 0.12772729 892   0.843812  0.3990
    ## distance_to_yeast25               0.1622222 0.12772729 892   1.270067  0.2044
    ## distance_to_yeast32               0.0977778 0.12772729 892   0.765520  0.4442
    ## distance_to_yeast41               0.0822222 0.12772729 892   0.643733  0.5199
    ## distance_to_yeast48              -0.0866667 0.12772729 892  -0.678529  0.4976
    ## distance_to_yeast55               0.0522222 0.12772729 892   0.408857  0.6827
    ## YeastEMM_F3                       0.1388889 0.12772729 892   1.087386  0.2772
    ## YeastEMM_F34                     -0.0400000 0.12772729 892  -0.313167  0.7542
    ## YeastEMM_F47                      0.1044444 0.12772729 892   0.817714  0.4137
    ## YeastEMM_F48                      0.0166667 0.12772729 892   0.130486  0.8962
    ## YeastEMM_F49                     -0.1788889 0.12772729 892  -1.400553  0.1617
    ## YeastEMM_F63                     -0.2100000 0.12772729 892  -1.644128  0.1005
    ## YeastEMM_F64                     -0.0466667 0.12772729 892  -0.365362  0.7149
    ## YeastEMM_F65                     -0.1200000 0.12772729 892  -0.939502  0.3477
    ## YeastEMM_F66                     -0.0133333 0.12772729 892  -0.104389  0.9169
    ## YeastEMM_F7                      -0.0988889 0.12772729 892  -0.774219  0.4390
    ## YeastEMM_F70                      0.2588889 0.12772729 892   2.026888  0.0430
    ## YeastEMM_F89                     -0.1311111 0.12772729 892  -1.026493  0.3049
    ## YeastSP_F14                      -1.7966667 0.12772729 892 -14.066428  0.0000
    ## YeastZAN_F3                      -0.0444444 0.12772729 892  -0.347964  0.7279
    ## YeastZAN_F4                      -0.0588889 0.12772729 892  -0.461052  0.6449
    ## distance_to_yeast17:YeastEMM_F3  -0.2688889 0.18063366 892  -1.488587  0.1369
    ## distance_to_yeast25:YeastEMM_F3  -0.2288889 0.18063366 892  -1.267144  0.2054
    ## distance_to_yeast32:YeastEMM_F3  -0.0488889 0.18063366 892  -0.270652  0.7867
    ## distance_to_yeast41:YeastEMM_F3  -0.2055556 0.18063366 892  -1.137969  0.2554
    ## distance_to_yeast48:YeastEMM_F3   0.0355556 0.18063366 892   0.196838  0.8440
    ## distance_to_yeast55:YeastEMM_F3  -0.0255556 0.18063366 892  -0.141477  0.8875
    ## distance_to_yeast17:YeastEMM_F34 -0.3711111 0.18063366 892  -2.054496  0.0402
    ## distance_to_yeast25:YeastEMM_F34 -0.2611111 0.18063366 892  -1.445529  0.1487
    ## distance_to_yeast32:YeastEMM_F34 -0.0700000 0.18063366 892  -0.387525  0.6985
    ## distance_to_yeast41:YeastEMM_F34 -0.1944444 0.18063366 892  -1.076457  0.2820
    ## distance_to_yeast48:YeastEMM_F34  0.2411111 0.18063366 892   1.334807  0.1823
    ## distance_to_yeast55:YeastEMM_F34 -0.1088889 0.18063366 892  -0.602816  0.5468
    ## distance_to_yeast17:YeastEMM_F47 -0.2466667 0.18063366 892  -1.365563  0.1724
    ## distance_to_yeast25:YeastEMM_F47 -0.4800000 0.18063366 892  -2.657312  0.0080
    ## distance_to_yeast32:YeastEMM_F47 -0.2266667 0.18063366 892  -1.254842  0.2099
    ## distance_to_yeast41:YeastEMM_F47 -0.3866667 0.18063366 892  -2.140612  0.0326
    ## distance_to_yeast48:YeastEMM_F47 -0.3011111 0.18063366 892  -1.666971  0.0959
    ## distance_to_yeast55:YeastEMM_F47 -0.4533333 0.18063366 892  -2.509684  0.0123
    ## distance_to_yeast17:YeastEMM_F48 -0.0211111 0.18063366 892  -0.116873  0.9070
    ## distance_to_yeast25:YeastEMM_F48  0.0622222 0.18063366 892   0.344466  0.7306
    ## distance_to_yeast32:YeastEMM_F48  0.1811111 0.18063366 892   1.002643  0.3163
    ## distance_to_yeast41:YeastEMM_F48  0.0122222 0.18063366 892   0.067663  0.9461
    ## distance_to_yeast48:YeastEMM_F48  0.0877778 0.18063366 892   0.485944  0.6271
    ## distance_to_yeast55:YeastEMM_F48 -0.0466667 0.18063366 892  -0.258350  0.7962
    ## distance_to_yeast17:YeastEMM_F49  0.1377778 0.18063366 892   0.762747  0.4458
    ## distance_to_yeast25:YeastEMM_F49  0.1100000 0.18063366 892   0.608967  0.5427
    ## distance_to_yeast32:YeastEMM_F49  0.2722222 0.18063366 892   1.507040  0.1322
    ## distance_to_yeast41:YeastEMM_F49  0.0355556 0.18063366 892   0.196838  0.8440
    ## distance_to_yeast48:YeastEMM_F49  0.3366667 0.18063366 892   1.863809  0.0627
    ## distance_to_yeast55:YeastEMM_F49  0.2344444 0.18063366 892   1.297900  0.1947
    ## distance_to_yeast17:YeastEMM_F63  0.1911111 0.18063366 892   1.058004  0.2903
    ## distance_to_yeast25:YeastEMM_F63  0.0022222 0.18063366 892   0.012302  0.9902
    ## distance_to_yeast32:YeastEMM_F63  0.1655556 0.18063366 892   0.916527  0.3596
    ## distance_to_yeast41:YeastEMM_F63 -0.0788889 0.18063366 892  -0.436734  0.6624
    ## distance_to_yeast48:YeastEMM_F63  0.1977778 0.18063366 892   1.094911  0.2739
    ## distance_to_yeast55:YeastEMM_F63  0.0511111 0.18063366 892   0.282955  0.7773
    ## distance_to_yeast17:YeastEMM_F64 -0.0233333 0.18063366 892  -0.129175  0.8972
    ## distance_to_yeast25:YeastEMM_F64 -0.0333333 0.18063366 892  -0.184536  0.8536
    ## distance_to_yeast32:YeastEMM_F64  0.0944444 0.18063366 892   0.522851  0.6012
    ## distance_to_yeast41:YeastEMM_F64  0.0300000 0.18063366 892   0.166082  0.8681
    ## distance_to_yeast48:YeastEMM_F64  0.2977778 0.18063366 892   1.648518  0.0996
    ## distance_to_yeast55:YeastEMM_F64 -0.0600000 0.18063366 892  -0.332164  0.7398
    ## distance_to_yeast17:YeastEMM_F65  0.1888889 0.18063366 892   1.045701  0.2960
    ## distance_to_yeast25:YeastEMM_F65  0.0233333 0.18063366 892   0.129175  0.8972
    ## distance_to_yeast32:YeastEMM_F65 -0.1422222 0.18063366 892  -0.787352  0.4313
    ## distance_to_yeast41:YeastEMM_F65  0.0911111 0.18063366 892   0.504397  0.6141
    ## distance_to_yeast48:YeastEMM_F65  0.2244444 0.18063366 892   1.242539  0.2144
    ## distance_to_yeast55:YeastEMM_F65  0.3644444 0.18063366 892   2.017589  0.0439
    ## distance_to_yeast17:YeastEMM_F66 -0.0500000 0.18063366 892  -0.276803  0.7820
    ## distance_to_yeast25:YeastEMM_F66 -0.0600000 0.18063366 892  -0.332164  0.7398
    ## distance_to_yeast32:YeastEMM_F66 -0.0911111 0.18063366 892  -0.504397  0.6141
    ## distance_to_yeast41:YeastEMM_F66  0.0066667 0.18063366 892   0.036907  0.9706
    ## distance_to_yeast48:YeastEMM_F66  0.1022222 0.18063366 892   0.565909  0.5716
    ## distance_to_yeast55:YeastEMM_F66  0.3233333 0.18063366 892   1.789995  0.0738
    ## distance_to_yeast17:YeastEMM_F7  -0.2211111 0.18063366 892  -1.224086  0.2212
    ## distance_to_yeast25:YeastEMM_F7  -0.1877778 0.18063366 892  -1.039550  0.2988
    ## distance_to_yeast32:YeastEMM_F7   0.0555556 0.18063366 892   0.307559  0.7585
    ## distance_to_yeast41:YeastEMM_F7  -0.1188889 0.18063366 892  -0.658177  0.5106
    ## distance_to_yeast48:YeastEMM_F7   0.1466667 0.18063366 892   0.811956  0.4170
    ## distance_to_yeast55:YeastEMM_F7   0.0377778 0.18063366 892   0.209140  0.8344
    ## distance_to_yeast17:YeastEMM_F70 -0.0611111 0.18063366 892  -0.338315  0.7352
    ## distance_to_yeast25:YeastEMM_F70 -0.1133333 0.18063366 892  -0.627421  0.5305
    ## distance_to_yeast32:YeastEMM_F70 -0.3600000 0.18063366 892  -1.992984  0.0466
    ## distance_to_yeast41:YeastEMM_F70 -0.6277778 0.18063366 892  -3.475420  0.0005
    ## distance_to_yeast48:YeastEMM_F70 -0.1811111 0.18063366 892  -1.002643  0.3163
    ## distance_to_yeast55:YeastEMM_F70 -0.1788889 0.18063366 892  -0.990341  0.3223
    ## distance_to_yeast17:YeastEMM_F89  0.1544444 0.18063366 892   0.855015  0.3928
    ## distance_to_yeast25:YeastEMM_F89 -0.0188889 0.18063366 892  -0.104570  0.9167
    ## distance_to_yeast32:YeastEMM_F89  0.0933333 0.18063366 892   0.516700  0.6055
    ## distance_to_yeast41:YeastEMM_F89  0.0077778 0.18063366 892   0.043058  0.9657
    ## distance_to_yeast48:YeastEMM_F89  0.1400000 0.18063366 892   0.775049  0.4385
    ## distance_to_yeast55:YeastEMM_F89  0.0855556 0.18063366 892   0.473641  0.6359
    ## distance_to_yeast17:YeastSP_F14   1.8533333 0.18063366 892  10.260177  0.0000
    ## distance_to_yeast25:YeastSP_F14   1.9455556 0.18063366 892  10.770725  0.0000
    ## distance_to_yeast32:YeastSP_F14   1.7277778 0.18063366 892   9.565093  0.0000
    ## distance_to_yeast41:YeastSP_F14   1.7555556 0.18063366 892   9.718873  0.0000
    ## distance_to_yeast48:YeastSP_F14   1.8711111 0.18063366 892  10.358596  0.0000
    ## distance_to_yeast55:YeastSP_F14   1.6822222 0.18063366 892   9.312894  0.0000
    ## distance_to_yeast17:YeastZAN_F3  -0.1744444 0.18063366 892  -0.965736  0.3344
    ## distance_to_yeast25:YeastZAN_F3  -0.1811111 0.18063366 892  -1.002643  0.3163
    ## distance_to_yeast32:YeastZAN_F3  -0.0144444 0.18063366 892  -0.079965  0.9363
    ## distance_to_yeast41:YeastZAN_F3  -0.2088889 0.18063366 892  -1.156423  0.2478
    ## distance_to_yeast48:YeastZAN_F3  -0.0433333 0.18063366 892  -0.239896  0.8105
    ## distance_to_yeast55:YeastZAN_F3  -0.0322222 0.18063366 892  -0.178384  0.8585
    ## distance_to_yeast17:YeastZAN_F4  -0.1211111 0.18063366 892  -0.670479  0.5027
    ## distance_to_yeast25:YeastZAN_F4  -0.1011111 0.18063366 892  -0.559758  0.5758
    ## distance_to_yeast32:YeastZAN_F4   0.0888889 0.18063366 892   0.492095  0.6228
    ## distance_to_yeast41:YeastZAN_F4   0.1155556 0.18063366 892   0.639723  0.5225
    ## distance_to_yeast48:YeastZAN_F4   0.2455556 0.18063366 892   1.359412  0.1744
    ## distance_to_yeast55:YeastZAN_F4   0.0488889 0.18063366 892   0.270652  0.7867
    ##  Correlation: 
    ##                                  (Intr) DAI6   DAI8   ds__17 ds__25 ds__32
    ## DAI6                             -0.104                                   
    ## DAI8                             -0.104  0.500                            
    ## distance_to_yeast17              -0.636  0.000  0.000                     
    ## distance_to_yeast25              -0.636  0.000  0.000  0.500              
    ## distance_to_yeast32              -0.636  0.000  0.000  0.500  0.500       
    ## distance_to_yeast41              -0.636  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast48              -0.636  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast55              -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F3                      -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F34                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F47                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F48                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F49                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F63                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F64                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F65                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F66                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F7                      -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F70                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F89                     -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastSP_F14                      -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastZAN_F3                      -0.636  0.000  0.000  0.500  0.500  0.500
    ## YeastZAN_F4                      -0.636  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast17:YeastEMM_F3   0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F3   0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F3   0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F3   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F3   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F3   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F34  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F34  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F34  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F34  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F34  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F34  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F47  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F47  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F47  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F47  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F47  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F47  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F48  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F48  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F48  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F48  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F48  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F48  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F49  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F49  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F49  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F49  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F49  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F49  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F63  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F63  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F63  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F63  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F63  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F63  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F64  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F64  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F64  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F64  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F64  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F64  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F65  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F65  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F65  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F65  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F65  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F65  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F66  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F66  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F66  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F66  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F66  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F66  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F7   0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F7   0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F7   0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F7   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F7   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F7   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F70  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F70  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F70  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F70  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F70  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F70  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F89  0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F89  0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F89  0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F89  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F89  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F89  0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastSP_F14   0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastSP_F14   0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastSP_F14   0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastSP_F14   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastSP_F14   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastSP_F14   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastZAN_F3   0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastZAN_F3   0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastZAN_F3   0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastZAN_F3   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastZAN_F3   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastZAN_F3   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastZAN_F4   0.450  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastZAN_F4   0.450  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastZAN_F4   0.450  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastZAN_F4   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastZAN_F4   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastZAN_F4   0.450  0.000  0.000 -0.354 -0.354 -0.354
    ##                                  ds__41 ds__48 ds__55 YsEMM_F3 YEMM_F34
    ## DAI6                                                                   
    ## DAI8                                                                   
    ## distance_to_yeast17                                                    
    ## distance_to_yeast25                                                    
    ## distance_to_yeast32                                                    
    ## distance_to_yeast41                                                    
    ## distance_to_yeast48               0.500                                
    ## distance_to_yeast55               0.500  0.500                         
    ## YeastEMM_F3                       0.500  0.500  0.500                  
    ## YeastEMM_F34                      0.500  0.500  0.500  0.500           
    ## YeastEMM_F47                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F48                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F49                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F63                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F64                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F65                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F66                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F7                       0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F70                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F89                      0.500  0.500  0.500  0.500    0.500  
    ## YeastSP_F14                       0.500  0.500  0.500  0.500    0.500  
    ## YeastZAN_F3                       0.500  0.500  0.500  0.500    0.500  
    ## YeastZAN_F4                       0.500  0.500  0.500  0.500    0.500  
    ## distance_to_yeast17:YeastEMM_F3  -0.354 -0.354 -0.354 -0.707   -0.354  
    ## distance_to_yeast25:YeastEMM_F3  -0.354 -0.354 -0.354 -0.707   -0.354  
    ## distance_to_yeast32:YeastEMM_F3  -0.354 -0.354 -0.354 -0.707   -0.354  
    ## distance_to_yeast41:YeastEMM_F3  -0.707 -0.354 -0.354 -0.707   -0.354  
    ## distance_to_yeast48:YeastEMM_F3  -0.354 -0.707 -0.354 -0.707   -0.354  
    ## distance_to_yeast55:YeastEMM_F3  -0.354 -0.354 -0.707 -0.707   -0.354  
    ## distance_to_yeast17:YeastEMM_F34 -0.354 -0.354 -0.354 -0.354   -0.707  
    ## distance_to_yeast25:YeastEMM_F34 -0.354 -0.354 -0.354 -0.354   -0.707  
    ## distance_to_yeast32:YeastEMM_F34 -0.354 -0.354 -0.354 -0.354   -0.707  
    ## distance_to_yeast41:YeastEMM_F34 -0.707 -0.354 -0.354 -0.354   -0.707  
    ## distance_to_yeast48:YeastEMM_F34 -0.354 -0.707 -0.354 -0.354   -0.707  
    ## distance_to_yeast55:YeastEMM_F34 -0.354 -0.354 -0.707 -0.354   -0.707  
    ## distance_to_yeast17:YeastEMM_F47 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F47 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F47 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F47 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F47 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F47 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F48 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F48 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F48 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F48 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F48 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F48 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F49 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F49 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F49 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F49 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F49 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F49 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F63 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F63 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F63 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F63 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F63 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F63 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F64 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F64 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F64 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F64 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F64 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F64 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F65 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F65 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F65 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F65 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F65 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F65 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F66 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F66 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F66 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F66 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F66 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F66 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F7  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F7  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F7  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F7  -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F7  -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F7  -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F70 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F70 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F70 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F70 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F70 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F70 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F89 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F89 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F89 -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F89 -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F89 -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F89 -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastSP_F14  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastSP_F14  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastSP_F14  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastSP_F14  -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastSP_F14  -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastSP_F14  -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastZAN_F3  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastZAN_F3  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastZAN_F3  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastZAN_F3  -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastZAN_F3  -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastZAN_F3  -0.354 -0.354 -0.707 -0.354   -0.354  
    ## distance_to_yeast17:YeastZAN_F4  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastZAN_F4  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastZAN_F4  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastZAN_F4  -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastZAN_F4  -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastZAN_F4  -0.354 -0.354 -0.707 -0.354   -0.354  
    ##                                  YEMM_F47 YEMM_F48 YEMM_F49 YEMM_F63 YEMM_F64
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
    ## YeastEMM_F48                      0.500                                      
    ## YeastEMM_F49                      0.500    0.500                             
    ## YeastEMM_F63                      0.500    0.500    0.500                    
    ## YeastEMM_F64                      0.500    0.500    0.500    0.500           
    ## YeastEMM_F65                      0.500    0.500    0.500    0.500    0.500  
    ## YeastEMM_F66                      0.500    0.500    0.500    0.500    0.500  
    ## YeastEMM_F7                       0.500    0.500    0.500    0.500    0.500  
    ## YeastEMM_F70                      0.500    0.500    0.500    0.500    0.500  
    ## YeastEMM_F89                      0.500    0.500    0.500    0.500    0.500  
    ## YeastSP_F14                       0.500    0.500    0.500    0.500    0.500  
    ## YeastZAN_F3                       0.500    0.500    0.500    0.500    0.500  
    ## YeastZAN_F4                       0.500    0.500    0.500    0.500    0.500  
    ## distance_to_yeast17:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F63 -0.354   -0.354   -0.354   -0.707   -0.354  
    ## distance_to_yeast25:YeastEMM_F63 -0.354   -0.354   -0.354   -0.707   -0.354  
    ## distance_to_yeast32:YeastEMM_F63 -0.354   -0.354   -0.354   -0.707   -0.354  
    ## distance_to_yeast41:YeastEMM_F63 -0.354   -0.354   -0.354   -0.707   -0.354  
    ## distance_to_yeast48:YeastEMM_F63 -0.354   -0.354   -0.354   -0.707   -0.354  
    ## distance_to_yeast55:YeastEMM_F63 -0.354   -0.354   -0.354   -0.707   -0.354  
    ## distance_to_yeast17:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.707  
    ## distance_to_yeast25:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.707  
    ## distance_to_yeast32:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.707  
    ## distance_to_yeast41:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.707  
    ## distance_to_yeast48:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.707  
    ## distance_to_yeast55:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.707  
    ## distance_to_yeast17:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F66 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F66 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F66 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F66 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F66 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F66 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F7  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F7  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F7  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F7  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F7  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F7  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast17:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast25:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast32:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast41:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast48:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354  
    ## distance_to_yeast55:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354  
    ##                                  YEMM_F65 YEMM_F66 YsEMM_F7 YEMM_F70 YEMM_F8
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
    ## YeastEMM_F63                                                                
    ## YeastEMM_F64                                                                
    ## YeastEMM_F65                                                                
    ## YeastEMM_F66                      0.500                                     
    ## YeastEMM_F7                       0.500    0.500                            
    ## YeastEMM_F70                      0.500    0.500    0.500                   
    ## YeastEMM_F89                      0.500    0.500    0.500    0.500          
    ## YeastSP_F14                       0.500    0.500    0.500    0.500    0.500 
    ## YeastZAN_F3                       0.500    0.500    0.500    0.500    0.500 
    ## YeastZAN_F4                       0.500    0.500    0.500    0.500    0.500 
    ## distance_to_yeast17:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F47 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F47 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F47 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F47 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F47 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F47 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F48 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F48 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F48 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F48 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F48 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F48 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F49 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F49 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F49 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F49 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F49 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F49 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F65 -0.707   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F65 -0.707   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F65 -0.707   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F65 -0.707   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F65 -0.707   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F65 -0.707   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F66 -0.354   -0.707   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F66 -0.354   -0.707   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F66 -0.354   -0.707   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F66 -0.354   -0.707   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F66 -0.354   -0.707   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F66 -0.354   -0.707   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F7  -0.354   -0.354   -0.707   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F7  -0.354   -0.354   -0.707   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F7  -0.354   -0.354   -0.707   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F7  -0.354   -0.354   -0.707   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F7  -0.354   -0.354   -0.707   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F7  -0.354   -0.354   -0.707   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F70 -0.354   -0.354   -0.354   -0.707   -0.354 
    ## distance_to_yeast25:YeastEMM_F70 -0.354   -0.354   -0.354   -0.707   -0.354 
    ## distance_to_yeast32:YeastEMM_F70 -0.354   -0.354   -0.354   -0.707   -0.354 
    ## distance_to_yeast41:YeastEMM_F70 -0.354   -0.354   -0.354   -0.707   -0.354 
    ## distance_to_yeast48:YeastEMM_F70 -0.354   -0.354   -0.354   -0.707   -0.354 
    ## distance_to_yeast55:YeastEMM_F70 -0.354   -0.354   -0.354   -0.707   -0.354 
    ## distance_to_yeast17:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.707 
    ## distance_to_yeast25:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.707 
    ## distance_to_yeast32:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.707 
    ## distance_to_yeast41:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.707 
    ## distance_to_yeast48:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.707 
    ## distance_to_yeast55:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354   -0.707 
    ## distance_to_yeast17:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastSP_F14  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast17:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast25:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast32:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast41:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast48:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354 
    ## distance_to_yeast55:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354   -0.354 
    ##                                  YSP_F1 YZAN_F3 YZAN_F4 ds__17:YEMM_F3
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
    ## YeastEMM_F63                                                          
    ## YeastEMM_F64                                                          
    ## YeastEMM_F65                                                          
    ## YeastEMM_F66                                                          
    ## YeastEMM_F7                                                           
    ## YeastEMM_F70                                                          
    ## YeastEMM_F89                                                          
    ## YeastSP_F14                                                           
    ## YeastZAN_F3                       0.500                               
    ## YeastZAN_F4                       0.500  0.500                        
    ## distance_to_yeast17:YeastEMM_F3  -0.354 -0.354  -0.354                
    ## distance_to_yeast25:YeastEMM_F3  -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast32:YeastEMM_F3  -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast41:YeastEMM_F3  -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast48:YeastEMM_F3  -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast55:YeastEMM_F3  -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast17:YeastEMM_F34 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F34 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F34 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F34 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F34 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F34 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F47 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F47 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F47 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F47 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F47 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F47 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F48 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F48 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F48 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F48 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F48 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F48 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F49 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F49 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F49 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F49 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F49 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F49 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F63 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F63 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F63 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F63 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F63 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F63 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F64 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F64 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F64 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F64 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F64 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F64 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F65 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F65 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F65 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F65 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F65 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F65 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F66 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F66 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F66 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F66 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F66 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F66 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F7  -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F7  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F7  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F7  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F7  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F7  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F70 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F70 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F70 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F70 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F70 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F70 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F89 -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F89 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F89 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F89 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F89 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F89 -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastSP_F14  -0.707 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastSP_F14  -0.707 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastSP_F14  -0.707 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastSP_F14  -0.707 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastSP_F14  -0.707 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastSP_F14  -0.707 -0.354  -0.354   0.250        
    ## distance_to_yeast17:YeastZAN_F3  -0.354 -0.707  -0.354   0.500        
    ## distance_to_yeast25:YeastZAN_F3  -0.354 -0.707  -0.354   0.250        
    ## distance_to_yeast32:YeastZAN_F3  -0.354 -0.707  -0.354   0.250        
    ## distance_to_yeast41:YeastZAN_F3  -0.354 -0.707  -0.354   0.250        
    ## distance_to_yeast48:YeastZAN_F3  -0.354 -0.707  -0.354   0.250        
    ## distance_to_yeast55:YeastZAN_F3  -0.354 -0.707  -0.354   0.250        
    ## distance_to_yeast17:YeastZAN_F4  -0.354 -0.354  -0.707   0.500        
    ## distance_to_yeast25:YeastZAN_F4  -0.354 -0.354  -0.707   0.250        
    ## distance_to_yeast32:YeastZAN_F4  -0.354 -0.354  -0.707   0.250        
    ## distance_to_yeast41:YeastZAN_F4  -0.354 -0.354  -0.707   0.250        
    ## distance_to_yeast48:YeastZAN_F4  -0.354 -0.354  -0.707   0.250        
    ## distance_to_yeast55:YeastZAN_F4  -0.354 -0.354  -0.707   0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3   0.500                                      
    ## distance_to_yeast41:YeastEMM_F3   0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F3   0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F3   0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F34  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F34  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F34  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F47  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F47  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
    ## YeastSP_F14                                                                  
    ## YeastZAN_F3                                                                  
    ## YeastZAN_F4                                                                  
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3   0.500                                      
    ## distance_to_yeast17:YeastEMM_F34  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F34  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F34  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F34  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F34  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F34  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F47  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F47  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F47  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F47  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F47  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F47  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F47  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F48  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F49  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ##                                  d__48:YEMM_F49 d__55:YEMM_F49 d__17:YEMM_F63
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63  0.500                                      
    ## distance_to_yeast41:YeastEMM_F63  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F63  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F63                                             
    ## distance_to_yeast25:YeastEMM_F63                                             
    ## distance_to_yeast32:YeastEMM_F63                                             
    ## distance_to_yeast41:YeastEMM_F63                                             
    ## distance_to_yeast48:YeastEMM_F63                                             
    ## distance_to_yeast55:YeastEMM_F63  0.500                                      
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ##                                  d__48:YEMM_F66 d__55:YEMM_F66 ds__17:YEMM_F7
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
    ##                                  ds__25:YEMM_F7 ds__32:YEMM_F7 ds__41:YEMM_F7
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7                                              
    ## distance_to_yeast25:YeastEMM_F7                                              
    ## distance_to_yeast32:YeastEMM_F7   0.500                                      
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F7   0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ##                                  ds__48:YEMM_F7 ds__55:YEMM_F7 d__17:YEMM_F70
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7                                              
    ## distance_to_yeast25:YeastEMM_F7                                              
    ## distance_to_yeast32:YeastEMM_F7                                              
    ## distance_to_yeast41:YeastEMM_F7                                              
    ## distance_to_yeast48:YeastEMM_F7                                              
    ## distance_to_yeast55:YeastEMM_F7   0.500                                      
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250        
    ##                                  d__25:YEMM_F70 d__32:YEMM_F70 d__41:YEMM_F70
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
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
    ## YeastEMM_F66                                                                 
    ## YeastEMM_F7                                                                  
    ## YeastEMM_F70                                                                 
    ## YeastEMM_F89                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7                                              
    ## distance_to_yeast25:YeastEMM_F7                                              
    ## distance_to_yeast32:YeastEMM_F7                                              
    ## distance_to_yeast41:YeastEMM_F7                                              
    ## distance_to_yeast48:YeastEMM_F7                                              
    ## distance_to_yeast55:YeastEMM_F7                                              
    ## distance_to_yeast17:YeastEMM_F70                                             
    ## distance_to_yeast25:YeastEMM_F70                                             
    ## distance_to_yeast32:YeastEMM_F70  0.500                                      
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F70  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ##                                  d__48:YEMM_F70 d__55:YEMM_F70 d__17:YEMM_F8
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
    ## YeastEMM_F63                                                                
    ## YeastEMM_F64                                                                
    ## YeastEMM_F65                                                                
    ## YeastEMM_F66                                                                
    ## YeastEMM_F7                                                                 
    ## YeastEMM_F70                                                                
    ## YeastEMM_F89                                                                
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
    ## distance_to_yeast17:YeastEMM_F7                                             
    ## distance_to_yeast25:YeastEMM_F7                                             
    ## distance_to_yeast32:YeastEMM_F7                                             
    ## distance_to_yeast41:YeastEMM_F7                                             
    ## distance_to_yeast48:YeastEMM_F7                                             
    ## distance_to_yeast55:YeastEMM_F7                                             
    ## distance_to_yeast17:YeastEMM_F70                                            
    ## distance_to_yeast25:YeastEMM_F70                                            
    ## distance_to_yeast32:YeastEMM_F70                                            
    ## distance_to_yeast41:YeastEMM_F70                                            
    ## distance_to_yeast48:YeastEMM_F70                                            
    ## distance_to_yeast55:YeastEMM_F70  0.500                                     
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250                      
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.500       
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500       
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.500       
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.500       
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.500       
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500       
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250       
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250       
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250       
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250       
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250       
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.500       
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250       
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250       
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250       
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.250       
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.250       
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.500       
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250       
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250       
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250       
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.250       
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.250       
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
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F7                                                               
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
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
    ## distance_to_yeast17:YeastEMM_F7                                           
    ## distance_to_yeast25:YeastEMM_F7                                           
    ## distance_to_yeast32:YeastEMM_F7                                           
    ## distance_to_yeast41:YeastEMM_F7                                           
    ## distance_to_yeast48:YeastEMM_F7                                           
    ## distance_to_yeast55:YeastEMM_F7                                           
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
    ## distance_to_yeast17:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastSP_F14   0.500         0.250         0.250       
    ## distance_to_yeast32:YeastSP_F14   0.250         0.500         0.250       
    ## distance_to_yeast41:YeastSP_F14   0.250         0.250         0.500       
    ## distance_to_yeast48:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast17:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F3   0.500         0.250         0.250       
    ## distance_to_yeast32:YeastZAN_F3   0.250         0.500         0.250       
    ## distance_to_yeast41:YeastZAN_F3   0.250         0.250         0.500       
    ## distance_to_yeast48:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast17:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F4   0.500         0.250         0.250       
    ## distance_to_yeast32:YeastZAN_F4   0.250         0.500         0.250       
    ## distance_to_yeast41:YeastZAN_F4   0.250         0.250         0.500       
    ## distance_to_yeast48:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastZAN_F4   0.250         0.250         0.250       
    ##                                  d__48:YEMM_F8 d__55:YEMM_F8 d__17:YS d__25:YS
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
    ## YeastEMM_F63                                                                  
    ## YeastEMM_F64                                                                  
    ## YeastEMM_F65                                                                  
    ## YeastEMM_F66                                                                  
    ## YeastEMM_F7                                                                   
    ## YeastEMM_F70                                                                  
    ## YeastEMM_F89                                                                  
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
    ## distance_to_yeast17:YeastEMM_F7                                               
    ## distance_to_yeast25:YeastEMM_F7                                               
    ## distance_to_yeast32:YeastEMM_F7                                               
    ## distance_to_yeast41:YeastEMM_F7                                               
    ## distance_to_yeast48:YeastEMM_F7                                               
    ## distance_to_yeast55:YeastEMM_F7                                               
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
    ## distance_to_yeast17:YeastSP_F14   0.250         0.250                         
    ## distance_to_yeast25:YeastSP_F14   0.250         0.250         0.500           
    ## distance_to_yeast32:YeastSP_F14   0.250         0.250         0.500    0.500  
    ## distance_to_yeast41:YeastSP_F14   0.250         0.250         0.500    0.500  
    ## distance_to_yeast48:YeastSP_F14   0.500         0.250         0.500    0.500  
    ## distance_to_yeast55:YeastSP_F14   0.250         0.500         0.500    0.500  
    ## distance_to_yeast17:YeastZAN_F3   0.250         0.250         0.500    0.250  
    ## distance_to_yeast25:YeastZAN_F3   0.250         0.250         0.250    0.500  
    ## distance_to_yeast32:YeastZAN_F3   0.250         0.250         0.250    0.250  
    ## distance_to_yeast41:YeastZAN_F3   0.250         0.250         0.250    0.250  
    ## distance_to_yeast48:YeastZAN_F3   0.500         0.250         0.250    0.250  
    ## distance_to_yeast55:YeastZAN_F3   0.250         0.500         0.250    0.250  
    ## distance_to_yeast17:YeastZAN_F4   0.250         0.250         0.500    0.250  
    ## distance_to_yeast25:YeastZAN_F4   0.250         0.250         0.250    0.500  
    ## distance_to_yeast32:YeastZAN_F4   0.250         0.250         0.250    0.250  
    ## distance_to_yeast41:YeastZAN_F4   0.250         0.250         0.250    0.250  
    ## distance_to_yeast48:YeastZAN_F4   0.500         0.250         0.250    0.250  
    ## distance_to_yeast55:YeastZAN_F4   0.250         0.500         0.250    0.250  
    ##                                  d__32:YS d__41:YS d__48:YS d__55:YS
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
    ## YeastEMM_F63                                                        
    ## YeastEMM_F64                                                        
    ## YeastEMM_F65                                                        
    ## YeastEMM_F66                                                        
    ## YeastEMM_F7                                                         
    ## YeastEMM_F70                                                        
    ## YeastEMM_F89                                                        
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
    ## distance_to_yeast17:YeastEMM_F7                                     
    ## distance_to_yeast25:YeastEMM_F7                                     
    ## distance_to_yeast32:YeastEMM_F7                                     
    ## distance_to_yeast41:YeastEMM_F7                                     
    ## distance_to_yeast48:YeastEMM_F7                                     
    ## distance_to_yeast55:YeastEMM_F7                                     
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
    ## distance_to_yeast17:YeastSP_F14                                     
    ## distance_to_yeast25:YeastSP_F14                                     
    ## distance_to_yeast32:YeastSP_F14                                     
    ## distance_to_yeast41:YeastSP_F14   0.500                             
    ## distance_to_yeast48:YeastSP_F14   0.500    0.500                    
    ## distance_to_yeast55:YeastSP_F14   0.500    0.500    0.500           
    ## distance_to_yeast17:YeastZAN_F3   0.250    0.250    0.250    0.250  
    ## distance_to_yeast25:YeastZAN_F3   0.250    0.250    0.250    0.250  
    ## distance_to_yeast32:YeastZAN_F3   0.500    0.250    0.250    0.250  
    ## distance_to_yeast41:YeastZAN_F3   0.250    0.500    0.250    0.250  
    ## distance_to_yeast48:YeastZAN_F3   0.250    0.250    0.500    0.250  
    ## distance_to_yeast55:YeastZAN_F3   0.250    0.250    0.250    0.500  
    ## distance_to_yeast17:YeastZAN_F4   0.250    0.250    0.250    0.250  
    ## distance_to_yeast25:YeastZAN_F4   0.250    0.250    0.250    0.250  
    ## distance_to_yeast32:YeastZAN_F4   0.500    0.250    0.250    0.250  
    ## distance_to_yeast41:YeastZAN_F4   0.250    0.500    0.250    0.250  
    ## distance_to_yeast48:YeastZAN_F4   0.250    0.250    0.500    0.250  
    ## distance_to_yeast55:YeastZAN_F4   0.250    0.250    0.250    0.500  
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
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F7                                                               
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
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
    ## distance_to_yeast17:YeastEMM_F7                                           
    ## distance_to_yeast25:YeastEMM_F7                                           
    ## distance_to_yeast32:YeastEMM_F7                                           
    ## distance_to_yeast41:YeastEMM_F7                                           
    ## distance_to_yeast48:YeastEMM_F7                                           
    ## distance_to_yeast55:YeastEMM_F7                                           
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
    ## distance_to_yeast17:YeastZAN_F4   0.500         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F4   0.250         0.500         0.250       
    ## distance_to_yeast32:YeastZAN_F4   0.250         0.250         0.500       
    ## distance_to_yeast41:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast48:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastZAN_F4   0.250         0.250         0.250       
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
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F7                                                               
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
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
    ## distance_to_yeast17:YeastEMM_F7                                           
    ## distance_to_yeast25:YeastEMM_F7                                           
    ## distance_to_yeast32:YeastEMM_F7                                           
    ## distance_to_yeast41:YeastEMM_F7                                           
    ## distance_to_yeast48:YeastEMM_F7                                           
    ## distance_to_yeast55:YeastEMM_F7                                           
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
    ## distance_to_yeast17:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast32:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast41:YeastZAN_F4   0.500         0.250         0.250       
    ## distance_to_yeast48:YeastZAN_F4   0.250         0.500         0.250       
    ## distance_to_yeast55:YeastZAN_F4   0.250         0.250         0.500       
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
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
    ## YeastEMM_F66                                                              
    ## YeastEMM_F7                                                               
    ## YeastEMM_F70                                                              
    ## YeastEMM_F89                                                              
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
    ## distance_to_yeast17:YeastEMM_F7                                           
    ## distance_to_yeast25:YeastEMM_F7                                           
    ## distance_to_yeast32:YeastEMM_F7                                           
    ## distance_to_yeast41:YeastEMM_F7                                           
    ## distance_to_yeast48:YeastEMM_F7                                           
    ## distance_to_yeast55:YeastEMM_F7                                           
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
    ##                                  d__41:YZAN_F4 d__48:YZAN_F4
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
    ## YeastEMM_F63                                                
    ## YeastEMM_F64                                                
    ## YeastEMM_F65                                                
    ## YeastEMM_F66                                                
    ## YeastEMM_F7                                                 
    ## YeastEMM_F70                                                
    ## YeastEMM_F89                                                
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
    ## distance_to_yeast17:YeastEMM_F7                             
    ## distance_to_yeast25:YeastEMM_F7                             
    ## distance_to_yeast32:YeastEMM_F7                             
    ## distance_to_yeast41:YeastEMM_F7                             
    ## distance_to_yeast48:YeastEMM_F7                             
    ## distance_to_yeast55:YeastEMM_F7                             
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
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -9.5098052 -0.4955909 -0.0177572  0.4758110  8.6693639 
    ## 
    ## Number of Observations: 1008
    ## Number of Groups: 3

``` r
anova(resultsB44) #for overall x variables
```

    ##                         numDF denDF  F-value p-value
    ## (Intercept)                 1   892 90.89249  <.0001
    ## DAI                         2   892 35.59007  <.0001
    ## distance_to_yeast           6   892  9.96241  <.0001
    ## Yeast                      15   892  7.49340  <.0001
    ## distance_to_yeast:Yeast    90   892  5.37046  <.0001

### Loop for running analysis for each day separately for EMM_B44

Days after inoculation (DAI) as a factor is always significantly
impacting the growth. Also, our plot will represent data for Day 6 thus
we want relevant stats and comparison on Day 6 to present in the plot.
So, loop was made for each Day data and removing DAI from the model and
keeping rest of it present.

``` r
# Extract all the unique days to be used for iteration in for loop
day <- unique(B44.no.1st.data$DAI)
# Initialize empty lists to store results
B44_Day <- list()
meansB44_Day <- list()
resultsB44_Day <- list()
model_summary_B44 <- list()
anova_table_B44 <- list()
Results_lsmeansB44 <- list()
filtered_comparisonsB44 <- list()

for (i in seq_along(day)) { #there are DAI 4, 6 and 8 that was extracted as unique days; sequence assigned will be 1 , 2 and 3 for each unique DAI
  
   # Filter data for the current day
  B44_Day[[i]] <- filter(B44.no.1st.data, DAI == day[i])
    # Fit the mixed effects model
  resultsB44_Day[[i]] <- lme(increase ~ distance_to_yeast * Yeast,
                            data = B44_Day[[i]],
                            random = ~1 | Replication,
                            na.action = na.omit) #Using mixed effect model here subsetting the data for the particular day

    model_summary_B44[[i]] <- summary(resultsB44_Day[[i]]) #store in the list value
  anova_table_B44[[i]] <- anova(resultsB44_Day[[i]]) #store in list value
   # Estimated marginal means
  meansB44_Day[[i]] <- emmeans(resultsB44_Day[[i]], ~ Yeast | distance_to_yeast)
    # Compact letter display with Tukey adjustment
  Results_lsmeansB44[[i]] <- multcomp::cld(meansB44_Day[[i]], alpha = 0.05,
                                          Letters = letters, reversed = TRUE,
                                          details = TRUE)
   # Print lsmeans results
  #print(Results_lsmeansB44[[i]])
  
  # Filter comparisons involving "Control"
  filtered_comparisonsB44[[i]] <- Results_lsmeansB44[[i]]$comparisons %>%
    filter(grepl("Control", contrast))
   # Print filtered contrasts
  print(filtered_comparisonsB44[[i]])
}
```

    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.85667 0.188 222   4.568  0.0009
    ##  Control - EMM_F63  0.18667 0.188 222   0.995  0.9998
    ##  Control - EMM_F49  0.14667 0.188 222   0.782  1.0000
    ##  Control - EMM_F64  0.12333 0.188 222   0.658  1.0000
    ##  Control - EMM_F65  0.08000 0.188 222   0.427  1.0000
    ##  Control - EMM_F89  0.07667 0.188 222   0.409  1.0000
    ##  Control - ZAN_F3   0.06667 0.188 222   0.356  1.0000
    ##  Control - EMM_F66  0.03333 0.188 222   0.178  1.0000
    ##  Control - EMM_F7   0.01667 0.188 222   0.089  1.0000
    ##  EMM_F34 - Control  0.04667 0.188 222   0.249  1.0000
    ##  ZAN_F4 - Control   0.09000 0.188 222   0.480  1.0000
    ##  EMM_F48 - Control  0.11000 0.188 222   0.587  1.0000
    ##  EMM_F47 - Control  0.19000 0.188 222   1.013  0.9997
    ##  EMM_F70 - Control  0.20667 0.188 222   1.102  0.9993
    ##  EMM_F3 - Control   0.25000 0.188 222   1.333  0.9942
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34  0.36333 0.188 222   1.937  0.8555
    ##  Control - EMM_F7   0.21333 0.188 222   1.138  0.9990
    ##  Control - ZAN_F3   0.19333 0.188 222   1.031  0.9997
    ##  Control - EMM_F3   0.19000 0.188 222   1.013  0.9997
    ##  Control - EMM_F64  0.16000 0.188 222   0.853  1.0000
    ##  Control - EMM_F47  0.15000 0.188 222   0.800  1.0000
    ##  Control - SP_F14   0.13000 0.188 222   0.693  1.0000
    ##  Control - EMM_F66  0.12667 0.188 222   0.675  1.0000
    ##  Control - ZAN_F4   0.11333 0.188 222   0.604  1.0000
    ##  Control - EMM_F49  0.02333 0.188 222   0.124  1.0000
    ##  Control - EMM_F63  0.02333 0.188 222   0.124  1.0000
    ##  EMM_F65 - Control  0.01333 0.188 222   0.071  1.0000
    ##  EMM_F89 - Control  0.02000 0.188 222   0.107  1.0000
    ##  EMM_F48 - Control  0.08333 0.188 222   0.444  1.0000
    ##  EMM_F70 - Control  0.19667 0.188 222   1.049  0.9996
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34  0.37333 0.188 222   1.991  0.8278
    ##  Control - EMM_F47  0.34000 0.188 222   1.813  0.9093
    ##  Control - ZAN_F3   0.33000 0.188 222   1.760  0.9276
    ##  Control - EMM_F7   0.31000 0.188 222   1.653  0.9563
    ##  Control - EMM_F63  0.24000 0.188 222   1.280  0.9963
    ##  Control - EMM_F65  0.16000 0.188 222   0.853  1.0000
    ##  Control - ZAN_F4   0.15333 0.188 222   0.818  1.0000
    ##  Control - EMM_F89  0.14667 0.188 222   0.782  1.0000
    ##  Control - EMM_F3   0.14000 0.188 222   0.747  1.0000
    ##  Control - EMM_F66  0.13667 0.188 222   0.729  1.0000
    ##  Control - EMM_F64  0.13667 0.188 222   0.729  1.0000
    ##  Control - EMM_F49  0.12000 0.188 222   0.640  1.0000
    ##  Control - SP_F14   0.08667 0.188 222   0.462  1.0000
    ##  EMM_F70 - Control  0.09000 0.188 222   0.480  1.0000
    ##  EMM_F48 - Control  0.21667 0.188 222   1.155  0.9988
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65  0.34000 0.188 222   1.813  0.9093
    ##  Control - EMM_F34  0.22667 0.188 222   1.209  0.9980
    ##  Control - EMM_F66  0.22000 0.188 222   1.173  0.9986
    ##  Control - SP_F14   0.15333 0.188 222   0.818  1.0000
    ##  Control - EMM_F63  0.13667 0.188 222   0.729  1.0000
    ##  Control - EMM_F47  0.08000 0.188 222   0.427  1.0000
    ##  Control - EMM_F89  0.08000 0.188 222   0.427  1.0000
    ##  Control - EMM_F70  0.06000 0.188 222   0.320  1.0000
    ##  Control - EMM_F7   0.05000 0.188 222   0.267  1.0000
    ##  Control - ZAN_F3   0.04667 0.188 222   0.249  1.0000
    ##  Control - EMM_F64  0.03000 0.188 222   0.160  1.0000
    ##  Control - ZAN_F4   0.00667 0.188 222   0.036  1.0000
    ##  EMM_F3 - Control   0.06667 0.188 222   0.356  1.0000
    ##  EMM_F49 - Control  0.08000 0.188 222   0.427  1.0000
    ##  EMM_F48 - Control  0.32333 0.188 222   1.724  0.9383
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F70  0.32667 0.188 222   1.742  0.9331
    ##  Control - EMM_F34  0.31000 0.188 222   1.653  0.9563
    ##  Control - EMM_F47  0.27667 0.188 222   1.475  0.9842
    ##  Control - EMM_F63  0.21333 0.188 222   1.138  0.9990
    ##  Control - ZAN_F3   0.20333 0.188 222   1.084  0.9994
    ##  Control - EMM_F7   0.13000 0.188 222   0.693  1.0000
    ##  Control - EMM_F3   0.08333 0.188 222   0.444  1.0000
    ##  Control - EMM_F48  0.08000 0.188 222   0.427  1.0000
    ##  Control - SP_F14   0.08000 0.188 222   0.427  1.0000
    ##  Control - EMM_F66  0.04667 0.188 222   0.249  1.0000
    ##  Control - EMM_F64  0.04667 0.188 222   0.249  1.0000
    ##  Control - EMM_F65  0.02667 0.188 222   0.142  1.0000
    ##  Control - EMM_F49  0.02333 0.188 222   0.124  1.0000
    ##  EMM_F89 - Control  0.02000 0.188 222   0.107  1.0000
    ##  ZAN_F4 - Control   0.19667 0.188 222   1.049  0.9996
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.08333 0.188 222   0.444  1.0000
    ##  Control - ZAN_F3   0.02667 0.188 222   0.142  1.0000
    ##  EMM_F48 - Control  0.02000 0.188 222   0.107  1.0000
    ##  EMM_F70 - Control  0.04000 0.188 222   0.213  1.0000
    ##  EMM_F49 - Control  0.09667 0.188 222   0.515  1.0000
    ##  EMM_F65 - Control  0.10667 0.188 222   0.569  1.0000
    ##  SP_F14 - Control   0.11000 0.188 222   0.587  1.0000
    ##  EMM_F89 - Control  0.11000 0.188 222   0.587  1.0000
    ##  EMM_F66 - Control  0.11667 0.188 222   0.622  1.0000
    ##  EMM_F34 - Control  0.12667 0.188 222   0.675  1.0000
    ##  EMM_F63 - Control  0.12667 0.188 222   0.675  1.0000
    ##  EMM_F7 - Control   0.13667 0.188 222   0.729  1.0000
    ##  EMM_F3 - Control   0.15667 0.188 222   0.835  1.0000
    ##  ZAN_F4 - Control   0.32000 0.188 222   1.706  0.9432
    ##  EMM_F64 - Control  0.32333 0.188 222   1.724  0.9383
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.34667 0.188 222   1.849  0.8955
    ##  Control - EMM_F63  0.22333 0.188 222   1.191  0.9983
    ##  Control - EMM_F64  0.18333 0.188 222   0.978  0.9998
    ##  Control - EMM_F34  0.18000 0.188 222   0.960  0.9999
    ##  Control - SP_F14   0.06000 0.188 222   0.320  1.0000
    ##  ZAN_F3 - Control   0.07667 0.188 222   0.409  1.0000
    ##  EMM_F48 - Control  0.08333 0.188 222   0.444  1.0000
    ##  EMM_F7 - Control   0.09333 0.188 222   0.498  1.0000
    ##  ZAN_F4 - Control   0.11000 0.188 222   0.587  1.0000
    ##  EMM_F89 - Control  0.11667 0.188 222   0.622  1.0000
    ##  EMM_F49 - Control  0.14333 0.188 222   0.764  1.0000
    ##  EMM_F3 - Control   0.18333 0.188 222   0.978  0.9998
    ##  EMM_F65 - Control  0.20000 0.188 222   1.067  0.9995
    ##  EMM_F70 - Control  0.28000 0.188 222   1.493  0.9823
    ##  EMM_F66 - Control  0.32667 0.188 222   1.742  0.9331
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   2.42667 0.285 222   8.502  <.0001
    ##  Control - EMM_F49  0.24333 0.285 222   0.853  1.0000
    ##  Control - EMM_F63  0.23667 0.285 222   0.829  1.0000
    ##  Control - EMM_F89  0.18000 0.285 222   0.631  1.0000
    ##  Control - EMM_F7   0.16667 0.285 222   0.584  1.0000
    ##  Control - EMM_F34  0.12000 0.285 222   0.420  1.0000
    ##  Control - EMM_F65  0.11667 0.285 222   0.409  1.0000
    ##  Control - EMM_F64  0.06333 0.285 222   0.222  1.0000
    ##  Control - ZAN_F4   0.06000 0.285 222   0.210  1.0000
    ##  Control - ZAN_F3   0.03333 0.285 222   0.117  1.0000
    ##  Control - EMM_F48  0.02000 0.285 222   0.070  1.0000
    ##  EMM_F66 - Control  0.04667 0.285 222   0.164  1.0000
    ##  EMM_F47 - Control  0.08333 0.285 222   0.292  1.0000
    ##  EMM_F3 - Control   0.16333 0.285 222   0.572  1.0000
    ##  EMM_F70 - Control  0.41333 0.285 222   1.448  0.9868
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34  0.39000 0.285 222   1.366  0.9926
    ##  Control - EMM_F7   0.30000 0.285 222   1.051  0.9996
    ##  Control - ZAN_F3   0.20000 0.285 222   0.701  1.0000
    ##  Control - ZAN_F4   0.09333 0.285 222   0.327  1.0000
    ##  Control - EMM_F47  0.05333 0.285 222   0.187  1.0000
    ##  Control - EMM_F3   0.00333 0.285 222   0.012  1.0000
    ##  EMM_F63 - Control  0.04333 0.285 222   0.152  1.0000
    ##  EMM_F64 - Control  0.05333 0.285 222   0.187  1.0000
    ##  EMM_F66 - Control  0.05667 0.285 222   0.199  1.0000
    ##  EMM_F48 - Control  0.06000 0.285 222   0.210  1.0000
    ##  EMM_F49 - Control  0.06000 0.285 222   0.210  1.0000
    ##  EMM_F89 - Control  0.10667 0.285 222   0.374  1.0000
    ##  EMM_F65 - Control  0.19000 0.285 222   0.666  1.0000
    ##  EMM_F70 - Control  0.25667 0.285 222   0.899  0.9999
    ##  SP_F14 - Control   0.27333 0.285 222   0.958  0.9999
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.40667 0.285 222   1.425  0.9887
    ##  Control - EMM_F7   0.28333 0.285 222   0.993  0.9998
    ##  Control - EMM_F34  0.25667 0.285 222   0.899  0.9999
    ##  Control - ZAN_F3   0.18000 0.285 222   0.631  1.0000
    ##  Control - ZAN_F4   0.12667 0.285 222   0.444  1.0000
    ##  Control - EMM_F89  0.11667 0.285 222   0.409  1.0000
    ##  Control - EMM_F63  0.11333 0.285 222   0.397  1.0000
    ##  Control - EMM_F64  0.03333 0.285 222   0.117  1.0000
    ##  EMM_F66 - Control  0.00333 0.285 222   0.012  1.0000
    ##  EMM_F49 - Control  0.02667 0.285 222   0.093  1.0000
    ##  EMM_F48 - Control  0.04333 0.285 222   0.152  1.0000
    ##  EMM_F3 - Control   0.06667 0.285 222   0.234  1.0000
    ##  EMM_F65 - Control  0.06667 0.285 222   0.234  1.0000
    ##  SP_F14 - Control   0.07333 0.285 222   0.257  1.0000
    ##  EMM_F70 - Control  0.19667 0.285 222   0.689  1.0000
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.17000 0.285 222   0.596  1.0000
    ##  Control - EMM_F65  0.15667 0.285 222   0.549  1.0000
    ##  Control - EMM_F70  0.15000 0.285 222   0.526  1.0000
    ##  Control - ZAN_F3   0.13333 0.285 222   0.467  1.0000
    ##  Control - EMM_F7   0.09333 0.285 222   0.327  1.0000
    ##  Control - EMM_F34  0.06333 0.285 222   0.222  1.0000
    ##  Control - SP_F14   0.06000 0.285 222   0.210  1.0000
    ##  Control - EMM_F66  0.05667 0.285 222   0.199  1.0000
    ##  EMM_F63 - Control  0.02333 0.285 222   0.082  1.0000
    ##  EMM_F89 - Control  0.05667 0.285 222   0.199  1.0000
    ##  ZAN_F4 - Control   0.06333 0.285 222   0.222  1.0000
    ##  EMM_F64 - Control  0.08000 0.285 222   0.280  1.0000
    ##  EMM_F49 - Control  0.08333 0.285 222   0.292  1.0000
    ##  EMM_F3 - Control   0.19667 0.285 222   0.689  1.0000
    ##  EMM_F48 - Control  0.22667 0.285 222   0.794  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F70  0.37667 0.285 222   1.320  0.9948
    ##  Control - ZAN_F3   0.35000 0.285 222   1.226  0.9976
    ##  Control - EMM_F63  0.31333 0.285 222   1.098  0.9993
    ##  Control - EMM_F47  0.26667 0.285 222   0.934  0.9999
    ##  Control - EMM_F34  0.22000 0.285 222   0.771  1.0000
    ##  Control - EMM_F7   0.21667 0.285 222   0.759  1.0000
    ##  Control - EMM_F49  0.21000 0.285 222   0.736  1.0000
    ##  Control - EMM_F89  0.10333 0.285 222   0.362  1.0000
    ##  Control - SP_F14   0.01667 0.285 222   0.058  1.0000
    ##  Control - EMM_F64  0.01333 0.285 222   0.047  1.0000
    ##  Control - EMM_F65  0.01333 0.285 222   0.047  1.0000
    ##  Control - ZAN_F4   0.01000 0.285 222   0.035  1.0000
    ##  EMM_F48 - Control  0.01000 0.285 222   0.035  1.0000
    ##  EMM_F66 - Control  0.02000 0.285 222   0.070  1.0000
    ##  EMM_F3 - Control   0.04667 0.285 222   0.164  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.17333 0.285 222   0.607  1.0000
    ##  Control - ZAN_F3   0.02000 0.285 222   0.070  1.0000
    ##  EMM_F63 - Control  0.06667 0.285 222   0.234  1.0000
    ##  EMM_F89 - Control  0.07000 0.285 222   0.245  1.0000
    ##  EMM_F70 - Control  0.09667 0.285 222   0.339  1.0000
    ##  EMM_F7 - Control   0.12333 0.285 222   0.432  1.0000
    ##  SP_F14 - Control   0.14000 0.285 222   0.491  1.0000
    ##  EMM_F66 - Control  0.14333 0.285 222   0.502  1.0000
    ##  EMM_F65 - Control  0.18667 0.285 222   0.654  1.0000
    ##  ZAN_F4 - Control   0.21333 0.285 222   0.747  1.0000
    ##  EMM_F48 - Control  0.24000 0.285 222   0.841  1.0000
    ##  EMM_F49 - Control  0.28333 0.285 222   0.993  0.9998
    ##  EMM_F34 - Control  0.29333 0.285 222   1.028  0.9997
    ##  EMM_F64 - Control  0.30333 0.285 222   1.063  0.9995
    ##  EMM_F3 - Control   0.34000 0.285 222   1.191  0.9983
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.34667 0.285 222   1.215  0.9979
    ##  Control - EMM_F63  0.22667 0.285 222   0.794  1.0000
    ##  Control - EMM_F48  0.15333 0.285 222   0.537  1.0000
    ##  Control - SP_F14   0.13667 0.285 222   0.479  1.0000
    ##  Control - EMM_F34  0.12333 0.285 222   0.432  1.0000
    ##  Control - EMM_F64  0.11667 0.285 222   0.409  1.0000
    ##  Control - ZAN_F3   0.08667 0.285 222   0.304  1.0000
    ##  Control - EMM_F7   0.06333 0.285 222   0.222  1.0000
    ##  Control - ZAN_F4   0.04333 0.285 222   0.152  1.0000
    ##  Control - EMM_F89  0.04000 0.285 222   0.140  1.0000
    ##  EMM_F70 - Control  0.05000 0.285 222   0.175  1.0000
    ##  EMM_F49 - Control  0.11000 0.285 222   0.385  1.0000
    ##  EMM_F3 - Control   0.17000 0.285 222   0.596  1.0000
    ##  EMM_F65 - Control  0.24333 0.285 222   0.853  1.0000
    ##  EMM_F66 - Control  0.31333 0.285 222   1.098  0.9993
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   2.10667 0.233 222   9.026  <.0001
    ##  Control - EMM_F63  0.20667 0.233 222   0.886  1.0000
    ##  Control - ZAN_F4   0.20667 0.233 222   0.886  1.0000
    ##  Control - EMM_F65  0.16333 0.233 222   0.700  1.0000
    ##  Control - EMM_F49  0.14667 0.233 222   0.628  1.0000
    ##  Control - EMM_F89  0.13667 0.233 222   0.586  1.0000
    ##  Control - EMM_F7   0.11333 0.233 222   0.486  1.0000
    ##  Control - EMM_F66  0.05333 0.233 222   0.229  1.0000
    ##  Control - EMM_F34  0.04667 0.233 222   0.200  1.0000
    ##  Control - EMM_F48  0.04000 0.233 222   0.171  1.0000
    ##  Control - ZAN_F3   0.03333 0.233 222   0.143  1.0000
    ##  EMM_F3 - Control   0.00333 0.233 222   0.014  1.0000
    ##  EMM_F47 - Control  0.04000 0.233 222   0.171  1.0000
    ##  EMM_F64 - Control  0.04667 0.233 222   0.200  1.0000
    ##  EMM_F70 - Control  0.15667 0.233 222   0.671  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34  0.48000 0.233 222   2.057  0.7900
    ##  Control - EMM_F7   0.44667 0.233 222   1.914  0.8669
    ##  Control - ZAN_F4   0.33333 0.233 222   1.428  0.9885
    ##  Control - ZAN_F3   0.26333 0.233 222   1.128  0.9991
    ##  Control - EMM_F47  0.22333 0.233 222   0.957  0.9999
    ##  Control - EMM_F3   0.19667 0.233 222   0.843  1.0000
    ##  Control - EMM_F49  0.16000 0.233 222   0.686  1.0000
    ##  Control - EMM_F48  0.15667 0.233 222   0.671  1.0000
    ##  Control - EMM_F66  0.12000 0.233 222   0.514  1.0000
    ##  Control - EMM_F64  0.10333 0.233 222   0.443  1.0000
    ##  Control - EMM_F63  0.07667 0.233 222   0.328  1.0000
    ##  Control - EMM_F89  0.05667 0.233 222   0.243  1.0000
    ##  EMM_F65 - Control  0.00333 0.233 222   0.014  1.0000
    ##  SP_F14 - Control   0.02667 0.233 222   0.114  1.0000
    ##  EMM_F70 - Control  0.14000 0.233 222   0.600  1.0000
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.38000 0.233 222   1.628  0.9615
    ##  Control - EMM_F34  0.27333 0.233 222   1.171  0.9986
    ##  Control - EMM_F63  0.27000 0.233 222   1.157  0.9988
    ##  Control - EMM_F7   0.26667 0.233 222   1.143  0.9989
    ##  Control - ZAN_F4   0.20000 0.233 222   0.857  1.0000
    ##  Control - EMM_F3   0.19667 0.233 222   0.843  1.0000
    ##  Control - EMM_F65  0.19667 0.233 222   0.843  1.0000
    ##  Control - EMM_F89  0.18667 0.233 222   0.800  1.0000
    ##  Control - ZAN_F3   0.16667 0.233 222   0.714  1.0000
    ##  Control - EMM_F49  0.11333 0.233 222   0.486  1.0000
    ##  Control - EMM_F66  0.08667 0.233 222   0.371  1.0000
    ##  Control - EMM_F64  0.07000 0.233 222   0.300  1.0000
    ##  Control - EMM_F48  0.02333 0.233 222   0.100  1.0000
    ##  EMM_F70 - Control  0.15000 0.233 222   0.643  1.0000
    ##  SP_F14 - Control   0.46000 0.233 222   1.971  0.8384
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65  0.29000 0.233 222   1.243  0.9973
    ##  Control - EMM_F47  0.11667 0.233 222   0.500  1.0000
    ##  Control - EMM_F70  0.09333 0.233 222   0.400  1.0000
    ##  Control - EMM_F89  0.09000 0.233 222   0.386  1.0000
    ##  Control - EMM_F34  0.04000 0.233 222   0.171  1.0000
    ##  Control - EMM_F66  0.03667 0.233 222   0.157  1.0000
    ##  Control - EMM_F63  0.02000 0.233 222   0.086  1.0000
    ##  ZAN_F3 - Control   0.00333 0.233 222   0.014  1.0000
    ##  SP_F14 - Control   0.00667 0.233 222   0.029  1.0000
    ##  EMM_F3 - Control   0.00667 0.233 222   0.029  1.0000
    ##  EMM_F7 - Control   0.01333 0.233 222   0.057  1.0000
    ##  ZAN_F4 - Control   0.03333 0.233 222   0.143  1.0000
    ##  EMM_F48 - Control  0.04333 0.233 222   0.186  1.0000
    ##  EMM_F64 - Control  0.09333 0.233 222   0.400  1.0000
    ##  EMM_F49 - Control  0.11667 0.233 222   0.500  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F70  0.40333 0.233 222   1.728  0.9372
    ##  Control - EMM_F63  0.34000 0.233 222   1.457  0.9860
    ##  Control - EMM_F7   0.30667 0.233 222   1.314  0.9950
    ##  Control - EMM_F47  0.30333 0.233 222   1.300  0.9956
    ##  Control - EMM_F89  0.28667 0.233 222   1.228  0.9976
    ##  Control - ZAN_F3   0.20667 0.233 222   0.886  1.0000
    ##  Control - EMM_F49  0.19667 0.233 222   0.843  1.0000
    ##  Control - EMM_F34  0.17333 0.233 222   0.743  1.0000
    ##  Control - EMM_F3   0.16333 0.233 222   0.700  1.0000
    ##  Control - EMM_F65  0.04667 0.233 222   0.200  1.0000
    ##  Control - SP_F14   0.02667 0.233 222   0.114  1.0000
    ##  Control - ZAN_F4   0.01667 0.233 222   0.071  1.0000
    ##  EMM_F66 - Control  0.00667 0.233 222   0.029  1.0000
    ##  EMM_F64 - Control  0.01000 0.233 222   0.043  1.0000
    ##  EMM_F48 - Control  0.15667 0.233 222   0.671  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.33333 0.233 222   1.428  0.9885
    ##  Control - EMM_F63  0.23000 0.233 222   0.985  0.9998
    ##  Control - ZAN_F3   0.21667 0.233 222   0.928  0.9999
    ##  Control - EMM_F89  0.15333 0.233 222   0.657  1.0000
    ##  Control - EMM_F7   0.11667 0.233 222   0.500  1.0000
    ##  Control - SP_F14   0.02667 0.233 222   0.114  1.0000
    ##  EMM_F66 - Control  0.00667 0.233 222   0.029  1.0000
    ##  EMM_F65 - Control  0.02000 0.233 222   0.086  1.0000
    ##  EMM_F3 - Control   0.02667 0.233 222   0.114  1.0000
    ##  ZAN_F4 - Control   0.02667 0.233 222   0.114  1.0000
    ##  EMM_F48 - Control  0.05333 0.233 222   0.229  1.0000
    ##  EMM_F49 - Control  0.09333 0.233 222   0.400  1.0000
    ##  EMM_F70 - Control  0.09667 0.233 222   0.414  1.0000
    ##  EMM_F64 - Control  0.12667 0.233 222   0.543  1.0000
    ##  EMM_F34 - Control  0.18333 0.233 222   0.786  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.35333 0.233 222   1.514  0.9799
    ##  Control - ZAN_F3   0.22000 0.233 222   0.943  0.9999
    ##  Control - EMM_F89  0.21333 0.233 222   0.914  0.9999
    ##  Control - EMM_F7   0.21333 0.233 222   0.914  0.9999
    ##  Control - SP_F14   0.14667 0.233 222   0.628  1.0000
    ##  Control - EMM_F34  0.14333 0.233 222   0.614  1.0000
    ##  Control - ZAN_F4   0.09667 0.233 222   0.414  1.0000
    ##  Control - EMM_F70  0.09000 0.233 222   0.386  1.0000
    ##  Control - EMM_F49  0.08667 0.233 222   0.371  1.0000
    ##  Control - EMM_F63  0.02667 0.233 222   0.114  1.0000
    ##  Control - EMM_F64  0.02000 0.233 222   0.086  1.0000
    ##  Control - EMM_F48  0.02000 0.233 222   0.086  1.0000
    ##  Control - EMM_F3   0.01333 0.233 222   0.057  1.0000
    ##  EMM_F66 - Control  0.29000 0.233 222   1.243  0.9973
    ##  EMM_F65 - Control  0.29000 0.233 222   1.243  0.9973
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates

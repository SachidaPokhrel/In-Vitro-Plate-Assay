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

## ***Hymenobacter polaris* EMM_B17**

Plot is generated using loop around the 4 different classes of yeast,
coming up with 4 plots as an output which will be combined in one plot.

``` r
#read data
B17 <- read.csv("CoCultureAssay/CoCultureAssayData/2024-07-17_PeaceAssay_B17.csv", na.strings = "na") 


#load cbb color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
str(B17)
```

    ## 'data.frame':    1536 obs. of  9 variables:
    ##  $ Bacteria               : chr  "B17" "B17" "B17" "B17" ...
    ##  $ Yeast                  : chr  "EMM_F3" "EMM_F3" "EMM_F3" "EMM_F3" ...
    ##  $ Class                  : chr  "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" ...
    ##  $ Replication            : int  1 1 1 1 1 1 1 1 2 2 ...
    ##  $ DAI                    : int  2 2 2 2 2 2 2 2 2 2 ...
    ##  $ distance_to_yeast      : num  0 11.4 17.4 24.9 31.8 ...
    ##  $ colony_diameter        : num  6.88 6.88 7.32 6.87 7.16 7.09 7.09 6.92 7.98 7.03 ...
    ##  $ colony_diameter_control: num  7.79 7.79 7.79 7.79 7.52 6.98 7.67 7.67 6.19 7.31 ...
    ##  $ increase               : num  0 0 0 0 0 0 0 0 0 0 ...

``` r
#round off of the distance to yeast so that it is visibly a categorical variable
B17$distance_to_yeast <- round(B17$distance_to_yeast, 0)
#change into categorical dataset
B17$Replication=as.factor(B17$Replication)
B17$Yeast=as.factor(B17$Yeast)
B17$DAI=as.factor(B17$DAI)
B17$distance_to_yeast=as.factor(B17$distance_to_yeast)

#remove the data for the distance zero
B17.no.contact <- B17[B17$distance_to_yeast != 0, ]
head(B17.no.contact)
```

    ##   Bacteria  Yeast           Class Replication DAI distance_to_yeast
    ## 2      B17 EMM_F3 Dothideomycetes           1   2                11
    ## 3      B17 EMM_F3 Dothideomycetes           1   2                17
    ## 4      B17 EMM_F3 Dothideomycetes           1   2                25
    ## 5      B17 EMM_F3 Dothideomycetes           1   2                32
    ## 6      B17 EMM_F3 Dothideomycetes           1   2                41
    ## 7      B17 EMM_F3 Dothideomycetes           1   2                48
    ##   colony_diameter colony_diameter_control increase
    ## 2            6.88                    7.79        0
    ## 3            7.32                    7.79        0
    ## 4            6.87                    7.79        0
    ## 5            7.16                    7.52        0
    ## 6            7.09                    6.98        0
    ## 7            7.09                    7.67        0

``` r
# remove data first data since it is a reference for increase in colony size and donot follow the assumption of linear model i.e. normal distributuion
B17.Day6 <- B17.no.contact[B17.no.contact$DAI == "6",]


plot_listB17 <- list()

# Get unique classes
unique_class <- unique(B17.Day6$Class[B17.Day6$Class != "Control"])

# Loop through each class
for (i in seq_along(unique_class)) {
  
  # Filter for Control and the current class
  data <- B17.Day6 %>%
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
  plot_listB17[[i]] <- p
}

combined_plotB17 <- ggarrange(plotlist = plot_listB17, ncol = 2, nrow = 2, common.legend = FALSE)
```

    ## Warning: Removed 7 rows containing non-finite outside the scale range
    ## (`stat_summary()`).
    ## Removed 7 rows containing non-finite outside the scale range
    ## (`stat_summary()`).

``` r
# Annotation for the title
final_plotB17 <- annotate_figure(combined_plotB17,
                                top = text_grob(
    expression("Impact on growth of"~italic("Hymenobacter polaris")~"EMM_B17 by Yeast"), face = "bold", color = "Blue2", size = 14, hjust = 0.5))

print(final_plotB17)
```

![](EMM_B17_files/figure-gfm/Plot%20for%20EMM_B17-1.png)<!-- -->

### Stats for *Hymenobacter polaris* EMM_B17

We are using linear mixed model. Our dependent variable or y is increase
(increase in colony diameter from 1st data) and independent variables
are different Yeast isolates, days after inoculation (DAI), and distance
to yeast which is the distance between the yeast and bacterial colony in
the plate. Each plate is replicated 3 times.

``` r
#filter data to remove 1st day data since the first data is taken as a base to measure the increase in colony size to rule out the variability that is caused by the drop inoculation. So, initially the increase in the colony diameter for 1st data for all colony is "0" that violates the assumption of normality, thus we remove that from analysis. This would be similar for all the bacterial isolates.
# Unique days
B17.no.1st.data <- B17.no.contact[B17.no.contact$DAI != "2",]
B17try <- lme(increase~DAI*distance_to_yeast*Yeast, data = B17.no.1st.data, random = ~1|Replication, na.action = na.omit)
anova(B17try)
```

    ##                             numDF denDF   F-value p-value
    ## (Intercept)                     1   649 1042.6016  <.0001
    ## DAI                             2   649  187.9019  <.0001
    ## distance_to_yeast               6   649    9.5336  <.0001
    ## Yeast                          15   649    9.0201  <.0001
    ## DAI:distance_to_yeast          12   649    0.8551  0.5932
    ## DAI:Yeast                      30   649    1.0345  0.4172
    ## distance_to_yeast:Yeast        90   649    1.4769  0.0045
    ## DAI:distance_to_yeast:Yeast   180   649    0.1830  1.0000

``` r
resultsB17=lme(increase~DAI+distance_to_yeast*Yeast, data = B17.no.1st.data, random = ~1|Replication, na.action = na.omit)
summary(resultsB17)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: B17.no.1st.data 
    ##        AIC      BIC    logLik
    ##   1651.797 2205.341 -709.8983
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept)  Residual
    ## StdDev:   0.0769858 0.4706221
    ## 
    ## Fixed effects:  increase ~ DAI + distance_to_yeast * Yeast 
    ##                                       Value Std.Error  DF   t-value p-value
    ## (Intercept)                       1.5623978 0.1644198 871  9.502492  0.0000
    ## DAI6                              0.4599392 0.0366935 871 12.534617  0.0000
    ## DAI8                              0.7762006 0.0366935 871 21.153617  0.0000
    ## distance_to_yeast17              -0.7800000 0.2218534 871 -3.515836  0.0005
    ## distance_to_yeast25              -0.7111111 0.2218534 871 -3.205320  0.0014
    ## distance_to_yeast32              -0.4877778 0.2218534 871 -2.198649  0.0282
    ## distance_to_yeast41              -0.3866667 0.2218534 871 -1.742893  0.0817
    ## distance_to_yeast48              -0.3322222 0.2218534 871 -1.497486  0.1346
    ## distance_to_yeast55               0.1255556 0.2218534 871  0.565939  0.5716
    ## YeastEMM_F3                      -0.6755556 0.2218534 871 -3.045054  0.0024
    ## YeastEMM_F34                     -0.8166667 0.2218534 871 -3.681110  0.0002
    ## YeastEMM_F47                     -0.5877778 0.2218534 871 -2.649397  0.0082
    ## YeastEMM_F48                     -0.4155556 0.2218534 871 -1.873109  0.0614
    ## YeastEMM_F49                     -1.4166667 0.2218534 871 -6.385599  0.0000
    ## YeastEMM_F5                      -0.8633333 0.2218534 871 -3.891459  0.0001
    ## YeastEMM_F63                     -0.4833333 0.2218534 871 -2.178616  0.0296
    ## YeastEMM_F64                     -0.5355556 0.2218534 871 -2.414007  0.0160
    ## YeastEMM_F65                     -0.2466667 0.2218534 871 -1.111845  0.2665
    ## YeastEMM_F7                      -0.2771202 0.2482507 871 -1.116292  0.2646
    ## YeastEMM_F70                     -0.4033333 0.2218534 871 -1.818018  0.0694
    ## YeastEMM_F89                     -0.7133333 0.2218534 871 -3.215337  0.0014
    ## YeastSP_F14                      -0.8888889 0.2218534 871 -4.006650  0.0001
    ## YeastZAN_F3                      -0.2177778 0.2218534 871 -0.981629  0.3266
    ## YeastZAN_F4                      -0.4277778 0.2218534 871 -1.928200  0.0542
    ## distance_to_yeast17:YeastEMM_F3   0.7644444 0.3137481 871  2.436491  0.0150
    ## distance_to_yeast25:YeastEMM_F3   1.0433333 0.3137481 871  3.325386  0.0009
    ## distance_to_yeast32:YeastEMM_F3   0.8088889 0.3137481 871  2.578148  0.0101
    ## distance_to_yeast41:YeastEMM_F3   0.7122222 0.3137481 871  2.270045  0.0234
    ## distance_to_yeast48:YeastEMM_F3   0.7444444 0.3137481 871  2.372746  0.0179
    ## distance_to_yeast55:YeastEMM_F3   0.0588889 0.3137481 871  0.187695  0.8512
    ## distance_to_yeast17:YeastEMM_F34  1.0655556 0.3137481 871  3.396214  0.0007
    ## distance_to_yeast25:YeastEMM_F34  0.7866667 0.3137481 871  2.507320  0.0123
    ## distance_to_yeast32:YeastEMM_F34  0.9400000 0.3137481 871  2.996035  0.0028
    ## distance_to_yeast41:YeastEMM_F34  1.0266667 0.3137481 871  3.272265  0.0011
    ## distance_to_yeast48:YeastEMM_F34  1.1255556 0.3137481 871  3.587450  0.0004
    ## distance_to_yeast55:YeastEMM_F34  0.4488889 0.3137481 871  1.430730  0.1529
    ## distance_to_yeast17:YeastEMM_F47  0.7600000 0.3137481 871  2.422326  0.0156
    ## distance_to_yeast25:YeastEMM_F47  0.8711111 0.3137481 871  2.776467  0.0056
    ## distance_to_yeast32:YeastEMM_F47  0.3711111 0.3137481 871  1.182832  0.2372
    ## distance_to_yeast41:YeastEMM_F47  0.3622222 0.3137481 871  1.154500  0.2486
    ## distance_to_yeast48:YeastEMM_F47  0.3122222 0.3137481 871  0.995137  0.3199
    ## distance_to_yeast55:YeastEMM_F47  0.0200000 0.3137481 871  0.063745  0.9492
    ## distance_to_yeast17:YeastEMM_F48  0.6588889 0.3137481 871  2.100057  0.0360
    ## distance_to_yeast25:YeastEMM_F48  0.8788889 0.3137481 871  2.801257  0.0052
    ## distance_to_yeast32:YeastEMM_F48  0.6488889 0.3137481 871  2.068185  0.0389
    ## distance_to_yeast41:YeastEMM_F48  0.4444444 0.3137481 871  1.416565  0.1570
    ## distance_to_yeast48:YeastEMM_F48  0.3588889 0.3137481 871  1.143876  0.2530
    ## distance_to_yeast55:YeastEMM_F48  0.2944444 0.3137481 871  0.938474  0.3483
    ## distance_to_yeast17:YeastEMM_F49  1.2944444 0.3137481 871  4.125745  0.0000
    ## distance_to_yeast25:YeastEMM_F49  1.4977778 0.3137481 871  4.773823  0.0000
    ## distance_to_yeast32:YeastEMM_F49  1.3177778 0.3137481 871  4.200115  0.0000
    ## distance_to_yeast41:YeastEMM_F49  1.2711111 0.3137481 871  4.051375  0.0001
    ## distance_to_yeast48:YeastEMM_F49  1.2488889 0.3137481 871  3.980547  0.0001
    ## distance_to_yeast55:YeastEMM_F49  0.6388889 0.3137481 871  2.036312  0.0420
    ## distance_to_yeast17:YeastEMM_F5   1.1100000 0.3137481 871  3.537871  0.0004
    ## distance_to_yeast25:YeastEMM_F5   1.2244444 0.3137481 871  3.902636  0.0001
    ## distance_to_yeast32:YeastEMM_F5   0.7077778 0.3137481 871  2.255879  0.0243
    ## distance_to_yeast41:YeastEMM_F5   0.7333333 0.3137481 871  2.337332  0.0196
    ## distance_to_yeast48:YeastEMM_F5   0.8677778 0.3137481 871  2.765843  0.0058
    ## distance_to_yeast55:YeastEMM_F5   0.6422222 0.3137481 871  2.046936  0.0410
    ## distance_to_yeast17:YeastEMM_F63  0.6800000 0.3137481 871  2.167344  0.0305
    ## distance_to_yeast25:YeastEMM_F63  0.7944444 0.3137481 871  2.532110  0.0115
    ## distance_to_yeast32:YeastEMM_F63  0.6977778 0.3137481 871  2.224007  0.0264
    ## distance_to_yeast41:YeastEMM_F63  0.4055556 0.3137481 871  1.292615  0.1965
    ## distance_to_yeast48:YeastEMM_F63  0.6333333 0.3137481 871  2.018605  0.0438
    ## distance_to_yeast55:YeastEMM_F63  0.3422222 0.3137481 871  1.090755  0.2757
    ## distance_to_yeast17:YeastEMM_F64  0.6988889 0.3137481 871  2.227548  0.0262
    ## distance_to_yeast25:YeastEMM_F64  0.6888889 0.3137481 871  2.195675  0.0284
    ## distance_to_yeast32:YeastEMM_F64  0.6666667 0.3137481 871  2.124847  0.0339
    ## distance_to_yeast41:YeastEMM_F64  0.2911111 0.3137481 871  0.927850  0.3537
    ## distance_to_yeast48:YeastEMM_F64  0.4233333 0.3137481 871  1.349278  0.1776
    ## distance_to_yeast55:YeastEMM_F64  0.1777778 0.3137481 871  0.566626  0.5711
    ## distance_to_yeast17:YeastEMM_F65  0.6022222 0.3137481 871  1.919445  0.0553
    ## distance_to_yeast25:YeastEMM_F65  0.4433333 0.3137481 871  1.413023  0.1580
    ## distance_to_yeast32:YeastEMM_F65  0.3022222 0.3137481 871  0.963264  0.3357
    ## distance_to_yeast41:YeastEMM_F65  0.4688889 0.3137481 871  1.494476  0.1354
    ## distance_to_yeast48:YeastEMM_F65  0.5844444 0.3137481 871  1.862783  0.0628
    ## distance_to_yeast55:YeastEMM_F65  0.3511111 0.3137481 871  1.119086  0.2634
    ## distance_to_yeast17:YeastEMM_F7   0.4600000 0.3507810 871  1.311360  0.1901
    ## distance_to_yeast25:YeastEMM_F7   0.7844444 0.3507810 871  2.236280  0.0256
    ## distance_to_yeast32:YeastEMM_F7   0.4477778 0.3507810 871  1.276517  0.2021
    ## distance_to_yeast41:YeastEMM_F7   0.0016667 0.3507810 871  0.004751  0.9962
    ## distance_to_yeast48:YeastEMM_F7  -0.0927778 0.3507810 871 -0.264489  0.7915
    ## distance_to_yeast55:YeastEMM_F7  -1.0472222 0.3507810 871 -2.985402  0.0029
    ## distance_to_yeast17:YeastEMM_F70  0.8777778 0.3137481 871  2.797715  0.0053
    ## distance_to_yeast25:YeastEMM_F70  0.8666667 0.3137481 871  2.762301  0.0059
    ## distance_to_yeast32:YeastEMM_F70  0.1577778 0.3137481 871  0.502880  0.6152
    ## distance_to_yeast41:YeastEMM_F70  0.0977778 0.3137481 871  0.311644  0.7554
    ## distance_to_yeast48:YeastEMM_F70  0.1133333 0.3137481 871  0.361224  0.7180
    ## distance_to_yeast55:YeastEMM_F70  0.1033333 0.3137481 871  0.329351  0.7420
    ## distance_to_yeast17:YeastEMM_F89  0.7611111 0.3137481 871  2.425867  0.0155
    ## distance_to_yeast25:YeastEMM_F89  1.0288889 0.3137481 871  3.279347  0.0011
    ## distance_to_yeast32:YeastEMM_F89  0.7100000 0.3137481 871  2.262962  0.0239
    ## distance_to_yeast41:YeastEMM_F89  0.5100000 0.3137481 871  1.625508  0.1044
    ## distance_to_yeast48:YeastEMM_F89  0.8122222 0.3137481 871  2.588772  0.0098
    ## distance_to_yeast55:YeastEMM_F89  0.6722222 0.3137481 871  2.142554  0.0324
    ## distance_to_yeast17:YeastSP_F14   0.6833333 0.3137481 871  2.177968  0.0297
    ## distance_to_yeast25:YeastSP_F14   0.4377778 0.3137481 871  1.395316  0.1633
    ## distance_to_yeast32:YeastSP_F14   0.1922222 0.3137481 871  0.612664  0.5403
    ## distance_to_yeast41:YeastSP_F14   0.1788889 0.3137481 871  0.570167  0.5687
    ## distance_to_yeast48:YeastSP_F14   0.2144444 0.3137481 871  0.683493  0.4945
    ## distance_to_yeast55:YeastSP_F14   0.0588889 0.3137481 871  0.187695  0.8512
    ## distance_to_yeast17:YeastZAN_F3   0.2866667 0.3137481 871  0.913684  0.3611
    ## distance_to_yeast25:YeastZAN_F3   0.3100000 0.3137481 871  0.988054  0.3234
    ## distance_to_yeast32:YeastZAN_F3   0.4488889 0.3137481 871  1.430730  0.1529
    ## distance_to_yeast41:YeastZAN_F3   0.4722222 0.3137481 871  1.505100  0.1327
    ## distance_to_yeast48:YeastZAN_F3   0.6055556 0.3137481 871  1.930070  0.0539
    ## distance_to_yeast55:YeastZAN_F3   0.3488889 0.3137481 871  1.112003  0.2664
    ## distance_to_yeast17:YeastZAN_F4   0.7588889 0.3137481 871  2.418784  0.0158
    ## distance_to_yeast25:YeastZAN_F4   0.9522222 0.3137481 871  3.034990  0.0025
    ## distance_to_yeast32:YeastZAN_F4   0.6088889 0.3137481 871  1.940694  0.0526
    ## distance_to_yeast41:YeastZAN_F4   0.5911111 0.3137481 871  1.884031  0.0599
    ## distance_to_yeast48:YeastZAN_F4   0.3688889 0.3137481 871  1.175749  0.2400
    ## distance_to_yeast55:YeastZAN_F4   0.0511111 0.3137481 871  0.162905  0.8706
    ##  Correlation: 
    ##                                  (Intr) DAI6   DAI8   ds__17 ds__25 ds__32
    ## DAI6                             -0.112                                   
    ## DAI8                             -0.112  0.500                            
    ## distance_to_yeast17              -0.675  0.000  0.000                     
    ## distance_to_yeast25              -0.675  0.000  0.000  0.500              
    ## distance_to_yeast32              -0.675  0.000  0.000  0.500  0.500       
    ## distance_to_yeast41              -0.675  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast48              -0.675  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast55              -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F3                      -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F34                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F47                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F48                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F49                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F5                      -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F63                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F64                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F65                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F7                      -0.603  0.000  0.000  0.447  0.447  0.447
    ## YeastEMM_F70                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F89                     -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastSP_F14                      -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastZAN_F3                      -0.675  0.000  0.000  0.500  0.500  0.500
    ## YeastZAN_F4                      -0.675  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast17:YeastEMM_F3   0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F3   0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F3   0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F3   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F3   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F3   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F34  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F34  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F34  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F34  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F34  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F34  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F47  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F47  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F47  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F47  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F47  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F47  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F48  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F48  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F48  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F48  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F48  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F48  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F49  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F49  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F49  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F49  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F49  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F49  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F5   0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F5   0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F5   0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F5   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F5   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F5   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F63  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F63  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F63  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F63  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F63  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F63  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F64  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F64  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F64  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F64  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F64  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F64  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F65  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F65  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F65  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F65  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F65  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F65  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F7   0.427  0.000  0.000 -0.632 -0.316 -0.316
    ## distance_to_yeast25:YeastEMM_F7   0.427  0.000  0.000 -0.316 -0.632 -0.316
    ## distance_to_yeast32:YeastEMM_F7   0.427  0.000  0.000 -0.316 -0.316 -0.632
    ## distance_to_yeast41:YeastEMM_F7   0.427  0.000  0.000 -0.316 -0.316 -0.316
    ## distance_to_yeast48:YeastEMM_F7   0.427  0.000  0.000 -0.316 -0.316 -0.316
    ## distance_to_yeast55:YeastEMM_F7   0.427  0.000  0.000 -0.316 -0.316 -0.316
    ## distance_to_yeast17:YeastEMM_F70  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F70  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F70  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F70  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F70  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F70  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F89  0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F89  0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F89  0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F89  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F89  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F89  0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastSP_F14   0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastSP_F14   0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastSP_F14   0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastSP_F14   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastSP_F14   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastSP_F14   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastZAN_F3   0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastZAN_F3   0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastZAN_F3   0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastZAN_F3   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastZAN_F3   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastZAN_F3   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastZAN_F4   0.477  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastZAN_F4   0.477  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastZAN_F4   0.477  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastZAN_F4   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastZAN_F4   0.477  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastZAN_F4   0.477  0.000  0.000 -0.354 -0.354 -0.354
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
    ## YeastEMM_F5                       0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F63                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F64                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F65                      0.500  0.500  0.500  0.500    0.500  
    ## YeastEMM_F7                       0.447  0.447  0.447  0.447    0.447  
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
    ## distance_to_yeast17:YeastEMM_F5  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F5  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F5  -0.354 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F5  -0.707 -0.354 -0.354 -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F5  -0.354 -0.707 -0.354 -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F5  -0.354 -0.354 -0.707 -0.354   -0.354  
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
    ## distance_to_yeast17:YeastEMM_F7  -0.316 -0.316 -0.316 -0.316   -0.316  
    ## distance_to_yeast25:YeastEMM_F7  -0.316 -0.316 -0.316 -0.316   -0.316  
    ## distance_to_yeast32:YeastEMM_F7  -0.316 -0.316 -0.316 -0.316   -0.316  
    ## distance_to_yeast41:YeastEMM_F7  -0.632 -0.316 -0.316 -0.316   -0.316  
    ## distance_to_yeast48:YeastEMM_F7  -0.316 -0.632 -0.316 -0.316   -0.316  
    ## distance_to_yeast55:YeastEMM_F7  -0.316 -0.316 -0.632 -0.316   -0.316  
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
    ## YeastEMM_F48                      0.500                                     
    ## YeastEMM_F49                      0.500    0.500                            
    ## YeastEMM_F5                       0.500    0.500    0.500                   
    ## YeastEMM_F63                      0.500    0.500    0.500    0.500          
    ## YeastEMM_F64                      0.500    0.500    0.500    0.500   0.500  
    ## YeastEMM_F65                      0.500    0.500    0.500    0.500   0.500  
    ## YeastEMM_F7                       0.447    0.447    0.447    0.447   0.447  
    ## YeastEMM_F70                      0.500    0.500    0.500    0.500   0.500  
    ## YeastEMM_F89                      0.500    0.500    0.500    0.500   0.500  
    ## YeastSP_F14                       0.500    0.500    0.500    0.500   0.500  
    ## YeastZAN_F3                       0.500    0.500    0.500    0.500   0.500  
    ## YeastZAN_F4                       0.500    0.500    0.500    0.500   0.500  
    ## distance_to_yeast17:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F34 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F47 -0.707   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F48 -0.354   -0.707   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F49 -0.354   -0.354   -0.707   -0.354  -0.354  
    ## distance_to_yeast17:YeastEMM_F5  -0.354   -0.354   -0.354   -0.707  -0.354  
    ## distance_to_yeast25:YeastEMM_F5  -0.354   -0.354   -0.354   -0.707  -0.354  
    ## distance_to_yeast32:YeastEMM_F5  -0.354   -0.354   -0.354   -0.707  -0.354  
    ## distance_to_yeast41:YeastEMM_F5  -0.354   -0.354   -0.354   -0.707  -0.354  
    ## distance_to_yeast48:YeastEMM_F5  -0.354   -0.354   -0.354   -0.707  -0.354  
    ## distance_to_yeast55:YeastEMM_F5  -0.354   -0.354   -0.354   -0.707  -0.354  
    ## distance_to_yeast17:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354  -0.707  
    ## distance_to_yeast25:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354  -0.707  
    ## distance_to_yeast32:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354  -0.707  
    ## distance_to_yeast41:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354  -0.707  
    ## distance_to_yeast48:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354  -0.707  
    ## distance_to_yeast55:YeastEMM_F63 -0.354   -0.354   -0.354   -0.354  -0.707  
    ## distance_to_yeast17:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F64 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F65 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastEMM_F7  -0.316   -0.316   -0.316   -0.316  -0.316  
    ## distance_to_yeast25:YeastEMM_F7  -0.316   -0.316   -0.316   -0.316  -0.316  
    ## distance_to_yeast32:YeastEMM_F7  -0.316   -0.316   -0.316   -0.316  -0.316  
    ## distance_to_yeast41:YeastEMM_F7  -0.316   -0.316   -0.316   -0.316  -0.316  
    ## distance_to_yeast48:YeastEMM_F7  -0.316   -0.316   -0.316   -0.316  -0.316  
    ## distance_to_yeast55:YeastEMM_F7  -0.316   -0.316   -0.316   -0.316  -0.316  
    ## distance_to_yeast17:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F70 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastEMM_F89 -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastSP_F14  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastSP_F14  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastSP_F14  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastSP_F14  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastSP_F14  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastSP_F14  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastZAN_F3  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast17:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast25:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast32:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast41:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast48:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354  -0.354  
    ## distance_to_yeast55:YeastZAN_F4  -0.354   -0.354   -0.354   -0.354  -0.354  
    ##                                  YEMM_F64 YEMM_F65 YsEMM_F7 YEMM_F70 YEMM_F8
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
    ## YeastEMM_F65                      0.500                                     
    ## YeastEMM_F7                       0.447    0.447                            
    ## YeastEMM_F70                      0.500    0.500    0.447                   
    ## YeastEMM_F89                      0.500    0.500    0.447    0.500          
    ## YeastSP_F14                       0.500    0.500    0.447    0.500    0.500 
    ## YeastZAN_F3                       0.500    0.500    0.447    0.500    0.500 
    ## YeastZAN_F4                       0.500    0.500    0.447    0.500    0.500 
    ## distance_to_yeast17:YeastEMM_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F34 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F34 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F34 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F34 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F34 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F34 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F47 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F47 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F47 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F47 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F47 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F47 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F48 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F48 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F48 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F48 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F48 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F48 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F49 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F49 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F49 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F49 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F49 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F49 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F5  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F5  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F5  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F5  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F5  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F5  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F63 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F63 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F63 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F63 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F63 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F63 -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F64 -0.707   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F64 -0.707   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F64 -0.707   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F64 -0.707   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F64 -0.707   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F64 -0.707   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F65 -0.354   -0.707   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastEMM_F65 -0.354   -0.707   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastEMM_F65 -0.354   -0.707   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastEMM_F65 -0.354   -0.707   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastEMM_F65 -0.354   -0.707   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastEMM_F65 -0.354   -0.707   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastEMM_F7  -0.316   -0.316   -0.707   -0.316   -0.316 
    ## distance_to_yeast25:YeastEMM_F7  -0.316   -0.316   -0.707   -0.316   -0.316 
    ## distance_to_yeast32:YeastEMM_F7  -0.316   -0.316   -0.707   -0.316   -0.316 
    ## distance_to_yeast41:YeastEMM_F7  -0.316   -0.316   -0.707   -0.316   -0.316 
    ## distance_to_yeast48:YeastEMM_F7  -0.316   -0.316   -0.707   -0.316   -0.316 
    ## distance_to_yeast55:YeastEMM_F7  -0.316   -0.316   -0.707   -0.316   -0.316 
    ## distance_to_yeast17:YeastEMM_F70 -0.354   -0.354   -0.316   -0.707   -0.354 
    ## distance_to_yeast25:YeastEMM_F70 -0.354   -0.354   -0.316   -0.707   -0.354 
    ## distance_to_yeast32:YeastEMM_F70 -0.354   -0.354   -0.316   -0.707   -0.354 
    ## distance_to_yeast41:YeastEMM_F70 -0.354   -0.354   -0.316   -0.707   -0.354 
    ## distance_to_yeast48:YeastEMM_F70 -0.354   -0.354   -0.316   -0.707   -0.354 
    ## distance_to_yeast55:YeastEMM_F70 -0.354   -0.354   -0.316   -0.707   -0.354 
    ## distance_to_yeast17:YeastEMM_F89 -0.354   -0.354   -0.316   -0.354   -0.707 
    ## distance_to_yeast25:YeastEMM_F89 -0.354   -0.354   -0.316   -0.354   -0.707 
    ## distance_to_yeast32:YeastEMM_F89 -0.354   -0.354   -0.316   -0.354   -0.707 
    ## distance_to_yeast41:YeastEMM_F89 -0.354   -0.354   -0.316   -0.354   -0.707 
    ## distance_to_yeast48:YeastEMM_F89 -0.354   -0.354   -0.316   -0.354   -0.707 
    ## distance_to_yeast55:YeastEMM_F89 -0.354   -0.354   -0.316   -0.354   -0.707 
    ## distance_to_yeast17:YeastSP_F14  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastSP_F14  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastSP_F14  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastSP_F14  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastSP_F14  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastSP_F14  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastZAN_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastZAN_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastZAN_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastZAN_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastZAN_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastZAN_F3  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast17:YeastZAN_F4  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast25:YeastZAN_F4  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast32:YeastZAN_F4  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast41:YeastZAN_F4  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast48:YeastZAN_F4  -0.354   -0.354   -0.316   -0.354   -0.354 
    ## distance_to_yeast55:YeastZAN_F4  -0.354   -0.354   -0.316   -0.354   -0.354 
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
    ## YeastEMM_F5                                                           
    ## YeastEMM_F63                                                          
    ## YeastEMM_F64                                                          
    ## YeastEMM_F65                                                          
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
    ## distance_to_yeast17:YeastEMM_F5  -0.354 -0.354  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F5  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F5  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F5  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F5  -0.354 -0.354  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F5  -0.354 -0.354  -0.354   0.250        
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
    ## distance_to_yeast17:YeastEMM_F7  -0.316 -0.316  -0.316   0.447        
    ## distance_to_yeast25:YeastEMM_F7  -0.316 -0.316  -0.316   0.224        
    ## distance_to_yeast32:YeastEMM_F7  -0.316 -0.316  -0.316   0.224        
    ## distance_to_yeast41:YeastEMM_F7  -0.316 -0.316  -0.316   0.224        
    ## distance_to_yeast48:YeastEMM_F7  -0.316 -0.316  -0.316   0.224        
    ## distance_to_yeast55:YeastEMM_F7  -0.316 -0.316  -0.316   0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.250          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.224          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.500          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast25:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast48:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.447          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.250          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.224          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.500          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast25:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast48:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.447          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.250          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.224          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.500          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast25:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast48:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.447          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.250          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.224          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.500          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast25:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast48:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.447          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F5   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F5   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.250          0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.224          0.224        
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
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250                      
    ## distance_to_yeast25:YeastEMM_F5   0.250          0.250          0.500       
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.250          0.500       
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.500       
    ## distance_to_yeast48:YeastEMM_F5   0.500          0.250          0.500       
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.500          0.500       
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.447       
    ## distance_to_yeast25:YeastEMM_F7   0.224          0.224          0.224       
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.224          0.224       
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.224       
    ## distance_to_yeast48:YeastEMM_F7   0.447          0.224          0.224       
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.447          0.224       
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
    ## distance_to_yeast17:YeastEMM_F5                                           
    ## distance_to_yeast25:YeastEMM_F5                                           
    ## distance_to_yeast32:YeastEMM_F5   0.500                                   
    ## distance_to_yeast41:YeastEMM_F5   0.500         0.500                     
    ## distance_to_yeast48:YeastEMM_F5   0.500         0.500         0.500       
    ## distance_to_yeast55:YeastEMM_F5   0.500         0.500         0.500       
    ## distance_to_yeast17:YeastEMM_F63  0.250         0.250         0.250       
    ## distance_to_yeast25:YeastEMM_F63  0.500         0.250         0.250       
    ## distance_to_yeast32:YeastEMM_F63  0.250         0.500         0.250       
    ## distance_to_yeast41:YeastEMM_F63  0.250         0.250         0.500       
    ## distance_to_yeast48:YeastEMM_F63  0.250         0.250         0.250       
    ## distance_to_yeast55:YeastEMM_F63  0.250         0.250         0.250       
    ## distance_to_yeast17:YeastEMM_F64  0.250         0.250         0.250       
    ## distance_to_yeast25:YeastEMM_F64  0.500         0.250         0.250       
    ## distance_to_yeast32:YeastEMM_F64  0.250         0.500         0.250       
    ## distance_to_yeast41:YeastEMM_F64  0.250         0.250         0.500       
    ## distance_to_yeast48:YeastEMM_F64  0.250         0.250         0.250       
    ## distance_to_yeast55:YeastEMM_F64  0.250         0.250         0.250       
    ## distance_to_yeast17:YeastEMM_F65  0.250         0.250         0.250       
    ## distance_to_yeast25:YeastEMM_F65  0.500         0.250         0.250       
    ## distance_to_yeast32:YeastEMM_F65  0.250         0.500         0.250       
    ## distance_to_yeast41:YeastEMM_F65  0.250         0.250         0.500       
    ## distance_to_yeast48:YeastEMM_F65  0.250         0.250         0.250       
    ## distance_to_yeast55:YeastEMM_F65  0.250         0.250         0.250       
    ## distance_to_yeast17:YeastEMM_F7   0.224         0.224         0.224       
    ## distance_to_yeast25:YeastEMM_F7   0.447         0.224         0.224       
    ## distance_to_yeast32:YeastEMM_F7   0.224         0.447         0.224       
    ## distance_to_yeast41:YeastEMM_F7   0.224         0.224         0.447       
    ## distance_to_yeast48:YeastEMM_F7   0.224         0.224         0.224       
    ## distance_to_yeast55:YeastEMM_F7   0.224         0.224         0.224       
    ## distance_to_yeast17:YeastEMM_F70  0.250         0.250         0.250       
    ## distance_to_yeast25:YeastEMM_F70  0.500         0.250         0.250       
    ## distance_to_yeast32:YeastEMM_F70  0.250         0.500         0.250       
    ## distance_to_yeast41:YeastEMM_F70  0.250         0.250         0.500       
    ## distance_to_yeast48:YeastEMM_F70  0.250         0.250         0.250       
    ## distance_to_yeast55:YeastEMM_F70  0.250         0.250         0.250       
    ## distance_to_yeast17:YeastEMM_F89  0.250         0.250         0.250       
    ## distance_to_yeast25:YeastEMM_F89  0.500         0.250         0.250       
    ## distance_to_yeast32:YeastEMM_F89  0.250         0.500         0.250       
    ## distance_to_yeast41:YeastEMM_F89  0.250         0.250         0.500       
    ## distance_to_yeast48:YeastEMM_F89  0.250         0.250         0.250       
    ## distance_to_yeast55:YeastEMM_F89  0.250         0.250         0.250       
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
    ## distance_to_yeast17:YeastEMM_F5                                            
    ## distance_to_yeast25:YeastEMM_F5                                            
    ## distance_to_yeast32:YeastEMM_F5                                            
    ## distance_to_yeast41:YeastEMM_F5                                            
    ## distance_to_yeast48:YeastEMM_F5                                            
    ## distance_to_yeast55:YeastEMM_F5   0.500                                    
    ## distance_to_yeast17:YeastEMM_F63  0.250         0.250                      
    ## distance_to_yeast25:YeastEMM_F63  0.250         0.250         0.500        
    ## distance_to_yeast32:YeastEMM_F63  0.250         0.250         0.500        
    ## distance_to_yeast41:YeastEMM_F63  0.250         0.250         0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.500         0.250         0.500        
    ## distance_to_yeast55:YeastEMM_F63  0.250         0.500         0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.250         0.250         0.500        
    ## distance_to_yeast25:YeastEMM_F64  0.250         0.250         0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250         0.250         0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.250         0.250         0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.500         0.250         0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250         0.500         0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.250         0.250         0.500        
    ## distance_to_yeast25:YeastEMM_F65  0.250         0.250         0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250         0.250         0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250         0.250         0.250        
    ## distance_to_yeast48:YeastEMM_F65  0.500         0.250         0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250         0.500         0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.224         0.224         0.447        
    ## distance_to_yeast25:YeastEMM_F7   0.224         0.224         0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224         0.224         0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224         0.224         0.224        
    ## distance_to_yeast48:YeastEMM_F7   0.447         0.224         0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224         0.447         0.224        
    ## distance_to_yeast17:YeastEMM_F70  0.250         0.250         0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250         0.250         0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250         0.250         0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250         0.250         0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500         0.250         0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250         0.500         0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.250         0.250         0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.250         0.250         0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250         0.250         0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.250         0.250         0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.500         0.250         0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250         0.500         0.250        
    ## distance_to_yeast17:YeastSP_F14   0.250         0.250         0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250         0.250         0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250         0.250         0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250         0.250         0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500         0.250         0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250         0.500         0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.250         0.250         0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.250         0.250         0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250         0.250         0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.250         0.250         0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.500         0.250         0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250         0.500         0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.250         0.250         0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.250         0.250         0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250         0.250         0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.250         0.250         0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.500         0.250         0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250         0.500         0.250        
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.224          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast25:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast48:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.447          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F65  0.500          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.500          0.250        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.224          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F65  0.250          0.250                       
    ## distance_to_yeast25:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast32:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F65  0.250          0.250          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.500          0.250          0.500        
    ## distance_to_yeast55:YeastEMM_F65  0.250          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast25:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast48:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.447          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F7   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F7   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.224          0.224        
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
    ##                                  d__48:YEMM_F65 d__55:YEMM_F65 ds__17:YEMM_F7
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
    ## distance_to_yeast17:YeastEMM_F7   0.224          0.224                       
    ## distance_to_yeast25:YeastEMM_F7   0.224          0.224          0.500        
    ## distance_to_yeast32:YeastEMM_F7   0.224          0.224          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.224          0.224          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.447          0.224          0.500        
    ## distance_to_yeast55:YeastEMM_F7   0.224          0.447          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.447        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.224        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.224        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.224        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.224        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.224        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.447        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.224        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.224        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.224        
    ## distance_to_yeast48:YeastEMM_F89  0.500          0.250          0.224        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.500          0.224        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.447        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.224        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.224        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.224        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.224        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.224        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.447        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.224        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.224        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.224        
    ## distance_to_yeast48:YeastZAN_F3   0.500          0.250          0.224        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.500          0.224        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.447        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.224        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.224        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.224        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.224        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7                                              
    ## distance_to_yeast25:YeastEMM_F7                                              
    ## distance_to_yeast32:YeastEMM_F7   0.500                                      
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F7   0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F70  0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F70  0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F70  0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F70  0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F70  0.224          0.224          0.224        
    ## distance_to_yeast17:YeastEMM_F89  0.224          0.224          0.224        
    ## distance_to_yeast25:YeastEMM_F89  0.447          0.224          0.224        
    ## distance_to_yeast32:YeastEMM_F89  0.224          0.447          0.224        
    ## distance_to_yeast41:YeastEMM_F89  0.224          0.224          0.447        
    ## distance_to_yeast48:YeastEMM_F89  0.224          0.224          0.224        
    ## distance_to_yeast55:YeastEMM_F89  0.224          0.224          0.224        
    ## distance_to_yeast17:YeastSP_F14   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastSP_F14   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastSP_F14   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastSP_F14   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastSP_F14   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastSP_F14   0.224          0.224          0.224        
    ## distance_to_yeast17:YeastZAN_F3   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastZAN_F3   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastZAN_F3   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastZAN_F3   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastZAN_F3   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastZAN_F3   0.224          0.224          0.224        
    ## distance_to_yeast17:YeastZAN_F4   0.224          0.224          0.224        
    ## distance_to_yeast25:YeastZAN_F4   0.447          0.224          0.224        
    ## distance_to_yeast32:YeastZAN_F4   0.224          0.447          0.224        
    ## distance_to_yeast41:YeastZAN_F4   0.224          0.224          0.447        
    ## distance_to_yeast48:YeastZAN_F4   0.224          0.224          0.224        
    ## distance_to_yeast55:YeastZAN_F4   0.224          0.224          0.224        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## distance_to_yeast17:YeastEMM_F7                                              
    ## distance_to_yeast25:YeastEMM_F7                                              
    ## distance_to_yeast32:YeastEMM_F7                                              
    ## distance_to_yeast41:YeastEMM_F7                                              
    ## distance_to_yeast48:YeastEMM_F7                                              
    ## distance_to_yeast55:YeastEMM_F7   0.500                                      
    ## distance_to_yeast17:YeastEMM_F70  0.224          0.224                       
    ## distance_to_yeast25:YeastEMM_F70  0.224          0.224          0.500        
    ## distance_to_yeast32:YeastEMM_F70  0.224          0.224          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.224          0.224          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.447          0.224          0.500        
    ## distance_to_yeast55:YeastEMM_F70  0.224          0.447          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.224          0.224          0.500        
    ## distance_to_yeast25:YeastEMM_F89  0.224          0.224          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.224          0.224          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.224          0.224          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.447          0.224          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.224          0.447          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.224          0.224          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.224          0.224          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.224          0.224          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.224          0.224          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.447          0.224          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.224          0.447          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.224          0.224          0.500        
    ## distance_to_yeast25:YeastZAN_F3   0.224          0.224          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.224          0.224          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.224          0.224          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.447          0.224          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.224          0.447          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.224          0.224          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.224          0.224          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.224          0.224          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.224          0.224          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.447          0.224          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.224          0.447          0.250        
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
    ## YeastEMM_F5                                                                  
    ## YeastEMM_F63                                                                 
    ## YeastEMM_F64                                                                 
    ## YeastEMM_F65                                                                 
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
    ## YeastEMM_F5                                                                 
    ## YeastEMM_F63                                                                
    ## YeastEMM_F64                                                                
    ## YeastEMM_F65                                                                
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
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
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
    ## YeastEMM_F5                                                                   
    ## YeastEMM_F63                                                                  
    ## YeastEMM_F64                                                                  
    ## YeastEMM_F65                                                                  
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
    ## YeastEMM_F5                                                         
    ## YeastEMM_F63                                                        
    ## YeastEMM_F64                                                        
    ## YeastEMM_F65                                                        
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
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
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
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
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
    ## YeastEMM_F5                                                               
    ## YeastEMM_F63                                                              
    ## YeastEMM_F64                                                              
    ## YeastEMM_F65                                                              
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
    ## YeastEMM_F5                                                 
    ## YeastEMM_F63                                                
    ## YeastEMM_F64                                                
    ## YeastEMM_F65                                                
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
    ##          Min           Q1          Med           Q3          Max 
    ## -4.887204324 -0.566252239 -0.002177437  0.604449129  3.337324837 
    ## 
    ## Number of Observations: 987
    ## Number of Groups: 3

``` r
#Anova(resultsB5)
anova(resultsB17)
```

    ##                         numDF denDF   F-value p-value
    ## (Intercept)                 1   871 1042.2003  <.0001
    ## DAI                         2   871  226.2931  <.0001
    ## distance_to_yeast           6   871   11.4815  <.0001
    ## Yeast                      15   871   10.8649  <.0001
    ## distance_to_yeast:Yeast    90   871    1.7787  <.0001

### Loop for running analysis for each day separately for EMM_B17

Days after inoculation (DAI) as a factor is always significantly
impacting the growth. Also, our plot will represent data for Day 6 thus
we want relevant stats and comparison on Day 6 to present in the plot.
So, loop was made for each Day data and removing DAI from the model and
keeping rest of it present.

``` r
day <- unique(B17.no.1st.data$DAI)
# Initialize empty lists to store results
B17_Day <- list()
meansB17_Day <- list()
resultsB17_Day <- list()
model_summary_B17 <- list()
anova_table_B17 <- list()
Results_lsmeansB17 <- list()
filtered_comparisonsB17 <- list()

for (i in seq_along(day)) {
    # Filter data for the current day
  B17_Day[[i]] <- filter(B17.no.1st.data, DAI == day[i])
    # Fit the mixed effects model
  resultsB17_Day[[i]] <- lme(increase ~ distance_to_yeast * Yeast,
                            data = B17_Day[[i]],
                            random = ~1 | Replication,
                            na.action = na.omit)
    model_summary_B17[[i]] <- summary(resultsB17_Day[[i]])
  anova_table_B17[[i]] <- anova(resultsB17_Day[[i]])
   # Estimated marginal means
  meansB17_Day[[i]] <- emmeans(resultsB17_Day[[i]], ~ Yeast | distance_to_yeast)
    # Compact letter display with Tukey adjustment
  Results_lsmeansB17[[i]] <- multcomp::cld(meansB17_Day[[i]], alpha = 0.05,
                                          Letters = letters, reversed = TRUE,
                                          details = TRUE)
   # Print lsmeans results
  #print(Results_lsmeansB17[[i]])
  
  # Filter comparisons involving "Control"
  filtered_comparisonsB17[[i]] <- Results_lsmeansB17[[i]]$comparisons %>%
    filter(grepl("Control", contrast))
   # Print filtered contrasts
  print(filtered_comparisonsB17[[i]])
}
```

    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49  1.53000 0.424 215   3.608  0.0322
    ##  Control - EMM_F5   0.69000 0.424 215   1.627  0.9617
    ##  Control - EMM_F64  0.59333 0.424 215   1.399  0.9906
    ##  Control - EMM_F89  0.57000 0.424 215   1.344  0.9937
    ##  Control - EMM_F34  0.39000 0.424 215   0.920  0.9999
    ##  Control - EMM_F3   0.31333 0.424 215   0.739  1.0000
    ##  Control - EMM_F7   0.29025 0.474 215   0.612  1.0000
    ##  Control - ZAN_F4   0.26667 0.424 215   0.629  1.0000
    ##  Control - EMM_F65  0.24333 0.424 215   0.574  1.0000
    ##  Control - EMM_F70  0.24000 0.424 215   0.566  1.0000
    ##  Control - EMM_F63  0.23000 0.424 215   0.542  1.0000
    ##  Control - EMM_F48  0.22333 0.424 215   0.527  1.0000
    ##  Control - SP_F14   0.19667 0.424 215   0.464  1.0000
    ##  Control - EMM_F47  0.14333 0.424 215   0.338  1.0000
    ##  Control - ZAN_F3   0.04667 0.424 215   0.110  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F49 - Control  0.08667 0.424 215   0.204  1.0000
    ##  EMM_F65 - Control  0.12667 0.424 215   0.299  1.0000
    ##  EMM_F89 - Control  0.17000 0.424 215   0.401  1.0000
    ##  ZAN_F3 - Control   0.17667 0.424 215   0.417  1.0000
    ##  ZAN_F4 - Control   0.22000 0.424 215   0.519  1.0000
    ##  EMM_F64 - Control  0.29667 0.424 215   0.700  1.0000
    ##  SP_F14 - Control   0.32000 0.424 215   0.755  1.0000
    ##  EMM_F70 - Control  0.32333 0.424 215   0.762  1.0000
    ##  EMM_F47 - Control  0.37000 0.424 215   0.872  1.0000
    ##  EMM_F48 - Control  0.37667 0.424 215   0.888  1.0000
    ##  EMM_F7 - Control   0.38975 0.474 215   0.822  1.0000
    ##  EMM_F63 - Control  0.41333 0.424 215   0.975  0.9998
    ##  EMM_F5 - Control   0.45667 0.424 215   1.077  0.9995
    ##  EMM_F3 - Control   0.53333 0.424 215   1.258  0.9969
    ##  EMM_F34 - Control  0.57667 0.424 215   1.360  0.9929
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  SP_F14 - Control   0.05000 0.424 215   0.118  1.0000
    ##  EMM_F65 - Control  0.25000 0.424 215   0.589  1.0000
    ##  ZAN_F3 - Control   0.30667 0.424 215   0.723  1.0000
    ##  EMM_F64 - Control  0.32000 0.424 215   0.755  1.0000
    ##  EMM_F34 - Control  0.44333 0.424 215   1.045  0.9996
    ##  EMM_F70 - Control  0.52000 0.424 215   1.226  0.9976
    ##  EMM_F49 - Control  0.57000 0.424 215   1.344  0.9937
    ##  EMM_F63 - Control  0.58333 0.424 215   1.375  0.9920
    ##  ZAN_F4 - Control   0.61000 0.424 215   1.438  0.9876
    ##  EMM_F3 - Control   0.64667 0.424 215   1.525  0.9785
    ##  EMM_F89 - Control  0.67000 0.424 215   1.580  0.9704
    ##  EMM_F47 - Control  0.67667 0.424 215   1.596  0.9677
    ##  EMM_F48 - Control  0.67667 0.424 215   1.596  0.9677
    ##  EMM_F5 - Control   0.76000 0.424 215   1.792  0.9167
    ##  EMM_F7 - Control   1.02141 0.474 215   2.154  0.7277
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.62333 0.424 215   1.470  0.9847
    ##  Control - EMM_F70  0.35000 0.424 215   0.825  1.0000
    ##  Control - EMM_F5   0.20333 0.424 215   0.479  1.0000
    ##  Control - EMM_F65  0.18667 0.424 215   0.440  1.0000
    ##  Control - EMM_F47  0.16333 0.424 215   0.385  1.0000
    ##  Control - EMM_F89  0.10333 0.424 215   0.244  1.0000
    ##  Control - EMM_F49  0.07667 0.424 215   0.181  1.0000
    ##  Control - ZAN_F3   0.04333 0.424 215   0.102  1.0000
    ##  EMM_F64 - Control  0.07000 0.424 215   0.165  1.0000
    ##  ZAN_F4 - Control   0.12667 0.424 215   0.299  1.0000
    ##  EMM_F3 - Control   0.13667 0.424 215   0.322  1.0000
    ##  EMM_F7 - Control   0.18641 0.474 215   0.393  1.0000
    ##  EMM_F63 - Control  0.19333 0.424 215   0.456  1.0000
    ##  EMM_F48 - Control  0.22333 0.424 215   0.527  1.0000
    ##  EMM_F34 - Control  0.23667 0.424 215   0.558  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.31333 0.424 215   0.739  1.0000
    ##  Control - EMM_F70  0.16333 0.424 215   0.385  1.0000
    ##  Control - EMM_F64  0.10000 0.424 215   0.236  1.0000
    ##  Control - EMM_F5   0.07333 0.424 215   0.173  1.0000
    ##  Control - EMM_F47  0.06333 0.424 215   0.149  1.0000
    ##  Control - EMM_F7   0.06025 0.474 215   0.127  1.0000
    ##  EMM_F49 - Control  0.00667 0.424 215   0.016  1.0000
    ##  EMM_F89 - Control  0.04000 0.424 215   0.094  1.0000
    ##  EMM_F48 - Control  0.10000 0.424 215   0.236  1.0000
    ##  EMM_F63 - Control  0.17667 0.424 215   0.417  1.0000
    ##  EMM_F65 - Control  0.22333 0.424 215   0.527  1.0000
    ##  EMM_F3 - Control   0.30333 0.424 215   0.715  1.0000
    ##  ZAN_F4 - Control   0.32667 0.424 215   0.770  1.0000
    ##  ZAN_F3 - Control   0.33000 0.424 215   0.778  1.0000
    ##  EMM_F34 - Control  0.37333 0.424 215   0.880  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.33333 0.424 215   0.786  1.0000
    ##  Control - EMM_F70  0.20333 0.424 215   0.479  1.0000
    ##  Control - EMM_F7   0.11859 0.474 215   0.250  1.0000
    ##  Control - EMM_F48  0.11333 0.424 215   0.267  1.0000
    ##  Control - EMM_F5   0.09000 0.424 215   0.212  1.0000
    ##  Control - ZAN_F4   0.04333 0.424 215   0.102  1.0000
    ##  EMM_F64 - Control  0.05667 0.424 215   0.134  1.0000
    ##  EMM_F49 - Control  0.06333 0.424 215   0.149  1.0000
    ##  EMM_F47 - Control  0.12333 0.424 215   0.291  1.0000
    ##  EMM_F89 - Control  0.14333 0.424 215   0.338  1.0000
    ##  EMM_F34 - Control  0.23333 0.424 215   0.550  1.0000
    ##  EMM_F63 - Control  0.23667 0.424 215   0.558  1.0000
    ##  EMM_F3 - Control   0.27000 0.424 215   0.637  1.0000
    ##  EMM_F65 - Control  0.38333 0.424 215   0.904  0.9999
    ##  ZAN_F3 - Control   0.40000 0.424 215   0.943  0.9999
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F7   1.33525 0.474 215   2.816  0.2683
    ##  Control - EMM_F49  0.71667 0.424 215   1.690  0.9474
    ##  Control - EMM_F3   0.40667 0.424 215   0.959  0.9999
    ##  Control - SP_F14   0.32000 0.424 215   0.755  1.0000
    ##  Control - ZAN_F4   0.27333 0.424 215   0.645  1.0000
    ##  Control - EMM_F34  0.21667 0.424 215   0.511  1.0000
    ##  Control - EMM_F70  0.19667 0.424 215   0.464  1.0000
    ##  Control - EMM_F47  0.16667 0.424 215   0.393  1.0000
    ##  Control - EMM_F48  0.01333 0.424 215   0.031  1.0000
    ##  EMM_F5 - Control   0.02667 0.424 215   0.063  1.0000
    ##  EMM_F65 - Control  0.07000 0.424 215   0.165  1.0000
    ##  EMM_F64 - Control  0.11000 0.424 215   0.259  1.0000
    ##  EMM_F63 - Control  0.19667 0.424 215   0.464  1.0000
    ##  EMM_F89 - Control  0.24000 0.424 215   0.566  1.0000
    ##  ZAN_F3 - Control   0.38000 0.424 215   0.896  0.9999
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49  1.43000 0.414 215   3.453  0.0522
    ##  Control - EMM_F5   1.09333 0.414 215   2.640  0.3770
    ##  Control - EMM_F34  1.03667 0.414 215   2.503  0.4734
    ##  Control - SP_F14   1.02667 0.414 215   2.479  0.4911
    ##  Control - EMM_F89  0.99333 0.414 215   2.398  0.5509
    ##  Control - EMM_F3   0.89333 0.414 215   2.157  0.7257
    ##  Control - EMM_F47  0.83000 0.414 215   2.004  0.8204
    ##  Control - EMM_F63  0.72333 0.414 215   1.747  0.9316
    ##  Control - EMM_F64  0.71333 0.414 215   1.722  0.9387
    ##  Control - EMM_F48  0.55667 0.414 215   1.344  0.9937
    ##  Control - ZAN_F4   0.52000 0.414 215   1.256  0.9969
    ##  Control - EMM_F70  0.51667 0.414 215   1.248  0.9971
    ##  Control - EMM_F7   0.48633 0.463 215   1.050  0.9996
    ##  Control - ZAN_F3   0.32667 0.414 215   0.789  1.0000
    ##  Control - EMM_F65  0.25333 0.414 215   0.612  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.44000 0.414 215   1.062  0.9995
    ##  Control - EMM_F49  0.31667 0.414 215   0.765  1.0000
    ##  Control - EMM_F3   0.23667 0.414 215   0.571  1.0000
    ##  Control - EMM_F63  0.15667 0.414 215   0.378  1.0000
    ##  Control - EMM_F47  0.13000 0.414 215   0.314  1.0000
    ##  Control - EMM_F89  0.07333 0.414 215   0.177  1.0000
    ##  Control - EMM_F64  0.05000 0.414 215   0.121  1.0000
    ##  Control - EMM_F7   0.04300 0.463 215   0.093  1.0000
    ##  Control - ZAN_F3   0.02333 0.414 215   0.056  1.0000
    ##  EMM_F48 - Control  0.11000 0.414 215   0.266  1.0000
    ##  EMM_F34 - Control  0.14000 0.414 215   0.338  1.0000
    ##  EMM_F5 - Control   0.14333 0.414 215   0.346  1.0000
    ##  ZAN_F4 - Control   0.42000 0.414 215   1.014  0.9997
    ##  EMM_F65 - Control  0.45000 0.414 215   1.087  0.9994
    ##  EMM_F70 - Control  0.52333 0.414 215   1.264  0.9967
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.29667 0.414 215   0.716  1.0000
    ##  Control - EMM_F34  0.02667 0.414 215   0.064  1.0000
    ##  EMM_F49 - Control  0.08000 0.414 215   0.193  1.0000
    ##  ZAN_F3 - Control   0.14333 0.414 215   0.346  1.0000
    ##  EMM_F64 - Control  0.16667 0.414 215   0.402  1.0000
    ##  EMM_F47 - Control  0.29667 0.414 215   0.716  1.0000
    ##  EMM_F89 - Control  0.30000 0.414 215   0.724  1.0000
    ##  EMM_F63 - Control  0.36333 0.414 215   0.877  1.0000
    ##  EMM_F5 - Control   0.37000 0.414 215   0.893  0.9999
    ##  EMM_F65 - Control  0.38333 0.414 215   0.926  0.9999
    ##  EMM_F3 - Control   0.40333 0.414 215   0.974  0.9998
    ##  EMM_F48 - Control  0.55667 0.414 215   1.344  0.9937
    ##  EMM_F7 - Control   0.60867 0.463 215   1.314  0.9951
    ##  EMM_F70 - Control  0.63000 0.414 215   1.521  0.9789
    ##  ZAN_F4 - Control   0.76667 0.414 215   1.851  0.8944
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.58333 0.414 215   1.408  0.9899
    ##  Control - EMM_F70  0.11333 0.414 215   0.274  1.0000
    ##  Control - EMM_F49  0.06667 0.414 215   0.161  1.0000
    ##  Control - EMM_F47  0.00000 0.414 215   0.000  1.0000
    ##  EMM_F5 - Control   0.00000 0.414 215   0.000  1.0000
    ##  EMM_F89 - Control  0.17333 0.414 215   0.419  1.0000
    ##  EMM_F34 - Control  0.17333 0.414 215   0.419  1.0000
    ##  EMM_F64 - Control  0.18333 0.414 215   0.443  1.0000
    ##  EMM_F48 - Control  0.25667 0.414 215   0.620  1.0000
    ##  EMM_F3 - Control   0.26333 0.414 215   0.636  1.0000
    ##  EMM_F63 - Control  0.29000 0.414 215   0.700  1.0000
    ##  EMM_F7 - Control   0.29367 0.463 215   0.634  1.0000
    ##  EMM_F65 - Control  0.32333 0.414 215   0.781  1.0000
    ##  ZAN_F3 - Control   0.42000 0.414 215   1.014  0.9997
    ##  ZAN_F4 - Control   0.42000 0.414 215   1.014  0.9997
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.73333 0.414 215   1.771  0.9240
    ##  Control - EMM_F70  0.39667 0.414 215   0.958  0.9999
    ##  Control - EMM_F64  0.35000 0.414 215   0.845  1.0000
    ##  Control - EMM_F7   0.26633 0.463 215   0.575  1.0000
    ##  Control - EMM_F63  0.21667 0.414 215   0.523  1.0000
    ##  Control - EMM_F47  0.15667 0.414 215   0.378  1.0000
    ##  Control - EMM_F5   0.14000 0.414 215   0.338  1.0000
    ##  Control - EMM_F49  0.13000 0.414 215   0.314  1.0000
    ##  Control - EMM_F89  0.11000 0.414 215   0.266  1.0000
    ##  Control - EMM_F3   0.01333 0.414 215   0.032  1.0000
    ##  EMM_F48 - Control  0.00333 0.414 215   0.008  1.0000
    ##  ZAN_F4 - Control   0.17667 0.414 215   0.427  1.0000
    ##  ZAN_F3 - Control   0.22000 0.414 215   0.531  1.0000
    ##  EMM_F65 - Control  0.25333 0.414 215   0.612  1.0000
    ##  EMM_F34 - Control  0.32667 0.414 215   0.789  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.64333 0.414 215   1.553  0.9745
    ##  Control - EMM_F47  0.42333 0.414 215   1.022  0.9997
    ##  Control - EMM_F70  0.38000 0.414 215   0.918  0.9999
    ##  Control - EMM_F7   0.34300 0.463 215   0.740  1.0000
    ##  Control - EMM_F49  0.18000 0.414 215   0.435  1.0000
    ##  Control - EMM_F64  0.04667 0.414 215   0.113  1.0000
    ##  EMM_F48 - Control  0.00667 0.414 215   0.016  1.0000
    ##  ZAN_F4 - Control   0.02000 0.414 215   0.048  1.0000
    ##  EMM_F89 - Control  0.06000 0.414 215   0.145  1.0000
    ##  EMM_F63 - Control  0.15333 0.414 215   0.370  1.0000
    ##  EMM_F3 - Control   0.18333 0.414 215   0.443  1.0000
    ##  EMM_F5 - Control   0.32000 0.414 215   0.773  1.0000
    ##  EMM_F65 - Control  0.33000 0.414 215   0.797  1.0000
    ##  ZAN_F3 - Control   0.41333 0.414 215   0.998  0.9998
    ##  EMM_F34 - Control  0.47667 0.414 215   1.151  0.9988
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F7   1.16800 0.463 215   2.521  0.4606
    ##  Control - SP_F14   1.07000 0.414 215   2.584  0.4157
    ##  Control - EMM_F47  0.83000 0.414 215   2.004  0.8204
    ##  Control - EMM_F3   0.63667 0.414 215   1.537  0.9768
    ##  Control - EMM_F49  0.61333 0.414 215   1.481  0.9836
    ##  Control - EMM_F64  0.50000 0.414 215   1.207  0.9980
    ##  Control - EMM_F34  0.38333 0.414 215   0.926  0.9999
    ##  Control - EMM_F70  0.38000 0.414 215   0.918  0.9999
    ##  Control - EMM_F63  0.31000 0.414 215   0.749  1.0000
    ##  Control - ZAN_F4   0.23667 0.414 215   0.571  1.0000
    ##  Control - EMM_F5   0.22667 0.414 215   0.547  1.0000
    ##  Control - EMM_F89  0.07000 0.414 215   0.169  1.0000
    ##  Control - EMM_F48  0.02667 0.414 215   0.064  1.0000
    ##  ZAN_F3 - Control   0.16333 0.414 215   0.394  1.0000
    ##  EMM_F65 - Control  0.28333 0.414 215   0.684  1.0000
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   1.44333 0.430 215   3.360  0.0687
    ##  Control - EMM_F49  1.29000 0.430 215   3.003  0.1766
    ##  Control - EMM_F34  1.02333 0.430 215   2.382  0.5630
    ##  Control - EMM_F3   0.82000 0.430 215   1.909  0.8691
    ##  Control - EMM_F5   0.80667 0.430 215   1.878  0.8831
    ##  Control - EMM_F47  0.79000 0.430 215   1.839  0.8992
    ##  Control - EMM_F89  0.57667 0.430 215   1.342  0.9938
    ##  Control - EMM_F63  0.49667 0.430 215   1.156  0.9988
    ##  Control - ZAN_F4   0.49667 0.430 215   1.156  0.9988
    ##  Control - EMM_F48  0.46667 0.430 215   1.086  0.9994
    ##  Control - EMM_F70  0.45333 0.430 215   1.055  0.9996
    ##  Control - EMM_F64  0.30000 0.430 215   0.698  1.0000
    ##  Control - ZAN_F3   0.28000 0.430 215   0.652  1.0000
    ##  Control - EMM_F65  0.24333 0.430 215   0.566  1.0000
    ##  Control - EMM_F7   0.02263 0.481 215   0.047  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.49667 0.430 215   1.156  0.9988
    ##  Control - EMM_F49  0.13667 0.430 215   0.318  1.0000
    ##  Control - EMM_F3   0.03000 0.430 215   0.070  1.0000
    ##  EMM_F34 - Control  0.03000 0.430 215   0.070  1.0000
    ##  EMM_F89 - Control  0.04667 0.430 215   0.109  1.0000
    ##  ZAN_F3 - Control   0.05333 0.430 215   0.124  1.0000
    ##  EMM_F5 - Control   0.14000 0.430 215   0.326  1.0000
    ##  EMM_F7 - Control   0.23403 0.481 215   0.487  1.0000
    ##  EMM_F48 - Control  0.24333 0.430 215   0.566  1.0000
    ##  EMM_F64 - Control  0.24333 0.430 215   0.566  1.0000
    ##  EMM_F47 - Control  0.27667 0.430 215   0.644  1.0000
    ##  EMM_F63 - Control  0.33333 0.430 215   0.776  1.0000
    ##  ZAN_F4 - Control   0.35333 0.430 215   0.823  1.0000
    ##  EMM_F65 - Control  0.49000 0.430 215   1.141  0.9990
    ##  EMM_F70 - Control  0.57667 0.430 215   1.342  0.9938
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   1.10667 0.430 215   2.576  0.4208
    ##  Control - EMM_F34  0.50667 0.430 215   1.180  0.9985
    ##  Control - EMM_F49  0.40667 0.430 215   0.947  0.9999
    ##  Control - ZAN_F3   0.17333 0.430 215   0.404  1.0000
    ##  Control - EMM_F47  0.12333 0.430 215   0.287  1.0000
    ##  Control - EMM_F7   0.07597 0.481 215   0.158  1.0000
    ##  Control - EMM_F5   0.04667 0.430 215   0.109  1.0000
    ##  Control - EMM_F65  0.04333 0.430 215   0.101  1.0000
    ##  Control - EMM_F64  0.02667 0.430 215   0.062  1.0000
    ##  Control - EMM_F89  0.02333 0.430 215   0.054  1.0000
    ##  Control - EMM_F63  0.01333 0.430 215   0.031  1.0000
    ##  EMM_F3 - Control   0.05333 0.430 215   0.124  1.0000
    ##  EMM_F48 - Control  0.15667 0.430 215   0.365  1.0000
    ##  ZAN_F4 - Control   0.19667 0.430 215   0.458  1.0000
    ##  EMM_F70 - Control  0.24000 0.430 215   0.559  1.0000
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.88333 0.430 215   2.056  0.7901
    ##  Control - EMM_F47  0.48667 0.430 215   1.133  0.9990
    ##  Control - EMM_F70  0.27333 0.430 215   0.636  1.0000
    ##  Control - EMM_F5   0.26333 0.430 215   0.613  1.0000
    ##  Control - EMM_F49  0.15333 0.430 215   0.357  1.0000
    ##  Control - EMM_F89  0.08000 0.430 215   0.186  1.0000
    ##  Control - EMM_F34  0.04000 0.430 215   0.093  1.0000
    ##  Control - ZAN_F4   0.00333 0.430 215   0.008  1.0000
    ##  EMM_F3 - Control   0.00000 0.430 215   0.000  1.0000
    ##  EMM_F65 - Control  0.03000 0.430 215   0.070  1.0000
    ##  EMM_F7 - Control   0.06403 0.481 215   0.133  1.0000
    ##  EMM_F64 - Control  0.14000 0.430 215   0.326  1.0000
    ##  EMM_F63 - Control  0.16000 0.430 215   0.372  1.0000
    ##  EMM_F48 - Control  0.22000 0.430 215   0.512  1.0000
    ##  ZAN_F3 - Control   0.31667 0.430 215   0.737  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   1.08333 0.430 215   2.522  0.4596
    ##  Control - EMM_F89  0.54000 0.430 215   1.257  0.9969
    ##  Control - EMM_F7   0.46763 0.481 215   0.973  0.9998
    ##  Control - EMM_F47  0.45667 0.430 215   1.063  0.9995
    ##  Control - EMM_F70  0.35667 0.430 215   0.830  1.0000
    ##  Control - EMM_F49  0.31333 0.430 215   0.729  1.0000
    ##  Control - EMM_F64  0.28333 0.430 215   0.660  1.0000
    ##  Control - EMM_F63  0.19333 0.430 215   0.450  1.0000
    ##  Control - EMM_F3   0.18000 0.430 215   0.419  1.0000
    ##  Control - EMM_F5   0.17667 0.430 215   0.411  1.0000
    ##  Control - EMM_F34  0.07000 0.430 215   0.163  1.0000
    ##  Control - EMM_F48  0.01667 0.430 215   0.039  1.0000
    ##  Control - ZAN_F4   0.01333 0.430 215   0.031  1.0000
    ##  EMM_F65 - Control  0.19000 0.430 215   0.442  1.0000
    ##  ZAN_F3 - Control   0.21333 0.430 215   0.497  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   1.04667 0.430 215   2.437  0.5224
    ##  Control - EMM_F7   0.61597 0.481 215   1.282  0.9962
    ##  Control - EMM_F47  0.52667 0.430 215   1.226  0.9976
    ##  Control - EMM_F49  0.38667 0.430 215   0.900  0.9999
    ##  Control - EMM_F64  0.34667 0.430 215   0.807  1.0000
    ##  Control - EMM_F70  0.28667 0.430 215   0.667  1.0000
    ##  Control - EMM_F3   0.24667 0.430 215   0.574  1.0000
    ##  Control - EMM_F5   0.21667 0.430 215   0.504  1.0000
    ##  Control - ZAN_F4   0.15333 0.430 215   0.357  1.0000
    ##  Control - EMM_F48  0.06333 0.430 215   0.147  1.0000
    ##  EMM_F63 - Control  0.06000 0.430 215   0.140  1.0000
    ##  EMM_F89 - Control  0.09333 0.430 215   0.217  1.0000
    ##  EMM_F34 - Control  0.21667 0.430 215   0.504  1.0000
    ##  EMM_F65 - Control  0.30000 0.430 215   0.698  1.0000
    ##  ZAN_F3 - Control   0.35000 0.430 215   0.815  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F7   1.43763 0.481 215   2.992  0.1815
    ##  Control - SP_F14   1.10000 0.430 215   2.561  0.4318
    ##  Control - EMM_F49  1.00333 0.430 215   2.336  0.5977
    ##  Control - EMM_F3   0.80667 0.430 215   1.878  0.8831
    ##  Control - EMM_F47  0.70667 0.430 215   1.645  0.9579
    ##  Control - EMM_F64  0.68333 0.430 215   1.591  0.9685
    ##  Control - ZAN_F4   0.62000 0.430 215   1.443  0.9872
    ##  Control - EMM_F34  0.50333 0.430 215   1.172  0.9986
    ##  Control - EMM_F5   0.46333 0.430 215   1.079  0.9995
    ##  Control - EMM_F48  0.32333 0.430 215   0.753  1.0000
    ##  Control - EMM_F70  0.32333 0.430 215   0.753  1.0000
    ##  Control - EMM_F63  0.31000 0.430 215   0.722  1.0000
    ##  Control - EMM_F89  0.29333 0.430 215   0.683  1.0000
    ##  Control - ZAN_F3   0.15000 0.430 215   0.349  1.0000
    ##  Control - EMM_F65  0.04000 0.430 215   0.093  1.0000
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates

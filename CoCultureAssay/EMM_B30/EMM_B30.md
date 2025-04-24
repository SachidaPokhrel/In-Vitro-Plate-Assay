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

## ***Curtobacterium flaccumfaciens* EMM_B30**

Plot is generated using loop around the 4 different classes of yeast,
coming up with 4 plots as an output which will be combined in one plot.

``` r
#read data
B30 <- read.csv("CoCultureAssay/CoCultureAssayData/2024-07-21_PeaceAssay_B30.csv", na.strings = "na") 


#load cbb color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
str(B30)
```

    ## 'data.frame':    1536 obs. of  9 variables:
    ##  $ Bacteria               : chr  "B30" "B30" "B30" "B30" ...
    ##  $ Yeast                  : chr  "EMM_F3" "EMM_F3" "EMM_F3" "EMM_F3" ...
    ##  $ Class                  : chr  "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" ...
    ##  $ Replication            : int  1 1 1 1 1 1 1 1 2 2 ...
    ##  $ DAI                    : int  2 2 2 2 2 2 2 2 2 2 ...
    ##  $ distance_to_yeast      : num  0 11.4 17.4 24.9 31.8 ...
    ##  $ colony_diameter        : num  7.42 7.42 8.06 8.06 8.22 7.99 8.11 6.61 7.05 7.83 ...
    ##  $ colony_diameter_control: num  8.01 8 7.69 7.41 6.92 7.72 6.82 7.08 7.71 7.54 ...
    ##  $ increase               : num  0 0 0 0 0 0 0 0 0 0 ...

``` r
#round off of the distance to yeast so that it is visibly a categorical variable
B30$distance_to_yeast <- round(B30$distance_to_yeast, 0)
#change into categorical dataset
B30$Replication=as.factor(B30$Replication)
B30$Yeast=as.factor(B30$Yeast)
B30$DAI=as.factor(B30$DAI)
B30$distance_to_yeast=as.factor(B30$distance_to_yeast)

#remove the data for the distance zero
B30.no.contact <- B30[B30$distance_to_yeast != 0, ]
head(B30.no.contact)
```

    ##   Bacteria  Yeast           Class Replication DAI distance_to_yeast
    ## 2      B30 EMM_F3 Dothideomycetes           1   2                11
    ## 3      B30 EMM_F3 Dothideomycetes           1   2                17
    ## 4      B30 EMM_F3 Dothideomycetes           1   2                25
    ## 5      B30 EMM_F3 Dothideomycetes           1   2                32
    ## 6      B30 EMM_F3 Dothideomycetes           1   2                41
    ## 7      B30 EMM_F3 Dothideomycetes           1   2                48
    ##   colony_diameter colony_diameter_control increase
    ## 2            7.42                    8.00        0
    ## 3            8.06                    7.69        0
    ## 4            8.06                    7.41        0
    ## 5            8.22                    6.92        0
    ## 6            7.99                    7.72        0
    ## 7            8.11                    6.82        0

``` r
# remove data first data since it is a reference for increase in colony size and donot follow the assumption of linear model i.e. normal distribution
B30.no.1st.data <-B30.no.contact[B30.no.contact$DAI != "2",]


B30.Day6 <- B30.no.contact[B30.no.contact$DAI == "6",]


plot_listB30 <- list()

# Get unique classes
unique_class <- unique(B30.Day6$Class[B30.Day6$Class != "Control"])

# Loop through each class
for (i in seq_along(unique_class)) {
  
  # Filter for Control and the current class
  data <- B30.Day6 %>%
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
  plot_listB30[[i]] <- p
}

combined_plotB30 <- ggarrange(plotlist = plot_listB30, ncol = 2, nrow = 2, common.legend = FALSE)
```

    ## Warning: Removed 7 rows containing non-finite outside the scale range
    ## (`stat_summary()`).
    ## Removed 7 rows containing non-finite outside the scale range
    ## (`stat_summary()`).
    ## Removed 7 rows containing non-finite outside the scale range
    ## (`stat_summary()`).
    ## Removed 7 rows containing non-finite outside the scale range
    ## (`stat_summary()`).

``` r
# Annotation for the title
final_plotB30 <- annotate_figure(combined_plotB30,
                               top = text_grob(
    expression("Impact on growth of"~italic("Curtobacterium flaccumfaciens")~"EMM_B30 by Yeast"),, color = "Blue2", face = "bold", size = 14, hjust = 0.5))

print(final_plotB30)
```

![](EMM_B30_files/figure-gfm/Plot%20for%20EMM_B30-1.png)<!-- -->

***After first bacteria everything is repeated in similar way for all
dataset***

### Stats for *Curtobacterium flaccumfaciens* EMM_B30

We are using linear mixed model. Our dependent variable or y is increase
(increase in colony diameter from 1st data) and independent variables
are different Yeast isolates, days after inoculation (DAI), and distance
to yeast which is the distance between the yeast and bacterial colony in
the plate. Each plate is replicated 3 times.

``` r
#filter data to remove 1st day data since the first data is taken as a base to measure the increase in colony size to rule out the variability that is caused by the drop inoculation. So, initially the increase in the colony diameter for 1st data for all colony is "0" that violates the assumption of normality, thus we remove that from analysis. This would be similar for all the bacterial isolates.
B30.no.1st.data <- B30.no.contact[B30.no.contact$DAI != "2",]
B30try <- lme(increase~DAI*distance_to_yeast*Yeast, data = B30.no.1st.data, random = ~1|Replication, na.action = na.omit)
anova(B30try)
```

    ##                             numDF denDF   F-value p-value
    ## (Intercept)                     1   614 1387.2730  <.0001
    ## DAI                             2   614   47.2606  <.0001
    ## distance_to_yeast               6   614    3.2216  0.0040
    ## Yeast                          15   614    6.2149  <.0001
    ## DAI:distance_to_yeast          12   614    0.1294  0.9998
    ## DAI:Yeast                      30   614    1.1796  0.2360
    ## distance_to_yeast:Yeast        90   614    1.5588  0.0015
    ## DAI:distance_to_yeast:Yeast   180   614    0.2569  1.0000

``` r
resultsB30=lme(increase~DAI+distance_to_yeast*Yeast, data = B30.no.1st.data, random = ~1|Replication, na.action = na.omit)
summary(resultsB30)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: B30.no.1st.data 
    ##        AIC      BIC    logLik
    ##   1517.143 2065.942 -642.5717
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##          (Intercept)  Residual
    ## StdDev: 4.051196e-05 0.4488257
    ## 
    ## Fixed effects:  increase ~ DAI + distance_to_yeast * Yeast 
    ##                                       Value Std.Error  DF   t-value p-value
    ## (Intercept)                       0.1411622 0.1509854 836  0.934939  0.3501
    ## DAI6                              0.2114083 0.0350495 836  6.031705  0.0000
    ## DAI8                              0.3784383 0.0359489 836 10.527115  0.0000
    ## distance_to_yeast17              -0.0688889 0.2115784 836 -0.325595  0.7448
    ## distance_to_yeast25              -0.0555556 0.2115784 836 -0.262577  0.7929
    ## distance_to_yeast32              -0.0400000 0.2115784 836 -0.189055  0.8501
    ## distance_to_yeast41               0.0188889 0.2115784 836  0.089276  0.9289
    ## distance_to_yeast48               0.0311111 0.2115784 836  0.147043  0.8831
    ## distance_to_yeast55               0.1200000 0.2115784 836  0.567166  0.5708
    ## YeastEMM_F3                       0.2688889 0.2115784 836  1.270871  0.2041
    ## YeastEMM_F34                     -0.1600499 0.2181059 836 -0.733817  0.4633
    ## YeastEMM_F48                      0.4949501 0.2181059 836  2.269310  0.0235
    ## YeastEMM_F49                      0.1631673 0.2262055 836  0.721323  0.4709
    ## YeastEMM_F5                       0.6011111 0.2115784 836  2.841079  0.0046
    ## YeastEMM_F63                      0.3733333 0.2115784 836  1.764515  0.0780
    ## YeastEMM_F64                      0.0588889 0.2115784 836  0.278331  0.7808
    ## YeastEMM_F65                      0.0637001 0.2181059 836  0.292060  0.7703
    ## YeastEMM_F66                      0.2611111 0.2115784 836  1.234110  0.2175
    ## YeastEMM_F7                       0.2088889 0.2115784 836  0.987288  0.3238
    ## YeastEMM_F70                     -0.0455556 0.2115784 836 -0.215313  0.8296
    ## YeastEMM_F89                      0.0174530 0.2262055 836  0.077156  0.9385
    ## YeastSP_F14                      -0.2477778 0.2115784 836 -1.171092  0.2419
    ## YeastZAN_F3                       0.8787001 0.2181059 836  4.028776  0.0001
    ## YeastZAN_F4                      -0.1477778 0.2115784 836 -0.698454  0.4851
    ## distance_to_yeast17:YeastEMM_F3  -0.1166667 0.2992171 836 -0.389906  0.6967
    ## distance_to_yeast25:YeastEMM_F3  -0.1144444 0.2992171 836 -0.382480  0.7022
    ## distance_to_yeast32:YeastEMM_F3  -0.2522222 0.2992171 836 -0.842941  0.3995
    ## distance_to_yeast41:YeastEMM_F3  -0.1511111 0.2992171 836 -0.505022  0.6137
    ## distance_to_yeast48:YeastEMM_F3  -0.1600000 0.2992171 836 -0.534729  0.5930
    ## distance_to_yeast55:YeastEMM_F3  -0.2333333 0.2992171 836 -0.779813  0.4357
    ## distance_to_yeast17:YeastEMM_F34  0.2613889 0.3084259 836  0.847493  0.3970
    ## distance_to_yeast25:YeastEMM_F34  0.0768056 0.3084259 836  0.249024  0.8034
    ## distance_to_yeast32:YeastEMM_F34  0.1575000 0.3084259 836  0.510657  0.6097
    ## distance_to_yeast41:YeastEMM_F34  0.2411111 0.3084259 836  0.781747  0.4346
    ## distance_to_yeast48:YeastEMM_F34  0.3101389 0.3084259 836  1.005554  0.3149
    ## distance_to_yeast55:YeastEMM_F34 -0.0487500 0.3084259 836 -0.158061  0.8744
    ## distance_to_yeast17:YeastEMM_F48  0.2176389 0.3084259 836  0.705644  0.4806
    ## distance_to_yeast25:YeastEMM_F48  0.0980556 0.3084259 836  0.317923  0.7506
    ## distance_to_yeast32:YeastEMM_F48 -0.0137500 0.3084259 836 -0.044581  0.9645
    ## distance_to_yeast41:YeastEMM_F48 -0.2463889 0.3084259 836 -0.798859  0.4246
    ## distance_to_yeast48:YeastEMM_F48 -0.4098611 0.3084259 836 -1.328880  0.1843
    ## distance_to_yeast55:YeastEMM_F48 -0.1250000 0.3084259 836 -0.405284  0.6854
    ## distance_to_yeast17:YeastEMM_F49  0.2260317 0.3198765 836  0.706622  0.4800
    ## distance_to_yeast25:YeastEMM_F49  0.0069841 0.3198765 836  0.021834  0.9826
    ## distance_to_yeast32:YeastEMM_F49 -0.0414286 0.3198765 836 -0.129514  0.8970
    ## distance_to_yeast41:YeastEMM_F49 -0.3088889 0.3198765 836 -0.965650  0.3345
    ## distance_to_yeast48:YeastEMM_F49 -0.0139683 0.3198765 836 -0.043668  0.9652
    ## distance_to_yeast55:YeastEMM_F49 -0.3528571 0.3198765 836 -1.103104  0.2703
    ## distance_to_yeast17:YeastEMM_F5  -0.0811111 0.2992171 836 -0.271078  0.7864
    ## distance_to_yeast25:YeastEMM_F5  -0.0644444 0.2992171 836 -0.215377  0.8295
    ## distance_to_yeast32:YeastEMM_F5  -0.3500000 0.2992171 836 -1.169719  0.2424
    ## distance_to_yeast41:YeastEMM_F5  -0.4566667 0.2992171 836 -1.526205  0.1273
    ## distance_to_yeast48:YeastEMM_F5  -0.3344444 0.2992171 836 -1.117732  0.2640
    ## distance_to_yeast55:YeastEMM_F5  -0.3011111 0.2992171 836 -1.006330  0.3145
    ## distance_to_yeast17:YeastEMM_F63  0.0522222 0.2992171 836  0.174530  0.8615
    ## distance_to_yeast25:YeastEMM_F63  0.0977778 0.2992171 836  0.326779  0.7439
    ## distance_to_yeast32:YeastEMM_F63  0.3311111 0.2992171 836  1.106592  0.2688
    ## distance_to_yeast41:YeastEMM_F63  0.1366667 0.2992171 836  0.456748  0.6480
    ## distance_to_yeast48:YeastEMM_F63 -0.0488889 0.2992171 836 -0.163389  0.8703
    ## distance_to_yeast55:YeastEMM_F63 -0.1544444 0.2992171 836 -0.516162  0.6059
    ## distance_to_yeast17:YeastEMM_F64  0.2822222 0.2992171 836  0.943202  0.3458
    ## distance_to_yeast25:YeastEMM_F64  0.4700000 0.2992171 836  1.570766  0.1166
    ## distance_to_yeast32:YeastEMM_F64  0.5900000 0.2992171 836  1.971812  0.0490
    ## distance_to_yeast41:YeastEMM_F64  0.3688889 0.2992171 836  1.232847  0.2180
    ## distance_to_yeast48:YeastEMM_F64  0.4233333 0.2992171 836  1.414803  0.1575
    ## distance_to_yeast55:YeastEMM_F64  0.9100000 0.2992171 836  3.041270  0.0024
    ## distance_to_yeast17:YeastEMM_F65  0.4788889 0.3084259 836  1.552687  0.1209
    ## distance_to_yeast25:YeastEMM_F65  0.4280556 0.3084259 836  1.387872  0.1655
    ## distance_to_yeast32:YeastEMM_F65 -0.1350000 0.3084259 836 -0.437706  0.6617
    ## distance_to_yeast41:YeastEMM_F65  0.2586111 0.3084259 836  0.838487  0.4020
    ## distance_to_yeast48:YeastEMM_F65  0.6238889 0.3084259 836  2.022816  0.0434
    ## distance_to_yeast55:YeastEMM_F65  0.2425000 0.3084259 836  0.786250  0.4319
    ## distance_to_yeast17:YeastEMM_F66 -0.0777778 0.2992171 836 -0.259938  0.7950
    ## distance_to_yeast25:YeastEMM_F66 -0.2722222 0.2992171 836 -0.909782  0.3632
    ## distance_to_yeast32:YeastEMM_F66 -0.0066667 0.2992171 836 -0.022280  0.9822
    ## distance_to_yeast41:YeastEMM_F66 -0.1622222 0.2992171 836 -0.542156  0.5879
    ## distance_to_yeast48:YeastEMM_F66  0.1633333 0.2992171 836  0.545869  0.5853
    ## distance_to_yeast55:YeastEMM_F66 -0.0566667 0.2992171 836 -0.189383  0.8498
    ## distance_to_yeast17:YeastEMM_F7   0.0855556 0.2992171 836  0.285931  0.7750
    ## distance_to_yeast25:YeastEMM_F7   0.3355556 0.2992171 836  1.121445  0.2624
    ## distance_to_yeast32:YeastEMM_F7   0.0944444 0.2992171 836  0.315639  0.7524
    ## distance_to_yeast41:YeastEMM_F7   0.2044444 0.2992171 836  0.683265  0.4946
    ## distance_to_yeast48:YeastEMM_F7   0.0533333 0.2992171 836  0.178243  0.8586
    ## distance_to_yeast55:YeastEMM_F7   0.1244444 0.2992171 836  0.415900  0.6776
    ## distance_to_yeast17:YeastEMM_F70  0.2788889 0.2992171 836  0.932062  0.3516
    ## distance_to_yeast25:YeastEMM_F70  0.2555556 0.2992171 836  0.854081  0.3933
    ## distance_to_yeast32:YeastEMM_F70  0.3688889 0.2992171 836  1.232847  0.2180
    ## distance_to_yeast41:YeastEMM_F70  0.1977778 0.2992171 836  0.660984  0.5088
    ## distance_to_yeast48:YeastEMM_F70  0.4222222 0.2992171 836  1.411090  0.1586
    ## distance_to_yeast55:YeastEMM_F70  0.2033333 0.2992171 836  0.679551  0.4970
    ## distance_to_yeast17:YeastEMM_F89  0.1331746 0.3198765 836  0.416331  0.6773
    ## distance_to_yeast25:YeastEMM_F89  0.0826984 0.3198765 836  0.258532  0.7961
    ## distance_to_yeast32:YeastEMM_F89  0.3285714 0.3198765 836  1.027182  0.3046
    ## distance_to_yeast41:YeastEMM_F89  0.3739683 0.3198765 836  1.169102  0.2427
    ## distance_to_yeast48:YeastEMM_F89  0.6688889 0.3198765 836  2.091085  0.0368
    ## distance_to_yeast55:YeastEMM_F89  0.0585714 0.3198765 836  0.183106  0.8548
    ## distance_to_yeast17:YeastSP_F14   0.4633333 0.2992171 836  1.548485  0.1219
    ## distance_to_yeast25:YeastSP_F14   0.5933333 0.2992171 836  1.982953  0.0477
    ## distance_to_yeast32:YeastSP_F14   0.4444444 0.2992171 836  1.485358  0.1378
    ## distance_to_yeast41:YeastSP_F14   0.8266667 0.2992171 836  2.762765  0.0059
    ## distance_to_yeast48:YeastSP_F14   0.6544444 0.2992171 836  2.187189  0.0290
    ## distance_to_yeast55:YeastSP_F14   0.7588889 0.2992171 836  2.536248  0.0114
    ## distance_to_yeast17:YeastZAN_F3  -0.7473611 0.3084259 836 -2.423146  0.0156
    ## distance_to_yeast25:YeastZAN_F3  -0.7081944 0.3084259 836 -2.296157  0.0219
    ## distance_to_yeast32:YeastZAN_F3  -0.3150000 0.3084259 836 -1.021315  0.3074
    ## distance_to_yeast41:YeastZAN_F3  -0.2938889 0.3084259 836 -0.952867  0.3409
    ## distance_to_yeast48:YeastZAN_F3   0.1938889 0.3084259 836  0.628640  0.5298
    ## distance_to_yeast55:YeastZAN_F3  -1.0737500 0.3084259 836 -3.481387  0.0005
    ## distance_to_yeast17:YeastZAN_F4   0.1588889 0.2992171 836  0.531015  0.5955
    ## distance_to_yeast25:YeastZAN_F4   0.2877778 0.2992171 836  0.961769  0.3364
    ## distance_to_yeast32:YeastZAN_F4   0.6022222 0.2992171 836  2.012660  0.0445
    ## distance_to_yeast41:YeastZAN_F4   0.3388889 0.2992171 836  1.132585  0.2577
    ## distance_to_yeast48:YeastZAN_F4   0.4544444 0.2992171 836  1.518778  0.1292
    ## distance_to_yeast55:YeastZAN_F4   0.4888889 0.2992171 836  1.633894  0.1027
    ##  Correlation: 
    ##                                  (Intr) DAI6   DAI8   ds__17 ds__25 ds__32
    ## DAI6                             -0.115                                   
    ## DAI8                             -0.116  0.478                            
    ## distance_to_yeast17              -0.701  0.000  0.000                     
    ## distance_to_yeast25              -0.701  0.000  0.000  0.500              
    ## distance_to_yeast32              -0.701  0.000  0.000  0.500  0.500       
    ## distance_to_yeast41              -0.701  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast48              -0.701  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast55              -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F3                      -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F34                     -0.681  0.000  0.011  0.485  0.485  0.485
    ## YeastEMM_F48                     -0.681  0.000  0.011  0.485  0.485  0.485
    ## YeastEMM_F49                     -0.657  0.011  0.011  0.468  0.468  0.468
    ## YeastEMM_F5                      -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F63                     -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F64                     -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F65                     -0.681  0.000  0.011  0.485  0.485  0.485
    ## YeastEMM_F66                     -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F7                      -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F70                     -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastEMM_F89                     -0.657  0.011  0.011  0.468  0.468  0.468
    ## YeastSP_F14                      -0.701  0.000  0.000  0.500  0.500  0.500
    ## YeastZAN_F3                      -0.681  0.000  0.011  0.485  0.485  0.485
    ## YeastZAN_F4                      -0.701  0.000  0.000  0.500  0.500  0.500
    ## distance_to_yeast17:YeastEMM_F3   0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F3   0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F3   0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F3   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F3   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F3   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F34  0.481  0.000  0.000 -0.686 -0.343 -0.343
    ## distance_to_yeast25:YeastEMM_F34  0.481  0.000  0.000 -0.343 -0.686 -0.343
    ## distance_to_yeast32:YeastEMM_F34  0.481  0.000  0.000 -0.343 -0.343 -0.686
    ## distance_to_yeast41:YeastEMM_F34  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast48:YeastEMM_F34  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast55:YeastEMM_F34  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast17:YeastEMM_F48  0.481  0.000  0.000 -0.686 -0.343 -0.343
    ## distance_to_yeast25:YeastEMM_F48  0.481  0.000  0.000 -0.343 -0.686 -0.343
    ## distance_to_yeast32:YeastEMM_F48  0.481  0.000  0.000 -0.343 -0.343 -0.686
    ## distance_to_yeast41:YeastEMM_F48  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast48:YeastEMM_F48  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast55:YeastEMM_F48  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast17:YeastEMM_F49  0.463  0.000  0.000 -0.661 -0.331 -0.331
    ## distance_to_yeast25:YeastEMM_F49  0.463  0.000  0.000 -0.331 -0.661 -0.331
    ## distance_to_yeast32:YeastEMM_F49  0.463  0.000  0.000 -0.331 -0.331 -0.661
    ## distance_to_yeast41:YeastEMM_F49  0.463  0.000  0.000 -0.331 -0.331 -0.331
    ## distance_to_yeast48:YeastEMM_F49  0.463  0.000  0.000 -0.331 -0.331 -0.331
    ## distance_to_yeast55:YeastEMM_F49  0.463  0.000  0.000 -0.331 -0.331 -0.331
    ## distance_to_yeast17:YeastEMM_F5   0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F5   0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F5   0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F5   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F5   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F5   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F63  0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F63  0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F63  0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F63  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F63  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F63  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F64  0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F64  0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F64  0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F64  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F64  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F64  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F65  0.481  0.000  0.000 -0.686 -0.343 -0.343
    ## distance_to_yeast25:YeastEMM_F65  0.481  0.000  0.000 -0.343 -0.686 -0.343
    ## distance_to_yeast32:YeastEMM_F65  0.481  0.000  0.000 -0.343 -0.343 -0.686
    ## distance_to_yeast41:YeastEMM_F65  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast48:YeastEMM_F65  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast55:YeastEMM_F65  0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast17:YeastEMM_F66  0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F66  0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F66  0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F66  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F66  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F66  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F7   0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F7   0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F7   0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F7   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F7   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F7   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F70  0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F70  0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastEMM_F70  0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastEMM_F70  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastEMM_F70  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastEMM_F70  0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F89  0.463  0.000  0.000 -0.661 -0.331 -0.331
    ## distance_to_yeast25:YeastEMM_F89  0.463  0.000  0.000 -0.331 -0.661 -0.331
    ## distance_to_yeast32:YeastEMM_F89  0.463  0.000  0.000 -0.331 -0.331 -0.661
    ## distance_to_yeast41:YeastEMM_F89  0.463  0.000  0.000 -0.331 -0.331 -0.331
    ## distance_to_yeast48:YeastEMM_F89  0.463  0.000  0.000 -0.331 -0.331 -0.331
    ## distance_to_yeast55:YeastEMM_F89  0.463  0.000  0.000 -0.331 -0.331 -0.331
    ## distance_to_yeast17:YeastSP_F14   0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastSP_F14   0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastSP_F14   0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastSP_F14   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastSP_F14   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastSP_F14   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastZAN_F3   0.481  0.000  0.000 -0.686 -0.343 -0.343
    ## distance_to_yeast25:YeastZAN_F3   0.481  0.000  0.000 -0.343 -0.686 -0.343
    ## distance_to_yeast32:YeastZAN_F3   0.481  0.000  0.000 -0.343 -0.343 -0.686
    ## distance_to_yeast41:YeastZAN_F3   0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast48:YeastZAN_F3   0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast55:YeastZAN_F3   0.481  0.000  0.000 -0.343 -0.343 -0.343
    ## distance_to_yeast17:YeastZAN_F4   0.495  0.000  0.000 -0.707 -0.354 -0.354
    ## distance_to_yeast25:YeastZAN_F4   0.495  0.000  0.000 -0.354 -0.707 -0.354
    ## distance_to_yeast32:YeastZAN_F4   0.495  0.000  0.000 -0.354 -0.354 -0.707
    ## distance_to_yeast41:YeastZAN_F4   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast48:YeastZAN_F4   0.495  0.000  0.000 -0.354 -0.354 -0.354
    ## distance_to_yeast55:YeastZAN_F4   0.495  0.000  0.000 -0.354 -0.354 -0.354
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
    ## YeastEMM_F34                      0.485  0.485  0.485  0.485           
    ## YeastEMM_F48                      0.485  0.485  0.485  0.485    0.471  
    ## YeastEMM_F49                      0.468  0.468  0.468  0.468    0.454  
    ## YeastEMM_F5                       0.500  0.500  0.500  0.500    0.485  
    ## YeastEMM_F63                      0.500  0.500  0.500  0.500    0.485  
    ## YeastEMM_F64                      0.500  0.500  0.500  0.500    0.485  
    ## YeastEMM_F65                      0.485  0.485  0.485  0.485    0.471  
    ## YeastEMM_F66                      0.500  0.500  0.500  0.500    0.485  
    ## YeastEMM_F7                       0.500  0.500  0.500  0.500    0.485  
    ## YeastEMM_F70                      0.500  0.500  0.500  0.500    0.485  
    ## YeastEMM_F89                      0.468  0.468  0.468  0.468    0.454  
    ## YeastSP_F14                       0.500  0.500  0.500  0.500    0.485  
    ## YeastZAN_F3                       0.485  0.485  0.485  0.485    0.471  
    ## YeastZAN_F4                       0.500  0.500  0.500  0.500    0.485  
    ## distance_to_yeast17:YeastEMM_F3  -0.354 -0.354 -0.354 -0.707   -0.343  
    ## distance_to_yeast25:YeastEMM_F3  -0.354 -0.354 -0.354 -0.707   -0.343  
    ## distance_to_yeast32:YeastEMM_F3  -0.354 -0.354 -0.354 -0.707   -0.343  
    ## distance_to_yeast41:YeastEMM_F3  -0.707 -0.354 -0.354 -0.707   -0.343  
    ## distance_to_yeast48:YeastEMM_F3  -0.354 -0.707 -0.354 -0.707   -0.343  
    ## distance_to_yeast55:YeastEMM_F3  -0.354 -0.354 -0.707 -0.707   -0.343  
    ## distance_to_yeast17:YeastEMM_F34 -0.343 -0.343 -0.343 -0.343   -0.707  
    ## distance_to_yeast25:YeastEMM_F34 -0.343 -0.343 -0.343 -0.343   -0.707  
    ## distance_to_yeast32:YeastEMM_F34 -0.343 -0.343 -0.343 -0.343   -0.707  
    ## distance_to_yeast41:YeastEMM_F34 -0.686 -0.343 -0.343 -0.343   -0.707  
    ## distance_to_yeast48:YeastEMM_F34 -0.343 -0.686 -0.343 -0.343   -0.707  
    ## distance_to_yeast55:YeastEMM_F34 -0.343 -0.343 -0.686 -0.343   -0.707  
    ## distance_to_yeast17:YeastEMM_F48 -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast25:YeastEMM_F48 -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast32:YeastEMM_F48 -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast41:YeastEMM_F48 -0.686 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast48:YeastEMM_F48 -0.343 -0.686 -0.343 -0.343   -0.333  
    ## distance_to_yeast55:YeastEMM_F48 -0.343 -0.343 -0.686 -0.343   -0.333  
    ## distance_to_yeast17:YeastEMM_F49 -0.331 -0.331 -0.331 -0.331   -0.321  
    ## distance_to_yeast25:YeastEMM_F49 -0.331 -0.331 -0.331 -0.331   -0.321  
    ## distance_to_yeast32:YeastEMM_F49 -0.331 -0.331 -0.331 -0.331   -0.321  
    ## distance_to_yeast41:YeastEMM_F49 -0.661 -0.331 -0.331 -0.331   -0.321  
    ## distance_to_yeast48:YeastEMM_F49 -0.331 -0.661 -0.331 -0.331   -0.321  
    ## distance_to_yeast55:YeastEMM_F49 -0.331 -0.331 -0.661 -0.331   -0.321  
    ## distance_to_yeast17:YeastEMM_F5  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast25:YeastEMM_F5  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast32:YeastEMM_F5  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast41:YeastEMM_F5  -0.707 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast48:YeastEMM_F5  -0.354 -0.707 -0.354 -0.354   -0.343  
    ## distance_to_yeast55:YeastEMM_F5  -0.354 -0.354 -0.707 -0.354   -0.343  
    ## distance_to_yeast17:YeastEMM_F63 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast25:YeastEMM_F63 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast32:YeastEMM_F63 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast41:YeastEMM_F63 -0.707 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast48:YeastEMM_F63 -0.354 -0.707 -0.354 -0.354   -0.343  
    ## distance_to_yeast55:YeastEMM_F63 -0.354 -0.354 -0.707 -0.354   -0.343  
    ## distance_to_yeast17:YeastEMM_F64 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast25:YeastEMM_F64 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast32:YeastEMM_F64 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast41:YeastEMM_F64 -0.707 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast48:YeastEMM_F64 -0.354 -0.707 -0.354 -0.354   -0.343  
    ## distance_to_yeast55:YeastEMM_F64 -0.354 -0.354 -0.707 -0.354   -0.343  
    ## distance_to_yeast17:YeastEMM_F65 -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast25:YeastEMM_F65 -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast32:YeastEMM_F65 -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast41:YeastEMM_F65 -0.686 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast48:YeastEMM_F65 -0.343 -0.686 -0.343 -0.343   -0.333  
    ## distance_to_yeast55:YeastEMM_F65 -0.343 -0.343 -0.686 -0.343   -0.333  
    ## distance_to_yeast17:YeastEMM_F66 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast25:YeastEMM_F66 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast32:YeastEMM_F66 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast41:YeastEMM_F66 -0.707 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast48:YeastEMM_F66 -0.354 -0.707 -0.354 -0.354   -0.343  
    ## distance_to_yeast55:YeastEMM_F66 -0.354 -0.354 -0.707 -0.354   -0.343  
    ## distance_to_yeast17:YeastEMM_F7  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast25:YeastEMM_F7  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast32:YeastEMM_F7  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast41:YeastEMM_F7  -0.707 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast48:YeastEMM_F7  -0.354 -0.707 -0.354 -0.354   -0.343  
    ## distance_to_yeast55:YeastEMM_F7  -0.354 -0.354 -0.707 -0.354   -0.343  
    ## distance_to_yeast17:YeastEMM_F70 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast25:YeastEMM_F70 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast32:YeastEMM_F70 -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast41:YeastEMM_F70 -0.707 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast48:YeastEMM_F70 -0.354 -0.707 -0.354 -0.354   -0.343  
    ## distance_to_yeast55:YeastEMM_F70 -0.354 -0.354 -0.707 -0.354   -0.343  
    ## distance_to_yeast17:YeastEMM_F89 -0.331 -0.331 -0.331 -0.331   -0.321  
    ## distance_to_yeast25:YeastEMM_F89 -0.331 -0.331 -0.331 -0.331   -0.321  
    ## distance_to_yeast32:YeastEMM_F89 -0.331 -0.331 -0.331 -0.331   -0.321  
    ## distance_to_yeast41:YeastEMM_F89 -0.661 -0.331 -0.331 -0.331   -0.321  
    ## distance_to_yeast48:YeastEMM_F89 -0.331 -0.661 -0.331 -0.331   -0.321  
    ## distance_to_yeast55:YeastEMM_F89 -0.331 -0.331 -0.661 -0.331   -0.321  
    ## distance_to_yeast17:YeastSP_F14  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast25:YeastSP_F14  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast32:YeastSP_F14  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast41:YeastSP_F14  -0.707 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast48:YeastSP_F14  -0.354 -0.707 -0.354 -0.354   -0.343  
    ## distance_to_yeast55:YeastSP_F14  -0.354 -0.354 -0.707 -0.354   -0.343  
    ## distance_to_yeast17:YeastZAN_F3  -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast25:YeastZAN_F3  -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast32:YeastZAN_F3  -0.343 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast41:YeastZAN_F3  -0.686 -0.343 -0.343 -0.343   -0.333  
    ## distance_to_yeast48:YeastZAN_F3  -0.343 -0.686 -0.343 -0.343   -0.333  
    ## distance_to_yeast55:YeastZAN_F3  -0.343 -0.343 -0.686 -0.343   -0.333  
    ## distance_to_yeast17:YeastZAN_F4  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast25:YeastZAN_F4  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast32:YeastZAN_F4  -0.354 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast41:YeastZAN_F4  -0.707 -0.354 -0.354 -0.354   -0.343  
    ## distance_to_yeast48:YeastZAN_F4  -0.354 -0.707 -0.354 -0.354   -0.343  
    ## distance_to_yeast55:YeastZAN_F4  -0.354 -0.354 -0.707 -0.354   -0.343  
    ##                                  YEMM_F48 YEMM_F49 YEMM_F5 YEMM_F63 YEMM_F64
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
    ## YeastEMM_F48                                                                
    ## YeastEMM_F49                      0.454                                     
    ## YeastEMM_F5                       0.485    0.468                            
    ## YeastEMM_F63                      0.485    0.468    0.500                   
    ## YeastEMM_F64                      0.485    0.468    0.500   0.500           
    ## YeastEMM_F65                      0.471    0.454    0.485   0.485    0.485  
    ## YeastEMM_F66                      0.485    0.468    0.500   0.500    0.500  
    ## YeastEMM_F7                       0.485    0.468    0.500   0.500    0.500  
    ## YeastEMM_F70                      0.485    0.468    0.500   0.500    0.500  
    ## YeastEMM_F89                      0.454    0.438    0.468   0.468    0.468  
    ## YeastSP_F14                       0.485    0.468    0.500   0.500    0.500  
    ## YeastZAN_F3                       0.471    0.454    0.485   0.485    0.485  
    ## YeastZAN_F4                       0.485    0.468    0.500   0.500    0.500  
    ## distance_to_yeast17:YeastEMM_F3  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F3  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F3  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F3  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F3  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F3  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F34 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast25:YeastEMM_F34 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast32:YeastEMM_F34 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast41:YeastEMM_F34 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast48:YeastEMM_F34 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast55:YeastEMM_F34 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast17:YeastEMM_F48 -0.707   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast25:YeastEMM_F48 -0.707   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast32:YeastEMM_F48 -0.707   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast41:YeastEMM_F48 -0.707   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast48:YeastEMM_F48 -0.707   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast55:YeastEMM_F48 -0.707   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast17:YeastEMM_F49 -0.321   -0.707   -0.331  -0.331   -0.331  
    ## distance_to_yeast25:YeastEMM_F49 -0.321   -0.707   -0.331  -0.331   -0.331  
    ## distance_to_yeast32:YeastEMM_F49 -0.321   -0.707   -0.331  -0.331   -0.331  
    ## distance_to_yeast41:YeastEMM_F49 -0.321   -0.707   -0.331  -0.331   -0.331  
    ## distance_to_yeast48:YeastEMM_F49 -0.321   -0.707   -0.331  -0.331   -0.331  
    ## distance_to_yeast55:YeastEMM_F49 -0.321   -0.707   -0.331  -0.331   -0.331  
    ## distance_to_yeast17:YeastEMM_F5  -0.343   -0.331   -0.707  -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F5  -0.343   -0.331   -0.707  -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F5  -0.343   -0.331   -0.707  -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F5  -0.343   -0.331   -0.707  -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F5  -0.343   -0.331   -0.707  -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F5  -0.343   -0.331   -0.707  -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F63 -0.343   -0.331   -0.354  -0.707   -0.354  
    ## distance_to_yeast25:YeastEMM_F63 -0.343   -0.331   -0.354  -0.707   -0.354  
    ## distance_to_yeast32:YeastEMM_F63 -0.343   -0.331   -0.354  -0.707   -0.354  
    ## distance_to_yeast41:YeastEMM_F63 -0.343   -0.331   -0.354  -0.707   -0.354  
    ## distance_to_yeast48:YeastEMM_F63 -0.343   -0.331   -0.354  -0.707   -0.354  
    ## distance_to_yeast55:YeastEMM_F63 -0.343   -0.331   -0.354  -0.707   -0.354  
    ## distance_to_yeast17:YeastEMM_F64 -0.343   -0.331   -0.354  -0.354   -0.707  
    ## distance_to_yeast25:YeastEMM_F64 -0.343   -0.331   -0.354  -0.354   -0.707  
    ## distance_to_yeast32:YeastEMM_F64 -0.343   -0.331   -0.354  -0.354   -0.707  
    ## distance_to_yeast41:YeastEMM_F64 -0.343   -0.331   -0.354  -0.354   -0.707  
    ## distance_to_yeast48:YeastEMM_F64 -0.343   -0.331   -0.354  -0.354   -0.707  
    ## distance_to_yeast55:YeastEMM_F64 -0.343   -0.331   -0.354  -0.354   -0.707  
    ## distance_to_yeast17:YeastEMM_F65 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast25:YeastEMM_F65 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast32:YeastEMM_F65 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast41:YeastEMM_F65 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast48:YeastEMM_F65 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast55:YeastEMM_F65 -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast17:YeastEMM_F66 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F66 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F66 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F66 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F66 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F66 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F7  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F7  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F7  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F7  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F7  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F7  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F70 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast25:YeastEMM_F70 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast32:YeastEMM_F70 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast41:YeastEMM_F70 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast48:YeastEMM_F70 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast55:YeastEMM_F70 -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast17:YeastEMM_F89 -0.321   -0.309   -0.331  -0.331   -0.331  
    ## distance_to_yeast25:YeastEMM_F89 -0.321   -0.309   -0.331  -0.331   -0.331  
    ## distance_to_yeast32:YeastEMM_F89 -0.321   -0.309   -0.331  -0.331   -0.331  
    ## distance_to_yeast41:YeastEMM_F89 -0.321   -0.309   -0.331  -0.331   -0.331  
    ## distance_to_yeast48:YeastEMM_F89 -0.321   -0.309   -0.331  -0.331   -0.331  
    ## distance_to_yeast55:YeastEMM_F89 -0.321   -0.309   -0.331  -0.331   -0.331  
    ## distance_to_yeast17:YeastSP_F14  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast25:YeastSP_F14  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast32:YeastSP_F14  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast41:YeastSP_F14  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast48:YeastSP_F14  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast55:YeastSP_F14  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast17:YeastZAN_F3  -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast25:YeastZAN_F3  -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast32:YeastZAN_F3  -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast41:YeastZAN_F3  -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast48:YeastZAN_F3  -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast55:YeastZAN_F3  -0.333   -0.321   -0.343  -0.343   -0.343  
    ## distance_to_yeast17:YeastZAN_F4  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast25:YeastZAN_F4  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast32:YeastZAN_F4  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast41:YeastZAN_F4  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast48:YeastZAN_F4  -0.343   -0.331   -0.354  -0.354   -0.354  
    ## distance_to_yeast55:YeastZAN_F4  -0.343   -0.331   -0.354  -0.354   -0.354  
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
    ## YeastEMM_F48                                                                
    ## YeastEMM_F49                                                                
    ## YeastEMM_F5                                                                 
    ## YeastEMM_F63                                                                
    ## YeastEMM_F64                                                                
    ## YeastEMM_F65                                                                
    ## YeastEMM_F66                      0.485                                     
    ## YeastEMM_F7                       0.485    0.500                            
    ## YeastEMM_F70                      0.485    0.500    0.500                   
    ## YeastEMM_F89                      0.454    0.468    0.468    0.468          
    ## YeastSP_F14                       0.485    0.500    0.500    0.500    0.468 
    ## YeastZAN_F3                       0.471    0.485    0.485    0.485    0.454 
    ## YeastZAN_F4                       0.485    0.500    0.500    0.500    0.468 
    ## distance_to_yeast17:YeastEMM_F3  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast25:YeastEMM_F3  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast32:YeastEMM_F3  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast41:YeastEMM_F3  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast48:YeastEMM_F3  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast55:YeastEMM_F3  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast17:YeastEMM_F34 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast25:YeastEMM_F34 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast32:YeastEMM_F34 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast41:YeastEMM_F34 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast48:YeastEMM_F34 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast55:YeastEMM_F34 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast17:YeastEMM_F48 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast25:YeastEMM_F48 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast32:YeastEMM_F48 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast41:YeastEMM_F48 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast48:YeastEMM_F48 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast55:YeastEMM_F48 -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast17:YeastEMM_F49 -0.321   -0.331   -0.331   -0.331   -0.309 
    ## distance_to_yeast25:YeastEMM_F49 -0.321   -0.331   -0.331   -0.331   -0.309 
    ## distance_to_yeast32:YeastEMM_F49 -0.321   -0.331   -0.331   -0.331   -0.309 
    ## distance_to_yeast41:YeastEMM_F49 -0.321   -0.331   -0.331   -0.331   -0.309 
    ## distance_to_yeast48:YeastEMM_F49 -0.321   -0.331   -0.331   -0.331   -0.309 
    ## distance_to_yeast55:YeastEMM_F49 -0.321   -0.331   -0.331   -0.331   -0.309 
    ## distance_to_yeast17:YeastEMM_F5  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast25:YeastEMM_F5  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast32:YeastEMM_F5  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast41:YeastEMM_F5  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast48:YeastEMM_F5  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast55:YeastEMM_F5  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast17:YeastEMM_F63 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast25:YeastEMM_F63 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast32:YeastEMM_F63 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast41:YeastEMM_F63 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast48:YeastEMM_F63 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast55:YeastEMM_F63 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast17:YeastEMM_F64 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast25:YeastEMM_F64 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast32:YeastEMM_F64 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast41:YeastEMM_F64 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast48:YeastEMM_F64 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast55:YeastEMM_F64 -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast17:YeastEMM_F65 -0.707   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast25:YeastEMM_F65 -0.707   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast32:YeastEMM_F65 -0.707   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast41:YeastEMM_F65 -0.707   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast48:YeastEMM_F65 -0.707   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast55:YeastEMM_F65 -0.707   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast17:YeastEMM_F66 -0.343   -0.707   -0.354   -0.354   -0.331 
    ## distance_to_yeast25:YeastEMM_F66 -0.343   -0.707   -0.354   -0.354   -0.331 
    ## distance_to_yeast32:YeastEMM_F66 -0.343   -0.707   -0.354   -0.354   -0.331 
    ## distance_to_yeast41:YeastEMM_F66 -0.343   -0.707   -0.354   -0.354   -0.331 
    ## distance_to_yeast48:YeastEMM_F66 -0.343   -0.707   -0.354   -0.354   -0.331 
    ## distance_to_yeast55:YeastEMM_F66 -0.343   -0.707   -0.354   -0.354   -0.331 
    ## distance_to_yeast17:YeastEMM_F7  -0.343   -0.354   -0.707   -0.354   -0.331 
    ## distance_to_yeast25:YeastEMM_F7  -0.343   -0.354   -0.707   -0.354   -0.331 
    ## distance_to_yeast32:YeastEMM_F7  -0.343   -0.354   -0.707   -0.354   -0.331 
    ## distance_to_yeast41:YeastEMM_F7  -0.343   -0.354   -0.707   -0.354   -0.331 
    ## distance_to_yeast48:YeastEMM_F7  -0.343   -0.354   -0.707   -0.354   -0.331 
    ## distance_to_yeast55:YeastEMM_F7  -0.343   -0.354   -0.707   -0.354   -0.331 
    ## distance_to_yeast17:YeastEMM_F70 -0.343   -0.354   -0.354   -0.707   -0.331 
    ## distance_to_yeast25:YeastEMM_F70 -0.343   -0.354   -0.354   -0.707   -0.331 
    ## distance_to_yeast32:YeastEMM_F70 -0.343   -0.354   -0.354   -0.707   -0.331 
    ## distance_to_yeast41:YeastEMM_F70 -0.343   -0.354   -0.354   -0.707   -0.331 
    ## distance_to_yeast48:YeastEMM_F70 -0.343   -0.354   -0.354   -0.707   -0.331 
    ## distance_to_yeast55:YeastEMM_F70 -0.343   -0.354   -0.354   -0.707   -0.331 
    ## distance_to_yeast17:YeastEMM_F89 -0.321   -0.331   -0.331   -0.331   -0.707 
    ## distance_to_yeast25:YeastEMM_F89 -0.321   -0.331   -0.331   -0.331   -0.707 
    ## distance_to_yeast32:YeastEMM_F89 -0.321   -0.331   -0.331   -0.331   -0.707 
    ## distance_to_yeast41:YeastEMM_F89 -0.321   -0.331   -0.331   -0.331   -0.707 
    ## distance_to_yeast48:YeastEMM_F89 -0.321   -0.331   -0.331   -0.331   -0.707 
    ## distance_to_yeast55:YeastEMM_F89 -0.321   -0.331   -0.331   -0.331   -0.707 
    ## distance_to_yeast17:YeastSP_F14  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast25:YeastSP_F14  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast32:YeastSP_F14  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast41:YeastSP_F14  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast48:YeastSP_F14  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast55:YeastSP_F14  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast17:YeastZAN_F3  -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast25:YeastZAN_F3  -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast32:YeastZAN_F3  -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast41:YeastZAN_F3  -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast48:YeastZAN_F3  -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast55:YeastZAN_F3  -0.333   -0.343   -0.343   -0.343   -0.321 
    ## distance_to_yeast17:YeastZAN_F4  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast25:YeastZAN_F4  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast32:YeastZAN_F4  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast41:YeastZAN_F4  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast48:YeastZAN_F4  -0.343   -0.354   -0.354   -0.354   -0.331 
    ## distance_to_yeast55:YeastZAN_F4  -0.343   -0.354   -0.354   -0.354   -0.331 
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
    ## YeastEMM_F48                                                          
    ## YeastEMM_F49                                                          
    ## YeastEMM_F5                                                           
    ## YeastEMM_F63                                                          
    ## YeastEMM_F64                                                          
    ## YeastEMM_F65                                                          
    ## YeastEMM_F66                                                          
    ## YeastEMM_F7                                                           
    ## YeastEMM_F70                                                          
    ## YeastEMM_F89                                                          
    ## YeastSP_F14                                                           
    ## YeastZAN_F3                       0.485                               
    ## YeastZAN_F4                       0.500  0.485                        
    ## distance_to_yeast17:YeastEMM_F3  -0.354 -0.343  -0.354                
    ## distance_to_yeast25:YeastEMM_F3  -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast32:YeastEMM_F3  -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast41:YeastEMM_F3  -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast48:YeastEMM_F3  -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast55:YeastEMM_F3  -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast17:YeastEMM_F34 -0.343 -0.333  -0.343   0.485        
    ## distance_to_yeast25:YeastEMM_F34 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast32:YeastEMM_F34 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast41:YeastEMM_F34 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast48:YeastEMM_F34 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast55:YeastEMM_F34 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast17:YeastEMM_F48 -0.343 -0.333  -0.343   0.485        
    ## distance_to_yeast25:YeastEMM_F48 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast32:YeastEMM_F48 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast41:YeastEMM_F48 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast48:YeastEMM_F48 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast55:YeastEMM_F48 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast17:YeastEMM_F49 -0.331 -0.321  -0.331   0.468        
    ## distance_to_yeast25:YeastEMM_F49 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast32:YeastEMM_F49 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast41:YeastEMM_F49 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast48:YeastEMM_F49 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast55:YeastEMM_F49 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast17:YeastEMM_F5  -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F5  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F5  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F5  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F5  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F5  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F63 -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F63 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F63 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F63 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F63 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F63 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F64 -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F64 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F64 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F64 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F64 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F64 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F65 -0.343 -0.333  -0.343   0.485        
    ## distance_to_yeast25:YeastEMM_F65 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast32:YeastEMM_F65 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast41:YeastEMM_F65 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast48:YeastEMM_F65 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast55:YeastEMM_F65 -0.343 -0.333  -0.343   0.243        
    ## distance_to_yeast17:YeastEMM_F66 -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F66 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F66 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F66 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F66 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F66 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F7  -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F7  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F7  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F7  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F7  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F7  -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F70 -0.354 -0.343  -0.354   0.500        
    ## distance_to_yeast25:YeastEMM_F70 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast32:YeastEMM_F70 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast41:YeastEMM_F70 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast48:YeastEMM_F70 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast55:YeastEMM_F70 -0.354 -0.343  -0.354   0.250        
    ## distance_to_yeast17:YeastEMM_F89 -0.331 -0.321  -0.331   0.468        
    ## distance_to_yeast25:YeastEMM_F89 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast32:YeastEMM_F89 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast41:YeastEMM_F89 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast48:YeastEMM_F89 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast55:YeastEMM_F89 -0.331 -0.321  -0.331   0.234        
    ## distance_to_yeast17:YeastSP_F14  -0.707 -0.343  -0.354   0.500        
    ## distance_to_yeast25:YeastSP_F14  -0.707 -0.343  -0.354   0.250        
    ## distance_to_yeast32:YeastSP_F14  -0.707 -0.343  -0.354   0.250        
    ## distance_to_yeast41:YeastSP_F14  -0.707 -0.343  -0.354   0.250        
    ## distance_to_yeast48:YeastSP_F14  -0.707 -0.343  -0.354   0.250        
    ## distance_to_yeast55:YeastSP_F14  -0.707 -0.343  -0.354   0.250        
    ## distance_to_yeast17:YeastZAN_F3  -0.343 -0.707  -0.343   0.485        
    ## distance_to_yeast25:YeastZAN_F3  -0.343 -0.707  -0.343   0.243        
    ## distance_to_yeast32:YeastZAN_F3  -0.343 -0.707  -0.343   0.243        
    ## distance_to_yeast41:YeastZAN_F3  -0.343 -0.707  -0.343   0.243        
    ## distance_to_yeast48:YeastZAN_F3  -0.343 -0.707  -0.343   0.243        
    ## distance_to_yeast55:YeastZAN_F3  -0.343 -0.707  -0.343   0.243        
    ## distance_to_yeast17:YeastZAN_F4  -0.354 -0.343  -0.707   0.500        
    ## distance_to_yeast25:YeastZAN_F4  -0.354 -0.343  -0.707   0.250        
    ## distance_to_yeast32:YeastZAN_F4  -0.354 -0.343  -0.707   0.250        
    ## distance_to_yeast41:YeastZAN_F4  -0.354 -0.343  -0.707   0.250        
    ## distance_to_yeast48:YeastZAN_F4  -0.354 -0.343  -0.707   0.250        
    ## distance_to_yeast55:YeastZAN_F4  -0.354 -0.343  -0.707   0.250        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F34  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F34  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F34  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F34  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F34  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F34  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F48  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F48  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F48  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F48  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F48  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F48  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F49  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F49  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F49  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F49  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F49  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F49  0.234          0.234          0.234        
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
    ## distance_to_yeast17:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F65  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F65  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F65  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F65  0.243          0.243          0.243        
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.243          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F34  0.243          0.243                       
    ## distance_to_yeast25:YeastEMM_F34  0.243          0.243          0.500        
    ## distance_to_yeast32:YeastEMM_F34  0.243          0.243          0.500        
    ## distance_to_yeast41:YeastEMM_F34  0.243          0.243          0.500        
    ## distance_to_yeast48:YeastEMM_F34  0.485          0.243          0.500        
    ## distance_to_yeast55:YeastEMM_F34  0.243          0.485          0.500        
    ## distance_to_yeast17:YeastEMM_F48  0.243          0.243          0.471        
    ## distance_to_yeast25:YeastEMM_F48  0.243          0.243          0.235        
    ## distance_to_yeast32:YeastEMM_F48  0.243          0.243          0.235        
    ## distance_to_yeast41:YeastEMM_F48  0.243          0.243          0.235        
    ## distance_to_yeast48:YeastEMM_F48  0.485          0.243          0.235        
    ## distance_to_yeast55:YeastEMM_F48  0.243          0.485          0.235        
    ## distance_to_yeast17:YeastEMM_F49  0.234          0.234          0.454        
    ## distance_to_yeast25:YeastEMM_F49  0.234          0.234          0.227        
    ## distance_to_yeast32:YeastEMM_F49  0.234          0.234          0.227        
    ## distance_to_yeast41:YeastEMM_F49  0.234          0.234          0.227        
    ## distance_to_yeast48:YeastEMM_F49  0.468          0.234          0.227        
    ## distance_to_yeast55:YeastEMM_F49  0.234          0.468          0.227        
    ## distance_to_yeast17:YeastEMM_F5   0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F5   0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F5   0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F5   0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F5   0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F5   0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F65  0.243          0.243          0.471        
    ## distance_to_yeast25:YeastEMM_F65  0.243          0.243          0.235        
    ## distance_to_yeast32:YeastEMM_F65  0.243          0.243          0.235        
    ## distance_to_yeast41:YeastEMM_F65  0.243          0.243          0.235        
    ## distance_to_yeast48:YeastEMM_F65  0.485          0.243          0.235        
    ## distance_to_yeast55:YeastEMM_F65  0.243          0.485          0.235        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.454        
    ## distance_to_yeast25:YeastEMM_F89  0.234          0.234          0.227        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.234          0.227        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.227        
    ## distance_to_yeast48:YeastEMM_F89  0.468          0.234          0.227        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.468          0.227        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.485        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.243        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.243        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.243        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.243        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.243        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.471        
    ## distance_to_yeast25:YeastZAN_F3   0.243          0.243          0.235        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.243          0.235        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.235        
    ## distance_to_yeast48:YeastZAN_F3   0.485          0.243          0.235        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.485          0.235        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.485        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.243        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.243        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.243        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.243        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F48  0.235          0.235          0.235        
    ## distance_to_yeast25:YeastEMM_F48  0.471          0.235          0.235        
    ## distance_to_yeast32:YeastEMM_F48  0.235          0.471          0.235        
    ## distance_to_yeast41:YeastEMM_F48  0.235          0.235          0.471        
    ## distance_to_yeast48:YeastEMM_F48  0.235          0.235          0.235        
    ## distance_to_yeast55:YeastEMM_F48  0.235          0.235          0.235        
    ## distance_to_yeast17:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast25:YeastEMM_F49  0.454          0.227          0.227        
    ## distance_to_yeast32:YeastEMM_F49  0.227          0.454          0.227        
    ## distance_to_yeast41:YeastEMM_F49  0.227          0.227          0.454        
    ## distance_to_yeast48:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast55:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast17:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F5   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F5   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F5   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F63  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F63  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F63  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F64  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F64  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F64  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast25:YeastEMM_F65  0.471          0.235          0.235        
    ## distance_to_yeast32:YeastEMM_F65  0.235          0.471          0.235        
    ## distance_to_yeast41:YeastEMM_F65  0.235          0.235          0.471        
    ## distance_to_yeast48:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast55:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast17:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F66  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F66  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F66  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F7   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F7   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F7   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F70  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F70  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F70  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast25:YeastEMM_F89  0.454          0.227          0.227        
    ## distance_to_yeast32:YeastEMM_F89  0.227          0.454          0.227        
    ## distance_to_yeast41:YeastEMM_F89  0.227          0.227          0.454        
    ## distance_to_yeast48:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast55:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast17:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastSP_F14   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastSP_F14   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastSP_F14   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast17:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast25:YeastZAN_F3   0.471          0.235          0.235        
    ## distance_to_yeast32:YeastZAN_F3   0.235          0.471          0.235        
    ## distance_to_yeast41:YeastZAN_F3   0.235          0.235          0.471        
    ## distance_to_yeast48:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast55:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast17:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F4   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F4   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F4   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F4   0.243          0.243          0.243        
    ##                                  d__48:YEMM_F34 d__55:YEMM_F34 d__17:YEMM_F48
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F48  0.235          0.235                       
    ## distance_to_yeast25:YeastEMM_F48  0.235          0.235          0.500        
    ## distance_to_yeast32:YeastEMM_F48  0.235          0.235          0.500        
    ## distance_to_yeast41:YeastEMM_F48  0.235          0.235          0.500        
    ## distance_to_yeast48:YeastEMM_F48  0.471          0.235          0.500        
    ## distance_to_yeast55:YeastEMM_F48  0.235          0.471          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.227          0.227          0.454        
    ## distance_to_yeast25:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast32:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast41:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast48:YeastEMM_F49  0.454          0.227          0.227        
    ## distance_to_yeast55:YeastEMM_F49  0.227          0.454          0.227        
    ## distance_to_yeast17:YeastEMM_F5   0.243          0.243          0.485        
    ## distance_to_yeast25:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast41:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast48:YeastEMM_F5   0.485          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F5   0.243          0.485          0.243        
    ## distance_to_yeast17:YeastEMM_F63  0.243          0.243          0.485        
    ## distance_to_yeast25:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast41:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast48:YeastEMM_F63  0.485          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F63  0.243          0.485          0.243        
    ## distance_to_yeast17:YeastEMM_F64  0.243          0.243          0.485        
    ## distance_to_yeast25:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast41:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast48:YeastEMM_F64  0.485          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F64  0.243          0.485          0.243        
    ## distance_to_yeast17:YeastEMM_F65  0.235          0.235          0.471        
    ## distance_to_yeast25:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast32:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast41:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast48:YeastEMM_F65  0.471          0.235          0.235        
    ## distance_to_yeast55:YeastEMM_F65  0.235          0.471          0.235        
    ## distance_to_yeast17:YeastEMM_F66  0.243          0.243          0.485        
    ## distance_to_yeast25:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast41:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast48:YeastEMM_F66  0.485          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F66  0.243          0.485          0.243        
    ## distance_to_yeast17:YeastEMM_F7   0.243          0.243          0.485        
    ## distance_to_yeast25:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast41:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast48:YeastEMM_F7   0.485          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F7   0.243          0.485          0.243        
    ## distance_to_yeast17:YeastEMM_F70  0.243          0.243          0.485        
    ## distance_to_yeast25:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast41:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast48:YeastEMM_F70  0.485          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F70  0.243          0.485          0.243        
    ## distance_to_yeast17:YeastEMM_F89  0.227          0.227          0.454        
    ## distance_to_yeast25:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast32:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast41:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast48:YeastEMM_F89  0.454          0.227          0.227        
    ## distance_to_yeast55:YeastEMM_F89  0.227          0.454          0.227        
    ## distance_to_yeast17:YeastSP_F14   0.243          0.243          0.485        
    ## distance_to_yeast25:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast32:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast41:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast48:YeastSP_F14   0.485          0.243          0.243        
    ## distance_to_yeast55:YeastSP_F14   0.243          0.485          0.243        
    ## distance_to_yeast17:YeastZAN_F3   0.235          0.235          0.471        
    ## distance_to_yeast25:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast32:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast41:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast48:YeastZAN_F3   0.471          0.235          0.235        
    ## distance_to_yeast55:YeastZAN_F3   0.235          0.471          0.235        
    ## distance_to_yeast17:YeastZAN_F4   0.243          0.243          0.485        
    ## distance_to_yeast25:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast41:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast48:YeastZAN_F4   0.485          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F4   0.243          0.485          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48  0.500                                      
    ## distance_to_yeast41:YeastEMM_F48  0.500          0.500                       
    ## distance_to_yeast48:YeastEMM_F48  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F48  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast25:YeastEMM_F49  0.454          0.227          0.227        
    ## distance_to_yeast32:YeastEMM_F49  0.227          0.454          0.227        
    ## distance_to_yeast41:YeastEMM_F49  0.227          0.227          0.454        
    ## distance_to_yeast48:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast55:YeastEMM_F49  0.227          0.227          0.227        
    ## distance_to_yeast17:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F5   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F5   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F5   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F5   0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F63  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F63  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F63  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F63  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F64  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F64  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F64  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F64  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast25:YeastEMM_F65  0.471          0.235          0.235        
    ## distance_to_yeast32:YeastEMM_F65  0.235          0.471          0.235        
    ## distance_to_yeast41:YeastEMM_F65  0.235          0.235          0.471        
    ## distance_to_yeast48:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast55:YeastEMM_F65  0.235          0.235          0.235        
    ## distance_to_yeast17:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F66  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F66  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F66  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F7   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F7   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F7   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F70  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F70  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F70  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast25:YeastEMM_F89  0.454          0.227          0.227        
    ## distance_to_yeast32:YeastEMM_F89  0.227          0.454          0.227        
    ## distance_to_yeast41:YeastEMM_F89  0.227          0.227          0.454        
    ## distance_to_yeast48:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast55:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast17:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastSP_F14   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastSP_F14   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastSP_F14   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast17:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast25:YeastZAN_F3   0.471          0.235          0.235        
    ## distance_to_yeast32:YeastZAN_F3   0.235          0.471          0.235        
    ## distance_to_yeast41:YeastZAN_F3   0.235          0.235          0.471        
    ## distance_to_yeast48:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast55:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast17:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F4   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F4   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F4   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F4   0.243          0.243          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F48                                             
    ## distance_to_yeast25:YeastEMM_F48                                             
    ## distance_to_yeast32:YeastEMM_F48                                             
    ## distance_to_yeast41:YeastEMM_F48                                             
    ## distance_to_yeast48:YeastEMM_F48                                             
    ## distance_to_yeast55:YeastEMM_F48  0.500                                      
    ## distance_to_yeast17:YeastEMM_F49  0.227          0.227                       
    ## distance_to_yeast25:YeastEMM_F49  0.227          0.227          0.500        
    ## distance_to_yeast32:YeastEMM_F49  0.227          0.227          0.500        
    ## distance_to_yeast41:YeastEMM_F49  0.227          0.227          0.500        
    ## distance_to_yeast48:YeastEMM_F49  0.454          0.227          0.500        
    ## distance_to_yeast55:YeastEMM_F49  0.227          0.454          0.500        
    ## distance_to_yeast17:YeastEMM_F5   0.243          0.243          0.468        
    ## distance_to_yeast25:YeastEMM_F5   0.243          0.243          0.234        
    ## distance_to_yeast32:YeastEMM_F5   0.243          0.243          0.234        
    ## distance_to_yeast41:YeastEMM_F5   0.243          0.243          0.234        
    ## distance_to_yeast48:YeastEMM_F5   0.485          0.243          0.234        
    ## distance_to_yeast55:YeastEMM_F5   0.243          0.485          0.234        
    ## distance_to_yeast17:YeastEMM_F63  0.243          0.243          0.468        
    ## distance_to_yeast25:YeastEMM_F63  0.243          0.243          0.234        
    ## distance_to_yeast32:YeastEMM_F63  0.243          0.243          0.234        
    ## distance_to_yeast41:YeastEMM_F63  0.243          0.243          0.234        
    ## distance_to_yeast48:YeastEMM_F63  0.485          0.243          0.234        
    ## distance_to_yeast55:YeastEMM_F63  0.243          0.485          0.234        
    ## distance_to_yeast17:YeastEMM_F64  0.243          0.243          0.468        
    ## distance_to_yeast25:YeastEMM_F64  0.243          0.243          0.234        
    ## distance_to_yeast32:YeastEMM_F64  0.243          0.243          0.234        
    ## distance_to_yeast41:YeastEMM_F64  0.243          0.243          0.234        
    ## distance_to_yeast48:YeastEMM_F64  0.485          0.243          0.234        
    ## distance_to_yeast55:YeastEMM_F64  0.243          0.485          0.234        
    ## distance_to_yeast17:YeastEMM_F65  0.235          0.235          0.454        
    ## distance_to_yeast25:YeastEMM_F65  0.235          0.235          0.227        
    ## distance_to_yeast32:YeastEMM_F65  0.235          0.235          0.227        
    ## distance_to_yeast41:YeastEMM_F65  0.235          0.235          0.227        
    ## distance_to_yeast48:YeastEMM_F65  0.471          0.235          0.227        
    ## distance_to_yeast55:YeastEMM_F65  0.235          0.471          0.227        
    ## distance_to_yeast17:YeastEMM_F66  0.243          0.243          0.468        
    ## distance_to_yeast25:YeastEMM_F66  0.243          0.243          0.234        
    ## distance_to_yeast32:YeastEMM_F66  0.243          0.243          0.234        
    ## distance_to_yeast41:YeastEMM_F66  0.243          0.243          0.234        
    ## distance_to_yeast48:YeastEMM_F66  0.485          0.243          0.234        
    ## distance_to_yeast55:YeastEMM_F66  0.243          0.485          0.234        
    ## distance_to_yeast17:YeastEMM_F7   0.243          0.243          0.468        
    ## distance_to_yeast25:YeastEMM_F7   0.243          0.243          0.234        
    ## distance_to_yeast32:YeastEMM_F7   0.243          0.243          0.234        
    ## distance_to_yeast41:YeastEMM_F7   0.243          0.243          0.234        
    ## distance_to_yeast48:YeastEMM_F7   0.485          0.243          0.234        
    ## distance_to_yeast55:YeastEMM_F7   0.243          0.485          0.234        
    ## distance_to_yeast17:YeastEMM_F70  0.243          0.243          0.468        
    ## distance_to_yeast25:YeastEMM_F70  0.243          0.243          0.234        
    ## distance_to_yeast32:YeastEMM_F70  0.243          0.243          0.234        
    ## distance_to_yeast41:YeastEMM_F70  0.243          0.243          0.234        
    ## distance_to_yeast48:YeastEMM_F70  0.485          0.243          0.234        
    ## distance_to_yeast55:YeastEMM_F70  0.243          0.485          0.234        
    ## distance_to_yeast17:YeastEMM_F89  0.227          0.227          0.437        
    ## distance_to_yeast25:YeastEMM_F89  0.227          0.227          0.219        
    ## distance_to_yeast32:YeastEMM_F89  0.227          0.227          0.219        
    ## distance_to_yeast41:YeastEMM_F89  0.227          0.227          0.219        
    ## distance_to_yeast48:YeastEMM_F89  0.454          0.227          0.219        
    ## distance_to_yeast55:YeastEMM_F89  0.227          0.454          0.219        
    ## distance_to_yeast17:YeastSP_F14   0.243          0.243          0.468        
    ## distance_to_yeast25:YeastSP_F14   0.243          0.243          0.234        
    ## distance_to_yeast32:YeastSP_F14   0.243          0.243          0.234        
    ## distance_to_yeast41:YeastSP_F14   0.243          0.243          0.234        
    ## distance_to_yeast48:YeastSP_F14   0.485          0.243          0.234        
    ## distance_to_yeast55:YeastSP_F14   0.243          0.485          0.234        
    ## distance_to_yeast17:YeastZAN_F3   0.235          0.235          0.454        
    ## distance_to_yeast25:YeastZAN_F3   0.235          0.235          0.227        
    ## distance_to_yeast32:YeastZAN_F3   0.235          0.235          0.227        
    ## distance_to_yeast41:YeastZAN_F3   0.235          0.235          0.227        
    ## distance_to_yeast48:YeastZAN_F3   0.471          0.235          0.227        
    ## distance_to_yeast55:YeastZAN_F3   0.235          0.471          0.227        
    ## distance_to_yeast17:YeastZAN_F4   0.243          0.243          0.468        
    ## distance_to_yeast25:YeastZAN_F4   0.243          0.243          0.234        
    ## distance_to_yeast32:YeastZAN_F4   0.243          0.243          0.234        
    ## distance_to_yeast41:YeastZAN_F4   0.243          0.243          0.234        
    ## distance_to_yeast48:YeastZAN_F4   0.485          0.243          0.234        
    ## distance_to_yeast55:YeastZAN_F4   0.243          0.485          0.234        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F5   0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F5   0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F5   0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F5   0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F5   0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F5   0.234          0.234          0.234        
    ## distance_to_yeast17:YeastEMM_F63  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F63  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F63  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F63  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F63  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F63  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastEMM_F64  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F64  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F64  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F64  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F64  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F64  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastEMM_F65  0.227          0.227          0.227        
    ## distance_to_yeast25:YeastEMM_F65  0.454          0.227          0.227        
    ## distance_to_yeast32:YeastEMM_F65  0.227          0.454          0.227        
    ## distance_to_yeast41:YeastEMM_F65  0.227          0.227          0.454        
    ## distance_to_yeast48:YeastEMM_F65  0.227          0.227          0.227        
    ## distance_to_yeast55:YeastEMM_F65  0.227          0.227          0.227        
    ## distance_to_yeast17:YeastEMM_F66  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F66  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F66  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F66  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F66  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F66  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastEMM_F7   0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F7   0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F7   0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F7   0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F7   0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F7   0.234          0.234          0.234        
    ## distance_to_yeast17:YeastEMM_F70  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F70  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F70  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F70  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F70  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F70  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastEMM_F89  0.219          0.219          0.219        
    ## distance_to_yeast25:YeastEMM_F89  0.437          0.219          0.219        
    ## distance_to_yeast32:YeastEMM_F89  0.219          0.437          0.219        
    ## distance_to_yeast41:YeastEMM_F89  0.219          0.219          0.437        
    ## distance_to_yeast48:YeastEMM_F89  0.219          0.219          0.219        
    ## distance_to_yeast55:YeastEMM_F89  0.219          0.219          0.219        
    ## distance_to_yeast17:YeastSP_F14   0.234          0.234          0.234        
    ## distance_to_yeast25:YeastSP_F14   0.468          0.234          0.234        
    ## distance_to_yeast32:YeastSP_F14   0.234          0.468          0.234        
    ## distance_to_yeast41:YeastSP_F14   0.234          0.234          0.468        
    ## distance_to_yeast48:YeastSP_F14   0.234          0.234          0.234        
    ## distance_to_yeast55:YeastSP_F14   0.234          0.234          0.234        
    ## distance_to_yeast17:YeastZAN_F3   0.227          0.227          0.227        
    ## distance_to_yeast25:YeastZAN_F3   0.454          0.227          0.227        
    ## distance_to_yeast32:YeastZAN_F3   0.227          0.454          0.227        
    ## distance_to_yeast41:YeastZAN_F3   0.227          0.227          0.454        
    ## distance_to_yeast48:YeastZAN_F3   0.227          0.227          0.227        
    ## distance_to_yeast55:YeastZAN_F3   0.227          0.227          0.227        
    ## distance_to_yeast17:YeastZAN_F4   0.234          0.234          0.234        
    ## distance_to_yeast25:YeastZAN_F4   0.468          0.234          0.234        
    ## distance_to_yeast32:YeastZAN_F4   0.234          0.468          0.234        
    ## distance_to_yeast41:YeastZAN_F4   0.234          0.234          0.468        
    ## distance_to_yeast48:YeastZAN_F4   0.234          0.234          0.234        
    ## distance_to_yeast55:YeastZAN_F4   0.234          0.234          0.234        
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
    ## YeastEMM_F48                                                                
    ## YeastEMM_F49                                                                
    ## YeastEMM_F5                                                                 
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
    ## distance_to_yeast17:YeastEMM_F5   0.234          0.234                      
    ## distance_to_yeast25:YeastEMM_F5   0.234          0.234          0.500       
    ## distance_to_yeast32:YeastEMM_F5   0.234          0.234          0.500       
    ## distance_to_yeast41:YeastEMM_F5   0.234          0.234          0.500       
    ## distance_to_yeast48:YeastEMM_F5   0.468          0.234          0.500       
    ## distance_to_yeast55:YeastEMM_F5   0.234          0.468          0.500       
    ## distance_to_yeast17:YeastEMM_F63  0.234          0.234          0.500       
    ## distance_to_yeast25:YeastEMM_F63  0.234          0.234          0.250       
    ## distance_to_yeast32:YeastEMM_F63  0.234          0.234          0.250       
    ## distance_to_yeast41:YeastEMM_F63  0.234          0.234          0.250       
    ## distance_to_yeast48:YeastEMM_F63  0.468          0.234          0.250       
    ## distance_to_yeast55:YeastEMM_F63  0.234          0.468          0.250       
    ## distance_to_yeast17:YeastEMM_F64  0.234          0.234          0.500       
    ## distance_to_yeast25:YeastEMM_F64  0.234          0.234          0.250       
    ## distance_to_yeast32:YeastEMM_F64  0.234          0.234          0.250       
    ## distance_to_yeast41:YeastEMM_F64  0.234          0.234          0.250       
    ## distance_to_yeast48:YeastEMM_F64  0.468          0.234          0.250       
    ## distance_to_yeast55:YeastEMM_F64  0.234          0.468          0.250       
    ## distance_to_yeast17:YeastEMM_F65  0.227          0.227          0.485       
    ## distance_to_yeast25:YeastEMM_F65  0.227          0.227          0.243       
    ## distance_to_yeast32:YeastEMM_F65  0.227          0.227          0.243       
    ## distance_to_yeast41:YeastEMM_F65  0.227          0.227          0.243       
    ## distance_to_yeast48:YeastEMM_F65  0.454          0.227          0.243       
    ## distance_to_yeast55:YeastEMM_F65  0.227          0.454          0.243       
    ## distance_to_yeast17:YeastEMM_F66  0.234          0.234          0.500       
    ## distance_to_yeast25:YeastEMM_F66  0.234          0.234          0.250       
    ## distance_to_yeast32:YeastEMM_F66  0.234          0.234          0.250       
    ## distance_to_yeast41:YeastEMM_F66  0.234          0.234          0.250       
    ## distance_to_yeast48:YeastEMM_F66  0.468          0.234          0.250       
    ## distance_to_yeast55:YeastEMM_F66  0.234          0.468          0.250       
    ## distance_to_yeast17:YeastEMM_F7   0.234          0.234          0.500       
    ## distance_to_yeast25:YeastEMM_F7   0.234          0.234          0.250       
    ## distance_to_yeast32:YeastEMM_F7   0.234          0.234          0.250       
    ## distance_to_yeast41:YeastEMM_F7   0.234          0.234          0.250       
    ## distance_to_yeast48:YeastEMM_F7   0.468          0.234          0.250       
    ## distance_to_yeast55:YeastEMM_F7   0.234          0.468          0.250       
    ## distance_to_yeast17:YeastEMM_F70  0.234          0.234          0.500       
    ## distance_to_yeast25:YeastEMM_F70  0.234          0.234          0.250       
    ## distance_to_yeast32:YeastEMM_F70  0.234          0.234          0.250       
    ## distance_to_yeast41:YeastEMM_F70  0.234          0.234          0.250       
    ## distance_to_yeast48:YeastEMM_F70  0.468          0.234          0.250       
    ## distance_to_yeast55:YeastEMM_F70  0.234          0.468          0.250       
    ## distance_to_yeast17:YeastEMM_F89  0.219          0.219          0.468       
    ## distance_to_yeast25:YeastEMM_F89  0.219          0.219          0.234       
    ## distance_to_yeast32:YeastEMM_F89  0.219          0.219          0.234       
    ## distance_to_yeast41:YeastEMM_F89  0.219          0.219          0.234       
    ## distance_to_yeast48:YeastEMM_F89  0.438          0.219          0.234       
    ## distance_to_yeast55:YeastEMM_F89  0.219          0.437          0.234       
    ## distance_to_yeast17:YeastSP_F14   0.234          0.234          0.500       
    ## distance_to_yeast25:YeastSP_F14   0.234          0.234          0.250       
    ## distance_to_yeast32:YeastSP_F14   0.234          0.234          0.250       
    ## distance_to_yeast41:YeastSP_F14   0.234          0.234          0.250       
    ## distance_to_yeast48:YeastSP_F14   0.468          0.234          0.250       
    ## distance_to_yeast55:YeastSP_F14   0.234          0.468          0.250       
    ## distance_to_yeast17:YeastZAN_F3   0.227          0.227          0.485       
    ## distance_to_yeast25:YeastZAN_F3   0.227          0.227          0.243       
    ## distance_to_yeast32:YeastZAN_F3   0.227          0.227          0.243       
    ## distance_to_yeast41:YeastZAN_F3   0.227          0.227          0.243       
    ## distance_to_yeast48:YeastZAN_F3   0.454          0.227          0.243       
    ## distance_to_yeast55:YeastZAN_F3   0.227          0.454          0.243       
    ## distance_to_yeast17:YeastZAN_F4   0.234          0.234          0.500       
    ## distance_to_yeast25:YeastZAN_F4   0.234          0.234          0.250       
    ## distance_to_yeast32:YeastZAN_F4   0.234          0.234          0.250       
    ## distance_to_yeast41:YeastZAN_F4   0.234          0.234          0.250       
    ## distance_to_yeast48:YeastZAN_F4   0.468          0.234          0.250       
    ## distance_to_yeast55:YeastZAN_F4   0.234          0.468          0.250       
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
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
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
    ## distance_to_yeast17:YeastEMM_F65  0.243         0.243         0.243       
    ## distance_to_yeast25:YeastEMM_F65  0.485         0.243         0.243       
    ## distance_to_yeast32:YeastEMM_F65  0.243         0.485         0.243       
    ## distance_to_yeast41:YeastEMM_F65  0.243         0.243         0.485       
    ## distance_to_yeast48:YeastEMM_F65  0.243         0.243         0.243       
    ## distance_to_yeast55:YeastEMM_F65  0.243         0.243         0.243       
    ## distance_to_yeast17:YeastEMM_F66  0.250         0.250         0.250       
    ## distance_to_yeast25:YeastEMM_F66  0.500         0.250         0.250       
    ## distance_to_yeast32:YeastEMM_F66  0.250         0.500         0.250       
    ## distance_to_yeast41:YeastEMM_F66  0.250         0.250         0.500       
    ## distance_to_yeast48:YeastEMM_F66  0.250         0.250         0.250       
    ## distance_to_yeast55:YeastEMM_F66  0.250         0.250         0.250       
    ## distance_to_yeast17:YeastEMM_F7   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastEMM_F7   0.500         0.250         0.250       
    ## distance_to_yeast32:YeastEMM_F7   0.250         0.500         0.250       
    ## distance_to_yeast41:YeastEMM_F7   0.250         0.250         0.500       
    ## distance_to_yeast48:YeastEMM_F7   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastEMM_F7   0.250         0.250         0.250       
    ## distance_to_yeast17:YeastEMM_F70  0.250         0.250         0.250       
    ## distance_to_yeast25:YeastEMM_F70  0.500         0.250         0.250       
    ## distance_to_yeast32:YeastEMM_F70  0.250         0.500         0.250       
    ## distance_to_yeast41:YeastEMM_F70  0.250         0.250         0.500       
    ## distance_to_yeast48:YeastEMM_F70  0.250         0.250         0.250       
    ## distance_to_yeast55:YeastEMM_F70  0.250         0.250         0.250       
    ## distance_to_yeast17:YeastEMM_F89  0.234         0.234         0.234       
    ## distance_to_yeast25:YeastEMM_F89  0.468         0.234         0.234       
    ## distance_to_yeast32:YeastEMM_F89  0.234         0.468         0.234       
    ## distance_to_yeast41:YeastEMM_F89  0.234         0.234         0.468       
    ## distance_to_yeast48:YeastEMM_F89  0.234         0.234         0.234       
    ## distance_to_yeast55:YeastEMM_F89  0.234         0.234         0.234       
    ## distance_to_yeast17:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastSP_F14   0.500         0.250         0.250       
    ## distance_to_yeast32:YeastSP_F14   0.250         0.500         0.250       
    ## distance_to_yeast41:YeastSP_F14   0.250         0.250         0.500       
    ## distance_to_yeast48:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast17:YeastZAN_F3   0.243         0.243         0.243       
    ## distance_to_yeast25:YeastZAN_F3   0.485         0.243         0.243       
    ## distance_to_yeast32:YeastZAN_F3   0.243         0.485         0.243       
    ## distance_to_yeast41:YeastZAN_F3   0.243         0.243         0.485       
    ## distance_to_yeast48:YeastZAN_F3   0.243         0.243         0.243       
    ## distance_to_yeast55:YeastZAN_F3   0.243         0.243         0.243       
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
    ## YeastEMM_F48                                                               
    ## YeastEMM_F49                                                               
    ## YeastEMM_F5                                                                
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
    ## distance_to_yeast17:YeastEMM_F65  0.243         0.243         0.485        
    ## distance_to_yeast25:YeastEMM_F65  0.243         0.243         0.243        
    ## distance_to_yeast32:YeastEMM_F65  0.243         0.243         0.243        
    ## distance_to_yeast41:YeastEMM_F65  0.243         0.243         0.243        
    ## distance_to_yeast48:YeastEMM_F65  0.485         0.243         0.243        
    ## distance_to_yeast55:YeastEMM_F65  0.243         0.485         0.243        
    ## distance_to_yeast17:YeastEMM_F66  0.250         0.250         0.500        
    ## distance_to_yeast25:YeastEMM_F66  0.250         0.250         0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250         0.250         0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.250         0.250         0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.500         0.250         0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250         0.500         0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.250         0.250         0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.250         0.250         0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250         0.250         0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.250         0.250         0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.500         0.250         0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250         0.500         0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.250         0.250         0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.250         0.250         0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250         0.250         0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.250         0.250         0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.500         0.250         0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250         0.500         0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.234         0.234         0.468        
    ## distance_to_yeast25:YeastEMM_F89  0.234         0.234         0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234         0.234         0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234         0.234         0.234        
    ## distance_to_yeast48:YeastEMM_F89  0.468         0.234         0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234         0.468         0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250         0.250         0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250         0.250         0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250         0.250         0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250         0.250         0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500         0.250         0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250         0.500         0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243         0.243         0.485        
    ## distance_to_yeast25:YeastZAN_F3   0.243         0.243         0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243         0.243         0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243         0.243         0.243        
    ## distance_to_yeast48:YeastZAN_F3   0.485         0.243         0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243         0.485         0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F65  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F65  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F65  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F65  0.243          0.243          0.243        
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.243          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F65  0.243          0.243          0.485        
    ## distance_to_yeast25:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast41:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast48:YeastEMM_F65  0.485          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F65  0.243          0.485          0.243        
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast25:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast48:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast25:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast48:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.485          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F65  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F65  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F65  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F65  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F65  0.243          0.243          0.243        
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.243          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F65  0.243          0.243                       
    ## distance_to_yeast25:YeastEMM_F65  0.243          0.243          0.500        
    ## distance_to_yeast32:YeastEMM_F65  0.243          0.243          0.500        
    ## distance_to_yeast41:YeastEMM_F65  0.243          0.243          0.500        
    ## distance_to_yeast48:YeastEMM_F65  0.485          0.243          0.500        
    ## distance_to_yeast55:YeastEMM_F65  0.243          0.485          0.500        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.485        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.243        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.243        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.243        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.250          0.243        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.500          0.243        
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.454        
    ## distance_to_yeast25:YeastEMM_F89  0.234          0.234          0.227        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.234          0.227        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.227        
    ## distance_to_yeast48:YeastEMM_F89  0.468          0.234          0.227        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.468          0.227        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.485        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.243        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.243        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.243        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.243        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.243        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.471        
    ## distance_to_yeast25:YeastZAN_F3   0.243          0.243          0.235        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.243          0.235        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.235        
    ## distance_to_yeast48:YeastZAN_F3   0.485          0.243          0.235        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.485          0.235        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.485        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.243        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.243        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.243        
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.243        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F66  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F66  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F66  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F66  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F7   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F7   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F7   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F7   0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast25:YeastEMM_F70  0.485          0.243          0.243        
    ## distance_to_yeast32:YeastEMM_F70  0.243          0.485          0.243        
    ## distance_to_yeast41:YeastEMM_F70  0.243          0.243          0.485        
    ## distance_to_yeast48:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast55:YeastEMM_F70  0.243          0.243          0.243        
    ## distance_to_yeast17:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast25:YeastEMM_F89  0.454          0.227          0.227        
    ## distance_to_yeast32:YeastEMM_F89  0.227          0.454          0.227        
    ## distance_to_yeast41:YeastEMM_F89  0.227          0.227          0.454        
    ## distance_to_yeast48:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast55:YeastEMM_F89  0.227          0.227          0.227        
    ## distance_to_yeast17:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastSP_F14   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastSP_F14   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastSP_F14   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastSP_F14   0.243          0.243          0.243        
    ## distance_to_yeast17:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast25:YeastZAN_F3   0.471          0.235          0.235        
    ## distance_to_yeast32:YeastZAN_F3   0.235          0.471          0.235        
    ## distance_to_yeast41:YeastZAN_F3   0.235          0.235          0.471        
    ## distance_to_yeast48:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast55:YeastZAN_F3   0.235          0.235          0.235        
    ## distance_to_yeast17:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F4   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F4   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F4   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F4   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F4   0.243          0.243          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F66  0.243          0.243                       
    ## distance_to_yeast25:YeastEMM_F66  0.243          0.243          0.500        
    ## distance_to_yeast32:YeastEMM_F66  0.243          0.243          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.243          0.243          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.485          0.243          0.500        
    ## distance_to_yeast55:YeastEMM_F66  0.243          0.485          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.243          0.243          0.500        
    ## distance_to_yeast25:YeastEMM_F7   0.243          0.243          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.243          0.243          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.243          0.243          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.485          0.243          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.243          0.485          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.243          0.243          0.500        
    ## distance_to_yeast25:YeastEMM_F70  0.243          0.243          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.243          0.243          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.243          0.243          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.485          0.243          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.243          0.485          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.227          0.227          0.468        
    ## distance_to_yeast25:YeastEMM_F89  0.227          0.227          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.227          0.227          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.227          0.227          0.234        
    ## distance_to_yeast48:YeastEMM_F89  0.454          0.227          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.227          0.454          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.243          0.243          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.243          0.243          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.243          0.243          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.243          0.243          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.485          0.243          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.243          0.485          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.235          0.235          0.485        
    ## distance_to_yeast25:YeastZAN_F3   0.235          0.235          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.235          0.235          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.235          0.235          0.243        
    ## distance_to_yeast48:YeastZAN_F3   0.471          0.235          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.235          0.471          0.243        
    ## distance_to_yeast17:YeastZAN_F4   0.243          0.243          0.500        
    ## distance_to_yeast25:YeastZAN_F4   0.243          0.243          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.243          0.243          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.243          0.243          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.485          0.243          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.243          0.485          0.250        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.243          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast25:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast48:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast25:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast48:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.485          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.243          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast25:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast48:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast25:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast48:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.485          0.243        
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
    ## YeastEMM_F48                                                                 
    ## YeastEMM_F49                                                                 
    ## YeastEMM_F5                                                                  
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast25:YeastEMM_F89  0.468          0.234          0.234        
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.468          0.234        
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.468        
    ## distance_to_yeast48:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.234          0.234        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast25:YeastZAN_F3   0.485          0.243          0.243        
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.485          0.243        
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.485        
    ## distance_to_yeast48:YeastZAN_F3   0.243          0.243          0.243        
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.243          0.243        
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
    ## YeastEMM_F48                                                                
    ## YeastEMM_F49                                                                
    ## YeastEMM_F5                                                                 
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
    ## distance_to_yeast17:YeastEMM_F89  0.234          0.234                      
    ## distance_to_yeast25:YeastEMM_F89  0.234          0.234          0.500       
    ## distance_to_yeast32:YeastEMM_F89  0.234          0.234          0.500       
    ## distance_to_yeast41:YeastEMM_F89  0.234          0.234          0.500       
    ## distance_to_yeast48:YeastEMM_F89  0.468          0.234          0.500       
    ## distance_to_yeast55:YeastEMM_F89  0.234          0.468          0.500       
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.468       
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.234       
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.234       
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.234       
    ## distance_to_yeast48:YeastSP_F14   0.500          0.250          0.234       
    ## distance_to_yeast55:YeastSP_F14   0.250          0.500          0.234       
    ## distance_to_yeast17:YeastZAN_F3   0.243          0.243          0.454       
    ## distance_to_yeast25:YeastZAN_F3   0.243          0.243          0.227       
    ## distance_to_yeast32:YeastZAN_F3   0.243          0.243          0.227       
    ## distance_to_yeast41:YeastZAN_F3   0.243          0.243          0.227       
    ## distance_to_yeast48:YeastZAN_F3   0.485          0.243          0.227       
    ## distance_to_yeast55:YeastZAN_F3   0.243          0.485          0.227       
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.468       
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.234       
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.234       
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.234       
    ## distance_to_yeast48:YeastZAN_F4   0.500          0.250          0.234       
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.500          0.234       
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
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
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
    ## distance_to_yeast17:YeastSP_F14   0.234         0.234         0.234       
    ## distance_to_yeast25:YeastSP_F14   0.468         0.234         0.234       
    ## distance_to_yeast32:YeastSP_F14   0.234         0.468         0.234       
    ## distance_to_yeast41:YeastSP_F14   0.234         0.234         0.468       
    ## distance_to_yeast48:YeastSP_F14   0.234         0.234         0.234       
    ## distance_to_yeast55:YeastSP_F14   0.234         0.234         0.234       
    ## distance_to_yeast17:YeastZAN_F3   0.227         0.227         0.227       
    ## distance_to_yeast25:YeastZAN_F3   0.454         0.227         0.227       
    ## distance_to_yeast32:YeastZAN_F3   0.227         0.454         0.227       
    ## distance_to_yeast41:YeastZAN_F3   0.227         0.227         0.454       
    ## distance_to_yeast48:YeastZAN_F3   0.227         0.227         0.227       
    ## distance_to_yeast55:YeastZAN_F3   0.227         0.227         0.227       
    ## distance_to_yeast17:YeastZAN_F4   0.234         0.234         0.234       
    ## distance_to_yeast25:YeastZAN_F4   0.468         0.234         0.234       
    ## distance_to_yeast32:YeastZAN_F4   0.234         0.468         0.234       
    ## distance_to_yeast41:YeastZAN_F4   0.234         0.234         0.468       
    ## distance_to_yeast48:YeastZAN_F4   0.234         0.234         0.234       
    ## distance_to_yeast55:YeastZAN_F4   0.234         0.234         0.234       
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
    ## YeastEMM_F48                                                                  
    ## YeastEMM_F49                                                                  
    ## YeastEMM_F5                                                                   
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
    ## distance_to_yeast17:YeastSP_F14   0.234         0.234                         
    ## distance_to_yeast25:YeastSP_F14   0.234         0.234         0.500           
    ## distance_to_yeast32:YeastSP_F14   0.234         0.234         0.500    0.500  
    ## distance_to_yeast41:YeastSP_F14   0.234         0.234         0.500    0.500  
    ## distance_to_yeast48:YeastSP_F14   0.468         0.234         0.500    0.500  
    ## distance_to_yeast55:YeastSP_F14   0.234         0.468         0.500    0.500  
    ## distance_to_yeast17:YeastZAN_F3   0.227         0.227         0.485    0.243  
    ## distance_to_yeast25:YeastZAN_F3   0.227         0.227         0.243    0.485  
    ## distance_to_yeast32:YeastZAN_F3   0.227         0.227         0.243    0.243  
    ## distance_to_yeast41:YeastZAN_F3   0.227         0.227         0.243    0.243  
    ## distance_to_yeast48:YeastZAN_F3   0.454         0.227         0.243    0.243  
    ## distance_to_yeast55:YeastZAN_F3   0.227         0.454         0.243    0.243  
    ## distance_to_yeast17:YeastZAN_F4   0.234         0.234         0.500    0.250  
    ## distance_to_yeast25:YeastZAN_F4   0.234         0.234         0.250    0.500  
    ## distance_to_yeast32:YeastZAN_F4   0.234         0.234         0.250    0.250  
    ## distance_to_yeast41:YeastZAN_F4   0.234         0.234         0.250    0.250  
    ## distance_to_yeast48:YeastZAN_F4   0.468         0.234         0.250    0.250  
    ## distance_to_yeast55:YeastZAN_F4   0.234         0.468         0.250    0.250  
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
    ## YeastEMM_F48                                                        
    ## YeastEMM_F49                                                        
    ## YeastEMM_F5                                                         
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
    ## distance_to_yeast17:YeastZAN_F3   0.243    0.243    0.243    0.243  
    ## distance_to_yeast25:YeastZAN_F3   0.243    0.243    0.243    0.243  
    ## distance_to_yeast32:YeastZAN_F3   0.485    0.243    0.243    0.243  
    ## distance_to_yeast41:YeastZAN_F3   0.243    0.485    0.243    0.243  
    ## distance_to_yeast48:YeastZAN_F3   0.243    0.243    0.485    0.243  
    ## distance_to_yeast55:YeastZAN_F3   0.243    0.243    0.243    0.485  
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
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
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
    ## distance_to_yeast17:YeastZAN_F4   0.485         0.243         0.243       
    ## distance_to_yeast25:YeastZAN_F4   0.243         0.485         0.243       
    ## distance_to_yeast32:YeastZAN_F4   0.243         0.243         0.485       
    ## distance_to_yeast41:YeastZAN_F4   0.243         0.243         0.243       
    ## distance_to_yeast48:YeastZAN_F4   0.243         0.243         0.243       
    ## distance_to_yeast55:YeastZAN_F4   0.243         0.243         0.243       
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
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
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
    ## distance_to_yeast17:YeastZAN_F4   0.243         0.243         0.243       
    ## distance_to_yeast25:YeastZAN_F4   0.243         0.243         0.243       
    ## distance_to_yeast32:YeastZAN_F4   0.243         0.243         0.243       
    ## distance_to_yeast41:YeastZAN_F4   0.485         0.243         0.243       
    ## distance_to_yeast48:YeastZAN_F4   0.243         0.485         0.243       
    ## distance_to_yeast55:YeastZAN_F4   0.243         0.243         0.485       
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
    ## YeastEMM_F48                                                              
    ## YeastEMM_F49                                                              
    ## YeastEMM_F5                                                               
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
    ## YeastEMM_F48                                                
    ## YeastEMM_F49                                                
    ## YeastEMM_F5                                                 
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
    ## -3.5816149 -0.5701631 -0.1114478  0.5024321  4.3978943 
    ## 
    ## Number of Observations: 952
    ## Number of Groups: 3

``` r
#Anova(resultsB5)
anova(resultsB30)
```

    ##                         numDF denDF   F-value p-value
    ## (Intercept)                 1   836 1663.2205  <.0001
    ## DAI                         2   836   56.6615  <.0001
    ## distance_to_yeast           6   836    3.8624   8e-04
    ## Yeast                      15   836    7.4512  <.0001
    ## distance_to_yeast:Yeast    90   836    1.8717  <.0001

### Loop for running analysis for each day separately for EMM_B30

Days after inoculation (DAI) as a factor is always significantly
impacting the growth. Also, our plot will represent data for Day 6 thus
we want relevant stats and comparison on Day 6 to present in the plot.
So, loop was made for each Day data and removing DAI from the model and
keeping rest of it present.

``` r
# Unique days
day <- unique(B30.no.1st.data$DAI)
# Initialize empty lists to store results
B30_Day <- list()
meansB30_Day <- list()
resultsB30_Day <- list()
model_summary_B30 <- list()
anova_table_B30 <- list()
Results_lsmeansB30 <- list()
filtered_comparisonsB30 <- list()

for (i in seq_along(day)) {
    # Filter data for the current day
  B30_Day[[i]] <- filter(B30.no.1st.data, DAI == day[i])
    # Fit the mixed effects model
  resultsB30_Day[[i]] <- lme(increase ~ distance_to_yeast * Yeast,
                            data = B30_Day[[i]],
                            random = ~1 | Replication,
                            na.action = na.omit)
    model_summary_B30[[i]] <- summary(resultsB30_Day[[i]])
  anova_table_B30[[i]] <- anova(resultsB30_Day[[i]])
   # Estimated marginal means
  meansB30_Day[[i]] <- emmeans(resultsB30_Day[[i]], ~ Yeast | distance_to_yeast)
    # Compact letter display with Tukey adjustment
  Results_lsmeansB30[[i]] <- multcomp::cld(meansB30_Day[[i]], alpha = 0.05,
                                          Letters = letters, reversed = TRUE,
                                          details = TRUE)
   # Print lsmeans results
  print(Results_lsmeansB30[[i]])
  
  # Filter comparisons involving "Control"
  filtered_comparisonsB30[[i]] <- Results_lsmeansB30[[i]]$comparisons %>%
    filter(grepl("Control", contrast))
   # Print filtered contrasts
  print(filtered_comparisonsB30[[i]])
}
```

    ## $emmeans
    ## distance_to_yeast = 11:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  ZAN_F3   1.01333 0.278  2   -0.185    2.211  a    
    ##  EMM_F48  0.63000 0.278  2   -0.568    1.828  a    
    ##  EMM_F63  0.59667 0.278  2   -0.601    1.795  a    
    ##  EMM_F5   0.51667 0.278  2   -0.681    1.715  a    
    ##  EMM_F66  0.51333 0.278  2   -0.685    1.711  a    
    ##  SP_F14   0.43000 0.278  2   -0.768    1.628  a    
    ##  EMM_F3   0.40333 0.278  2   -0.795    1.601  a    
    ##  Control  0.31000 0.278  2   -0.888    1.508  a    
    ##  EMM_F64  0.27333 0.278  2   -0.925    1.471  a    
    ##  EMM_F49  0.26333 0.278  2   -0.935    1.461  a    
    ##  EMM_F65  0.25333 0.278  2   -0.945    1.451  a    
    ##  EMM_F70  0.23000 0.278  2   -0.968    1.428  a    
    ##  EMM_F7   0.15333 0.278  2   -1.045    1.351  a    
    ##  EMM_F89  0.14000 0.278  2   -1.058    1.338  a    
    ##  EMM_F34  0.00000 0.278  2   -1.198    1.198  a    
    ##  ZAN_F4  -0.20667 0.278  2   -1.405    0.991  a    
    ## 
    ## distance_to_yeast = 17:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  EMM_F48  0.61667 0.278  2   -0.581    1.815  a    
    ##  EMM_F49  0.58667 0.278  2   -0.611    1.785  a    
    ##  EMM_F64  0.49000 0.278  2   -0.708    1.688  a    
    ##  EMM_F63  0.44333 0.278  2   -0.755    1.641  a    
    ##  EMM_F70  0.44333 0.278  2   -0.755    1.641  a    
    ##  EMM_F5   0.43333 0.278  2   -0.765    1.631  a    
    ##  EMM_F65  0.42667 0.278  2   -0.771    1.625  a    
    ##  EMM_F89  0.40000 0.278  2   -0.798    1.598  a    
    ##  EMM_F3   0.37333 0.278  2   -0.825    1.571  a    
    ##  EMM_F34  0.32333 0.278  2   -0.875    1.521  a    
    ##  EMM_F66  0.29333 0.278  2   -0.905    1.491  a    
    ##  ZAN_F3   0.28667 0.278  2   -0.911    1.485  a    
    ##  EMM_F7   0.28000 0.278  2   -0.918    1.478  a    
    ##  SP_F14   0.12667 0.278  2   -1.071    1.325  a    
    ##  Control  0.05667 0.278  2   -1.141    1.255  a    
    ##  ZAN_F4  -0.24667 0.278  2   -1.445    0.951  a    
    ## 
    ## distance_to_yeast = 25:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  EMM_F48  0.79000 0.278  2   -0.408    1.988  a    
    ##  EMM_F5   0.67667 0.278  2   -0.521    1.875  a    
    ##  EMM_F65  0.58000 0.278  2   -0.618    1.778  a    
    ##  EMM_F64  0.57667 0.278  2   -0.621    1.775  a    
    ##  EMM_F70  0.56000 0.278  2   -0.638    1.758  a    
    ##  ZAN_F3   0.50000 0.278  2   -0.698    1.698  a    
    ##  EMM_F7   0.47333 0.278  2   -0.725    1.671  a    
    ##  EMM_F49  0.43333 0.278  2   -0.765    1.631  a    
    ##  EMM_F63  0.43000 0.278  2   -0.768    1.628  a    
    ##  EMM_F3   0.34667 0.278  2   -0.851    1.545  a    
    ##  EMM_F89  0.32333 0.278  2   -0.875    1.521  a    
    ##  SP_F14   0.26667 0.278  2   -0.931    1.465  a    
    ##  EMM_F66  0.13667 0.278  2   -1.061    1.335  a    
    ##  Control  0.10333 0.278  2   -1.095    1.301  a    
    ##  EMM_F34 -0.07333 0.278  2   -1.271    1.125  a    
    ##  ZAN_F4  -0.12333 0.278  2   -1.321    1.075  a    
    ## 
    ## distance_to_yeast = 32:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  ZAN_F3   0.80000 0.278  2   -0.398    1.998  a    
    ##  EMM_F63  0.73000 0.278  2   -0.468    1.928  a    
    ##  EMM_F70  0.68667 0.278  2   -0.511    1.885  a    
    ##  EMM_F48  0.58667 0.278  2   -0.611    1.785  a    
    ##  EMM_F64  0.56667 0.278  2   -0.631    1.765  a    
    ##  EMM_F89  0.53000 0.278  2   -0.668    1.728  a    
    ##  EMM_F66  0.40333 0.278  2   -0.795    1.601  a    
    ##  ZAN_F4   0.39667 0.278  2   -0.801    1.595  a    
    ##  EMM_F49  0.37333 0.278  2   -0.825    1.571  a    
    ##  EMM_F5   0.29333 0.278  2   -0.905    1.491  a    
    ##  EMM_F7   0.27333 0.278  2   -0.925    1.471  a    
    ##  SP_F14   0.19333 0.278  2   -1.005    1.391  a    
    ##  Control  0.10333 0.278  2   -1.095    1.301  a    
    ##  EMM_F3   0.04000 0.278  2   -1.158    1.238  a    
    ##  EMM_F34 -0.00667 0.278  2   -1.205    1.191  a    
    ##  EMM_F65 -0.05333 0.278  2   -1.251    1.145  a    
    ## 
    ## distance_to_yeast = 41:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  SP_F14   0.84000 0.278  2   -0.358    2.038  a    
    ##  ZAN_F3   0.78333 0.278  2   -0.415    1.981  a    
    ##  EMM_F89  0.64667 0.278  2   -0.551    1.845  a    
    ##  EMM_F7   0.50333 0.278  2   -0.695    1.701  a    
    ##  EMM_F64  0.48667 0.278  2   -0.711    1.685  a    
    ##  EMM_F48  0.47667 0.278  2   -0.721    1.675  a    
    ##  EMM_F70  0.47667 0.278  2   -0.721    1.675  a    
    ##  EMM_F66  0.40000 0.278  2   -0.798    1.598  a    
    ##  EMM_F63  0.39667 0.278  2   -0.801    1.595  a    
    ##  EMM_F65  0.38667 0.278  2   -0.811    1.585  a    
    ##  EMM_F3   0.26667 0.278  2   -0.931    1.465  a    
    ##  EMM_F49  0.26333 0.278  2   -0.935    1.461  a    
    ##  Control  0.21667 0.278  2   -0.981    1.415  a    
    ##  EMM_F5   0.16333 0.278  2   -1.035    1.361  a    
    ##  ZAN_F4   0.09667 0.278  2   -1.101    1.295  a    
    ##  EMM_F34  0.04000 0.278  2   -1.158    1.238  a    
    ## 
    ## distance_to_yeast = 48:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  ZAN_F3   1.37000 0.278  2    0.172    2.568  a    
    ##  EMM_F89  0.95000 0.278  2   -0.248    2.148  a    
    ##  EMM_F49  0.77000 0.278  2   -0.428    1.968  a    
    ##  EMM_F65  0.75333 0.278  2   -0.445    1.951  a    
    ##  EMM_F70  0.70000 0.278  2   -0.498    1.898  a    
    ##  EMM_F66  0.60667 0.278  2   -0.591    1.805  a    
    ##  EMM_F64  0.57000 0.278  2   -0.628    1.768  a    
    ##  SP_F14   0.55000 0.278  2   -0.648    1.748  a    
    ##  EMM_F7   0.50333 0.278  2   -0.695    1.701  a    
    ##  EMM_F63  0.45667 0.278  2   -0.741    1.655  a    
    ##  Control  0.43000 0.278  2   -0.768    1.628  a    
    ##  ZAN_F4   0.28000 0.278  2   -0.918    1.478  a    
    ##  EMM_F48  0.25333 0.278  2   -0.945    1.451  a    
    ##  EMM_F5   0.23000 0.278  2   -0.968    1.428  a    
    ##  EMM_F3   0.17000 0.278  2   -1.028    1.368  a    
    ##  EMM_F34  0.03000 0.278  2   -1.168    1.228  a    
    ## 
    ## distance_to_yeast = 55:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  EMM_F64  1.38333 0.278  2    0.185    2.581  a    
    ##  SP_F14   0.68333 0.278  2   -0.515    1.881  ab   
    ##  EMM_F89  0.57667 0.278  2   -0.621    1.775  ab   
    ##  EMM_F7   0.56667 0.278  2   -0.631    1.765  ab   
    ##  Control  0.55000 0.278  2   -0.648    1.748  ab   
    ##  EMM_F48  0.50000 0.278  2   -0.698    1.698  ab   
    ##  EMM_F65  0.48333 0.278  2   -0.715    1.681  ab   
    ##  EMM_F63  0.45333 0.278  2   -0.745    1.651  ab   
    ##  ZAN_F3   0.41667 0.278  2   -0.781    1.615  ab   
    ##  ZAN_F4   0.40000 0.278  2   -0.798    1.598  ab   
    ##  EMM_F70  0.38000 0.278  2   -0.818    1.578  ab   
    ##  EMM_F5   0.32667 0.278  2   -0.871    1.525  ab   
    ##  EMM_F49  0.28667 0.278  2   -0.911    1.485  ab   
    ##  EMM_F66  0.22667 0.278  2   -0.971    1.425  ab   
    ##  EMM_F3   0.21333 0.278  2   -0.985    1.411  ab   
    ##  EMM_F34 -0.10000 0.278  2   -1.298    1.098   b   
    ## 
    ## Degrees-of-freedom method: containment 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F34 - ZAN_F4   0.20667 0.394 222   0.525  1.0000
    ##  EMM_F89 - ZAN_F4   0.34667 0.394 222   0.880  1.0000
    ##  EMM_F89 - EMM_F34  0.14000 0.394 222   0.356  1.0000
    ##  EMM_F7 - ZAN_F4    0.36000 0.394 222   0.914  0.9999
    ##  EMM_F7 - EMM_F34   0.15333 0.394 222   0.389  1.0000
    ##  EMM_F7 - EMM_F89   0.01333 0.394 222   0.034  1.0000
    ##  EMM_F70 - ZAN_F4   0.43667 0.394 222   1.109  0.9992
    ##  EMM_F70 - EMM_F34  0.23000 0.394 222   0.584  1.0000
    ##  EMM_F70 - EMM_F89  0.09000 0.394 222   0.229  1.0000
    ##  EMM_F70 - EMM_F7   0.07667 0.394 222   0.195  1.0000
    ##  EMM_F65 - ZAN_F4   0.46000 0.394 222   1.168  0.9986
    ##  EMM_F65 - EMM_F34  0.25333 0.394 222   0.643  1.0000
    ##  EMM_F65 - EMM_F89  0.11333 0.394 222   0.288  1.0000
    ##  EMM_F65 - EMM_F7   0.10000 0.394 222   0.254  1.0000
    ##  EMM_F65 - EMM_F70  0.02333 0.394 222   0.059  1.0000
    ##  EMM_F49 - ZAN_F4   0.47000 0.394 222   1.194  0.9983
    ##  EMM_F49 - EMM_F34  0.26333 0.394 222   0.669  1.0000
    ##  EMM_F49 - EMM_F89  0.12333 0.394 222   0.313  1.0000
    ##  EMM_F49 - EMM_F7   0.11000 0.394 222   0.279  1.0000
    ##  EMM_F49 - EMM_F70  0.03333 0.394 222   0.085  1.0000
    ##  EMM_F49 - EMM_F65  0.01000 0.394 222   0.025  1.0000
    ##  EMM_F64 - ZAN_F4   0.48000 0.394 222   1.219  0.9978
    ##  EMM_F64 - EMM_F34  0.27333 0.394 222   0.694  1.0000
    ##  EMM_F64 - EMM_F89  0.13333 0.394 222   0.339  1.0000
    ##  EMM_F64 - EMM_F7   0.12000 0.394 222   0.305  1.0000
    ##  EMM_F64 - EMM_F70  0.04333 0.394 222   0.110  1.0000
    ##  EMM_F64 - EMM_F65  0.02000 0.394 222   0.051  1.0000
    ##  EMM_F64 - EMM_F49  0.01000 0.394 222   0.025  1.0000
    ##  Control - ZAN_F4   0.51667 0.394 222   1.312  0.9951
    ##  Control - EMM_F34  0.31000 0.394 222   0.787  1.0000
    ##  Control - EMM_F89  0.17000 0.394 222   0.432  1.0000
    ##  Control - EMM_F7   0.15667 0.394 222   0.398  1.0000
    ##  Control - EMM_F70  0.08000 0.394 222   0.203  1.0000
    ##  Control - EMM_F65  0.05667 0.394 222   0.144  1.0000
    ##  Control - EMM_F49  0.04667 0.394 222   0.119  1.0000
    ##  Control - EMM_F64  0.03667 0.394 222   0.093  1.0000
    ##  EMM_F3 - ZAN_F4    0.61000 0.394 222   1.549  0.9752
    ##  EMM_F3 - EMM_F34   0.40333 0.394 222   1.024  0.9997
    ##  EMM_F3 - EMM_F89   0.26333 0.394 222   0.669  1.0000
    ##  EMM_F3 - EMM_F7    0.25000 0.394 222   0.635  1.0000
    ##  EMM_F3 - EMM_F70   0.17333 0.394 222   0.440  1.0000
    ##  EMM_F3 - EMM_F65   0.15000 0.394 222   0.381  1.0000
    ##  EMM_F3 - EMM_F49   0.14000 0.394 222   0.356  1.0000
    ##  EMM_F3 - EMM_F64   0.13000 0.394 222   0.330  1.0000
    ##  EMM_F3 - Control   0.09333 0.394 222   0.237  1.0000
    ##  SP_F14 - ZAN_F4    0.63667 0.394 222   1.617  0.9638
    ##  SP_F14 - EMM_F34   0.43000 0.394 222   1.092  0.9994
    ##  SP_F14 - EMM_F89   0.29000 0.394 222   0.736  1.0000
    ##  SP_F14 - EMM_F7    0.27667 0.394 222   0.703  1.0000
    ##  SP_F14 - EMM_F70   0.20000 0.394 222   0.508  1.0000
    ##  SP_F14 - EMM_F65   0.17667 0.394 222   0.449  1.0000
    ##  SP_F14 - EMM_F49   0.16667 0.394 222   0.423  1.0000
    ##  SP_F14 - EMM_F64   0.15667 0.394 222   0.398  1.0000
    ##  SP_F14 - Control   0.12000 0.394 222   0.305  1.0000
    ##  SP_F14 - EMM_F3    0.02667 0.394 222   0.068  1.0000
    ##  EMM_F66 - ZAN_F4   0.72000 0.394 222   1.828  0.9035
    ##  EMM_F66 - EMM_F34  0.51333 0.394 222   1.304  0.9954
    ##  EMM_F66 - EMM_F89  0.37333 0.394 222   0.948  0.9999
    ##  EMM_F66 - EMM_F7   0.36000 0.394 222   0.914  0.9999
    ##  EMM_F66 - EMM_F70  0.28333 0.394 222   0.720  1.0000
    ##  EMM_F66 - EMM_F65  0.26000 0.394 222   0.660  1.0000
    ##  EMM_F66 - EMM_F49  0.25000 0.394 222   0.635  1.0000
    ##  EMM_F66 - EMM_F64  0.24000 0.394 222   0.609  1.0000
    ##  EMM_F66 - Control  0.20333 0.394 222   0.516  1.0000
    ##  EMM_F66 - EMM_F3   0.11000 0.394 222   0.279  1.0000
    ##  EMM_F66 - SP_F14   0.08333 0.394 222   0.212  1.0000
    ##  EMM_F5 - ZAN_F4    0.72333 0.394 222   1.837  0.9002
    ##  EMM_F5 - EMM_F34   0.51667 0.394 222   1.312  0.9951
    ##  EMM_F5 - EMM_F89   0.37667 0.394 222   0.957  0.9999
    ##  EMM_F5 - EMM_F7    0.36333 0.394 222   0.923  0.9999
    ##  EMM_F5 - EMM_F70   0.28667 0.394 222   0.728  1.0000
    ##  EMM_F5 - EMM_F65   0.26333 0.394 222   0.669  1.0000
    ##  EMM_F5 - EMM_F49   0.25333 0.394 222   0.643  1.0000
    ##  EMM_F5 - EMM_F64   0.24333 0.394 222   0.618  1.0000
    ##  EMM_F5 - Control   0.20667 0.394 222   0.525  1.0000
    ##  EMM_F5 - EMM_F3    0.11333 0.394 222   0.288  1.0000
    ##  EMM_F5 - SP_F14    0.08667 0.394 222   0.220  1.0000
    ##  EMM_F5 - EMM_F66   0.00333 0.394 222   0.008  1.0000
    ##  EMM_F63 - ZAN_F4   0.80333 0.394 222   2.040  0.7999
    ##  EMM_F63 - EMM_F34  0.59667 0.394 222   1.515  0.9797
    ##  EMM_F63 - EMM_F89  0.45667 0.394 222   1.160  0.9987
    ##  EMM_F63 - EMM_F7   0.44333 0.394 222   1.126  0.9991
    ##  EMM_F63 - EMM_F70  0.36667 0.394 222   0.931  0.9999
    ##  EMM_F63 - EMM_F65  0.34333 0.394 222   0.872  1.0000
    ##  EMM_F63 - EMM_F49  0.33333 0.394 222   0.846  1.0000
    ##  EMM_F63 - EMM_F64  0.32333 0.394 222   0.821  1.0000
    ##  EMM_F63 - Control  0.28667 0.394 222   0.728  1.0000
    ##  EMM_F63 - EMM_F3   0.19333 0.394 222   0.491  1.0000
    ##  EMM_F63 - SP_F14   0.16667 0.394 222   0.423  1.0000
    ##  EMM_F63 - EMM_F66  0.08333 0.394 222   0.212  1.0000
    ##  EMM_F63 - EMM_F5   0.08000 0.394 222   0.203  1.0000
    ##  EMM_F48 - ZAN_F4   0.83667 0.394 222   2.125  0.7473
    ##  EMM_F48 - EMM_F34  0.63000 0.394 222   1.600  0.9670
    ##  EMM_F48 - EMM_F89  0.49000 0.394 222   1.244  0.9972
    ##  EMM_F48 - EMM_F7   0.47667 0.394 222   1.210  0.9980
    ##  EMM_F48 - EMM_F70  0.40000 0.394 222   1.016  0.9997
    ##  EMM_F48 - EMM_F65  0.37667 0.394 222   0.957  0.9999
    ##  EMM_F48 - EMM_F49  0.36667 0.394 222   0.931  0.9999
    ##  EMM_F48 - EMM_F64  0.35667 0.394 222   0.906  0.9999
    ##  EMM_F48 - Control  0.32000 0.394 222   0.813  1.0000
    ##  EMM_F48 - EMM_F3   0.22667 0.394 222   0.576  1.0000
    ##  EMM_F48 - SP_F14   0.20000 0.394 222   0.508  1.0000
    ##  EMM_F48 - EMM_F66  0.11667 0.394 222   0.296  1.0000
    ##  EMM_F48 - EMM_F5   0.11333 0.394 222   0.288  1.0000
    ##  EMM_F48 - EMM_F63  0.03333 0.394 222   0.085  1.0000
    ##  ZAN_F3 - ZAN_F4    1.22000 0.394 222   3.098  0.1396
    ##  ZAN_F3 - EMM_F34   1.01333 0.394 222   2.573  0.4227
    ##  ZAN_F3 - EMM_F89   0.87333 0.394 222   2.218  0.6837
    ##  ZAN_F3 - EMM_F7    0.86000 0.394 222   2.184  0.7074
    ##  ZAN_F3 - EMM_F70   0.78333 0.394 222   1.989  0.8287
    ##  ZAN_F3 - EMM_F65   0.76000 0.394 222   1.930  0.8592
    ##  ZAN_F3 - EMM_F49   0.75000 0.394 222   1.905  0.8712
    ##  ZAN_F3 - EMM_F64   0.74000 0.394 222   1.879  0.8826
    ##  ZAN_F3 - Control   0.70333 0.394 222   1.786  0.9189
    ##  ZAN_F3 - EMM_F3    0.61000 0.394 222   1.549  0.9752
    ##  ZAN_F3 - SP_F14    0.58333 0.394 222   1.481  0.9836
    ##  ZAN_F3 - EMM_F66   0.50000 0.394 222   1.270  0.9966
    ##  ZAN_F3 - EMM_F5    0.49667 0.394 222   1.261  0.9968
    ##  ZAN_F3 - EMM_F63   0.41667 0.394 222   1.058  0.9996
    ##  ZAN_F3 - EMM_F48   0.38333 0.394 222   0.973  0.9998
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F4   0.30333 0.394 222   0.770  1.0000
    ##  SP_F14 - ZAN_F4    0.37333 0.394 222   0.948  0.9999
    ##  SP_F14 - Control   0.07000 0.394 222   0.178  1.0000
    ##  EMM_F7 - ZAN_F4    0.52667 0.394 222   1.337  0.9940
    ##  EMM_F7 - Control   0.22333 0.394 222   0.567  1.0000
    ##  EMM_F7 - SP_F14    0.15333 0.394 222   0.389  1.0000
    ##  ZAN_F3 - ZAN_F4    0.53333 0.394 222   1.354  0.9932
    ##  ZAN_F3 - Control   0.23000 0.394 222   0.584  1.0000
    ##  ZAN_F3 - SP_F14    0.16000 0.394 222   0.406  1.0000
    ##  ZAN_F3 - EMM_F7    0.00667 0.394 222   0.017  1.0000
    ##  EMM_F66 - ZAN_F4   0.54000 0.394 222   1.371  0.9923
    ##  EMM_F66 - Control  0.23667 0.394 222   0.601  1.0000
    ##  EMM_F66 - SP_F14   0.16667 0.394 222   0.423  1.0000
    ##  EMM_F66 - EMM_F7   0.01333 0.394 222   0.034  1.0000
    ##  EMM_F66 - ZAN_F3   0.00667 0.394 222   0.017  1.0000
    ##  EMM_F34 - ZAN_F4   0.57000 0.394 222   1.448  0.9869
    ##  EMM_F34 - Control  0.26667 0.394 222   0.677  1.0000
    ##  EMM_F34 - SP_F14   0.19667 0.394 222   0.499  1.0000
    ##  EMM_F34 - EMM_F7   0.04333 0.394 222   0.110  1.0000
    ##  EMM_F34 - ZAN_F3   0.03667 0.394 222   0.093  1.0000
    ##  EMM_F34 - EMM_F66  0.03000 0.394 222   0.076  1.0000
    ##  EMM_F3 - ZAN_F4    0.62000 0.394 222   1.574  0.9713
    ##  EMM_F3 - Control   0.31667 0.394 222   0.804  1.0000
    ##  EMM_F3 - SP_F14    0.24667 0.394 222   0.626  1.0000
    ##  EMM_F3 - EMM_F7    0.09333 0.394 222   0.237  1.0000
    ##  EMM_F3 - ZAN_F3    0.08667 0.394 222   0.220  1.0000
    ##  EMM_F3 - EMM_F66   0.08000 0.394 222   0.203  1.0000
    ##  EMM_F3 - EMM_F34   0.05000 0.394 222   0.127  1.0000
    ##  EMM_F89 - ZAN_F4   0.64667 0.394 222   1.642  0.9586
    ##  EMM_F89 - Control  0.34333 0.394 222   0.872  1.0000
    ##  EMM_F89 - SP_F14   0.27333 0.394 222   0.694  1.0000
    ##  EMM_F89 - EMM_F7   0.12000 0.394 222   0.305  1.0000
    ##  EMM_F89 - ZAN_F3   0.11333 0.394 222   0.288  1.0000
    ##  EMM_F89 - EMM_F66  0.10667 0.394 222   0.271  1.0000
    ##  EMM_F89 - EMM_F34  0.07667 0.394 222   0.195  1.0000
    ##  EMM_F89 - EMM_F3   0.02667 0.394 222   0.068  1.0000
    ##  EMM_F65 - ZAN_F4   0.67333 0.394 222   1.710  0.9423
    ##  EMM_F65 - Control  0.37000 0.394 222   0.940  0.9999
    ##  EMM_F65 - SP_F14   0.30000 0.394 222   0.762  1.0000
    ##  EMM_F65 - EMM_F7   0.14667 0.394 222   0.372  1.0000
    ##  EMM_F65 - ZAN_F3   0.14000 0.394 222   0.356  1.0000
    ##  EMM_F65 - EMM_F66  0.13333 0.394 222   0.339  1.0000
    ##  EMM_F65 - EMM_F34  0.10333 0.394 222   0.262  1.0000
    ##  EMM_F65 - EMM_F3   0.05333 0.394 222   0.135  1.0000
    ##  EMM_F65 - EMM_F89  0.02667 0.394 222   0.068  1.0000
    ##  EMM_F5 - ZAN_F4    0.68000 0.394 222   1.727  0.9375
    ##  EMM_F5 - Control   0.37667 0.394 222   0.957  0.9999
    ##  EMM_F5 - SP_F14    0.30667 0.394 222   0.779  1.0000
    ##  EMM_F5 - EMM_F7    0.15333 0.394 222   0.389  1.0000
    ##  EMM_F5 - ZAN_F3    0.14667 0.394 222   0.372  1.0000
    ##  EMM_F5 - EMM_F66   0.14000 0.394 222   0.356  1.0000
    ##  EMM_F5 - EMM_F34   0.11000 0.394 222   0.279  1.0000
    ##  EMM_F5 - EMM_F3    0.06000 0.394 222   0.152  1.0000
    ##  EMM_F5 - EMM_F89   0.03333 0.394 222   0.085  1.0000
    ##  EMM_F5 - EMM_F65   0.00667 0.394 222   0.017  1.0000
    ##  EMM_F70 - ZAN_F4   0.69000 0.394 222   1.752  0.9299
    ##  EMM_F70 - Control  0.38667 0.394 222   0.982  0.9998
    ##  EMM_F70 - SP_F14   0.31667 0.394 222   0.804  1.0000
    ##  EMM_F70 - EMM_F7   0.16333 0.394 222   0.415  1.0000
    ##  EMM_F70 - ZAN_F3   0.15667 0.394 222   0.398  1.0000
    ##  EMM_F70 - EMM_F66  0.15000 0.394 222   0.381  1.0000
    ##  EMM_F70 - EMM_F34  0.12000 0.394 222   0.305  1.0000
    ##  EMM_F70 - EMM_F3   0.07000 0.394 222   0.178  1.0000
    ##  EMM_F70 - EMM_F89  0.04333 0.394 222   0.110  1.0000
    ##  EMM_F70 - EMM_F65  0.01667 0.394 222   0.042  1.0000
    ##  EMM_F70 - EMM_F5   0.01000 0.394 222   0.025  1.0000
    ##  EMM_F63 - ZAN_F4   0.69000 0.394 222   1.752  0.9299
    ##  EMM_F63 - Control  0.38667 0.394 222   0.982  0.9998
    ##  EMM_F63 - SP_F14   0.31667 0.394 222   0.804  1.0000
    ##  EMM_F63 - EMM_F7   0.16333 0.394 222   0.415  1.0000
    ##  EMM_F63 - ZAN_F3   0.15667 0.394 222   0.398  1.0000
    ##  EMM_F63 - EMM_F66  0.15000 0.394 222   0.381  1.0000
    ##  EMM_F63 - EMM_F34  0.12000 0.394 222   0.305  1.0000
    ##  EMM_F63 - EMM_F3   0.07000 0.394 222   0.178  1.0000
    ##  EMM_F63 - EMM_F89  0.04333 0.394 222   0.110  1.0000
    ##  EMM_F63 - EMM_F65  0.01667 0.394 222   0.042  1.0000
    ##  EMM_F63 - EMM_F5   0.01000 0.394 222   0.025  1.0000
    ##  EMM_F63 - EMM_F70  0.00000 0.394 222   0.000  1.0000
    ##  EMM_F64 - ZAN_F4   0.73667 0.394 222   1.871  0.8863
    ##  EMM_F64 - Control  0.43333 0.394 222   1.100  0.9993
    ##  EMM_F64 - SP_F14   0.36333 0.394 222   0.923  0.9999
    ##  EMM_F64 - EMM_F7   0.21000 0.394 222   0.533  1.0000
    ##  EMM_F64 - ZAN_F3   0.20333 0.394 222   0.516  1.0000
    ##  EMM_F64 - EMM_F66  0.19667 0.394 222   0.499  1.0000
    ##  EMM_F64 - EMM_F34  0.16667 0.394 222   0.423  1.0000
    ##  EMM_F64 - EMM_F3   0.11667 0.394 222   0.296  1.0000
    ##  EMM_F64 - EMM_F89  0.09000 0.394 222   0.229  1.0000
    ##  EMM_F64 - EMM_F65  0.06333 0.394 222   0.161  1.0000
    ##  EMM_F64 - EMM_F5   0.05667 0.394 222   0.144  1.0000
    ##  EMM_F64 - EMM_F70  0.04667 0.394 222   0.119  1.0000
    ##  EMM_F64 - EMM_F63  0.04667 0.394 222   0.119  1.0000
    ##  EMM_F49 - ZAN_F4   0.83333 0.394 222   2.116  0.7528
    ##  EMM_F49 - Control  0.53000 0.394 222   1.346  0.9936
    ##  EMM_F49 - SP_F14   0.46000 0.394 222   1.168  0.9986
    ##  EMM_F49 - EMM_F7   0.30667 0.394 222   0.779  1.0000
    ##  EMM_F49 - ZAN_F3   0.30000 0.394 222   0.762  1.0000
    ##  EMM_F49 - EMM_F66  0.29333 0.394 222   0.745  1.0000
    ##  EMM_F49 - EMM_F34  0.26333 0.394 222   0.669  1.0000
    ##  EMM_F49 - EMM_F3   0.21333 0.394 222   0.542  1.0000
    ##  EMM_F49 - EMM_F89  0.18667 0.394 222   0.474  1.0000
    ##  EMM_F49 - EMM_F65  0.16000 0.394 222   0.406  1.0000
    ##  EMM_F49 - EMM_F5   0.15333 0.394 222   0.389  1.0000
    ##  EMM_F49 - EMM_F70  0.14333 0.394 222   0.364  1.0000
    ##  EMM_F49 - EMM_F63  0.14333 0.394 222   0.364  1.0000
    ##  EMM_F49 - EMM_F64  0.09667 0.394 222   0.245  1.0000
    ##  EMM_F48 - ZAN_F4   0.86333 0.394 222   2.192  0.7015
    ##  EMM_F48 - Control  0.56000 0.394 222   1.422  0.9889
    ##  EMM_F48 - SP_F14   0.49000 0.394 222   1.244  0.9972
    ##  EMM_F48 - EMM_F7   0.33667 0.394 222   0.855  1.0000
    ##  EMM_F48 - ZAN_F3   0.33000 0.394 222   0.838  1.0000
    ##  EMM_F48 - EMM_F66  0.32333 0.394 222   0.821  1.0000
    ##  EMM_F48 - EMM_F34  0.29333 0.394 222   0.745  1.0000
    ##  EMM_F48 - EMM_F3   0.24333 0.394 222   0.618  1.0000
    ##  EMM_F48 - EMM_F89  0.21667 0.394 222   0.550  1.0000
    ##  EMM_F48 - EMM_F65  0.19000 0.394 222   0.483  1.0000
    ##  EMM_F48 - EMM_F5   0.18333 0.394 222   0.466  1.0000
    ##  EMM_F48 - EMM_F70  0.17333 0.394 222   0.440  1.0000
    ##  EMM_F48 - EMM_F63  0.17333 0.394 222   0.440  1.0000
    ##  EMM_F48 - EMM_F64  0.12667 0.394 222   0.322  1.0000
    ##  EMM_F48 - EMM_F49  0.03000 0.394 222   0.076  1.0000
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F34 - ZAN_F4   0.05000 0.394 222   0.127  1.0000
    ##  Control - ZAN_F4   0.22667 0.394 222   0.576  1.0000
    ##  Control - EMM_F34  0.17667 0.394 222   0.449  1.0000
    ##  EMM_F66 - ZAN_F4   0.26000 0.394 222   0.660  1.0000
    ##  EMM_F66 - EMM_F34  0.21000 0.394 222   0.533  1.0000
    ##  EMM_F66 - Control  0.03333 0.394 222   0.085  1.0000
    ##  SP_F14 - ZAN_F4    0.39000 0.394 222   0.990  0.9998
    ##  SP_F14 - EMM_F34   0.34000 0.394 222   0.863  1.0000
    ##  SP_F14 - Control   0.16333 0.394 222   0.415  1.0000
    ##  SP_F14 - EMM_F66   0.13000 0.394 222   0.330  1.0000
    ##  EMM_F89 - ZAN_F4   0.44667 0.394 222   1.134  0.9990
    ##  EMM_F89 - EMM_F34  0.39667 0.394 222   1.007  0.9998
    ##  EMM_F89 - Control  0.22000 0.394 222   0.559  1.0000
    ##  EMM_F89 - EMM_F66  0.18667 0.394 222   0.474  1.0000
    ##  EMM_F89 - SP_F14   0.05667 0.394 222   0.144  1.0000
    ##  EMM_F3 - ZAN_F4    0.47000 0.394 222   1.194  0.9983
    ##  EMM_F3 - EMM_F34   0.42000 0.394 222   1.067  0.9995
    ##  EMM_F3 - Control   0.24333 0.394 222   0.618  1.0000
    ##  EMM_F3 - EMM_F66   0.21000 0.394 222   0.533  1.0000
    ##  EMM_F3 - SP_F14    0.08000 0.394 222   0.203  1.0000
    ##  EMM_F3 - EMM_F89   0.02333 0.394 222   0.059  1.0000
    ##  EMM_F63 - ZAN_F4   0.55333 0.394 222   1.405  0.9902
    ##  EMM_F63 - EMM_F34  0.50333 0.394 222   1.278  0.9963
    ##  EMM_F63 - Control  0.32667 0.394 222   0.830  1.0000
    ##  EMM_F63 - EMM_F66  0.29333 0.394 222   0.745  1.0000
    ##  EMM_F63 - SP_F14   0.16333 0.394 222   0.415  1.0000
    ##  EMM_F63 - EMM_F89  0.10667 0.394 222   0.271  1.0000
    ##  EMM_F63 - EMM_F3   0.08333 0.394 222   0.212  1.0000
    ##  EMM_F49 - ZAN_F4   0.55667 0.394 222   1.414  0.9896
    ##  EMM_F49 - EMM_F34  0.50667 0.394 222   1.287  0.9960
    ##  EMM_F49 - Control  0.33000 0.394 222   0.838  1.0000
    ##  EMM_F49 - EMM_F66  0.29667 0.394 222   0.753  1.0000
    ##  EMM_F49 - SP_F14   0.16667 0.394 222   0.423  1.0000
    ##  EMM_F49 - EMM_F89  0.11000 0.394 222   0.279  1.0000
    ##  EMM_F49 - EMM_F3   0.08667 0.394 222   0.220  1.0000
    ##  EMM_F49 - EMM_F63  0.00333 0.394 222   0.008  1.0000
    ##  EMM_F7 - ZAN_F4    0.59667 0.394 222   1.515  0.9797
    ##  EMM_F7 - EMM_F34   0.54667 0.394 222   1.388  0.9913
    ##  EMM_F7 - Control   0.37000 0.394 222   0.940  0.9999
    ##  EMM_F7 - EMM_F66   0.33667 0.394 222   0.855  1.0000
    ##  EMM_F7 - SP_F14    0.20667 0.394 222   0.525  1.0000
    ##  EMM_F7 - EMM_F89   0.15000 0.394 222   0.381  1.0000
    ##  EMM_F7 - EMM_F3    0.12667 0.394 222   0.322  1.0000
    ##  EMM_F7 - EMM_F63   0.04333 0.394 222   0.110  1.0000
    ##  EMM_F7 - EMM_F49   0.04000 0.394 222   0.102  1.0000
    ##  ZAN_F3 - ZAN_F4    0.62333 0.394 222   1.583  0.9699
    ##  ZAN_F3 - EMM_F34   0.57333 0.394 222   1.456  0.9861
    ##  ZAN_F3 - Control   0.39667 0.394 222   1.007  0.9998
    ##  ZAN_F3 - EMM_F66   0.36333 0.394 222   0.923  0.9999
    ##  ZAN_F3 - SP_F14    0.23333 0.394 222   0.593  1.0000
    ##  ZAN_F3 - EMM_F89   0.17667 0.394 222   0.449  1.0000
    ##  ZAN_F3 - EMM_F3    0.15333 0.394 222   0.389  1.0000
    ##  ZAN_F3 - EMM_F63   0.07000 0.394 222   0.178  1.0000
    ##  ZAN_F3 - EMM_F49   0.06667 0.394 222   0.169  1.0000
    ##  ZAN_F3 - EMM_F7    0.02667 0.394 222   0.068  1.0000
    ##  EMM_F70 - ZAN_F4   0.68333 0.394 222   1.735  0.9351
    ##  EMM_F70 - EMM_F34  0.63333 0.394 222   1.608  0.9654
    ##  EMM_F70 - Control  0.45667 0.394 222   1.160  0.9987
    ##  EMM_F70 - EMM_F66  0.42333 0.394 222   1.075  0.9995
    ##  EMM_F70 - SP_F14   0.29333 0.394 222   0.745  1.0000
    ##  EMM_F70 - EMM_F89  0.23667 0.394 222   0.601  1.0000
    ##  EMM_F70 - EMM_F3   0.21333 0.394 222   0.542  1.0000
    ##  EMM_F70 - EMM_F63  0.13000 0.394 222   0.330  1.0000
    ##  EMM_F70 - EMM_F49  0.12667 0.394 222   0.322  1.0000
    ##  EMM_F70 - EMM_F7   0.08667 0.394 222   0.220  1.0000
    ##  EMM_F70 - ZAN_F3   0.06000 0.394 222   0.152  1.0000
    ##  EMM_F64 - ZAN_F4   0.70000 0.394 222   1.778  0.9217
    ##  EMM_F64 - EMM_F34  0.65000 0.394 222   1.651  0.9568
    ##  EMM_F64 - Control  0.47333 0.394 222   1.202  0.9981
    ##  EMM_F64 - EMM_F66  0.44000 0.394 222   1.117  0.9992
    ##  EMM_F64 - SP_F14   0.31000 0.394 222   0.787  1.0000
    ##  EMM_F64 - EMM_F89  0.25333 0.394 222   0.643  1.0000
    ##  EMM_F64 - EMM_F3   0.23000 0.394 222   0.584  1.0000
    ##  EMM_F64 - EMM_F63  0.14667 0.394 222   0.372  1.0000
    ##  EMM_F64 - EMM_F49  0.14333 0.394 222   0.364  1.0000
    ##  EMM_F64 - EMM_F7   0.10333 0.394 222   0.262  1.0000
    ##  EMM_F64 - ZAN_F3   0.07667 0.394 222   0.195  1.0000
    ##  EMM_F64 - EMM_F70  0.01667 0.394 222   0.042  1.0000
    ##  EMM_F65 - ZAN_F4   0.70333 0.394 222   1.786  0.9189
    ##  EMM_F65 - EMM_F34  0.65333 0.394 222   1.659  0.9549
    ##  EMM_F65 - Control  0.47667 0.394 222   1.210  0.9980
    ##  EMM_F65 - EMM_F66  0.44333 0.394 222   1.126  0.9991
    ##  EMM_F65 - SP_F14   0.31333 0.394 222   0.796  1.0000
    ##  EMM_F65 - EMM_F89  0.25667 0.394 222   0.652  1.0000
    ##  EMM_F65 - EMM_F3   0.23333 0.394 222   0.593  1.0000
    ##  EMM_F65 - EMM_F63  0.15000 0.394 222   0.381  1.0000
    ##  EMM_F65 - EMM_F49  0.14667 0.394 222   0.372  1.0000
    ##  EMM_F65 - EMM_F7   0.10667 0.394 222   0.271  1.0000
    ##  EMM_F65 - ZAN_F3   0.08000 0.394 222   0.203  1.0000
    ##  EMM_F65 - EMM_F70  0.02000 0.394 222   0.051  1.0000
    ##  EMM_F65 - EMM_F64  0.00333 0.394 222   0.008  1.0000
    ##  EMM_F5 - ZAN_F4    0.80000 0.394 222   2.032  0.8049
    ##  EMM_F5 - EMM_F34   0.75000 0.394 222   1.905  0.8712
    ##  EMM_F5 - Control   0.57333 0.394 222   1.456  0.9861
    ##  EMM_F5 - EMM_F66   0.54000 0.394 222   1.371  0.9923
    ##  EMM_F5 - SP_F14    0.41000 0.394 222   1.041  0.9996
    ##  EMM_F5 - EMM_F89   0.35333 0.394 222   0.897  0.9999
    ##  EMM_F5 - EMM_F3    0.33000 0.394 222   0.838  1.0000
    ##  EMM_F5 - EMM_F63   0.24667 0.394 222   0.626  1.0000
    ##  EMM_F5 - EMM_F49   0.24333 0.394 222   0.618  1.0000
    ##  EMM_F5 - EMM_F7    0.20333 0.394 222   0.516  1.0000
    ##  EMM_F5 - ZAN_F3    0.17667 0.394 222   0.449  1.0000
    ##  EMM_F5 - EMM_F70   0.11667 0.394 222   0.296  1.0000
    ##  EMM_F5 - EMM_F64   0.10000 0.394 222   0.254  1.0000
    ##  EMM_F5 - EMM_F65   0.09667 0.394 222   0.245  1.0000
    ##  EMM_F48 - ZAN_F4   0.91333 0.394 222   2.319  0.6098
    ##  EMM_F48 - EMM_F34  0.86333 0.394 222   2.192  0.7015
    ##  EMM_F48 - Control  0.68667 0.394 222   1.744  0.9325
    ##  EMM_F48 - EMM_F66  0.65333 0.394 222   1.659  0.9549
    ##  EMM_F48 - SP_F14   0.52333 0.394 222   1.329  0.9944
    ##  EMM_F48 - EMM_F89  0.46667 0.394 222   1.185  0.9984
    ##  EMM_F48 - EMM_F3   0.44333 0.394 222   1.126  0.9991
    ##  EMM_F48 - EMM_F63  0.36000 0.394 222   0.914  0.9999
    ##  EMM_F48 - EMM_F49  0.35667 0.394 222   0.906  0.9999
    ##  EMM_F48 - EMM_F7   0.31667 0.394 222   0.804  1.0000
    ##  EMM_F48 - ZAN_F3   0.29000 0.394 222   0.736  1.0000
    ##  EMM_F48 - EMM_F70  0.23000 0.394 222   0.584  1.0000
    ##  EMM_F48 - EMM_F64  0.21333 0.394 222   0.542  1.0000
    ##  EMM_F48 - EMM_F65  0.21000 0.394 222   0.533  1.0000
    ##  EMM_F48 - EMM_F5   0.11333 0.394 222   0.288  1.0000
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F34 - EMM_F65  0.04667 0.394 222   0.119  1.0000
    ##  EMM_F3 - EMM_F65   0.09333 0.394 222   0.237  1.0000
    ##  EMM_F3 - EMM_F34   0.04667 0.394 222   0.119  1.0000
    ##  Control - EMM_F65  0.15667 0.394 222   0.398  1.0000
    ##  Control - EMM_F34  0.11000 0.394 222   0.279  1.0000
    ##  Control - EMM_F3   0.06333 0.394 222   0.161  1.0000
    ##  SP_F14 - EMM_F65   0.24667 0.394 222   0.626  1.0000
    ##  SP_F14 - EMM_F34   0.20000 0.394 222   0.508  1.0000
    ##  SP_F14 - EMM_F3    0.15333 0.394 222   0.389  1.0000
    ##  SP_F14 - Control   0.09000 0.394 222   0.229  1.0000
    ##  EMM_F7 - EMM_F65   0.32667 0.394 222   0.830  1.0000
    ##  EMM_F7 - EMM_F34   0.28000 0.394 222   0.711  1.0000
    ##  EMM_F7 - EMM_F3    0.23333 0.394 222   0.593  1.0000
    ##  EMM_F7 - Control   0.17000 0.394 222   0.432  1.0000
    ##  EMM_F7 - SP_F14    0.08000 0.394 222   0.203  1.0000
    ##  EMM_F5 - EMM_F65   0.34667 0.394 222   0.880  1.0000
    ##  EMM_F5 - EMM_F34   0.30000 0.394 222   0.762  1.0000
    ##  EMM_F5 - EMM_F3    0.25333 0.394 222   0.643  1.0000
    ##  EMM_F5 - Control   0.19000 0.394 222   0.483  1.0000
    ##  EMM_F5 - SP_F14    0.10000 0.394 222   0.254  1.0000
    ##  EMM_F5 - EMM_F7    0.02000 0.394 222   0.051  1.0000
    ##  EMM_F49 - EMM_F65  0.42667 0.394 222   1.084  0.9994
    ##  EMM_F49 - EMM_F34  0.38000 0.394 222   0.965  0.9999
    ##  EMM_F49 - EMM_F3   0.33333 0.394 222   0.846  1.0000
    ##  EMM_F49 - Control  0.27000 0.394 222   0.686  1.0000
    ##  EMM_F49 - SP_F14   0.18000 0.394 222   0.457  1.0000
    ##  EMM_F49 - EMM_F7   0.10000 0.394 222   0.254  1.0000
    ##  EMM_F49 - EMM_F5   0.08000 0.394 222   0.203  1.0000
    ##  ZAN_F4 - EMM_F65   0.45000 0.394 222   1.143  0.9989
    ##  ZAN_F4 - EMM_F34   0.40333 0.394 222   1.024  0.9997
    ##  ZAN_F4 - EMM_F3    0.35667 0.394 222   0.906  0.9999
    ##  ZAN_F4 - Control   0.29333 0.394 222   0.745  1.0000
    ##  ZAN_F4 - SP_F14    0.20333 0.394 222   0.516  1.0000
    ##  ZAN_F4 - EMM_F7    0.12333 0.394 222   0.313  1.0000
    ##  ZAN_F4 - EMM_F5    0.10333 0.394 222   0.262  1.0000
    ##  ZAN_F4 - EMM_F49   0.02333 0.394 222   0.059  1.0000
    ##  EMM_F66 - EMM_F65  0.45667 0.394 222   1.160  0.9987
    ##  EMM_F66 - EMM_F34  0.41000 0.394 222   1.041  0.9996
    ##  EMM_F66 - EMM_F3   0.36333 0.394 222   0.923  0.9999
    ##  EMM_F66 - Control  0.30000 0.394 222   0.762  1.0000
    ##  EMM_F66 - SP_F14   0.21000 0.394 222   0.533  1.0000
    ##  EMM_F66 - EMM_F7   0.13000 0.394 222   0.330  1.0000
    ##  EMM_F66 - EMM_F5   0.11000 0.394 222   0.279  1.0000
    ##  EMM_F66 - EMM_F49  0.03000 0.394 222   0.076  1.0000
    ##  EMM_F66 - ZAN_F4   0.00667 0.394 222   0.017  1.0000
    ##  EMM_F89 - EMM_F65  0.58333 0.394 222   1.481  0.9836
    ##  EMM_F89 - EMM_F34  0.53667 0.394 222   1.363  0.9928
    ##  EMM_F89 - EMM_F3   0.49000 0.394 222   1.244  0.9972
    ##  EMM_F89 - Control  0.42667 0.394 222   1.084  0.9994
    ##  EMM_F89 - SP_F14   0.33667 0.394 222   0.855  1.0000
    ##  EMM_F89 - EMM_F7   0.25667 0.394 222   0.652  1.0000
    ##  EMM_F89 - EMM_F5   0.23667 0.394 222   0.601  1.0000
    ##  EMM_F89 - EMM_F49  0.15667 0.394 222   0.398  1.0000
    ##  EMM_F89 - ZAN_F4   0.13333 0.394 222   0.339  1.0000
    ##  EMM_F89 - EMM_F66  0.12667 0.394 222   0.322  1.0000
    ##  EMM_F64 - EMM_F65  0.62000 0.394 222   1.574  0.9713
    ##  EMM_F64 - EMM_F34  0.57333 0.394 222   1.456  0.9861
    ##  EMM_F64 - EMM_F3   0.52667 0.394 222   1.337  0.9940
    ##  EMM_F64 - Control  0.46333 0.394 222   1.177  0.9985
    ##  EMM_F64 - SP_F14   0.37333 0.394 222   0.948  0.9999
    ##  EMM_F64 - EMM_F7   0.29333 0.394 222   0.745  1.0000
    ##  EMM_F64 - EMM_F5   0.27333 0.394 222   0.694  1.0000
    ##  EMM_F64 - EMM_F49  0.19333 0.394 222   0.491  1.0000
    ##  EMM_F64 - ZAN_F4   0.17000 0.394 222   0.432  1.0000
    ##  EMM_F64 - EMM_F66  0.16333 0.394 222   0.415  1.0000
    ##  EMM_F64 - EMM_F89  0.03667 0.394 222   0.093  1.0000
    ##  EMM_F48 - EMM_F65  0.64000 0.394 222   1.625  0.9621
    ##  EMM_F48 - EMM_F34  0.59333 0.394 222   1.507  0.9808
    ##  EMM_F48 - EMM_F3   0.54667 0.394 222   1.388  0.9913
    ##  EMM_F48 - Control  0.48333 0.394 222   1.227  0.9976
    ##  EMM_F48 - SP_F14   0.39333 0.394 222   0.999  0.9998
    ##  EMM_F48 - EMM_F7   0.31333 0.394 222   0.796  1.0000
    ##  EMM_F48 - EMM_F5   0.29333 0.394 222   0.745  1.0000
    ##  EMM_F48 - EMM_F49  0.21333 0.394 222   0.542  1.0000
    ##  EMM_F48 - ZAN_F4   0.19000 0.394 222   0.483  1.0000
    ##  EMM_F48 - EMM_F66  0.18333 0.394 222   0.466  1.0000
    ##  EMM_F48 - EMM_F89  0.05667 0.394 222   0.144  1.0000
    ##  EMM_F48 - EMM_F64  0.02000 0.394 222   0.051  1.0000
    ##  EMM_F70 - EMM_F65  0.74000 0.394 222   1.879  0.8826
    ##  EMM_F70 - EMM_F34  0.69333 0.394 222   1.761  0.9273
    ##  EMM_F70 - EMM_F3   0.64667 0.394 222   1.642  0.9586
    ##  EMM_F70 - Control  0.58333 0.394 222   1.481  0.9836
    ##  EMM_F70 - SP_F14   0.49333 0.394 222   1.253  0.9970
    ##  EMM_F70 - EMM_F7   0.41333 0.394 222   1.050  0.9996
    ##  EMM_F70 - EMM_F5   0.39333 0.394 222   0.999  0.9998
    ##  EMM_F70 - EMM_F49  0.31333 0.394 222   0.796  1.0000
    ##  EMM_F70 - ZAN_F4   0.29000 0.394 222   0.736  1.0000
    ##  EMM_F70 - EMM_F66  0.28333 0.394 222   0.720  1.0000
    ##  EMM_F70 - EMM_F89  0.15667 0.394 222   0.398  1.0000
    ##  EMM_F70 - EMM_F64  0.12000 0.394 222   0.305  1.0000
    ##  EMM_F70 - EMM_F48  0.10000 0.394 222   0.254  1.0000
    ##  EMM_F63 - EMM_F65  0.78333 0.394 222   1.989  0.8287
    ##  EMM_F63 - EMM_F34  0.73667 0.394 222   1.871  0.8863
    ##  EMM_F63 - EMM_F3   0.69000 0.394 222   1.752  0.9299
    ##  EMM_F63 - Control  0.62667 0.394 222   1.591  0.9685
    ##  EMM_F63 - SP_F14   0.53667 0.394 222   1.363  0.9928
    ##  EMM_F63 - EMM_F7   0.45667 0.394 222   1.160  0.9987
    ##  EMM_F63 - EMM_F5   0.43667 0.394 222   1.109  0.9992
    ##  EMM_F63 - EMM_F49  0.35667 0.394 222   0.906  0.9999
    ##  EMM_F63 - ZAN_F4   0.33333 0.394 222   0.846  1.0000
    ##  EMM_F63 - EMM_F66  0.32667 0.394 222   0.830  1.0000
    ##  EMM_F63 - EMM_F89  0.20000 0.394 222   0.508  1.0000
    ##  EMM_F63 - EMM_F64  0.16333 0.394 222   0.415  1.0000
    ##  EMM_F63 - EMM_F48  0.14333 0.394 222   0.364  1.0000
    ##  EMM_F63 - EMM_F70  0.04333 0.394 222   0.110  1.0000
    ##  ZAN_F3 - EMM_F65   0.85333 0.394 222   2.167  0.7190
    ##  ZAN_F3 - EMM_F34   0.80667 0.394 222   2.049  0.7949
    ##  ZAN_F3 - EMM_F3    0.76000 0.394 222   1.930  0.8592
    ##  ZAN_F3 - Control   0.69667 0.394 222   1.769  0.9246
    ##  ZAN_F3 - SP_F14    0.60667 0.394 222   1.541  0.9764
    ##  ZAN_F3 - EMM_F7    0.52667 0.394 222   1.337  0.9940
    ##  ZAN_F3 - EMM_F5    0.50667 0.394 222   1.287  0.9960
    ##  ZAN_F3 - EMM_F49   0.42667 0.394 222   1.084  0.9994
    ##  ZAN_F3 - ZAN_F4    0.40333 0.394 222   1.024  0.9997
    ##  ZAN_F3 - EMM_F66   0.39667 0.394 222   1.007  0.9998
    ##  ZAN_F3 - EMM_F89   0.27000 0.394 222   0.686  1.0000
    ##  ZAN_F3 - EMM_F64   0.23333 0.394 222   0.593  1.0000
    ##  ZAN_F3 - EMM_F48   0.21333 0.394 222   0.542  1.0000
    ##  ZAN_F3 - EMM_F70   0.11333 0.394 222   0.288  1.0000
    ##  ZAN_F3 - EMM_F63   0.07000 0.394 222   0.178  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  ZAN_F4 - EMM_F34   0.05667 0.394 222   0.144  1.0000
    ##  EMM_F5 - EMM_F34   0.12333 0.394 222   0.313  1.0000
    ##  EMM_F5 - ZAN_F4    0.06667 0.394 222   0.169  1.0000
    ##  Control - EMM_F34  0.17667 0.394 222   0.449  1.0000
    ##  Control - ZAN_F4   0.12000 0.394 222   0.305  1.0000
    ##  Control - EMM_F5   0.05333 0.394 222   0.135  1.0000
    ##  EMM_F49 - EMM_F34  0.22333 0.394 222   0.567  1.0000
    ##  EMM_F49 - ZAN_F4   0.16667 0.394 222   0.423  1.0000
    ##  EMM_F49 - EMM_F5   0.10000 0.394 222   0.254  1.0000
    ##  EMM_F49 - Control  0.04667 0.394 222   0.119  1.0000
    ##  EMM_F3 - EMM_F34   0.22667 0.394 222   0.576  1.0000
    ##  EMM_F3 - ZAN_F4    0.17000 0.394 222   0.432  1.0000
    ##  EMM_F3 - EMM_F5    0.10333 0.394 222   0.262  1.0000
    ##  EMM_F3 - Control   0.05000 0.394 222   0.127  1.0000
    ##  EMM_F3 - EMM_F49   0.00333 0.394 222   0.008  1.0000
    ##  EMM_F65 - EMM_F34  0.34667 0.394 222   0.880  1.0000
    ##  EMM_F65 - ZAN_F4   0.29000 0.394 222   0.736  1.0000
    ##  EMM_F65 - EMM_F5   0.22333 0.394 222   0.567  1.0000
    ##  EMM_F65 - Control  0.17000 0.394 222   0.432  1.0000
    ##  EMM_F65 - EMM_F49  0.12333 0.394 222   0.313  1.0000
    ##  EMM_F65 - EMM_F3   0.12000 0.394 222   0.305  1.0000
    ##  EMM_F63 - EMM_F34  0.35667 0.394 222   0.906  0.9999
    ##  EMM_F63 - ZAN_F4   0.30000 0.394 222   0.762  1.0000
    ##  EMM_F63 - EMM_F5   0.23333 0.394 222   0.593  1.0000
    ##  EMM_F63 - Control  0.18000 0.394 222   0.457  1.0000
    ##  EMM_F63 - EMM_F49  0.13333 0.394 222   0.339  1.0000
    ##  EMM_F63 - EMM_F3   0.13000 0.394 222   0.330  1.0000
    ##  EMM_F63 - EMM_F65  0.01000 0.394 222   0.025  1.0000
    ##  EMM_F66 - EMM_F34  0.36000 0.394 222   0.914  0.9999
    ##  EMM_F66 - ZAN_F4   0.30333 0.394 222   0.770  1.0000
    ##  EMM_F66 - EMM_F5   0.23667 0.394 222   0.601  1.0000
    ##  EMM_F66 - Control  0.18333 0.394 222   0.466  1.0000
    ##  EMM_F66 - EMM_F49  0.13667 0.394 222   0.347  1.0000
    ##  EMM_F66 - EMM_F3   0.13333 0.394 222   0.339  1.0000
    ##  EMM_F66 - EMM_F65  0.01333 0.394 222   0.034  1.0000
    ##  EMM_F66 - EMM_F63  0.00333 0.394 222   0.008  1.0000
    ##  EMM_F70 - EMM_F34  0.43667 0.394 222   1.109  0.9992
    ##  EMM_F70 - ZAN_F4   0.38000 0.394 222   0.965  0.9999
    ##  EMM_F70 - EMM_F5   0.31333 0.394 222   0.796  1.0000
    ##  EMM_F70 - Control  0.26000 0.394 222   0.660  1.0000
    ##  EMM_F70 - EMM_F49  0.21333 0.394 222   0.542  1.0000
    ##  EMM_F70 - EMM_F3   0.21000 0.394 222   0.533  1.0000
    ##  EMM_F70 - EMM_F65  0.09000 0.394 222   0.229  1.0000
    ##  EMM_F70 - EMM_F63  0.08000 0.394 222   0.203  1.0000
    ##  EMM_F70 - EMM_F66  0.07667 0.394 222   0.195  1.0000
    ##  EMM_F48 - EMM_F34  0.43667 0.394 222   1.109  0.9992
    ##  EMM_F48 - ZAN_F4   0.38000 0.394 222   0.965  0.9999
    ##  EMM_F48 - EMM_F5   0.31333 0.394 222   0.796  1.0000
    ##  EMM_F48 - Control  0.26000 0.394 222   0.660  1.0000
    ##  EMM_F48 - EMM_F49  0.21333 0.394 222   0.542  1.0000
    ##  EMM_F48 - EMM_F3   0.21000 0.394 222   0.533  1.0000
    ##  EMM_F48 - EMM_F65  0.09000 0.394 222   0.229  1.0000
    ##  EMM_F48 - EMM_F63  0.08000 0.394 222   0.203  1.0000
    ##  EMM_F48 - EMM_F66  0.07667 0.394 222   0.195  1.0000
    ##  EMM_F48 - EMM_F70  0.00000 0.394 222   0.000  1.0000
    ##  EMM_F64 - EMM_F34  0.44667 0.394 222   1.134  0.9990
    ##  EMM_F64 - ZAN_F4   0.39000 0.394 222   0.990  0.9998
    ##  EMM_F64 - EMM_F5   0.32333 0.394 222   0.821  1.0000
    ##  EMM_F64 - Control  0.27000 0.394 222   0.686  1.0000
    ##  EMM_F64 - EMM_F49  0.22333 0.394 222   0.567  1.0000
    ##  EMM_F64 - EMM_F3   0.22000 0.394 222   0.559  1.0000
    ##  EMM_F64 - EMM_F65  0.10000 0.394 222   0.254  1.0000
    ##  EMM_F64 - EMM_F63  0.09000 0.394 222   0.229  1.0000
    ##  EMM_F64 - EMM_F66  0.08667 0.394 222   0.220  1.0000
    ##  EMM_F64 - EMM_F70  0.01000 0.394 222   0.025  1.0000
    ##  EMM_F64 - EMM_F48  0.01000 0.394 222   0.025  1.0000
    ##  EMM_F7 - EMM_F34   0.46333 0.394 222   1.177  0.9985
    ##  EMM_F7 - ZAN_F4    0.40667 0.394 222   1.033  0.9997
    ##  EMM_F7 - EMM_F5    0.34000 0.394 222   0.863  1.0000
    ##  EMM_F7 - Control   0.28667 0.394 222   0.728  1.0000
    ##  EMM_F7 - EMM_F49   0.24000 0.394 222   0.609  1.0000
    ##  EMM_F7 - EMM_F3    0.23667 0.394 222   0.601  1.0000
    ##  EMM_F7 - EMM_F65   0.11667 0.394 222   0.296  1.0000
    ##  EMM_F7 - EMM_F63   0.10667 0.394 222   0.271  1.0000
    ##  EMM_F7 - EMM_F66   0.10333 0.394 222   0.262  1.0000
    ##  EMM_F7 - EMM_F70   0.02667 0.394 222   0.068  1.0000
    ##  EMM_F7 - EMM_F48   0.02667 0.394 222   0.068  1.0000
    ##  EMM_F7 - EMM_F64   0.01667 0.394 222   0.042  1.0000
    ##  EMM_F89 - EMM_F34  0.60667 0.394 222   1.541  0.9764
    ##  EMM_F89 - ZAN_F4   0.55000 0.394 222   1.397  0.9907
    ##  EMM_F89 - EMM_F5   0.48333 0.394 222   1.227  0.9976
    ##  EMM_F89 - Control  0.43000 0.394 222   1.092  0.9994
    ##  EMM_F89 - EMM_F49  0.38333 0.394 222   0.973  0.9998
    ##  EMM_F89 - EMM_F3   0.38000 0.394 222   0.965  0.9999
    ##  EMM_F89 - EMM_F65  0.26000 0.394 222   0.660  1.0000
    ##  EMM_F89 - EMM_F63  0.25000 0.394 222   0.635  1.0000
    ##  EMM_F89 - EMM_F66  0.24667 0.394 222   0.626  1.0000
    ##  EMM_F89 - EMM_F70  0.17000 0.394 222   0.432  1.0000
    ##  EMM_F89 - EMM_F48  0.17000 0.394 222   0.432  1.0000
    ##  EMM_F89 - EMM_F64  0.16000 0.394 222   0.406  1.0000
    ##  EMM_F89 - EMM_F7   0.14333 0.394 222   0.364  1.0000
    ##  ZAN_F3 - EMM_F34   0.74333 0.394 222   1.888  0.8789
    ##  ZAN_F3 - ZAN_F4    0.68667 0.394 222   1.744  0.9325
    ##  ZAN_F3 - EMM_F5    0.62000 0.394 222   1.574  0.9713
    ##  ZAN_F3 - Control   0.56667 0.394 222   1.439  0.9876
    ##  ZAN_F3 - EMM_F49   0.52000 0.394 222   1.321  0.9948
    ##  ZAN_F3 - EMM_F3    0.51667 0.394 222   1.312  0.9951
    ##  ZAN_F3 - EMM_F65   0.39667 0.394 222   1.007  0.9998
    ##  ZAN_F3 - EMM_F63   0.38667 0.394 222   0.982  0.9998
    ##  ZAN_F3 - EMM_F66   0.38333 0.394 222   0.973  0.9998
    ##  ZAN_F3 - EMM_F70   0.30667 0.394 222   0.779  1.0000
    ##  ZAN_F3 - EMM_F48   0.30667 0.394 222   0.779  1.0000
    ##  ZAN_F3 - EMM_F64   0.29667 0.394 222   0.753  1.0000
    ##  ZAN_F3 - EMM_F7    0.28000 0.394 222   0.711  1.0000
    ##  ZAN_F3 - EMM_F89   0.13667 0.394 222   0.347  1.0000
    ##  SP_F14 - EMM_F34   0.80000 0.394 222   2.032  0.8049
    ##  SP_F14 - ZAN_F4    0.74333 0.394 222   1.888  0.8789
    ##  SP_F14 - EMM_F5    0.67667 0.394 222   1.718  0.9399
    ##  SP_F14 - Control   0.62333 0.394 222   1.583  0.9699
    ##  SP_F14 - EMM_F49   0.57667 0.394 222   1.464  0.9853
    ##  SP_F14 - EMM_F3    0.57333 0.394 222   1.456  0.9861
    ##  SP_F14 - EMM_F65   0.45333 0.394 222   1.151  0.9988
    ##  SP_F14 - EMM_F63   0.44333 0.394 222   1.126  0.9991
    ##  SP_F14 - EMM_F66   0.44000 0.394 222   1.117  0.9992
    ##  SP_F14 - EMM_F70   0.36333 0.394 222   0.923  0.9999
    ##  SP_F14 - EMM_F48   0.36333 0.394 222   0.923  0.9999
    ##  SP_F14 - EMM_F64   0.35333 0.394 222   0.897  0.9999
    ##  SP_F14 - EMM_F7    0.33667 0.394 222   0.855  1.0000
    ##  SP_F14 - EMM_F89   0.19333 0.394 222   0.491  1.0000
    ##  SP_F14 - ZAN_F3    0.05667 0.394 222   0.144  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F3 - EMM_F34   0.14000 0.394 222   0.356  1.0000
    ##  EMM_F5 - EMM_F34   0.20000 0.394 222   0.508  1.0000
    ##  EMM_F5 - EMM_F3    0.06000 0.394 222   0.152  1.0000
    ##  EMM_F48 - EMM_F34  0.22333 0.394 222   0.567  1.0000
    ##  EMM_F48 - EMM_F3   0.08333 0.394 222   0.212  1.0000
    ##  EMM_F48 - EMM_F5   0.02333 0.394 222   0.059  1.0000
    ##  ZAN_F4 - EMM_F34   0.25000 0.394 222   0.635  1.0000
    ##  ZAN_F4 - EMM_F3    0.11000 0.394 222   0.279  1.0000
    ##  ZAN_F4 - EMM_F5    0.05000 0.394 222   0.127  1.0000
    ##  ZAN_F4 - EMM_F48   0.02667 0.394 222   0.068  1.0000
    ##  Control - EMM_F34  0.40000 0.394 222   1.016  0.9997
    ##  Control - EMM_F3   0.26000 0.394 222   0.660  1.0000
    ##  Control - EMM_F5   0.20000 0.394 222   0.508  1.0000
    ##  Control - EMM_F48  0.17667 0.394 222   0.449  1.0000
    ##  Control - ZAN_F4   0.15000 0.394 222   0.381  1.0000
    ##  EMM_F63 - EMM_F34  0.42667 0.394 222   1.084  0.9994
    ##  EMM_F63 - EMM_F3   0.28667 0.394 222   0.728  1.0000
    ##  EMM_F63 - EMM_F5   0.22667 0.394 222   0.576  1.0000
    ##  EMM_F63 - EMM_F48  0.20333 0.394 222   0.516  1.0000
    ##  EMM_F63 - ZAN_F4   0.17667 0.394 222   0.449  1.0000
    ##  EMM_F63 - Control  0.02667 0.394 222   0.068  1.0000
    ##  EMM_F7 - EMM_F34   0.47333 0.394 222   1.202  0.9981
    ##  EMM_F7 - EMM_F3    0.33333 0.394 222   0.846  1.0000
    ##  EMM_F7 - EMM_F5    0.27333 0.394 222   0.694  1.0000
    ##  EMM_F7 - EMM_F48   0.25000 0.394 222   0.635  1.0000
    ##  EMM_F7 - ZAN_F4    0.22333 0.394 222   0.567  1.0000
    ##  EMM_F7 - Control   0.07333 0.394 222   0.186  1.0000
    ##  EMM_F7 - EMM_F63   0.04667 0.394 222   0.119  1.0000
    ##  SP_F14 - EMM_F34   0.52000 0.394 222   1.321  0.9948
    ##  SP_F14 - EMM_F3    0.38000 0.394 222   0.965  0.9999
    ##  SP_F14 - EMM_F5    0.32000 0.394 222   0.813  1.0000
    ##  SP_F14 - EMM_F48   0.29667 0.394 222   0.753  1.0000
    ##  SP_F14 - ZAN_F4    0.27000 0.394 222   0.686  1.0000
    ##  SP_F14 - Control   0.12000 0.394 222   0.305  1.0000
    ##  SP_F14 - EMM_F63   0.09333 0.394 222   0.237  1.0000
    ##  SP_F14 - EMM_F7    0.04667 0.394 222   0.119  1.0000
    ##  EMM_F64 - EMM_F34  0.54000 0.394 222   1.371  0.9923
    ##  EMM_F64 - EMM_F3   0.40000 0.394 222   1.016  0.9997
    ##  EMM_F64 - EMM_F5   0.34000 0.394 222   0.863  1.0000
    ##  EMM_F64 - EMM_F48  0.31667 0.394 222   0.804  1.0000
    ##  EMM_F64 - ZAN_F4   0.29000 0.394 222   0.736  1.0000
    ##  EMM_F64 - Control  0.14000 0.394 222   0.356  1.0000
    ##  EMM_F64 - EMM_F63  0.11333 0.394 222   0.288  1.0000
    ##  EMM_F64 - EMM_F7   0.06667 0.394 222   0.169  1.0000
    ##  EMM_F64 - SP_F14   0.02000 0.394 222   0.051  1.0000
    ##  EMM_F66 - EMM_F34  0.57667 0.394 222   1.464  0.9853
    ##  EMM_F66 - EMM_F3   0.43667 0.394 222   1.109  0.9992
    ##  EMM_F66 - EMM_F5   0.37667 0.394 222   0.957  0.9999
    ##  EMM_F66 - EMM_F48  0.35333 0.394 222   0.897  0.9999
    ##  EMM_F66 - ZAN_F4   0.32667 0.394 222   0.830  1.0000
    ##  EMM_F66 - Control  0.17667 0.394 222   0.449  1.0000
    ##  EMM_F66 - EMM_F63  0.15000 0.394 222   0.381  1.0000
    ##  EMM_F66 - EMM_F7   0.10333 0.394 222   0.262  1.0000
    ##  EMM_F66 - SP_F14   0.05667 0.394 222   0.144  1.0000
    ##  EMM_F66 - EMM_F64  0.03667 0.394 222   0.093  1.0000
    ##  EMM_F70 - EMM_F34  0.67000 0.394 222   1.701  0.9445
    ##  EMM_F70 - EMM_F3   0.53000 0.394 222   1.346  0.9936
    ##  EMM_F70 - EMM_F5   0.47000 0.394 222   1.194  0.9983
    ##  EMM_F70 - EMM_F48  0.44667 0.394 222   1.134  0.9990
    ##  EMM_F70 - ZAN_F4   0.42000 0.394 222   1.067  0.9995
    ##  EMM_F70 - Control  0.27000 0.394 222   0.686  1.0000
    ##  EMM_F70 - EMM_F63  0.24333 0.394 222   0.618  1.0000
    ##  EMM_F70 - EMM_F7   0.19667 0.394 222   0.499  1.0000
    ##  EMM_F70 - SP_F14   0.15000 0.394 222   0.381  1.0000
    ##  EMM_F70 - EMM_F64  0.13000 0.394 222   0.330  1.0000
    ##  EMM_F70 - EMM_F66  0.09333 0.394 222   0.237  1.0000
    ##  EMM_F65 - EMM_F34  0.72333 0.394 222   1.837  0.9002
    ##  EMM_F65 - EMM_F3   0.58333 0.394 222   1.481  0.9836
    ##  EMM_F65 - EMM_F5   0.52333 0.394 222   1.329  0.9944
    ##  EMM_F65 - EMM_F48  0.50000 0.394 222   1.270  0.9966
    ##  EMM_F65 - ZAN_F4   0.47333 0.394 222   1.202  0.9981
    ##  EMM_F65 - Control  0.32333 0.394 222   0.821  1.0000
    ##  EMM_F65 - EMM_F63  0.29667 0.394 222   0.753  1.0000
    ##  EMM_F65 - EMM_F7   0.25000 0.394 222   0.635  1.0000
    ##  EMM_F65 - SP_F14   0.20333 0.394 222   0.516  1.0000
    ##  EMM_F65 - EMM_F64  0.18333 0.394 222   0.466  1.0000
    ##  EMM_F65 - EMM_F66  0.14667 0.394 222   0.372  1.0000
    ##  EMM_F65 - EMM_F70  0.05333 0.394 222   0.135  1.0000
    ##  EMM_F49 - EMM_F34  0.74000 0.394 222   1.879  0.8826
    ##  EMM_F49 - EMM_F3   0.60000 0.394 222   1.524  0.9787
    ##  EMM_F49 - EMM_F5   0.54000 0.394 222   1.371  0.9923
    ##  EMM_F49 - EMM_F48  0.51667 0.394 222   1.312  0.9951
    ##  EMM_F49 - ZAN_F4   0.49000 0.394 222   1.244  0.9972
    ##  EMM_F49 - Control  0.34000 0.394 222   0.863  1.0000
    ##  EMM_F49 - EMM_F63  0.31333 0.394 222   0.796  1.0000
    ##  EMM_F49 - EMM_F7   0.26667 0.394 222   0.677  1.0000
    ##  EMM_F49 - SP_F14   0.22000 0.394 222   0.559  1.0000
    ##  EMM_F49 - EMM_F64  0.20000 0.394 222   0.508  1.0000
    ##  EMM_F49 - EMM_F66  0.16333 0.394 222   0.415  1.0000
    ##  EMM_F49 - EMM_F70  0.07000 0.394 222   0.178  1.0000
    ##  EMM_F49 - EMM_F65  0.01667 0.394 222   0.042  1.0000
    ##  EMM_F89 - EMM_F34  0.92000 0.394 222   2.336  0.5972
    ##  EMM_F89 - EMM_F3   0.78000 0.394 222   1.981  0.8332
    ##  EMM_F89 - EMM_F5   0.72000 0.394 222   1.828  0.9035
    ##  EMM_F89 - EMM_F48  0.69667 0.394 222   1.769  0.9246
    ##  EMM_F89 - ZAN_F4   0.67000 0.394 222   1.701  0.9445
    ##  EMM_F89 - Control  0.52000 0.394 222   1.321  0.9948
    ##  EMM_F89 - EMM_F63  0.49333 0.394 222   1.253  0.9970
    ##  EMM_F89 - EMM_F7   0.44667 0.394 222   1.134  0.9990
    ##  EMM_F89 - SP_F14   0.40000 0.394 222   1.016  0.9997
    ##  EMM_F89 - EMM_F64  0.38000 0.394 222   0.965  0.9999
    ##  EMM_F89 - EMM_F66  0.34333 0.394 222   0.872  1.0000
    ##  EMM_F89 - EMM_F70  0.25000 0.394 222   0.635  1.0000
    ##  EMM_F89 - EMM_F65  0.19667 0.394 222   0.499  1.0000
    ##  EMM_F89 - EMM_F49  0.18000 0.394 222   0.457  1.0000
    ##  ZAN_F3 - EMM_F34   1.34000 0.394 222   3.403  0.0604
    ##  ZAN_F3 - EMM_F3    1.20000 0.394 222   3.047  0.1584
    ##  ZAN_F3 - EMM_F5    1.14000 0.394 222   2.895  0.2261
    ##  ZAN_F3 - EMM_F48   1.11667 0.394 222   2.836  0.2571
    ##  ZAN_F3 - ZAN_F4    1.09000 0.394 222   2.768  0.2956
    ##  ZAN_F3 - Control   0.94000 0.394 222   2.387  0.5593
    ##  ZAN_F3 - EMM_F63   0.91333 0.394 222   2.319  0.6098
    ##  ZAN_F3 - EMM_F7    0.86667 0.394 222   2.201  0.6956
    ##  ZAN_F3 - SP_F14    0.82000 0.394 222   2.082  0.7743
    ##  ZAN_F3 - EMM_F64   0.80000 0.394 222   2.032  0.8049
    ##  ZAN_F3 - EMM_F66   0.76333 0.394 222   1.938  0.8550
    ##  ZAN_F3 - EMM_F70   0.67000 0.394 222   1.701  0.9445
    ##  ZAN_F3 - EMM_F65   0.61667 0.394 222   1.566  0.9726
    ##  ZAN_F3 - EMM_F49   0.60000 0.394 222   1.524  0.9787
    ##  ZAN_F3 - EMM_F89   0.42000 0.394 222   1.067  0.9995
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F3 - EMM_F34   0.31333 0.394 222   0.796  1.0000
    ##  EMM_F66 - EMM_F34  0.32667 0.394 222   0.830  1.0000
    ##  EMM_F66 - EMM_F3   0.01333 0.394 222   0.034  1.0000
    ##  EMM_F49 - EMM_F34  0.38667 0.394 222   0.982  0.9998
    ##  EMM_F49 - EMM_F3   0.07333 0.394 222   0.186  1.0000
    ##  EMM_F49 - EMM_F66  0.06000 0.394 222   0.152  1.0000
    ##  EMM_F5 - EMM_F34   0.42667 0.394 222   1.084  0.9994
    ##  EMM_F5 - EMM_F3    0.11333 0.394 222   0.288  1.0000
    ##  EMM_F5 - EMM_F66   0.10000 0.394 222   0.254  1.0000
    ##  EMM_F5 - EMM_F49   0.04000 0.394 222   0.102  1.0000
    ##  EMM_F70 - EMM_F34  0.48000 0.394 222   1.219  0.9978
    ##  EMM_F70 - EMM_F3   0.16667 0.394 222   0.423  1.0000
    ##  EMM_F70 - EMM_F66  0.15333 0.394 222   0.389  1.0000
    ##  EMM_F70 - EMM_F49  0.09333 0.394 222   0.237  1.0000
    ##  EMM_F70 - EMM_F5   0.05333 0.394 222   0.135  1.0000
    ##  ZAN_F4 - EMM_F34   0.50000 0.394 222   1.270  0.9966
    ##  ZAN_F4 - EMM_F3    0.18667 0.394 222   0.474  1.0000
    ##  ZAN_F4 - EMM_F66   0.17333 0.394 222   0.440  1.0000
    ##  ZAN_F4 - EMM_F49   0.11333 0.394 222   0.288  1.0000
    ##  ZAN_F4 - EMM_F5    0.07333 0.394 222   0.186  1.0000
    ##  ZAN_F4 - EMM_F70   0.02000 0.394 222   0.051  1.0000
    ##  ZAN_F3 - EMM_F34   0.51667 0.394 222   1.312  0.9951
    ##  ZAN_F3 - EMM_F3    0.20333 0.394 222   0.516  1.0000
    ##  ZAN_F3 - EMM_F66   0.19000 0.394 222   0.483  1.0000
    ##  ZAN_F3 - EMM_F49   0.13000 0.394 222   0.330  1.0000
    ##  ZAN_F3 - EMM_F5    0.09000 0.394 222   0.229  1.0000
    ##  ZAN_F3 - EMM_F70   0.03667 0.394 222   0.093  1.0000
    ##  ZAN_F3 - ZAN_F4    0.01667 0.394 222   0.042  1.0000
    ##  EMM_F63 - EMM_F34  0.55333 0.394 222   1.405  0.9902
    ##  EMM_F63 - EMM_F3   0.24000 0.394 222   0.609  1.0000
    ##  EMM_F63 - EMM_F66  0.22667 0.394 222   0.576  1.0000
    ##  EMM_F63 - EMM_F49  0.16667 0.394 222   0.423  1.0000
    ##  EMM_F63 - EMM_F5   0.12667 0.394 222   0.322  1.0000
    ##  EMM_F63 - EMM_F70  0.07333 0.394 222   0.186  1.0000
    ##  EMM_F63 - ZAN_F4   0.05333 0.394 222   0.135  1.0000
    ##  EMM_F63 - ZAN_F3   0.03667 0.394 222   0.093  1.0000
    ##  EMM_F65 - EMM_F34  0.58333 0.394 222   1.481  0.9836
    ##  EMM_F65 - EMM_F3   0.27000 0.394 222   0.686  1.0000
    ##  EMM_F65 - EMM_F66  0.25667 0.394 222   0.652  1.0000
    ##  EMM_F65 - EMM_F49  0.19667 0.394 222   0.499  1.0000
    ##  EMM_F65 - EMM_F5   0.15667 0.394 222   0.398  1.0000
    ##  EMM_F65 - EMM_F70  0.10333 0.394 222   0.262  1.0000
    ##  EMM_F65 - ZAN_F4   0.08333 0.394 222   0.212  1.0000
    ##  EMM_F65 - ZAN_F3   0.06667 0.394 222   0.169  1.0000
    ##  EMM_F65 - EMM_F63  0.03000 0.394 222   0.076  1.0000
    ##  EMM_F48 - EMM_F34  0.60000 0.394 222   1.524  0.9787
    ##  EMM_F48 - EMM_F3   0.28667 0.394 222   0.728  1.0000
    ##  EMM_F48 - EMM_F66  0.27333 0.394 222   0.694  1.0000
    ##  EMM_F48 - EMM_F49  0.21333 0.394 222   0.542  1.0000
    ##  EMM_F48 - EMM_F5   0.17333 0.394 222   0.440  1.0000
    ##  EMM_F48 - EMM_F70  0.12000 0.394 222   0.305  1.0000
    ##  EMM_F48 - ZAN_F4   0.10000 0.394 222   0.254  1.0000
    ##  EMM_F48 - ZAN_F3   0.08333 0.394 222   0.212  1.0000
    ##  EMM_F48 - EMM_F63  0.04667 0.394 222   0.119  1.0000
    ##  EMM_F48 - EMM_F65  0.01667 0.394 222   0.042  1.0000
    ##  Control - EMM_F34  0.65000 0.394 222   1.651  0.9568
    ##  Control - EMM_F3   0.33667 0.394 222   0.855  1.0000
    ##  Control - EMM_F66  0.32333 0.394 222   0.821  1.0000
    ##  Control - EMM_F49  0.26333 0.394 222   0.669  1.0000
    ##  Control - EMM_F5   0.22333 0.394 222   0.567  1.0000
    ##  Control - EMM_F70  0.17000 0.394 222   0.432  1.0000
    ##  Control - ZAN_F4   0.15000 0.394 222   0.381  1.0000
    ##  Control - ZAN_F3   0.13333 0.394 222   0.339  1.0000
    ##  Control - EMM_F63  0.09667 0.394 222   0.245  1.0000
    ##  Control - EMM_F65  0.06667 0.394 222   0.169  1.0000
    ##  Control - EMM_F48  0.05000 0.394 222   0.127  1.0000
    ##  EMM_F7 - EMM_F34   0.66667 0.394 222   1.693  0.9467
    ##  EMM_F7 - EMM_F3    0.35333 0.394 222   0.897  0.9999
    ##  EMM_F7 - EMM_F66   0.34000 0.394 222   0.863  1.0000
    ##  EMM_F7 - EMM_F49   0.28000 0.394 222   0.711  1.0000
    ##  EMM_F7 - EMM_F5    0.24000 0.394 222   0.609  1.0000
    ##  EMM_F7 - EMM_F70   0.18667 0.394 222   0.474  1.0000
    ##  EMM_F7 - ZAN_F4    0.16667 0.394 222   0.423  1.0000
    ##  EMM_F7 - ZAN_F3    0.15000 0.394 222   0.381  1.0000
    ##  EMM_F7 - EMM_F63   0.11333 0.394 222   0.288  1.0000
    ##  EMM_F7 - EMM_F65   0.08333 0.394 222   0.212  1.0000
    ##  EMM_F7 - EMM_F48   0.06667 0.394 222   0.169  1.0000
    ##  EMM_F7 - Control   0.01667 0.394 222   0.042  1.0000
    ##  EMM_F89 - EMM_F34  0.67667 0.394 222   1.718  0.9399
    ##  EMM_F89 - EMM_F3   0.36333 0.394 222   0.923  0.9999
    ##  EMM_F89 - EMM_F66  0.35000 0.394 222   0.889  0.9999
    ##  EMM_F89 - EMM_F49  0.29000 0.394 222   0.736  1.0000
    ##  EMM_F89 - EMM_F5   0.25000 0.394 222   0.635  1.0000
    ##  EMM_F89 - EMM_F70  0.19667 0.394 222   0.499  1.0000
    ##  EMM_F89 - ZAN_F4   0.17667 0.394 222   0.449  1.0000
    ##  EMM_F89 - ZAN_F3   0.16000 0.394 222   0.406  1.0000
    ##  EMM_F89 - EMM_F63  0.12333 0.394 222   0.313  1.0000
    ##  EMM_F89 - EMM_F65  0.09333 0.394 222   0.237  1.0000
    ##  EMM_F89 - EMM_F48  0.07667 0.394 222   0.195  1.0000
    ##  EMM_F89 - Control  0.02667 0.394 222   0.068  1.0000
    ##  EMM_F89 - EMM_F7   0.01000 0.394 222   0.025  1.0000
    ##  SP_F14 - EMM_F34   0.78333 0.394 222   1.989  0.8287
    ##  SP_F14 - EMM_F3    0.47000 0.394 222   1.194  0.9983
    ##  SP_F14 - EMM_F66   0.45667 0.394 222   1.160  0.9987
    ##  SP_F14 - EMM_F49   0.39667 0.394 222   1.007  0.9998
    ##  SP_F14 - EMM_F5    0.35667 0.394 222   0.906  0.9999
    ##  SP_F14 - EMM_F70   0.30333 0.394 222   0.770  1.0000
    ##  SP_F14 - ZAN_F4    0.28333 0.394 222   0.720  1.0000
    ##  SP_F14 - ZAN_F3    0.26667 0.394 222   0.677  1.0000
    ##  SP_F14 - EMM_F63   0.23000 0.394 222   0.584  1.0000
    ##  SP_F14 - EMM_F65   0.20000 0.394 222   0.508  1.0000
    ##  SP_F14 - EMM_F48   0.18333 0.394 222   0.466  1.0000
    ##  SP_F14 - Control   0.13333 0.394 222   0.339  1.0000
    ##  SP_F14 - EMM_F7    0.11667 0.394 222   0.296  1.0000
    ##  SP_F14 - EMM_F89   0.10667 0.394 222   0.271  1.0000
    ##  EMM_F64 - EMM_F34  1.48333 0.394 222   3.767  0.0189
    ##  EMM_F64 - EMM_F3   1.17000 0.394 222   2.971  0.1901
    ##  EMM_F64 - EMM_F66  1.15667 0.394 222   2.937  0.2056
    ##  EMM_F64 - EMM_F49  1.09667 0.394 222   2.785  0.2856
    ##  EMM_F64 - EMM_F5   1.05667 0.394 222   2.683  0.3481
    ##  EMM_F64 - EMM_F70  1.00333 0.394 222   2.548  0.4408
    ##  EMM_F64 - ZAN_F4   0.98333 0.394 222   2.497  0.4776
    ##  EMM_F64 - ZAN_F3   0.96667 0.394 222   2.455  0.5088
    ##  EMM_F64 - EMM_F63  0.93000 0.394 222   2.362  0.5783
    ##  EMM_F64 - EMM_F65  0.90000 0.394 222   2.286  0.6348
    ##  EMM_F64 - EMM_F48  0.88333 0.394 222   2.243  0.6655
    ##  EMM_F64 - Control  0.83333 0.394 222   2.116  0.7528
    ##  EMM_F64 - EMM_F7   0.81667 0.394 222   2.074  0.7795
    ##  EMM_F64 - EMM_F89  0.80667 0.394 222   2.049  0.7949
    ##  EMM_F64 - SP_F14   0.70000 0.394 222   1.778  0.9217
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F4    0.5167 0.394 222   1.312  0.9951
    ##  Control - EMM_F34   0.3100 0.394 222   0.787  1.0000
    ##  Control - EMM_F89   0.1700 0.394 222   0.432  1.0000
    ##  Control - EMM_F7    0.1567 0.394 222   0.398  1.0000
    ##  Control - EMM_F70   0.0800 0.394 222   0.203  1.0000
    ##  Control - EMM_F65   0.0567 0.394 222   0.144  1.0000
    ##  Control - EMM_F49   0.0467 0.394 222   0.119  1.0000
    ##  Control - EMM_F64   0.0367 0.394 222   0.093  1.0000
    ##  EMM_F3 - Control    0.0933 0.394 222   0.237  1.0000
    ##  SP_F14 - Control    0.1200 0.394 222   0.305  1.0000
    ##  EMM_F66 - Control   0.2033 0.394 222   0.516  1.0000
    ##  EMM_F5 - Control    0.2067 0.394 222   0.525  1.0000
    ##  EMM_F63 - Control   0.2867 0.394 222   0.728  1.0000
    ##  EMM_F48 - Control   0.3200 0.394 222   0.813  1.0000
    ##  ZAN_F3 - Control    0.7033 0.394 222   1.786  0.9189
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F4    0.3033 0.394 222   0.770  1.0000
    ##  SP_F14 - Control    0.0700 0.394 222   0.178  1.0000
    ##  EMM_F7 - Control    0.2233 0.394 222   0.567  1.0000
    ##  ZAN_F3 - Control    0.2300 0.394 222   0.584  1.0000
    ##  EMM_F66 - Control   0.2367 0.394 222   0.601  1.0000
    ##  EMM_F34 - Control   0.2667 0.394 222   0.677  1.0000
    ##  EMM_F3 - Control    0.3167 0.394 222   0.804  1.0000
    ##  EMM_F89 - Control   0.3433 0.394 222   0.872  1.0000
    ##  EMM_F65 - Control   0.3700 0.394 222   0.940  0.9999
    ##  EMM_F5 - Control    0.3767 0.394 222   0.957  0.9999
    ##  EMM_F70 - Control   0.3867 0.394 222   0.982  0.9998
    ##  EMM_F63 - Control   0.3867 0.394 222   0.982  0.9998
    ##  EMM_F64 - Control   0.4333 0.394 222   1.100  0.9993
    ##  EMM_F49 - Control   0.5300 0.394 222   1.346  0.9936
    ##  EMM_F48 - Control   0.5600 0.394 222   1.422  0.9889
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F4    0.2267 0.394 222   0.576  1.0000
    ##  Control - EMM_F34   0.1767 0.394 222   0.449  1.0000
    ##  EMM_F66 - Control   0.0333 0.394 222   0.085  1.0000
    ##  SP_F14 - Control    0.1633 0.394 222   0.415  1.0000
    ##  EMM_F89 - Control   0.2200 0.394 222   0.559  1.0000
    ##  EMM_F3 - Control    0.2433 0.394 222   0.618  1.0000
    ##  EMM_F63 - Control   0.3267 0.394 222   0.830  1.0000
    ##  EMM_F49 - Control   0.3300 0.394 222   0.838  1.0000
    ##  EMM_F7 - Control    0.3700 0.394 222   0.940  0.9999
    ##  ZAN_F3 - Control    0.3967 0.394 222   1.007  0.9998
    ##  EMM_F70 - Control   0.4567 0.394 222   1.160  0.9987
    ##  EMM_F64 - Control   0.4733 0.394 222   1.202  0.9981
    ##  EMM_F65 - Control   0.4767 0.394 222   1.210  0.9980
    ##  EMM_F5 - Control    0.5733 0.394 222   1.456  0.9861
    ##  EMM_F48 - Control   0.6867 0.394 222   1.744  0.9325
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65   0.1567 0.394 222   0.398  1.0000
    ##  Control - EMM_F34   0.1100 0.394 222   0.279  1.0000
    ##  Control - EMM_F3    0.0633 0.394 222   0.161  1.0000
    ##  SP_F14 - Control    0.0900 0.394 222   0.229  1.0000
    ##  EMM_F7 - Control    0.1700 0.394 222   0.432  1.0000
    ##  EMM_F5 - Control    0.1900 0.394 222   0.483  1.0000
    ##  EMM_F49 - Control   0.2700 0.394 222   0.686  1.0000
    ##  ZAN_F4 - Control    0.2933 0.394 222   0.745  1.0000
    ##  EMM_F66 - Control   0.3000 0.394 222   0.762  1.0000
    ##  EMM_F89 - Control   0.4267 0.394 222   1.084  0.9994
    ##  EMM_F64 - Control   0.4633 0.394 222   1.177  0.9985
    ##  EMM_F48 - Control   0.4833 0.394 222   1.227  0.9976
    ##  EMM_F70 - Control   0.5833 0.394 222   1.481  0.9836
    ##  EMM_F63 - Control   0.6267 0.394 222   1.591  0.9685
    ##  ZAN_F3 - Control    0.6967 0.394 222   1.769  0.9246
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34   0.1767 0.394 222   0.449  1.0000
    ##  Control - ZAN_F4    0.1200 0.394 222   0.305  1.0000
    ##  Control - EMM_F5    0.0533 0.394 222   0.135  1.0000
    ##  EMM_F49 - Control   0.0467 0.394 222   0.119  1.0000
    ##  EMM_F3 - Control    0.0500 0.394 222   0.127  1.0000
    ##  EMM_F65 - Control   0.1700 0.394 222   0.432  1.0000
    ##  EMM_F63 - Control   0.1800 0.394 222   0.457  1.0000
    ##  EMM_F66 - Control   0.1833 0.394 222   0.466  1.0000
    ##  EMM_F70 - Control   0.2600 0.394 222   0.660  1.0000
    ##  EMM_F48 - Control   0.2600 0.394 222   0.660  1.0000
    ##  EMM_F64 - Control   0.2700 0.394 222   0.686  1.0000
    ##  EMM_F7 - Control    0.2867 0.394 222   0.728  1.0000
    ##  EMM_F89 - Control   0.4300 0.394 222   1.092  0.9994
    ##  ZAN_F3 - Control    0.5667 0.394 222   1.439  0.9876
    ##  SP_F14 - Control    0.6233 0.394 222   1.583  0.9699
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34   0.4000 0.394 222   1.016  0.9997
    ##  Control - EMM_F3    0.2600 0.394 222   0.660  1.0000
    ##  Control - EMM_F5    0.2000 0.394 222   0.508  1.0000
    ##  Control - EMM_F48   0.1767 0.394 222   0.449  1.0000
    ##  Control - ZAN_F4    0.1500 0.394 222   0.381  1.0000
    ##  EMM_F63 - Control   0.0267 0.394 222   0.068  1.0000
    ##  EMM_F7 - Control    0.0733 0.394 222   0.186  1.0000
    ##  SP_F14 - Control    0.1200 0.394 222   0.305  1.0000
    ##  EMM_F64 - Control   0.1400 0.394 222   0.356  1.0000
    ##  EMM_F66 - Control   0.1767 0.394 222   0.449  1.0000
    ##  EMM_F70 - Control   0.2700 0.394 222   0.686  1.0000
    ##  EMM_F65 - Control   0.3233 0.394 222   0.821  1.0000
    ##  EMM_F49 - Control   0.3400 0.394 222   0.863  1.0000
    ##  EMM_F89 - Control   0.5200 0.394 222   1.321  0.9948
    ##  ZAN_F3 - Control    0.9400 0.394 222   2.387  0.5593
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34   0.6500 0.394 222   1.651  0.9568
    ##  Control - EMM_F3    0.3367 0.394 222   0.855  1.0000
    ##  Control - EMM_F66   0.3233 0.394 222   0.821  1.0000
    ##  Control - EMM_F49   0.2633 0.394 222   0.669  1.0000
    ##  Control - EMM_F5    0.2233 0.394 222   0.567  1.0000
    ##  Control - EMM_F70   0.1700 0.394 222   0.432  1.0000
    ##  Control - ZAN_F4    0.1500 0.394 222   0.381  1.0000
    ##  Control - ZAN_F3    0.1333 0.394 222   0.339  1.0000
    ##  Control - EMM_F63   0.0967 0.394 222   0.245  1.0000
    ##  Control - EMM_F65   0.0667 0.394 222   0.169  1.0000
    ##  Control - EMM_F48   0.0500 0.394 222   0.127  1.0000
    ##  EMM_F7 - Control    0.0167 0.394 222   0.042  1.0000
    ##  EMM_F89 - Control   0.0267 0.394 222   0.068  1.0000
    ##  SP_F14 - Control    0.1333 0.394 222   0.339  1.0000
    ##  EMM_F64 - Control   0.8333 0.394 222   2.116  0.7528
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## $emmeans
    ## distance_to_yeast = 11:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  ZAN_F3   1.22000 0.286  2 -0.00876    2.449  a    
    ##  EMM_F5   1.06333 0.286  2 -0.16542    2.292  ab   
    ##  EMM_F63  0.91333 0.286  2 -0.31542    2.142  ab   
    ##  EMM_F48  0.78667 0.286  2 -0.44209    2.015  ab   
    ##  EMM_F3   0.63000 0.286  2 -0.59876    1.859  ab   
    ##  EMM_F65  0.63000 0.286  2 -0.59876    1.859  ab   
    ##  EMM_F7   0.57333 0.286  2 -0.65542    1.802  ab   
    ##  EMM_F66  0.56333 0.286  2 -0.66542    1.792  ab   
    ##  EMM_F49  0.54000 0.350  2 -0.96492    2.045  ab   
    ##  EMM_F64  0.43333 0.286  2 -0.79542    1.662  ab   
    ##  EMM_F89  0.42000 0.350  2 -1.08492    1.925  ab   
    ##  ZAN_F4   0.39333 0.286  2 -0.83542    1.622  ab   
    ##  EMM_F70  0.15667 0.286  2 -1.07209    1.385  ab   
    ##  Control  0.08333 0.286  2 -1.14542    1.312  ab   
    ##  EMM_F34  0.00333 0.286  2 -1.22542    1.232  ab   
    ##  SP_F14  -0.23333 0.286  2 -1.46209    0.995   b   
    ## 
    ## distance_to_yeast = 17:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  EMM_F65  1.21000 0.286  2 -0.01876    2.439  a    
    ##  EMM_F48  1.07333 0.286  2 -0.15542    2.302  a    
    ##  EMM_F5   0.88333 0.286  2 -0.34542    2.112  a    
    ##  EMM_F63  0.88000 0.286  2 -0.34876    2.109  a    
    ##  EMM_F7   0.66667 0.286  2 -0.56209    1.895  a    
    ##  EMM_F64  0.60333 0.286  2 -0.62542    1.832  a    
    ##  EMM_F66  0.49333 0.286  2 -0.73542    1.722  a    
    ##  ZAN_F4   0.48333 0.286  2 -0.74542    1.712  a    
    ##  EMM_F70  0.48000 0.286  2 -0.74876    1.709  a    
    ##  EMM_F49  0.47500 0.350  2 -1.02992    1.980  a    
    ##  SP_F14   0.45667 0.286  2 -0.77209    1.685  a    
    ##  ZAN_F3   0.44333 0.286  2 -0.78542    1.672  a    
    ##  EMM_F3   0.36333 0.286  2 -0.86542    1.592  a    
    ##  Control  0.29667 0.286  2 -0.93209    1.525  a    
    ##  EMM_F34  0.24000 0.286  2 -0.98876    1.469  a    
    ##  EMM_F89  0.17500 0.350  2 -1.32992    1.680  a    
    ## 
    ## distance_to_yeast = 25:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  EMM_F48  0.97000 0.286  2 -0.25876    2.199  a    
    ##  EMM_F7   0.89667 0.286  2 -0.33209    2.125  a    
    ##  EMM_F63  0.88333 0.286  2 -0.34542    2.112  a    
    ##  EMM_F65  0.83667 0.286  2 -0.39209    2.065  a    
    ##  EMM_F64  0.83333 0.286  2 -0.39542    2.062  a    
    ##  EMM_F5   0.79000 0.286  2 -0.43876    2.019  a    
    ##  SP_F14   0.71667 0.286  2 -0.51209    1.945  a    
    ##  ZAN_F4   0.69667 0.286  2 -0.53209    1.925  a    
    ##  EMM_F3   0.45000 0.286  2 -0.77876    1.679  a    
    ##  ZAN_F3   0.44000 0.286  2 -0.78876    1.669  a    
    ##  EMM_F70  0.42667 0.286  2 -0.80209    1.655  a    
    ##  EMM_F89  0.32000 0.350  2 -1.18492    1.825  a    
    ##  EMM_F49  0.31500 0.350  2 -1.18992    1.820  a    
    ##  Control  0.30000 0.286  2 -0.92876    1.529  a    
    ##  EMM_F66  0.23333 0.286  2 -0.99542    1.462  a    
    ##  EMM_F34  0.22667 0.286  2 -1.00209    1.455  a    
    ## 
    ## distance_to_yeast = 32:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  EMM_F63  1.02000 0.286  2 -0.20876    2.249  a    
    ##  EMM_F64  0.98333 0.286  2 -0.24542    2.212  a    
    ##  ZAN_F3   0.94000 0.286  2 -0.28876    2.169  a    
    ##  EMM_F48  0.82000 0.286  2 -0.40876    2.049  a    
    ##  ZAN_F4   0.81333 0.286  2 -0.41542    2.042  a    
    ##  SP_F14   0.55333 0.286  2 -0.67542    1.782  a    
    ##  EMM_F70  0.54667 0.286  2 -0.68209    1.775  a    
    ##  EMM_F66  0.54333 0.286  2 -0.68542    1.772  a    
    ##  EMM_F5   0.53667 0.286  2 -0.69209    1.765  a    
    ##  EMM_F7   0.53000 0.286  2 -0.69876    1.759  a    
    ##  EMM_F34  0.43333 0.286  2 -0.79542    1.662  a    
    ##  EMM_F3   0.40000 0.286  2 -0.82876    1.629  a    
    ##  EMM_F49  0.39500 0.350  2 -1.10992    1.900  a    
    ##  EMM_F89  0.36000 0.350  2 -1.14492    1.865  a    
    ##  Control  0.31667 0.286  2 -0.91209    1.545  a    
    ##  EMM_F65  0.22667 0.286  2 -1.00209    1.455  a    
    ## 
    ## distance_to_yeast = 41:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  ZAN_F3   1.07000 0.286  2 -0.15876    2.299  a    
    ##  EMM_F63  1.01667 0.286  2 -0.21209    2.245  a    
    ##  SP_F14   1.00667 0.286  2 -0.22209    2.235  a    
    ##  EMM_F65  0.85000 0.286  2 -0.37876    2.079  a    
    ##  EMM_F7   0.74667 0.286  2 -0.48209    1.975  a    
    ##  EMM_F64  0.74333 0.286  2 -0.48542    1.972  a    
    ##  EMM_F89  0.68000 0.350  2 -0.82492    2.185  a    
    ##  ZAN_F4   0.55667 0.286  2 -0.67209    1.785  a    
    ##  EMM_F34  0.54000 0.286  2 -0.68876    1.769  a    
    ##  EMM_F5   0.52000 0.286  2 -0.70876    1.749  a    
    ##  EMM_F3   0.50333 0.286  2 -0.72542    1.732  a    
    ##  EMM_F48  0.46000 0.286  2 -0.76876    1.689  a    
    ##  Control  0.42000 0.286  2 -0.80876    1.649  a    
    ##  EMM_F66  0.39667 0.286  2 -0.83209    1.625  a    
    ##  EMM_F70  0.31667 0.286  2 -0.91209    1.545  a    
    ##  EMM_F49  0.04500 0.350  2 -1.45992    1.550  a    
    ## 
    ## distance_to_yeast = 48:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  ZAN_F3   1.47667 0.286  2  0.24791    2.705  a    
    ##  EMM_F89  1.24500 0.350  2 -0.25992    2.750  a    
    ##  EMM_F65  1.09333 0.286  2 -0.13542    2.322  a    
    ##  EMM_F64  0.91667 0.286  2 -0.31209    2.145  a    
    ##  EMM_F66  0.83000 0.286  2 -0.39876    2.059  a    
    ##  SP_F14   0.80667 0.286  2 -0.42209    2.035  a    
    ##  ZAN_F4   0.78000 0.286  2 -0.44876    2.009  a    
    ##  EMM_F5   0.71333 0.286  2 -0.51542    1.942  a    
    ##  EMM_F34  0.65667 0.286  2 -0.57209    1.885  a    
    ##  EMM_F63  0.60333 0.286  2 -0.62542    1.832  a    
    ##  EMM_F7   0.60000 0.286  2 -0.62876    1.829  a    
    ##  EMM_F3   0.54333 0.286  2 -0.68542    1.772  a    
    ##  EMM_F70  0.52333 0.286  2 -0.70542    1.752  a    
    ##  Control  0.36667 0.286  2 -0.86209    1.595  a    
    ##  EMM_F48  0.33000 0.286  2 -0.89876    1.559  a    
    ##  EMM_F49  0.12500 0.350  2 -1.37992    1.630  a    
    ## 
    ## distance_to_yeast = 55:
    ##  Yeast     emmean    SE df lower.CL upper.CL .group
    ##  EMM_F64  1.36000 0.286  2  0.13124    2.589  a    
    ##  SP_F14   1.00333 0.286  2 -0.22542    2.232  a    
    ##  EMM_F7   0.89333 0.286  2 -0.33542    2.122  a    
    ##  EMM_F65  0.86667 0.286  2 -0.36209    2.095  a    
    ##  EMM_F5   0.85667 0.286  2 -0.37209    2.085  a    
    ##  ZAN_F4   0.84667 0.286  2 -0.38209    2.075  a    
    ##  EMM_F66  0.76000 0.286  2 -0.46876    1.989  a    
    ##  EMM_F70  0.71000 0.286  2 -0.51876    1.939  a    
    ##  EMM_F48  0.69333 0.286  2 -0.53542    1.922  a    
    ##  EMM_F63  0.68333 0.286  2 -0.54542    1.912  a    
    ##  Control  0.52333 0.286  2 -0.70542    1.752  a    
    ##  EMM_F3   0.49000 0.286  2 -0.73876    1.719  a    
    ##  EMM_F89  0.27000 0.350  2 -1.23492    1.775  a    
    ##  ZAN_F3   0.27000 0.286  2 -0.95876    1.499  a    
    ##  EMM_F34  0.24667 0.286  2 -0.98209    1.475  a    
    ##  EMM_F49  0.12500 0.350  2 -1.37992    1.630  a    
    ## 
    ## Degrees-of-freedom method: containment 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F34 - SP_F14   0.23667 0.404 208   0.586  1.0000
    ##  Control - SP_F14   0.31667 0.404 208   0.784  1.0000
    ##  Control - EMM_F34  0.08000 0.404 208   0.198  1.0000
    ##  EMM_F70 - SP_F14   0.39000 0.404 208   0.966  0.9999
    ##  EMM_F70 - EMM_F34  0.15333 0.404 208   0.380  1.0000
    ##  EMM_F70 - Control  0.07333 0.404 208   0.182  1.0000
    ##  ZAN_F4 - SP_F14    0.62667 0.404 208   1.552  0.9747
    ##  ZAN_F4 - EMM_F34   0.39000 0.404 208   0.966  0.9999
    ##  ZAN_F4 - Control   0.31000 0.404 208   0.768  1.0000
    ##  ZAN_F4 - EMM_F70   0.23667 0.404 208   0.586  1.0000
    ##  EMM_F89 - SP_F14   0.65333 0.452 208   1.447  0.9868
    ##  EMM_F89 - EMM_F34  0.41667 0.452 208   0.923  0.9999
    ##  EMM_F89 - Control  0.33667 0.452 208   0.746  1.0000
    ##  EMM_F89 - EMM_F70  0.26333 0.452 208   0.583  1.0000
    ##  EMM_F89 - ZAN_F4   0.02667 0.452 208   0.059  1.0000
    ##  EMM_F64 - SP_F14   0.66667 0.404 208   1.651  0.9567
    ##  EMM_F64 - EMM_F34  0.43000 0.404 208   1.065  0.9995
    ##  EMM_F64 - Control  0.35000 0.404 208   0.867  1.0000
    ##  EMM_F64 - EMM_F70  0.27667 0.404 208   0.685  1.0000
    ##  EMM_F64 - ZAN_F4   0.04000 0.404 208   0.099  1.0000
    ##  EMM_F64 - EMM_F89  0.01333 0.452 208   0.030  1.0000
    ##  EMM_F49 - SP_F14   0.77333 0.452 208   1.713  0.9414
    ##  EMM_F49 - EMM_F34  0.53667 0.452 208   1.189  0.9983
    ##  EMM_F49 - Control  0.45667 0.452 208   1.011  0.9997
    ##  EMM_F49 - EMM_F70  0.38333 0.452 208   0.849  1.0000
    ##  EMM_F49 - ZAN_F4   0.14667 0.452 208   0.325  1.0000
    ##  EMM_F49 - EMM_F89  0.12000 0.495 208   0.243  1.0000
    ##  EMM_F49 - EMM_F64  0.10667 0.452 208   0.236  1.0000
    ##  EMM_F66 - SP_F14   0.79667 0.404 208   1.973  0.8374
    ##  EMM_F66 - EMM_F34  0.56000 0.404 208   1.387  0.9914
    ##  EMM_F66 - Control  0.48000 0.404 208   1.188  0.9983
    ##  EMM_F66 - EMM_F70  0.40667 0.404 208   1.007  0.9998
    ##  EMM_F66 - ZAN_F4   0.17000 0.404 208   0.421  1.0000
    ##  EMM_F66 - EMM_F89  0.14333 0.452 208   0.317  1.0000
    ##  EMM_F66 - EMM_F64  0.13000 0.404 208   0.322  1.0000
    ##  EMM_F66 - EMM_F49  0.02333 0.452 208   0.052  1.0000
    ##  EMM_F7 - SP_F14    0.80667 0.404 208   1.997  0.8241
    ##  EMM_F7 - EMM_F34   0.57000 0.404 208   1.411  0.9897
    ##  EMM_F7 - Control   0.49000 0.404 208   1.213  0.9979
    ##  EMM_F7 - EMM_F70   0.41667 0.404 208   1.032  0.9997
    ##  EMM_F7 - ZAN_F4    0.18000 0.404 208   0.446  1.0000
    ##  EMM_F7 - EMM_F89   0.15333 0.452 208   0.340  1.0000
    ##  EMM_F7 - EMM_F64   0.14000 0.404 208   0.347  1.0000
    ##  EMM_F7 - EMM_F49   0.03333 0.452 208   0.074  1.0000
    ##  EMM_F7 - EMM_F66   0.01000 0.404 208   0.025  1.0000
    ##  EMM_F65 - SP_F14   0.86333 0.404 208   2.138  0.7386
    ##  EMM_F65 - EMM_F34  0.62667 0.404 208   1.552  0.9747
    ##  EMM_F65 - Control  0.54667 0.404 208   1.354  0.9932
    ##  EMM_F65 - EMM_F70  0.47333 0.404 208   1.172  0.9986
    ##  EMM_F65 - ZAN_F4   0.23667 0.404 208   0.586  1.0000
    ##  EMM_F65 - EMM_F89  0.21000 0.452 208   0.465  1.0000
    ##  EMM_F65 - EMM_F64  0.19667 0.404 208   0.487  1.0000
    ##  EMM_F65 - EMM_F49  0.09000 0.452 208   0.199  1.0000
    ##  EMM_F65 - EMM_F66  0.06667 0.404 208   0.165  1.0000
    ##  EMM_F65 - EMM_F7   0.05667 0.404 208   0.140  1.0000
    ##  EMM_F3 - SP_F14    0.86333 0.404 208   2.138  0.7386
    ##  EMM_F3 - EMM_F34   0.62667 0.404 208   1.552  0.9747
    ##  EMM_F3 - Control   0.54667 0.404 208   1.354  0.9932
    ##  EMM_F3 - EMM_F70   0.47333 0.404 208   1.172  0.9986
    ##  EMM_F3 - ZAN_F4    0.23667 0.404 208   0.586  1.0000
    ##  EMM_F3 - EMM_F89   0.21000 0.452 208   0.465  1.0000
    ##  EMM_F3 - EMM_F64   0.19667 0.404 208   0.487  1.0000
    ##  EMM_F3 - EMM_F49   0.09000 0.452 208   0.199  1.0000
    ##  EMM_F3 - EMM_F66   0.06667 0.404 208   0.165  1.0000
    ##  EMM_F3 - EMM_F7    0.05667 0.404 208   0.140  1.0000
    ##  EMM_F3 - EMM_F65   0.00000 0.404 208   0.000  1.0000
    ##  EMM_F48 - SP_F14   1.02000 0.404 208   2.526  0.4572
    ##  EMM_F48 - EMM_F34  0.78333 0.404 208   1.940  0.8543
    ##  EMM_F48 - Control  0.70333 0.404 208   1.741  0.9331
    ##  EMM_F48 - EMM_F70  0.63000 0.404 208   1.560  0.9735
    ##  EMM_F48 - ZAN_F4   0.39333 0.404 208   0.974  0.9998
    ##  EMM_F48 - EMM_F89  0.36667 0.452 208   0.812  1.0000
    ##  EMM_F48 - EMM_F64  0.35333 0.404 208   0.875  1.0000
    ##  EMM_F48 - EMM_F49  0.24667 0.452 208   0.546  1.0000
    ##  EMM_F48 - EMM_F66  0.22333 0.404 208   0.553  1.0000
    ##  EMM_F48 - EMM_F7   0.21333 0.404 208   0.528  1.0000
    ##  EMM_F48 - EMM_F65  0.15667 0.404 208   0.388  1.0000
    ##  EMM_F48 - EMM_F3   0.15667 0.404 208   0.388  1.0000
    ##  EMM_F63 - SP_F14   1.14667 0.404 208   2.839  0.2558
    ##  EMM_F63 - EMM_F34  0.91000 0.404 208   2.253  0.6583
    ##  EMM_F63 - Control  0.83000 0.404 208   2.055  0.7908
    ##  EMM_F63 - EMM_F70  0.75667 0.404 208   1.874  0.8849
    ##  EMM_F63 - ZAN_F4   0.52000 0.404 208   1.288  0.9960
    ##  EMM_F63 - EMM_F89  0.49333 0.452 208   1.093  0.9994
    ##  EMM_F63 - EMM_F64  0.48000 0.404 208   1.188  0.9983
    ##  EMM_F63 - EMM_F49  0.37333 0.452 208   0.827  1.0000
    ##  EMM_F63 - EMM_F66  0.35000 0.404 208   0.867  1.0000
    ##  EMM_F63 - EMM_F7   0.34000 0.404 208   0.842  1.0000
    ##  EMM_F63 - EMM_F65  0.28333 0.404 208   0.702  1.0000
    ##  EMM_F63 - EMM_F3   0.28333 0.404 208   0.702  1.0000
    ##  EMM_F63 - EMM_F48  0.12667 0.404 208   0.314  1.0000
    ##  EMM_F5 - SP_F14    1.29667 0.404 208   3.211  0.1046
    ##  EMM_F5 - EMM_F34   1.06000 0.404 208   2.625  0.3876
    ##  EMM_F5 - Control   0.98000 0.404 208   2.427  0.5301
    ##  EMM_F5 - EMM_F70   0.90667 0.404 208   2.245  0.6643
    ##  EMM_F5 - ZAN_F4    0.67000 0.404 208   1.659  0.9548
    ##  EMM_F5 - EMM_F89   0.64333 0.452 208   1.425  0.9887
    ##  EMM_F5 - EMM_F64   0.63000 0.404 208   1.560  0.9735
    ##  EMM_F5 - EMM_F49   0.52333 0.452 208   1.159  0.9987
    ##  EMM_F5 - EMM_F66   0.50000 0.404 208   1.238  0.9974
    ##  EMM_F5 - EMM_F7    0.49000 0.404 208   1.213  0.9979
    ##  EMM_F5 - EMM_F65   0.43333 0.404 208   1.073  0.9995
    ##  EMM_F5 - EMM_F3    0.43333 0.404 208   1.073  0.9995
    ##  EMM_F5 - EMM_F48   0.27667 0.404 208   0.685  1.0000
    ##  EMM_F5 - EMM_F63   0.15000 0.404 208   0.371  1.0000
    ##  ZAN_F3 - SP_F14    1.45333 0.404 208   3.598  0.0334
    ##  ZAN_F3 - EMM_F34   1.21667 0.404 208   3.012  0.1730
    ##  ZAN_F3 - Control   1.13667 0.404 208   2.814  0.2694
    ##  ZAN_F3 - EMM_F70   1.06333 0.404 208   2.633  0.3820
    ##  ZAN_F3 - ZAN_F4    0.82667 0.404 208   2.047  0.7957
    ##  ZAN_F3 - EMM_F89   0.80000 0.452 208   1.772  0.9235
    ##  ZAN_F3 - EMM_F64   0.78667 0.404 208   1.948  0.8502
    ##  ZAN_F3 - EMM_F49   0.68000 0.452 208   1.506  0.9808
    ##  ZAN_F3 - EMM_F66   0.65667 0.404 208   1.626  0.9619
    ##  ZAN_F3 - EMM_F7    0.64667 0.404 208   1.601  0.9666
    ##  ZAN_F3 - EMM_F65   0.59000 0.404 208   1.461  0.9856
    ##  ZAN_F3 - EMM_F3    0.59000 0.404 208   1.461  0.9856
    ##  ZAN_F3 - EMM_F48   0.43333 0.404 208   1.073  0.9995
    ##  ZAN_F3 - EMM_F63   0.30667 0.404 208   0.759  1.0000
    ##  ZAN_F3 - EMM_F5    0.15667 0.404 208   0.388  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F34 - EMM_F89  0.06500 0.452 208   0.144  1.0000
    ##  Control - EMM_F89  0.12167 0.452 208   0.269  1.0000
    ##  Control - EMM_F34  0.05667 0.404 208   0.140  1.0000
    ##  EMM_F3 - EMM_F89   0.18833 0.452 208   0.417  1.0000
    ##  EMM_F3 - EMM_F34   0.12333 0.404 208   0.305  1.0000
    ##  EMM_F3 - Control   0.06667 0.404 208   0.165  1.0000
    ##  ZAN_F3 - EMM_F89   0.26833 0.452 208   0.594  1.0000
    ##  ZAN_F3 - EMM_F34   0.20333 0.404 208   0.503  1.0000
    ##  ZAN_F3 - Control   0.14667 0.404 208   0.363  1.0000
    ##  ZAN_F3 - EMM_F3    0.08000 0.404 208   0.198  1.0000
    ##  SP_F14 - EMM_F89   0.28167 0.452 208   0.624  1.0000
    ##  SP_F14 - EMM_F34   0.21667 0.404 208   0.536  1.0000
    ##  SP_F14 - Control   0.16000 0.404 208   0.396  1.0000
    ##  SP_F14 - EMM_F3    0.09333 0.404 208   0.231  1.0000
    ##  SP_F14 - ZAN_F3    0.01333 0.404 208   0.033  1.0000
    ##  EMM_F49 - EMM_F89  0.30000 0.495 208   0.607  1.0000
    ##  EMM_F49 - EMM_F34  0.23500 0.452 208   0.520  1.0000
    ##  EMM_F49 - Control  0.17833 0.452 208   0.395  1.0000
    ##  EMM_F49 - EMM_F3   0.11167 0.452 208   0.247  1.0000
    ##  EMM_F49 - ZAN_F3   0.03167 0.452 208   0.070  1.0000
    ##  EMM_F49 - SP_F14   0.01833 0.452 208   0.041  1.0000
    ##  EMM_F70 - EMM_F89  0.30500 0.452 208   0.675  1.0000
    ##  EMM_F70 - EMM_F34  0.24000 0.404 208   0.594  1.0000
    ##  EMM_F70 - Control  0.18333 0.404 208   0.454  1.0000
    ##  EMM_F70 - EMM_F3   0.11667 0.404 208   0.289  1.0000
    ##  EMM_F70 - ZAN_F3   0.03667 0.404 208   0.091  1.0000
    ##  EMM_F70 - SP_F14   0.02333 0.404 208   0.058  1.0000
    ##  EMM_F70 - EMM_F49  0.00500 0.452 208   0.011  1.0000
    ##  ZAN_F4 - EMM_F89   0.30833 0.452 208   0.683  1.0000
    ##  ZAN_F4 - EMM_F34   0.24333 0.404 208   0.602  1.0000
    ##  ZAN_F4 - Control   0.18667 0.404 208   0.462  1.0000
    ##  ZAN_F4 - EMM_F3    0.12000 0.404 208   0.297  1.0000
    ##  ZAN_F4 - ZAN_F3    0.04000 0.404 208   0.099  1.0000
    ##  ZAN_F4 - SP_F14    0.02667 0.404 208   0.066  1.0000
    ##  ZAN_F4 - EMM_F49   0.00833 0.452 208   0.018  1.0000
    ##  ZAN_F4 - EMM_F70   0.00333 0.404 208   0.008  1.0000
    ##  EMM_F66 - EMM_F89  0.31833 0.452 208   0.705  1.0000
    ##  EMM_F66 - EMM_F34  0.25333 0.404 208   0.627  1.0000
    ##  EMM_F66 - Control  0.19667 0.404 208   0.487  1.0000
    ##  EMM_F66 - EMM_F3   0.13000 0.404 208   0.322  1.0000
    ##  EMM_F66 - ZAN_F3   0.05000 0.404 208   0.124  1.0000
    ##  EMM_F66 - SP_F14   0.03667 0.404 208   0.091  1.0000
    ##  EMM_F66 - EMM_F49  0.01833 0.452 208   0.041  1.0000
    ##  EMM_F66 - EMM_F70  0.01333 0.404 208   0.033  1.0000
    ##  EMM_F66 - ZAN_F4   0.01000 0.404 208   0.025  1.0000
    ##  EMM_F64 - EMM_F89  0.42833 0.452 208   0.949  0.9999
    ##  EMM_F64 - EMM_F34  0.36333 0.404 208   0.900  0.9999
    ##  EMM_F64 - Control  0.30667 0.404 208   0.759  1.0000
    ##  EMM_F64 - EMM_F3   0.24000 0.404 208   0.594  1.0000
    ##  EMM_F64 - ZAN_F3   0.16000 0.404 208   0.396  1.0000
    ##  EMM_F64 - SP_F14   0.14667 0.404 208   0.363  1.0000
    ##  EMM_F64 - EMM_F49  0.12833 0.452 208   0.284  1.0000
    ##  EMM_F64 - EMM_F70  0.12333 0.404 208   0.305  1.0000
    ##  EMM_F64 - ZAN_F4   0.12000 0.404 208   0.297  1.0000
    ##  EMM_F64 - EMM_F66  0.11000 0.404 208   0.272  1.0000
    ##  EMM_F7 - EMM_F89   0.49167 0.452 208   1.089  0.9994
    ##  EMM_F7 - EMM_F34   0.42667 0.404 208   1.056  0.9996
    ##  EMM_F7 - Control   0.37000 0.404 208   0.916  0.9999
    ##  EMM_F7 - EMM_F3    0.30333 0.404 208   0.751  1.0000
    ##  EMM_F7 - ZAN_F3    0.22333 0.404 208   0.553  1.0000
    ##  EMM_F7 - SP_F14    0.21000 0.404 208   0.520  1.0000
    ##  EMM_F7 - EMM_F49   0.19167 0.452 208   0.424  1.0000
    ##  EMM_F7 - EMM_F70   0.18667 0.404 208   0.462  1.0000
    ##  EMM_F7 - ZAN_F4    0.18333 0.404 208   0.454  1.0000
    ##  EMM_F7 - EMM_F66   0.17333 0.404 208   0.429  1.0000
    ##  EMM_F7 - EMM_F64   0.06333 0.404 208   0.157  1.0000
    ##  EMM_F63 - EMM_F89  0.70500 0.452 208   1.561  0.9733
    ##  EMM_F63 - EMM_F34  0.64000 0.404 208   1.585  0.9695
    ##  EMM_F63 - Control  0.58333 0.404 208   1.444  0.9871
    ##  EMM_F63 - EMM_F3   0.51667 0.404 208   1.279  0.9962
    ##  EMM_F63 - ZAN_F3   0.43667 0.404 208   1.081  0.9994
    ##  EMM_F63 - SP_F14   0.42333 0.404 208   1.048  0.9996
    ##  EMM_F63 - EMM_F49  0.40500 0.452 208   0.897  0.9999
    ##  EMM_F63 - EMM_F70  0.40000 0.404 208   0.990  0.9998
    ##  EMM_F63 - ZAN_F4   0.39667 0.404 208   0.982  0.9998
    ##  EMM_F63 - EMM_F66  0.38667 0.404 208   0.957  0.9999
    ##  EMM_F63 - EMM_F64  0.27667 0.404 208   0.685  1.0000
    ##  EMM_F63 - EMM_F7   0.21333 0.404 208   0.528  1.0000
    ##  EMM_F5 - EMM_F89   0.70833 0.452 208   1.569  0.9721
    ##  EMM_F5 - EMM_F34   0.64333 0.404 208   1.593  0.9681
    ##  EMM_F5 - Control   0.58667 0.404 208   1.453  0.9863
    ##  EMM_F5 - EMM_F3    0.52000 0.404 208   1.288  0.9960
    ##  EMM_F5 - ZAN_F3    0.44000 0.404 208   1.089  0.9994
    ##  EMM_F5 - SP_F14    0.42667 0.404 208   1.056  0.9996
    ##  EMM_F5 - EMM_F49   0.40833 0.452 208   0.904  0.9999
    ##  EMM_F5 - EMM_F70   0.40333 0.404 208   0.999  0.9998
    ##  EMM_F5 - ZAN_F4    0.40000 0.404 208   0.990  0.9998
    ##  EMM_F5 - EMM_F66   0.39000 0.404 208   0.966  0.9999
    ##  EMM_F5 - EMM_F64   0.28000 0.404 208   0.693  1.0000
    ##  EMM_F5 - EMM_F7    0.21667 0.404 208   0.536  1.0000
    ##  EMM_F5 - EMM_F63   0.00333 0.404 208   0.008  1.0000
    ##  EMM_F48 - EMM_F89  0.89833 0.452 208   1.989  0.8284
    ##  EMM_F48 - EMM_F34  0.83333 0.404 208   2.063  0.7858
    ##  EMM_F48 - Control  0.77667 0.404 208   1.923  0.8623
    ##  EMM_F48 - EMM_F3   0.71000 0.404 208   1.758  0.9280
    ##  EMM_F48 - ZAN_F3   0.63000 0.404 208   1.560  0.9735
    ##  EMM_F48 - SP_F14   0.61667 0.404 208   1.527  0.9782
    ##  EMM_F48 - EMM_F49  0.59833 0.452 208   1.325  0.9946
    ##  EMM_F48 - EMM_F70  0.59333 0.404 208   1.469  0.9848
    ##  EMM_F48 - ZAN_F4   0.59000 0.404 208   1.461  0.9856
    ##  EMM_F48 - EMM_F66  0.58000 0.404 208   1.436  0.9878
    ##  EMM_F48 - EMM_F64  0.47000 0.404 208   1.164  0.9987
    ##  EMM_F48 - EMM_F7   0.40667 0.404 208   1.007  0.9998
    ##  EMM_F48 - EMM_F63  0.19333 0.404 208   0.479  1.0000
    ##  EMM_F48 - EMM_F5   0.19000 0.404 208   0.470  1.0000
    ##  EMM_F65 - EMM_F89  1.03500 0.452 208   2.292  0.6300
    ##  EMM_F65 - EMM_F34  0.97000 0.404 208   2.402  0.5485
    ##  EMM_F65 - Control  0.91333 0.404 208   2.261  0.6524
    ##  EMM_F65 - EMM_F3   0.84667 0.404 208   2.096  0.7653
    ##  EMM_F65 - ZAN_F3   0.76667 0.404 208   1.898  0.8739
    ##  EMM_F65 - SP_F14   0.75333 0.404 208   1.865  0.8884
    ##  EMM_F65 - EMM_F49  0.73500 0.452 208   1.628  0.9615
    ##  EMM_F65 - EMM_F70  0.73000 0.404 208   1.807  0.9111
    ##  EMM_F65 - ZAN_F4   0.72667 0.404 208   1.799  0.9141
    ##  EMM_F65 - EMM_F66  0.71667 0.404 208   1.774  0.9226
    ##  EMM_F65 - EMM_F64  0.60667 0.404 208   1.502  0.9812
    ##  EMM_F65 - EMM_F7   0.54333 0.404 208   1.345  0.9936
    ##  EMM_F65 - EMM_F63  0.33000 0.404 208   0.817  1.0000
    ##  EMM_F65 - EMM_F5   0.32667 0.404 208   0.809  1.0000
    ##  EMM_F65 - EMM_F48  0.13667 0.404 208   0.338  1.0000
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F66 - EMM_F34  0.00667 0.404 208   0.017  1.0000
    ##  Control - EMM_F34  0.07333 0.404 208   0.182  1.0000
    ##  Control - EMM_F66  0.06667 0.404 208   0.165  1.0000
    ##  EMM_F49 - EMM_F34  0.08833 0.452 208   0.196  1.0000
    ##  EMM_F49 - EMM_F66  0.08167 0.452 208   0.181  1.0000
    ##  EMM_F49 - Control  0.01500 0.452 208   0.033  1.0000
    ##  EMM_F89 - EMM_F34  0.09333 0.452 208   0.207  1.0000
    ##  EMM_F89 - EMM_F66  0.08667 0.452 208   0.192  1.0000
    ##  EMM_F89 - Control  0.02000 0.452 208   0.044  1.0000
    ##  EMM_F89 - EMM_F49  0.00500 0.495 208   0.010  1.0000
    ##  EMM_F70 - EMM_F34  0.20000 0.404 208   0.495  1.0000
    ##  EMM_F70 - EMM_F66  0.19333 0.404 208   0.479  1.0000
    ##  EMM_F70 - Control  0.12667 0.404 208   0.314  1.0000
    ##  EMM_F70 - EMM_F49  0.11167 0.452 208   0.247  1.0000
    ##  EMM_F70 - EMM_F89  0.10667 0.452 208   0.236  1.0000
    ##  ZAN_F3 - EMM_F34   0.21333 0.404 208   0.528  1.0000
    ##  ZAN_F3 - EMM_F66   0.20667 0.404 208   0.512  1.0000
    ##  ZAN_F3 - Control   0.14000 0.404 208   0.347  1.0000
    ##  ZAN_F3 - EMM_F49   0.12500 0.452 208   0.277  1.0000
    ##  ZAN_F3 - EMM_F89   0.12000 0.452 208   0.266  1.0000
    ##  ZAN_F3 - EMM_F70   0.01333 0.404 208   0.033  1.0000
    ##  EMM_F3 - EMM_F34   0.22333 0.404 208   0.553  1.0000
    ##  EMM_F3 - EMM_F66   0.21667 0.404 208   0.536  1.0000
    ##  EMM_F3 - Control   0.15000 0.404 208   0.371  1.0000
    ##  EMM_F3 - EMM_F49   0.13500 0.452 208   0.299  1.0000
    ##  EMM_F3 - EMM_F89   0.13000 0.452 208   0.288  1.0000
    ##  EMM_F3 - EMM_F70   0.02333 0.404 208   0.058  1.0000
    ##  EMM_F3 - ZAN_F3    0.01000 0.404 208   0.025  1.0000
    ##  ZAN_F4 - EMM_F34   0.47000 0.404 208   1.164  0.9987
    ##  ZAN_F4 - EMM_F66   0.46333 0.404 208   1.147  0.9989
    ##  ZAN_F4 - Control   0.39667 0.404 208   0.982  0.9998
    ##  ZAN_F4 - EMM_F49   0.38167 0.452 208   0.845  1.0000
    ##  ZAN_F4 - EMM_F89   0.37667 0.452 208   0.834  1.0000
    ##  ZAN_F4 - EMM_F70   0.27000 0.404 208   0.669  1.0000
    ##  ZAN_F4 - ZAN_F3    0.25667 0.404 208   0.636  1.0000
    ##  ZAN_F4 - EMM_F3    0.24667 0.404 208   0.611  1.0000
    ##  SP_F14 - EMM_F34   0.49000 0.404 208   1.213  0.9979
    ##  SP_F14 - EMM_F66   0.48333 0.404 208   1.197  0.9982
    ##  SP_F14 - Control   0.41667 0.404 208   1.032  0.9997
    ##  SP_F14 - EMM_F49   0.40167 0.452 208   0.890  0.9999
    ##  SP_F14 - EMM_F89   0.39667 0.452 208   0.878  1.0000
    ##  SP_F14 - EMM_F70   0.29000 0.404 208   0.718  1.0000
    ##  SP_F14 - ZAN_F3    0.27667 0.404 208   0.685  1.0000
    ##  SP_F14 - EMM_F3    0.26667 0.404 208   0.660  1.0000
    ##  SP_F14 - ZAN_F4    0.02000 0.404 208   0.050  1.0000
    ##  EMM_F5 - EMM_F34   0.56333 0.404 208   1.395  0.9908
    ##  EMM_F5 - EMM_F66   0.55667 0.404 208   1.378  0.9919
    ##  EMM_F5 - Control   0.49000 0.404 208   1.213  0.9979
    ##  EMM_F5 - EMM_F49   0.47500 0.452 208   1.052  0.9996
    ##  EMM_F5 - EMM_F89   0.47000 0.452 208   1.041  0.9996
    ##  EMM_F5 - EMM_F70   0.36333 0.404 208   0.900  0.9999
    ##  EMM_F5 - ZAN_F3    0.35000 0.404 208   0.867  1.0000
    ##  EMM_F5 - EMM_F3    0.34000 0.404 208   0.842  1.0000
    ##  EMM_F5 - ZAN_F4    0.09333 0.404 208   0.231  1.0000
    ##  EMM_F5 - SP_F14    0.07333 0.404 208   0.182  1.0000
    ##  EMM_F64 - EMM_F34  0.60667 0.404 208   1.502  0.9812
    ##  EMM_F64 - EMM_F66  0.60000 0.404 208   1.486  0.9831
    ##  EMM_F64 - Control  0.53333 0.404 208   1.321  0.9948
    ##  EMM_F64 - EMM_F49  0.51833 0.452 208   1.148  0.9989
    ##  EMM_F64 - EMM_F89  0.51333 0.452 208   1.137  0.9990
    ##  EMM_F64 - EMM_F70  0.40667 0.404 208   1.007  0.9998
    ##  EMM_F64 - ZAN_F3   0.39333 0.404 208   0.974  0.9998
    ##  EMM_F64 - EMM_F3   0.38333 0.404 208   0.949  0.9999
    ##  EMM_F64 - ZAN_F4   0.13667 0.404 208   0.338  1.0000
    ##  EMM_F64 - SP_F14   0.11667 0.404 208   0.289  1.0000
    ##  EMM_F64 - EMM_F5   0.04333 0.404 208   0.107  1.0000
    ##  EMM_F65 - EMM_F34  0.61000 0.404 208   1.510  0.9802
    ##  EMM_F65 - EMM_F66  0.60333 0.404 208   1.494  0.9822
    ##  EMM_F65 - Control  0.53667 0.404 208   1.329  0.9944
    ##  EMM_F65 - EMM_F49  0.52167 0.452 208   1.155  0.9988
    ##  EMM_F65 - EMM_F89  0.51667 0.452 208   1.144  0.9989
    ##  EMM_F65 - EMM_F70  0.41000 0.404 208   1.015  0.9997
    ##  EMM_F65 - ZAN_F3   0.39667 0.404 208   0.982  0.9998
    ##  EMM_F65 - EMM_F3   0.38667 0.404 208   0.957  0.9999
    ##  EMM_F65 - ZAN_F4   0.14000 0.404 208   0.347  1.0000
    ##  EMM_F65 - SP_F14   0.12000 0.404 208   0.297  1.0000
    ##  EMM_F65 - EMM_F5   0.04667 0.404 208   0.116  1.0000
    ##  EMM_F65 - EMM_F64  0.00333 0.404 208   0.008  1.0000
    ##  EMM_F63 - EMM_F34  0.65667 0.404 208   1.626  0.9619
    ##  EMM_F63 - EMM_F66  0.65000 0.404 208   1.609  0.9651
    ##  EMM_F63 - Control  0.58333 0.404 208   1.444  0.9871
    ##  EMM_F63 - EMM_F49  0.56833 0.452 208   1.259  0.9968
    ##  EMM_F63 - EMM_F89  0.56333 0.452 208   1.248  0.9971
    ##  EMM_F63 - EMM_F70  0.45667 0.404 208   1.131  0.9990
    ##  EMM_F63 - ZAN_F3   0.44333 0.404 208   1.098  0.9993
    ##  EMM_F63 - EMM_F3   0.43333 0.404 208   1.073  0.9995
    ##  EMM_F63 - ZAN_F4   0.18667 0.404 208   0.462  1.0000
    ##  EMM_F63 - SP_F14   0.16667 0.404 208   0.413  1.0000
    ##  EMM_F63 - EMM_F5   0.09333 0.404 208   0.231  1.0000
    ##  EMM_F63 - EMM_F64  0.05000 0.404 208   0.124  1.0000
    ##  EMM_F63 - EMM_F65  0.04667 0.404 208   0.116  1.0000
    ##  EMM_F7 - EMM_F34   0.67000 0.404 208   1.659  0.9548
    ##  EMM_F7 - EMM_F66   0.66333 0.404 208   1.642  0.9584
    ##  EMM_F7 - Control   0.59667 0.404 208   1.477  0.9839
    ##  EMM_F7 - EMM_F49   0.58167 0.452 208   1.288  0.9960
    ##  EMM_F7 - EMM_F89   0.57667 0.452 208   1.277  0.9963
    ##  EMM_F7 - EMM_F70   0.47000 0.404 208   1.164  0.9987
    ##  EMM_F7 - ZAN_F3    0.45667 0.404 208   1.131  0.9990
    ##  EMM_F7 - EMM_F3    0.44667 0.404 208   1.106  0.9993
    ##  EMM_F7 - ZAN_F4    0.20000 0.404 208   0.495  1.0000
    ##  EMM_F7 - SP_F14    0.18000 0.404 208   0.446  1.0000
    ##  EMM_F7 - EMM_F5    0.10667 0.404 208   0.264  1.0000
    ##  EMM_F7 - EMM_F64   0.06333 0.404 208   0.157  1.0000
    ##  EMM_F7 - EMM_F65   0.06000 0.404 208   0.149  1.0000
    ##  EMM_F7 - EMM_F63   0.01333 0.404 208   0.033  1.0000
    ##  EMM_F48 - EMM_F34  0.74333 0.404 208   1.841  0.8985
    ##  EMM_F48 - EMM_F66  0.73667 0.404 208   1.824  0.9050
    ##  EMM_F48 - Control  0.67000 0.404 208   1.659  0.9548
    ##  EMM_F48 - EMM_F49  0.65500 0.452 208   1.451  0.9865
    ##  EMM_F48 - EMM_F89  0.65000 0.452 208   1.440  0.9875
    ##  EMM_F48 - EMM_F70  0.54333 0.404 208   1.345  0.9936
    ##  EMM_F48 - ZAN_F3   0.53000 0.404 208   1.312  0.9951
    ##  EMM_F48 - EMM_F3   0.52000 0.404 208   1.288  0.9960
    ##  EMM_F48 - ZAN_F4   0.27333 0.404 208   0.677  1.0000
    ##  EMM_F48 - SP_F14   0.25333 0.404 208   0.627  1.0000
    ##  EMM_F48 - EMM_F5   0.18000 0.404 208   0.446  1.0000
    ##  EMM_F48 - EMM_F64  0.13667 0.404 208   0.338  1.0000
    ##  EMM_F48 - EMM_F65  0.13333 0.404 208   0.330  1.0000
    ##  EMM_F48 - EMM_F63  0.08667 0.404 208   0.215  1.0000
    ##  EMM_F48 - EMM_F7   0.07333 0.404 208   0.182  1.0000
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65  0.09000 0.404 208   0.223  1.0000
    ##  EMM_F89 - EMM_F65  0.13333 0.452 208   0.295  1.0000
    ##  EMM_F89 - Control  0.04333 0.452 208   0.096  1.0000
    ##  EMM_F49 - EMM_F65  0.16833 0.452 208   0.373  1.0000
    ##  EMM_F49 - Control  0.07833 0.452 208   0.173  1.0000
    ##  EMM_F49 - EMM_F89  0.03500 0.495 208   0.071  1.0000
    ##  EMM_F3 - EMM_F65   0.17333 0.404 208   0.429  1.0000
    ##  EMM_F3 - Control   0.08333 0.404 208   0.206  1.0000
    ##  EMM_F3 - EMM_F89   0.04000 0.452 208   0.089  1.0000
    ##  EMM_F3 - EMM_F49   0.00500 0.452 208   0.011  1.0000
    ##  EMM_F34 - EMM_F65  0.20667 0.404 208   0.512  1.0000
    ##  EMM_F34 - Control  0.11667 0.404 208   0.289  1.0000
    ##  EMM_F34 - EMM_F89  0.07333 0.452 208   0.162  1.0000
    ##  EMM_F34 - EMM_F49  0.03833 0.452 208   0.085  1.0000
    ##  EMM_F34 - EMM_F3   0.03333 0.404 208   0.083  1.0000
    ##  EMM_F7 - EMM_F65   0.30333 0.404 208   0.751  1.0000
    ##  EMM_F7 - Control   0.21333 0.404 208   0.528  1.0000
    ##  EMM_F7 - EMM_F89   0.17000 0.452 208   0.376  1.0000
    ##  EMM_F7 - EMM_F49   0.13500 0.452 208   0.299  1.0000
    ##  EMM_F7 - EMM_F3    0.13000 0.404 208   0.322  1.0000
    ##  EMM_F7 - EMM_F34   0.09667 0.404 208   0.239  1.0000
    ##  EMM_F5 - EMM_F65   0.31000 0.404 208   0.768  1.0000
    ##  EMM_F5 - Control   0.22000 0.404 208   0.545  1.0000
    ##  EMM_F5 - EMM_F89   0.17667 0.452 208   0.391  1.0000
    ##  EMM_F5 - EMM_F49   0.14167 0.452 208   0.314  1.0000
    ##  EMM_F5 - EMM_F3    0.13667 0.404 208   0.338  1.0000
    ##  EMM_F5 - EMM_F34   0.10333 0.404 208   0.256  1.0000
    ##  EMM_F5 - EMM_F7    0.00667 0.404 208   0.017  1.0000
    ##  EMM_F66 - EMM_F65  0.31667 0.404 208   0.784  1.0000
    ##  EMM_F66 - Control  0.22667 0.404 208   0.561  1.0000
    ##  EMM_F66 - EMM_F89  0.18333 0.452 208   0.406  1.0000
    ##  EMM_F66 - EMM_F49  0.14833 0.452 208   0.329  1.0000
    ##  EMM_F66 - EMM_F3   0.14333 0.404 208   0.355  1.0000
    ##  EMM_F66 - EMM_F34  0.11000 0.404 208   0.272  1.0000
    ##  EMM_F66 - EMM_F7   0.01333 0.404 208   0.033  1.0000
    ##  EMM_F66 - EMM_F5   0.00667 0.404 208   0.017  1.0000
    ##  EMM_F70 - EMM_F65  0.32000 0.404 208   0.792  1.0000
    ##  EMM_F70 - Control  0.23000 0.404 208   0.569  1.0000
    ##  EMM_F70 - EMM_F89  0.18667 0.452 208   0.413  1.0000
    ##  EMM_F70 - EMM_F49  0.15167 0.452 208   0.336  1.0000
    ##  EMM_F70 - EMM_F3   0.14667 0.404 208   0.363  1.0000
    ##  EMM_F70 - EMM_F34  0.11333 0.404 208   0.281  1.0000
    ##  EMM_F70 - EMM_F7   0.01667 0.404 208   0.041  1.0000
    ##  EMM_F70 - EMM_F5   0.01000 0.404 208   0.025  1.0000
    ##  EMM_F70 - EMM_F66  0.00333 0.404 208   0.008  1.0000
    ##  SP_F14 - EMM_F65   0.32667 0.404 208   0.809  1.0000
    ##  SP_F14 - Control   0.23667 0.404 208   0.586  1.0000
    ##  SP_F14 - EMM_F89   0.19333 0.452 208   0.428  1.0000
    ##  SP_F14 - EMM_F49   0.15833 0.452 208   0.351  1.0000
    ##  SP_F14 - EMM_F3    0.15333 0.404 208   0.380  1.0000
    ##  SP_F14 - EMM_F34   0.12000 0.404 208   0.297  1.0000
    ##  SP_F14 - EMM_F7    0.02333 0.404 208   0.058  1.0000
    ##  SP_F14 - EMM_F5    0.01667 0.404 208   0.041  1.0000
    ##  SP_F14 - EMM_F66   0.01000 0.404 208   0.025  1.0000
    ##  SP_F14 - EMM_F70   0.00667 0.404 208   0.017  1.0000
    ##  ZAN_F4 - EMM_F65   0.58667 0.404 208   1.453  0.9863
    ##  ZAN_F4 - Control   0.49667 0.404 208   1.230  0.9976
    ##  ZAN_F4 - EMM_F89   0.45333 0.452 208   1.004  0.9998
    ##  ZAN_F4 - EMM_F49   0.41833 0.452 208   0.926  0.9999
    ##  ZAN_F4 - EMM_F3    0.41333 0.404 208   1.023  0.9997
    ##  ZAN_F4 - EMM_F34   0.38000 0.404 208   0.941  0.9999
    ##  ZAN_F4 - EMM_F7    0.28333 0.404 208   0.702  1.0000
    ##  ZAN_F4 - EMM_F5    0.27667 0.404 208   0.685  1.0000
    ##  ZAN_F4 - EMM_F66   0.27000 0.404 208   0.669  1.0000
    ##  ZAN_F4 - EMM_F70   0.26667 0.404 208   0.660  1.0000
    ##  ZAN_F4 - SP_F14    0.26000 0.404 208   0.644  1.0000
    ##  EMM_F48 - EMM_F65  0.59333 0.404 208   1.469  0.9848
    ##  EMM_F48 - Control  0.50333 0.404 208   1.246  0.9972
    ##  EMM_F48 - EMM_F89  0.46000 0.452 208   1.019  0.9997
    ##  EMM_F48 - EMM_F49  0.42500 0.452 208   0.941  0.9999
    ##  EMM_F48 - EMM_F3   0.42000 0.404 208   1.040  0.9996
    ##  EMM_F48 - EMM_F34  0.38667 0.404 208   0.957  0.9999
    ##  EMM_F48 - EMM_F7   0.29000 0.404 208   0.718  1.0000
    ##  EMM_F48 - EMM_F5   0.28333 0.404 208   0.702  1.0000
    ##  EMM_F48 - EMM_F66  0.27667 0.404 208   0.685  1.0000
    ##  EMM_F48 - EMM_F70  0.27333 0.404 208   0.677  1.0000
    ##  EMM_F48 - SP_F14   0.26667 0.404 208   0.660  1.0000
    ##  EMM_F48 - ZAN_F4   0.00667 0.404 208   0.017  1.0000
    ##  ZAN_F3 - EMM_F65   0.71333 0.404 208   1.766  0.9253
    ##  ZAN_F3 - Control   0.62333 0.404 208   1.543  0.9759
    ##  ZAN_F3 - EMM_F89   0.58000 0.452 208   1.284  0.9961
    ##  ZAN_F3 - EMM_F49   0.54500 0.452 208   1.207  0.9980
    ##  ZAN_F3 - EMM_F3    0.54000 0.404 208   1.337  0.9940
    ##  ZAN_F3 - EMM_F34   0.50667 0.404 208   1.255  0.9970
    ##  ZAN_F3 - EMM_F7    0.41000 0.404 208   1.015  0.9997
    ##  ZAN_F3 - EMM_F5    0.40333 0.404 208   0.999  0.9998
    ##  ZAN_F3 - EMM_F66   0.39667 0.404 208   0.982  0.9998
    ##  ZAN_F3 - EMM_F70   0.39333 0.404 208   0.974  0.9998
    ##  ZAN_F3 - SP_F14    0.38667 0.404 208   0.957  0.9999
    ##  ZAN_F3 - ZAN_F4    0.12667 0.404 208   0.314  1.0000
    ##  ZAN_F3 - EMM_F48   0.12000 0.404 208   0.297  1.0000
    ##  EMM_F64 - EMM_F65  0.75667 0.404 208   1.874  0.8849
    ##  EMM_F64 - Control  0.66667 0.404 208   1.651  0.9567
    ##  EMM_F64 - EMM_F89  0.62333 0.452 208   1.380  0.9917
    ##  EMM_F64 - EMM_F49  0.58833 0.452 208   1.303  0.9954
    ##  EMM_F64 - EMM_F3   0.58333 0.404 208   1.444  0.9871
    ##  EMM_F64 - EMM_F34  0.55000 0.404 208   1.362  0.9928
    ##  EMM_F64 - EMM_F7   0.45333 0.404 208   1.122  0.9991
    ##  EMM_F64 - EMM_F5   0.44667 0.404 208   1.106  0.9993
    ##  EMM_F64 - EMM_F66  0.44000 0.404 208   1.089  0.9994
    ##  EMM_F64 - EMM_F70  0.43667 0.404 208   1.081  0.9994
    ##  EMM_F64 - SP_F14   0.43000 0.404 208   1.065  0.9995
    ##  EMM_F64 - ZAN_F4   0.17000 0.404 208   0.421  1.0000
    ##  EMM_F64 - EMM_F48  0.16333 0.404 208   0.404  1.0000
    ##  EMM_F64 - ZAN_F3   0.04333 0.404 208   0.107  1.0000
    ##  EMM_F63 - EMM_F65  0.79333 0.404 208   1.964  0.8417
    ##  EMM_F63 - Control  0.70333 0.404 208   1.741  0.9331
    ##  EMM_F63 - EMM_F89  0.66000 0.452 208   1.462  0.9855
    ##  EMM_F63 - EMM_F49  0.62500 0.452 208   1.384  0.9915
    ##  EMM_F63 - EMM_F3   0.62000 0.404 208   1.535  0.9771
    ##  EMM_F63 - EMM_F34  0.58667 0.404 208   1.453  0.9863
    ##  EMM_F63 - EMM_F7   0.49000 0.404 208   1.213  0.9979
    ##  EMM_F63 - EMM_F5   0.48333 0.404 208   1.197  0.9982
    ##  EMM_F63 - EMM_F66  0.47667 0.404 208   1.180  0.9984
    ##  EMM_F63 - EMM_F70  0.47333 0.404 208   1.172  0.9986
    ##  EMM_F63 - SP_F14   0.46667 0.404 208   1.155  0.9988
    ##  EMM_F63 - ZAN_F4   0.20667 0.404 208   0.512  1.0000
    ##  EMM_F63 - EMM_F48  0.20000 0.404 208   0.495  1.0000
    ##  EMM_F63 - ZAN_F3   0.08000 0.404 208   0.198  1.0000
    ##  EMM_F63 - EMM_F64  0.03667 0.404 208   0.091  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F70 - EMM_F49  0.27167 0.452 208   0.602  1.0000
    ##  EMM_F66 - EMM_F49  0.35167 0.452 208   0.779  1.0000
    ##  EMM_F66 - EMM_F70  0.08000 0.404 208   0.198  1.0000
    ##  Control - EMM_F49  0.37500 0.452 208   0.830  1.0000
    ##  Control - EMM_F70  0.10333 0.404 208   0.256  1.0000
    ##  Control - EMM_F66  0.02333 0.404 208   0.058  1.0000
    ##  EMM_F48 - EMM_F49  0.41500 0.452 208   0.919  0.9999
    ##  EMM_F48 - EMM_F70  0.14333 0.404 208   0.355  1.0000
    ##  EMM_F48 - EMM_F66  0.06333 0.404 208   0.157  1.0000
    ##  EMM_F48 - Control  0.04000 0.404 208   0.099  1.0000
    ##  EMM_F3 - EMM_F49   0.45833 0.452 208   1.015  0.9997
    ##  EMM_F3 - EMM_F70   0.18667 0.404 208   0.462  1.0000
    ##  EMM_F3 - EMM_F66   0.10667 0.404 208   0.264  1.0000
    ##  EMM_F3 - Control   0.08333 0.404 208   0.206  1.0000
    ##  EMM_F3 - EMM_F48   0.04333 0.404 208   0.107  1.0000
    ##  EMM_F5 - EMM_F49   0.47500 0.452 208   1.052  0.9996
    ##  EMM_F5 - EMM_F70   0.20333 0.404 208   0.503  1.0000
    ##  EMM_F5 - EMM_F66   0.12333 0.404 208   0.305  1.0000
    ##  EMM_F5 - Control   0.10000 0.404 208   0.248  1.0000
    ##  EMM_F5 - EMM_F48   0.06000 0.404 208   0.149  1.0000
    ##  EMM_F5 - EMM_F3    0.01667 0.404 208   0.041  1.0000
    ##  EMM_F34 - EMM_F49  0.49500 0.452 208   1.096  0.9993
    ##  EMM_F34 - EMM_F70  0.22333 0.404 208   0.553  1.0000
    ##  EMM_F34 - EMM_F66  0.14333 0.404 208   0.355  1.0000
    ##  EMM_F34 - Control  0.12000 0.404 208   0.297  1.0000
    ##  EMM_F34 - EMM_F48  0.08000 0.404 208   0.198  1.0000
    ##  EMM_F34 - EMM_F3   0.03667 0.404 208   0.091  1.0000
    ##  EMM_F34 - EMM_F5   0.02000 0.404 208   0.050  1.0000
    ##  ZAN_F4 - EMM_F49   0.51167 0.452 208   1.133  0.9990
    ##  ZAN_F4 - EMM_F70   0.24000 0.404 208   0.594  1.0000
    ##  ZAN_F4 - EMM_F66   0.16000 0.404 208   0.396  1.0000
    ##  ZAN_F4 - Control   0.13667 0.404 208   0.338  1.0000
    ##  ZAN_F4 - EMM_F48   0.09667 0.404 208   0.239  1.0000
    ##  ZAN_F4 - EMM_F3    0.05333 0.404 208   0.132  1.0000
    ##  ZAN_F4 - EMM_F5    0.03667 0.404 208   0.091  1.0000
    ##  ZAN_F4 - EMM_F34   0.01667 0.404 208   0.041  1.0000
    ##  EMM_F89 - EMM_F49  0.63500 0.495 208   1.284  0.9961
    ##  EMM_F89 - EMM_F70  0.36333 0.452 208   0.805  1.0000
    ##  EMM_F89 - EMM_F66  0.28333 0.452 208   0.627  1.0000
    ##  EMM_F89 - Control  0.26000 0.452 208   0.576  1.0000
    ##  EMM_F89 - EMM_F48  0.22000 0.452 208   0.487  1.0000
    ##  EMM_F89 - EMM_F3   0.17667 0.452 208   0.391  1.0000
    ##  EMM_F89 - EMM_F5   0.16000 0.452 208   0.354  1.0000
    ##  EMM_F89 - EMM_F34  0.14000 0.452 208   0.310  1.0000
    ##  EMM_F89 - ZAN_F4   0.12333 0.452 208   0.273  1.0000
    ##  EMM_F64 - EMM_F49  0.69833 0.452 208   1.547  0.9755
    ##  EMM_F64 - EMM_F70  0.42667 0.404 208   1.056  0.9996
    ##  EMM_F64 - EMM_F66  0.34667 0.404 208   0.858  1.0000
    ##  EMM_F64 - Control  0.32333 0.404 208   0.801  1.0000
    ##  EMM_F64 - EMM_F48  0.28333 0.404 208   0.702  1.0000
    ##  EMM_F64 - EMM_F3   0.24000 0.404 208   0.594  1.0000
    ##  EMM_F64 - EMM_F5   0.22333 0.404 208   0.553  1.0000
    ##  EMM_F64 - EMM_F34  0.20333 0.404 208   0.503  1.0000
    ##  EMM_F64 - ZAN_F4   0.18667 0.404 208   0.462  1.0000
    ##  EMM_F64 - EMM_F89  0.06333 0.452 208   0.140  1.0000
    ##  EMM_F7 - EMM_F49   0.70167 0.452 208   1.554  0.9744
    ##  EMM_F7 - EMM_F70   0.43000 0.404 208   1.065  0.9995
    ##  EMM_F7 - EMM_F66   0.35000 0.404 208   0.867  1.0000
    ##  EMM_F7 - Control   0.32667 0.404 208   0.809  1.0000
    ##  EMM_F7 - EMM_F48   0.28667 0.404 208   0.710  1.0000
    ##  EMM_F7 - EMM_F3    0.24333 0.404 208   0.602  1.0000
    ##  EMM_F7 - EMM_F5    0.22667 0.404 208   0.561  1.0000
    ##  EMM_F7 - EMM_F34   0.20667 0.404 208   0.512  1.0000
    ##  EMM_F7 - ZAN_F4    0.19000 0.404 208   0.470  1.0000
    ##  EMM_F7 - EMM_F89   0.06667 0.452 208   0.148  1.0000
    ##  EMM_F7 - EMM_F64   0.00333 0.404 208   0.008  1.0000
    ##  EMM_F65 - EMM_F49  0.80500 0.452 208   1.783  0.9198
    ##  EMM_F65 - EMM_F70  0.53333 0.404 208   1.321  0.9948
    ##  EMM_F65 - EMM_F66  0.45333 0.404 208   1.122  0.9991
    ##  EMM_F65 - Control  0.43000 0.404 208   1.065  0.9995
    ##  EMM_F65 - EMM_F48  0.39000 0.404 208   0.966  0.9999
    ##  EMM_F65 - EMM_F3   0.34667 0.404 208   0.858  1.0000
    ##  EMM_F65 - EMM_F5   0.33000 0.404 208   0.817  1.0000
    ##  EMM_F65 - EMM_F34  0.31000 0.404 208   0.768  1.0000
    ##  EMM_F65 - ZAN_F4   0.29333 0.404 208   0.726  1.0000
    ##  EMM_F65 - EMM_F89  0.17000 0.452 208   0.376  1.0000
    ##  EMM_F65 - EMM_F64  0.10667 0.404 208   0.264  1.0000
    ##  EMM_F65 - EMM_F7   0.10333 0.404 208   0.256  1.0000
    ##  SP_F14 - EMM_F49   0.96167 0.452 208   2.130  0.7438
    ##  SP_F14 - EMM_F70   0.69000 0.404 208   1.708  0.9425
    ##  SP_F14 - EMM_F66   0.61000 0.404 208   1.510  0.9802
    ##  SP_F14 - Control   0.58667 0.404 208   1.453  0.9863
    ##  SP_F14 - EMM_F48   0.54667 0.404 208   1.354  0.9932
    ##  SP_F14 - EMM_F3    0.50333 0.404 208   1.246  0.9972
    ##  SP_F14 - EMM_F5    0.48667 0.404 208   1.205  0.9980
    ##  SP_F14 - EMM_F34   0.46667 0.404 208   1.155  0.9988
    ##  SP_F14 - ZAN_F4    0.45000 0.404 208   1.114  0.9992
    ##  SP_F14 - EMM_F89   0.32667 0.452 208   0.723  1.0000
    ##  SP_F14 - EMM_F64   0.26333 0.404 208   0.652  1.0000
    ##  SP_F14 - EMM_F7    0.26000 0.404 208   0.644  1.0000
    ##  SP_F14 - EMM_F65   0.15667 0.404 208   0.388  1.0000
    ##  EMM_F63 - EMM_F49  0.97167 0.452 208   2.152  0.7291
    ##  EMM_F63 - EMM_F70  0.70000 0.404 208   1.733  0.9355
    ##  EMM_F63 - EMM_F66  0.62000 0.404 208   1.535  0.9771
    ##  EMM_F63 - Control  0.59667 0.404 208   1.477  0.9839
    ##  EMM_F63 - EMM_F48  0.55667 0.404 208   1.378  0.9919
    ##  EMM_F63 - EMM_F3   0.51333 0.404 208   1.271  0.9965
    ##  EMM_F63 - EMM_F5   0.49667 0.404 208   1.230  0.9976
    ##  EMM_F63 - EMM_F34  0.47667 0.404 208   1.180  0.9984
    ##  EMM_F63 - ZAN_F4   0.46000 0.404 208   1.139  0.9990
    ##  EMM_F63 - EMM_F89  0.33667 0.452 208   0.746  1.0000
    ##  EMM_F63 - EMM_F64  0.27333 0.404 208   0.677  1.0000
    ##  EMM_F63 - EMM_F7   0.27000 0.404 208   0.669  1.0000
    ##  EMM_F63 - EMM_F65  0.16667 0.404 208   0.413  1.0000
    ##  EMM_F63 - SP_F14   0.01000 0.404 208   0.025  1.0000
    ##  ZAN_F3 - EMM_F49   1.02500 0.452 208   2.270  0.6461
    ##  ZAN_F3 - EMM_F70   0.75333 0.404 208   1.865  0.8884
    ##  ZAN_F3 - EMM_F66   0.67333 0.404 208   1.667  0.9529
    ##  ZAN_F3 - Control   0.65000 0.404 208   1.609  0.9651
    ##  ZAN_F3 - EMM_F48   0.61000 0.404 208   1.510  0.9802
    ##  ZAN_F3 - EMM_F3    0.56667 0.404 208   1.403  0.9903
    ##  ZAN_F3 - EMM_F5    0.55000 0.404 208   1.362  0.9928
    ##  ZAN_F3 - EMM_F34   0.53000 0.404 208   1.312  0.9951
    ##  ZAN_F3 - ZAN_F4    0.51333 0.404 208   1.271  0.9965
    ##  ZAN_F3 - EMM_F89   0.39000 0.452 208   0.864  1.0000
    ##  ZAN_F3 - EMM_F64   0.32667 0.404 208   0.809  1.0000
    ##  ZAN_F3 - EMM_F7    0.32333 0.404 208   0.801  1.0000
    ##  ZAN_F3 - EMM_F65   0.22000 0.404 208   0.545  1.0000
    ##  ZAN_F3 - SP_F14    0.06333 0.404 208   0.157  1.0000
    ##  ZAN_F3 - EMM_F63   0.05333 0.404 208   0.132  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F48 - EMM_F49  0.20500 0.452 208   0.454  1.0000
    ##  Control - EMM_F49  0.24167 0.452 208   0.535  1.0000
    ##  Control - EMM_F48  0.03667 0.404 208   0.091  1.0000
    ##  EMM_F70 - EMM_F49  0.39833 0.452 208   0.882  1.0000
    ##  EMM_F70 - EMM_F48  0.19333 0.404 208   0.479  1.0000
    ##  EMM_F70 - Control  0.15667 0.404 208   0.388  1.0000
    ##  EMM_F3 - EMM_F49   0.41833 0.452 208   0.926  0.9999
    ##  EMM_F3 - EMM_F48   0.21333 0.404 208   0.528  1.0000
    ##  EMM_F3 - Control   0.17667 0.404 208   0.437  1.0000
    ##  EMM_F3 - EMM_F70   0.02000 0.404 208   0.050  1.0000
    ##  EMM_F7 - EMM_F49   0.47500 0.452 208   1.052  0.9996
    ##  EMM_F7 - EMM_F48   0.27000 0.404 208   0.669  1.0000
    ##  EMM_F7 - Control   0.23333 0.404 208   0.578  1.0000
    ##  EMM_F7 - EMM_F70   0.07667 0.404 208   0.190  1.0000
    ##  EMM_F7 - EMM_F3    0.05667 0.404 208   0.140  1.0000
    ##  EMM_F63 - EMM_F49  0.47833 0.452 208   1.059  0.9996
    ##  EMM_F63 - EMM_F48  0.27333 0.404 208   0.677  1.0000
    ##  EMM_F63 - Control  0.23667 0.404 208   0.586  1.0000
    ##  EMM_F63 - EMM_F70  0.08000 0.404 208   0.198  1.0000
    ##  EMM_F63 - EMM_F3   0.06000 0.404 208   0.149  1.0000
    ##  EMM_F63 - EMM_F7   0.00333 0.404 208   0.008  1.0000
    ##  EMM_F34 - EMM_F49  0.53167 0.452 208   1.177  0.9985
    ##  EMM_F34 - EMM_F48  0.32667 0.404 208   0.809  1.0000
    ##  EMM_F34 - Control  0.29000 0.404 208   0.718  1.0000
    ##  EMM_F34 - EMM_F70  0.13333 0.404 208   0.330  1.0000
    ##  EMM_F34 - EMM_F3   0.11333 0.404 208   0.281  1.0000
    ##  EMM_F34 - EMM_F7   0.05667 0.404 208   0.140  1.0000
    ##  EMM_F34 - EMM_F63  0.05333 0.404 208   0.132  1.0000
    ##  EMM_F5 - EMM_F49   0.58833 0.452 208   1.303  0.9954
    ##  EMM_F5 - EMM_F48   0.38333 0.404 208   0.949  0.9999
    ##  EMM_F5 - Control   0.34667 0.404 208   0.858  1.0000
    ##  EMM_F5 - EMM_F70   0.19000 0.404 208   0.470  1.0000
    ##  EMM_F5 - EMM_F3    0.17000 0.404 208   0.421  1.0000
    ##  EMM_F5 - EMM_F7    0.11333 0.404 208   0.281  1.0000
    ##  EMM_F5 - EMM_F63   0.11000 0.404 208   0.272  1.0000
    ##  EMM_F5 - EMM_F34   0.05667 0.404 208   0.140  1.0000
    ##  ZAN_F4 - EMM_F49   0.65500 0.452 208   1.451  0.9865
    ##  ZAN_F4 - EMM_F48   0.45000 0.404 208   1.114  0.9992
    ##  ZAN_F4 - Control   0.41333 0.404 208   1.023  0.9997
    ##  ZAN_F4 - EMM_F70   0.25667 0.404 208   0.636  1.0000
    ##  ZAN_F4 - EMM_F3    0.23667 0.404 208   0.586  1.0000
    ##  ZAN_F4 - EMM_F7    0.18000 0.404 208   0.446  1.0000
    ##  ZAN_F4 - EMM_F63   0.17667 0.404 208   0.437  1.0000
    ##  ZAN_F4 - EMM_F34   0.12333 0.404 208   0.305  1.0000
    ##  ZAN_F4 - EMM_F5    0.06667 0.404 208   0.165  1.0000
    ##  SP_F14 - EMM_F49   0.68167 0.452 208   1.510  0.9803
    ##  SP_F14 - EMM_F48   0.47667 0.404 208   1.180  0.9984
    ##  SP_F14 - Control   0.44000 0.404 208   1.089  0.9994
    ##  SP_F14 - EMM_F70   0.28333 0.404 208   0.702  1.0000
    ##  SP_F14 - EMM_F3    0.26333 0.404 208   0.652  1.0000
    ##  SP_F14 - EMM_F7    0.20667 0.404 208   0.512  1.0000
    ##  SP_F14 - EMM_F63   0.20333 0.404 208   0.503  1.0000
    ##  SP_F14 - EMM_F34   0.15000 0.404 208   0.371  1.0000
    ##  SP_F14 - EMM_F5    0.09333 0.404 208   0.231  1.0000
    ##  SP_F14 - ZAN_F4    0.02667 0.404 208   0.066  1.0000
    ##  EMM_F66 - EMM_F49  0.70500 0.452 208   1.561  0.9733
    ##  EMM_F66 - EMM_F48  0.50000 0.404 208   1.238  0.9974
    ##  EMM_F66 - Control  0.46333 0.404 208   1.147  0.9989
    ##  EMM_F66 - EMM_F70  0.30667 0.404 208   0.759  1.0000
    ##  EMM_F66 - EMM_F3   0.28667 0.404 208   0.710  1.0000
    ##  EMM_F66 - EMM_F7   0.23000 0.404 208   0.569  1.0000
    ##  EMM_F66 - EMM_F63  0.22667 0.404 208   0.561  1.0000
    ##  EMM_F66 - EMM_F34  0.17333 0.404 208   0.429  1.0000
    ##  EMM_F66 - EMM_F5   0.11667 0.404 208   0.289  1.0000
    ##  EMM_F66 - ZAN_F4   0.05000 0.404 208   0.124  1.0000
    ##  EMM_F66 - SP_F14   0.02333 0.404 208   0.058  1.0000
    ##  EMM_F64 - EMM_F49  0.79167 0.452 208   1.753  0.9295
    ##  EMM_F64 - EMM_F48  0.58667 0.404 208   1.453  0.9863
    ##  EMM_F64 - Control  0.55000 0.404 208   1.362  0.9928
    ##  EMM_F64 - EMM_F70  0.39333 0.404 208   0.974  0.9998
    ##  EMM_F64 - EMM_F3   0.37333 0.404 208   0.924  0.9999
    ##  EMM_F64 - EMM_F7   0.31667 0.404 208   0.784  1.0000
    ##  EMM_F64 - EMM_F63  0.31333 0.404 208   0.776  1.0000
    ##  EMM_F64 - EMM_F34  0.26000 0.404 208   0.644  1.0000
    ##  EMM_F64 - EMM_F5   0.20333 0.404 208   0.503  1.0000
    ##  EMM_F64 - ZAN_F4   0.13667 0.404 208   0.338  1.0000
    ##  EMM_F64 - SP_F14   0.11000 0.404 208   0.272  1.0000
    ##  EMM_F64 - EMM_F66  0.08667 0.404 208   0.215  1.0000
    ##  EMM_F65 - EMM_F49  0.96833 0.452 208   2.144  0.7341
    ##  EMM_F65 - EMM_F48  0.76333 0.404 208   1.890  0.8776
    ##  EMM_F65 - Control  0.72667 0.404 208   1.799  0.9141
    ##  EMM_F65 - EMM_F70  0.57000 0.404 208   1.411  0.9897
    ##  EMM_F65 - EMM_F3   0.55000 0.404 208   1.362  0.9928
    ##  EMM_F65 - EMM_F7   0.49333 0.404 208   1.222  0.9977
    ##  EMM_F65 - EMM_F63  0.49000 0.404 208   1.213  0.9979
    ##  EMM_F65 - EMM_F34  0.43667 0.404 208   1.081  0.9994
    ##  EMM_F65 - EMM_F5   0.38000 0.404 208   0.941  0.9999
    ##  EMM_F65 - ZAN_F4   0.31333 0.404 208   0.776  1.0000
    ##  EMM_F65 - SP_F14   0.28667 0.404 208   0.710  1.0000
    ##  EMM_F65 - EMM_F66  0.26333 0.404 208   0.652  1.0000
    ##  EMM_F65 - EMM_F64  0.17667 0.404 208   0.437  1.0000
    ##  EMM_F89 - EMM_F49  1.12000 0.495 208   2.264  0.6503
    ##  EMM_F89 - EMM_F48  0.91500 0.452 208   2.026  0.8077
    ##  EMM_F89 - Control  0.87833 0.452 208   1.945  0.8515
    ##  EMM_F89 - EMM_F70  0.72167 0.452 208   1.598  0.9671
    ##  EMM_F89 - EMM_F3   0.70167 0.452 208   1.554  0.9744
    ##  EMM_F89 - EMM_F7   0.64500 0.452 208   1.428  0.9884
    ##  EMM_F89 - EMM_F63  0.64167 0.452 208   1.421  0.9890
    ##  EMM_F89 - EMM_F34  0.58833 0.452 208   1.303  0.9954
    ##  EMM_F89 - EMM_F5   0.53167 0.452 208   1.177  0.9985
    ##  EMM_F89 - ZAN_F4   0.46500 0.452 208   1.030  0.9997
    ##  EMM_F89 - SP_F14   0.43833 0.452 208   0.971  0.9998
    ##  EMM_F89 - EMM_F66  0.41500 0.452 208   0.919  0.9999
    ##  EMM_F89 - EMM_F64  0.32833 0.452 208   0.727  1.0000
    ##  EMM_F89 - EMM_F65  0.15167 0.452 208   0.336  1.0000
    ##  ZAN_F3 - EMM_F49   1.35167 0.452 208   2.993  0.1810
    ##  ZAN_F3 - EMM_F48   1.14667 0.404 208   2.839  0.2558
    ##  ZAN_F3 - Control   1.11000 0.404 208   2.748  0.3078
    ##  ZAN_F3 - EMM_F70   0.95333 0.404 208   2.360  0.5793
    ##  ZAN_F3 - EMM_F3    0.93333 0.404 208   2.311  0.6161
    ##  ZAN_F3 - EMM_F7    0.87667 0.404 208   2.171  0.7164
    ##  ZAN_F3 - EMM_F63   0.87333 0.404 208   2.162  0.7220
    ##  ZAN_F3 - EMM_F34   0.82000 0.404 208   2.030  0.8054
    ##  ZAN_F3 - EMM_F5    0.76333 0.404 208   1.890  0.8776
    ##  ZAN_F3 - ZAN_F4    0.69667 0.404 208   1.725  0.9379
    ##  ZAN_F3 - SP_F14    0.67000 0.404 208   1.659  0.9548
    ##  ZAN_F3 - EMM_F66   0.64667 0.404 208   1.601  0.9666
    ##  ZAN_F3 - EMM_F64   0.56000 0.404 208   1.387  0.9914
    ##  ZAN_F3 - EMM_F65   0.38333 0.404 208   0.949  0.9999
    ##  ZAN_F3 - EMM_F89   0.23167 0.452 208   0.513  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F34 - EMM_F49  0.12167 0.452 208   0.269  1.0000
    ##  ZAN_F3 - EMM_F49   0.14500 0.452 208   0.321  1.0000
    ##  ZAN_F3 - EMM_F34   0.02333 0.404 208   0.058  1.0000
    ##  EMM_F89 - EMM_F49  0.14500 0.495 208   0.293  1.0000
    ##  EMM_F89 - EMM_F34  0.02333 0.452 208   0.052  1.0000
    ##  EMM_F89 - ZAN_F3   0.00000 0.452 208   0.000  1.0000
    ##  EMM_F3 - EMM_F49   0.36500 0.452 208   0.808  1.0000
    ##  EMM_F3 - EMM_F34   0.24333 0.404 208   0.602  1.0000
    ##  EMM_F3 - ZAN_F3    0.22000 0.404 208   0.545  1.0000
    ##  EMM_F3 - EMM_F89   0.22000 0.452 208   0.487  1.0000
    ##  Control - EMM_F49  0.39833 0.452 208   0.882  1.0000
    ##  Control - EMM_F34  0.27667 0.404 208   0.685  1.0000
    ##  Control - ZAN_F3   0.25333 0.404 208   0.627  1.0000
    ##  Control - EMM_F89  0.25333 0.452 208   0.561  1.0000
    ##  Control - EMM_F3   0.03333 0.404 208   0.083  1.0000
    ##  EMM_F63 - EMM_F49  0.55833 0.452 208   1.236  0.9974
    ##  EMM_F63 - EMM_F34  0.43667 0.404 208   1.081  0.9994
    ##  EMM_F63 - ZAN_F3   0.41333 0.404 208   1.023  0.9997
    ##  EMM_F63 - EMM_F89  0.41333 0.452 208   0.915  0.9999
    ##  EMM_F63 - EMM_F3   0.19333 0.404 208   0.479  1.0000
    ##  EMM_F63 - Control  0.16000 0.404 208   0.396  1.0000
    ##  EMM_F48 - EMM_F49  0.56833 0.452 208   1.259  0.9968
    ##  EMM_F48 - EMM_F34  0.44667 0.404 208   1.106  0.9993
    ##  EMM_F48 - ZAN_F3   0.42333 0.404 208   1.048  0.9996
    ##  EMM_F48 - EMM_F89  0.42333 0.452 208   0.938  0.9999
    ##  EMM_F48 - EMM_F3   0.20333 0.404 208   0.503  1.0000
    ##  EMM_F48 - Control  0.17000 0.404 208   0.421  1.0000
    ##  EMM_F48 - EMM_F63  0.01000 0.404 208   0.025  1.0000
    ##  EMM_F70 - EMM_F49  0.58500 0.452 208   1.296  0.9957
    ##  EMM_F70 - EMM_F34  0.46333 0.404 208   1.147  0.9989
    ##  EMM_F70 - ZAN_F3   0.44000 0.404 208   1.089  0.9994
    ##  EMM_F70 - EMM_F89  0.44000 0.452 208   0.974  0.9998
    ##  EMM_F70 - EMM_F3   0.22000 0.404 208   0.545  1.0000
    ##  EMM_F70 - Control  0.18667 0.404 208   0.462  1.0000
    ##  EMM_F70 - EMM_F63  0.02667 0.404 208   0.066  1.0000
    ##  EMM_F70 - EMM_F48  0.01667 0.404 208   0.041  1.0000
    ##  EMM_F66 - EMM_F49  0.63500 0.452 208   1.406  0.9900
    ##  EMM_F66 - EMM_F34  0.51333 0.404 208   1.271  0.9965
    ##  EMM_F66 - ZAN_F3   0.49000 0.404 208   1.213  0.9979
    ##  EMM_F66 - EMM_F89  0.49000 0.452 208   1.085  0.9994
    ##  EMM_F66 - EMM_F3   0.27000 0.404 208   0.669  1.0000
    ##  EMM_F66 - Control  0.23667 0.404 208   0.586  1.0000
    ##  EMM_F66 - EMM_F63  0.07667 0.404 208   0.190  1.0000
    ##  EMM_F66 - EMM_F48  0.06667 0.404 208   0.165  1.0000
    ##  EMM_F66 - EMM_F70  0.05000 0.404 208   0.124  1.0000
    ##  ZAN_F4 - EMM_F49   0.72167 0.452 208   1.598  0.9671
    ##  ZAN_F4 - EMM_F34   0.60000 0.404 208   1.486  0.9831
    ##  ZAN_F4 - ZAN_F3    0.57667 0.404 208   1.428  0.9884
    ##  ZAN_F4 - EMM_F89   0.57667 0.452 208   1.277  0.9963
    ##  ZAN_F4 - EMM_F3    0.35667 0.404 208   0.883  1.0000
    ##  ZAN_F4 - Control   0.32333 0.404 208   0.801  1.0000
    ##  ZAN_F4 - EMM_F63   0.16333 0.404 208   0.404  1.0000
    ##  ZAN_F4 - EMM_F48   0.15333 0.404 208   0.380  1.0000
    ##  ZAN_F4 - EMM_F70   0.13667 0.404 208   0.338  1.0000
    ##  ZAN_F4 - EMM_F66   0.08667 0.404 208   0.215  1.0000
    ##  EMM_F5 - EMM_F49   0.73167 0.452 208   1.620  0.9630
    ##  EMM_F5 - EMM_F34   0.61000 0.404 208   1.510  0.9802
    ##  EMM_F5 - ZAN_F3    0.58667 0.404 208   1.453  0.9863
    ##  EMM_F5 - EMM_F89   0.58667 0.452 208   1.299  0.9956
    ##  EMM_F5 - EMM_F3    0.36667 0.404 208   0.908  0.9999
    ##  EMM_F5 - Control   0.33333 0.404 208   0.825  1.0000
    ##  EMM_F5 - EMM_F63   0.17333 0.404 208   0.429  1.0000
    ##  EMM_F5 - EMM_F48   0.16333 0.404 208   0.404  1.0000
    ##  EMM_F5 - EMM_F70   0.14667 0.404 208   0.363  1.0000
    ##  EMM_F5 - EMM_F66   0.09667 0.404 208   0.239  1.0000
    ##  EMM_F5 - ZAN_F4    0.01000 0.404 208   0.025  1.0000
    ##  EMM_F65 - EMM_F49  0.74167 0.452 208   1.643  0.9584
    ##  EMM_F65 - EMM_F34  0.62000 0.404 208   1.535  0.9771
    ##  EMM_F65 - ZAN_F3   0.59667 0.404 208   1.477  0.9839
    ##  EMM_F65 - EMM_F89  0.59667 0.452 208   1.321  0.9947
    ##  EMM_F65 - EMM_F3   0.37667 0.404 208   0.933  0.9999
    ##  EMM_F65 - Control  0.34333 0.404 208   0.850  1.0000
    ##  EMM_F65 - EMM_F63  0.18333 0.404 208   0.454  1.0000
    ##  EMM_F65 - EMM_F48  0.17333 0.404 208   0.429  1.0000
    ##  EMM_F65 - EMM_F70  0.15667 0.404 208   0.388  1.0000
    ##  EMM_F65 - EMM_F66  0.10667 0.404 208   0.264  1.0000
    ##  EMM_F65 - ZAN_F4   0.02000 0.404 208   0.050  1.0000
    ##  EMM_F65 - EMM_F5   0.01000 0.404 208   0.025  1.0000
    ##  EMM_F7 - EMM_F49   0.76833 0.452 208   1.702  0.9443
    ##  EMM_F7 - EMM_F34   0.64667 0.404 208   1.601  0.9666
    ##  EMM_F7 - ZAN_F3    0.62333 0.404 208   1.543  0.9759
    ##  EMM_F7 - EMM_F89   0.62333 0.452 208   1.380  0.9917
    ##  EMM_F7 - EMM_F3    0.40333 0.404 208   0.999  0.9998
    ##  EMM_F7 - Control   0.37000 0.404 208   0.916  0.9999
    ##  EMM_F7 - EMM_F63   0.21000 0.404 208   0.520  1.0000
    ##  EMM_F7 - EMM_F48   0.20000 0.404 208   0.495  1.0000
    ##  EMM_F7 - EMM_F70   0.18333 0.404 208   0.454  1.0000
    ##  EMM_F7 - EMM_F66   0.13333 0.404 208   0.330  1.0000
    ##  EMM_F7 - ZAN_F4    0.04667 0.404 208   0.116  1.0000
    ##  EMM_F7 - EMM_F5    0.03667 0.404 208   0.091  1.0000
    ##  EMM_F7 - EMM_F65   0.02667 0.404 208   0.066  1.0000
    ##  SP_F14 - EMM_F49   0.87833 0.452 208   1.945  0.8515
    ##  SP_F14 - EMM_F34   0.75667 0.404 208   1.874  0.8849
    ##  SP_F14 - ZAN_F3    0.73333 0.404 208   1.816  0.9081
    ##  SP_F14 - EMM_F89   0.73333 0.452 208   1.624  0.9622
    ##  SP_F14 - EMM_F3    0.51333 0.404 208   1.271  0.9965
    ##  SP_F14 - Control   0.48000 0.404 208   1.188  0.9983
    ##  SP_F14 - EMM_F63   0.32000 0.404 208   0.792  1.0000
    ##  SP_F14 - EMM_F48   0.31000 0.404 208   0.768  1.0000
    ##  SP_F14 - EMM_F70   0.29333 0.404 208   0.726  1.0000
    ##  SP_F14 - EMM_F66   0.24333 0.404 208   0.602  1.0000
    ##  SP_F14 - ZAN_F4    0.15667 0.404 208   0.388  1.0000
    ##  SP_F14 - EMM_F5    0.14667 0.404 208   0.363  1.0000
    ##  SP_F14 - EMM_F65   0.13667 0.404 208   0.338  1.0000
    ##  SP_F14 - EMM_F7    0.11000 0.404 208   0.272  1.0000
    ##  EMM_F64 - EMM_F49  1.23500 0.452 208   2.735  0.3160
    ##  EMM_F64 - EMM_F34  1.11333 0.404 208   2.757  0.3029
    ##  EMM_F64 - ZAN_F3   1.09000 0.404 208   2.699  0.3386
    ##  EMM_F64 - EMM_F89  1.09000 0.452 208   2.414  0.5394
    ##  EMM_F64 - EMM_F3   0.87000 0.404 208   2.154  0.7276
    ##  EMM_F64 - Control  0.83667 0.404 208   2.072  0.7808
    ##  EMM_F64 - EMM_F63  0.67667 0.404 208   1.675  0.9509
    ##  EMM_F64 - EMM_F48  0.66667 0.404 208   1.651  0.9567
    ##  EMM_F64 - EMM_F70  0.65000 0.404 208   1.609  0.9651
    ##  EMM_F64 - EMM_F66  0.60000 0.404 208   1.486  0.9831
    ##  EMM_F64 - ZAN_F4   0.51333 0.404 208   1.271  0.9965
    ##  EMM_F64 - EMM_F5   0.50333 0.404 208   1.246  0.9972
    ##  EMM_F64 - EMM_F65  0.49333 0.404 208   1.222  0.9977
    ##  EMM_F64 - EMM_F7   0.46667 0.404 208   1.155  0.9988
    ##  EMM_F64 - SP_F14   0.35667 0.404 208   0.883  1.0000
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14    0.3167 0.404 208   0.784  1.0000
    ##  Control - EMM_F34   0.0800 0.404 208   0.198  1.0000
    ##  EMM_F70 - Control   0.0733 0.404 208   0.182  1.0000
    ##  ZAN_F4 - Control    0.3100 0.404 208   0.768  1.0000
    ##  EMM_F89 - Control   0.3367 0.452 208   0.746  1.0000
    ##  EMM_F64 - Control   0.3500 0.404 208   0.867  1.0000
    ##  EMM_F49 - Control   0.4567 0.452 208   1.011  0.9997
    ##  EMM_F66 - Control   0.4800 0.404 208   1.188  0.9983
    ##  EMM_F7 - Control    0.4900 0.404 208   1.213  0.9979
    ##  EMM_F65 - Control   0.5467 0.404 208   1.354  0.9932
    ##  EMM_F3 - Control    0.5467 0.404 208   1.354  0.9932
    ##  EMM_F48 - Control   0.7033 0.404 208   1.741  0.9331
    ##  EMM_F63 - Control   0.8300 0.404 208   2.055  0.7908
    ##  EMM_F5 - Control    0.9800 0.404 208   2.427  0.5301
    ##  ZAN_F3 - Control    1.1367 0.404 208   2.814  0.2694
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F89   0.1217 0.452 208   0.269  1.0000
    ##  Control - EMM_F34   0.0567 0.404 208   0.140  1.0000
    ##  EMM_F3 - Control    0.0667 0.404 208   0.165  1.0000
    ##  ZAN_F3 - Control    0.1467 0.404 208   0.363  1.0000
    ##  SP_F14 - Control    0.1600 0.404 208   0.396  1.0000
    ##  EMM_F49 - Control   0.1783 0.452 208   0.395  1.0000
    ##  EMM_F70 - Control   0.1833 0.404 208   0.454  1.0000
    ##  ZAN_F4 - Control    0.1867 0.404 208   0.462  1.0000
    ##  EMM_F66 - Control   0.1967 0.404 208   0.487  1.0000
    ##  EMM_F64 - Control   0.3067 0.404 208   0.759  1.0000
    ##  EMM_F7 - Control    0.3700 0.404 208   0.916  0.9999
    ##  EMM_F63 - Control   0.5833 0.404 208   1.444  0.9871
    ##  EMM_F5 - Control    0.5867 0.404 208   1.453  0.9863
    ##  EMM_F48 - Control   0.7767 0.404 208   1.923  0.8623
    ##  EMM_F65 - Control   0.9133 0.404 208   2.261  0.6524
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34   0.0733 0.404 208   0.182  1.0000
    ##  Control - EMM_F66   0.0667 0.404 208   0.165  1.0000
    ##  EMM_F49 - Control   0.0150 0.452 208   0.033  1.0000
    ##  EMM_F89 - Control   0.0200 0.452 208   0.044  1.0000
    ##  EMM_F70 - Control   0.1267 0.404 208   0.314  1.0000
    ##  ZAN_F3 - Control    0.1400 0.404 208   0.347  1.0000
    ##  EMM_F3 - Control    0.1500 0.404 208   0.371  1.0000
    ##  ZAN_F4 - Control    0.3967 0.404 208   0.982  0.9998
    ##  SP_F14 - Control    0.4167 0.404 208   1.032  0.9997
    ##  EMM_F5 - Control    0.4900 0.404 208   1.213  0.9979
    ##  EMM_F64 - Control   0.5333 0.404 208   1.321  0.9948
    ##  EMM_F65 - Control   0.5367 0.404 208   1.329  0.9944
    ##  EMM_F63 - Control   0.5833 0.404 208   1.444  0.9871
    ##  EMM_F7 - Control    0.5967 0.404 208   1.477  0.9839
    ##  EMM_F48 - Control   0.6700 0.404 208   1.659  0.9548
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65   0.0900 0.404 208   0.223  1.0000
    ##  EMM_F89 - Control   0.0433 0.452 208   0.096  1.0000
    ##  EMM_F49 - Control   0.0783 0.452 208   0.173  1.0000
    ##  EMM_F3 - Control    0.0833 0.404 208   0.206  1.0000
    ##  EMM_F34 - Control   0.1167 0.404 208   0.289  1.0000
    ##  EMM_F7 - Control    0.2133 0.404 208   0.528  1.0000
    ##  EMM_F5 - Control    0.2200 0.404 208   0.545  1.0000
    ##  EMM_F66 - Control   0.2267 0.404 208   0.561  1.0000
    ##  EMM_F70 - Control   0.2300 0.404 208   0.569  1.0000
    ##  SP_F14 - Control    0.2367 0.404 208   0.586  1.0000
    ##  ZAN_F4 - Control    0.4967 0.404 208   1.230  0.9976
    ##  EMM_F48 - Control   0.5033 0.404 208   1.246  0.9972
    ##  ZAN_F3 - Control    0.6233 0.404 208   1.543  0.9759
    ##  EMM_F64 - Control   0.6667 0.404 208   1.651  0.9567
    ##  EMM_F63 - Control   0.7033 0.404 208   1.741  0.9331
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49   0.3750 0.452 208   0.830  1.0000
    ##  Control - EMM_F70   0.1033 0.404 208   0.256  1.0000
    ##  Control - EMM_F66   0.0233 0.404 208   0.058  1.0000
    ##  EMM_F48 - Control   0.0400 0.404 208   0.099  1.0000
    ##  EMM_F3 - Control    0.0833 0.404 208   0.206  1.0000
    ##  EMM_F5 - Control    0.1000 0.404 208   0.248  1.0000
    ##  EMM_F34 - Control   0.1200 0.404 208   0.297  1.0000
    ##  ZAN_F4 - Control    0.1367 0.404 208   0.338  1.0000
    ##  EMM_F89 - Control   0.2600 0.452 208   0.576  1.0000
    ##  EMM_F64 - Control   0.3233 0.404 208   0.801  1.0000
    ##  EMM_F7 - Control    0.3267 0.404 208   0.809  1.0000
    ##  EMM_F65 - Control   0.4300 0.404 208   1.065  0.9995
    ##  SP_F14 - Control    0.5867 0.404 208   1.453  0.9863
    ##  EMM_F63 - Control   0.5967 0.404 208   1.477  0.9839
    ##  ZAN_F3 - Control    0.6500 0.404 208   1.609  0.9651
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49   0.2417 0.452 208   0.535  1.0000
    ##  Control - EMM_F48   0.0367 0.404 208   0.091  1.0000
    ##  EMM_F70 - Control   0.1567 0.404 208   0.388  1.0000
    ##  EMM_F3 - Control    0.1767 0.404 208   0.437  1.0000
    ##  EMM_F7 - Control    0.2333 0.404 208   0.578  1.0000
    ##  EMM_F63 - Control   0.2367 0.404 208   0.586  1.0000
    ##  EMM_F34 - Control   0.2900 0.404 208   0.718  1.0000
    ##  EMM_F5 - Control    0.3467 0.404 208   0.858  1.0000
    ##  ZAN_F4 - Control    0.4133 0.404 208   1.023  0.9997
    ##  SP_F14 - Control    0.4400 0.404 208   1.089  0.9994
    ##  EMM_F66 - Control   0.4633 0.404 208   1.147  0.9989
    ##  EMM_F64 - Control   0.5500 0.404 208   1.362  0.9928
    ##  EMM_F65 - Control   0.7267 0.404 208   1.799  0.9141
    ##  EMM_F89 - Control   0.8783 0.452 208   1.945  0.8515
    ##  ZAN_F3 - Control    1.1100 0.404 208   2.748  0.3078
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49   0.3983 0.452 208   0.882  1.0000
    ##  Control - EMM_F34   0.2767 0.404 208   0.685  1.0000
    ##  Control - ZAN_F3    0.2533 0.404 208   0.627  1.0000
    ##  Control - EMM_F89   0.2533 0.452 208   0.561  1.0000
    ##  Control - EMM_F3    0.0333 0.404 208   0.083  1.0000
    ##  EMM_F63 - Control   0.1600 0.404 208   0.396  1.0000
    ##  EMM_F48 - Control   0.1700 0.404 208   0.421  1.0000
    ##  EMM_F70 - Control   0.1867 0.404 208   0.462  1.0000
    ##  EMM_F66 - Control   0.2367 0.404 208   0.586  1.0000
    ##  ZAN_F4 - Control    0.3233 0.404 208   0.801  1.0000
    ##  EMM_F5 - Control    0.3333 0.404 208   0.825  1.0000
    ##  EMM_F65 - Control   0.3433 0.404 208   0.850  1.0000
    ##  EMM_F7 - Control    0.3700 0.404 208   0.916  0.9999
    ##  SP_F14 - Control    0.4800 0.404 208   1.188  0.9983
    ##  EMM_F64 - Control   0.8367 0.404 208   2.072  0.7808
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## $emmeans
    ## distance_to_yeast = 11:
    ##  Yeast    emmean    SE df lower.CL upper.CL .group
    ##  ZAN_F3   1.3980 0.352  2  -0.1176     2.91  a    
    ##  EMM_F5   1.2367 0.288  2  -0.0027     2.48  a    
    ##  EMM_F48  1.1330 0.352  2  -0.3826     2.65  a    
    ##  EMM_F7   0.9133 0.288  2  -0.3260     2.15  a    
    ##  EMM_F3   0.7867 0.288  2  -0.4527     2.03  a    
    ##  EMM_F49  0.7290 0.352  2  -0.7865     2.24  a    
    ##  EMM_F66  0.7200 0.288  2  -0.5194     1.96  a    
    ##  EMM_F63  0.6233 0.288  2  -0.6160     1.86  a    
    ##  Control  0.6200 0.288  2  -0.6194     1.86  a    
    ##  EMM_F34  0.5880 0.352  2  -0.9276     2.10  a    
    ##  EMM_F89  0.5330 0.352  2  -0.9826     2.05  a    
    ##  EMM_F70  0.4900 0.288  2  -0.7494     1.73  a    
    ##  EMM_F64  0.4833 0.288  2  -0.7560     1.72  a    
    ##  ZAN_F4   0.3833 0.288  2  -0.8560     1.62  a    
    ##  EMM_F65  0.2080 0.352  2  -1.3076     1.72  a    
    ##  SP_F14   0.0733 0.288  2  -1.1660     1.31  a    
    ## 
    ## distance_to_yeast = 17:
    ##  Yeast    emmean    SE df lower.CL upper.CL .group
    ##  EMM_F48  1.3180 0.352  2  -0.1976     2.83  a    
    ##  EMM_F5   1.0500 0.288  2  -0.1894     2.29  a    
    ##  SP_F14   0.8700 0.288  2  -0.3694     2.11  a    
    ##  EMM_F49  0.8590 0.352  2  -0.6565     2.37  a    
    ##  EMM_F63  0.7600 0.288  2  -0.4794     2.00  a    
    ##  EMM_F7   0.7433 0.288  2  -0.4960     1.98  a    
    ##  EMM_F64  0.7367 0.288  2  -0.5027     1.98  a    
    ##  EMM_F65  0.7180 0.352  2  -0.7976     2.23  a    
    ##  EMM_F89  0.6130 0.352  2  -0.9026     2.13  a    
    ##  ZAN_F4   0.6033 0.288  2  -0.6360     1.84  a    
    ##  EMM_F70  0.5833 0.288  2  -0.6560     1.82  a    
    ##  EMM_F66  0.5700 0.288  2  -0.6694     1.81  a    
    ##  EMM_F3   0.5267 0.288  2  -0.7127     1.77  a    
    ##  EMM_F34  0.5180 0.352  2  -0.9976     2.03  a    
    ##  Control  0.4533 0.288  2  -0.7860     1.69  a    
    ##  ZAN_F3   0.3880 0.352  2  -1.1276     1.90  a    
    ## 
    ## distance_to_yeast = 25:
    ##  Yeast    emmean    SE df lower.CL upper.CL .group
    ##  EMM_F7   1.1100 0.288  2  -0.1294     2.35  a    
    ##  EMM_F64  1.0233 0.288  2  -0.2160     2.26  a    
    ##  EMM_F5   0.9900 0.288  2  -0.2494     2.23  a    
    ##  EMM_F63  0.9467 0.288  2  -0.2927     2.19  a    
    ##  SP_F14   0.9000 0.288  2  -0.3394     2.14  a    
    ##  EMM_F65  0.8980 0.352  2  -0.6176     2.41  a    
    ##  EMM_F48  0.7880 0.352  2  -0.7276     2.30  a    
    ##  ZAN_F4   0.6933 0.288  2  -0.5460     1.93  a    
    ##  EMM_F49  0.5290 0.352  2  -0.9865     2.04  a    
    ##  EMM_F3   0.5133 0.288  2  -0.7260     1.75  a    
    ##  EMM_F70  0.4900 0.288  2  -0.7494     1.73  a    
    ##  EMM_F89  0.4530 0.352  2  -1.0626     1.97  a    
    ##  EMM_F34  0.4480 0.352  2  -1.0676     1.96  a    
    ##  Control  0.4433 0.288  2  -0.7960     1.68  a    
    ##  EMM_F66  0.4433 0.288  2  -0.7960     1.68  a    
    ##  ZAN_F3   0.2830 0.352  2  -1.2326     1.80  a    
    ## 
    ## distance_to_yeast = 32:
    ##  Yeast    emmean    SE df lower.CL upper.CL .group
    ##  EMM_F64  1.2900 0.288  2   0.0506     2.53  a    
    ##  EMM_F63  1.2567 0.288  2   0.0173     2.50  a    
    ##  ZAN_F4   1.0467 0.288  2  -0.1927     2.29  a    
    ##  EMM_F89  1.0180 0.352  2  -0.4976     2.53  a    
    ##  EMM_F7   1.0000 0.288  2  -0.2394     2.24  a    
    ##  EMM_F48  0.9330 0.352  2  -0.5826     2.45  a    
    ##  EMM_F5   0.8167 0.288  2  -0.4227     2.06  a    
    ##  SP_F14   0.7367 0.288  2  -0.5027     1.98  a    
    ##  ZAN_F3   0.7180 0.352  2  -0.7976     2.23  a    
    ##  EMM_F66  0.7100 0.288  2  -0.5294     1.95  a    
    ##  EMM_F70  0.6300 0.288  2  -0.6094     1.87  a    
    ##  EMM_F65  0.5730 0.352  2  -0.9426     2.09  a    
    ##  EMM_F3   0.5033 0.288  2  -0.7360     1.74  a    
    ##  Control  0.4733 0.288  2  -0.7660     1.71  a    
    ##  EMM_F49  0.4240 0.352  2  -1.0915     1.94  a    
    ##  EMM_F34  0.4230 0.352  2  -1.0926     1.94  a    
    ## 
    ## distance_to_yeast = 41:
    ##  Yeast    emmean    SE df lower.CL upper.CL .group
    ##  EMM_F63  1.1867 0.288  2  -0.0527     2.43  a    
    ##  EMM_F64  1.1233 0.288  2  -0.1160     2.36  a    
    ##  EMM_F7   1.0600 0.288  2  -0.1794     2.30  a    
    ##  ZAN_F4   0.9900 0.288  2  -0.2494     2.23  a    
    ##  SP_F14   0.9600 0.288  2  -0.2794     2.20  a    
    ##  EMM_F48  0.9430 0.352  2  -0.5726     2.46  a    
    ##  EMM_F89  0.8880 0.352  2  -0.6276     2.40  a    
    ##  ZAN_F3   0.8680 0.352  2  -0.6476     2.38  a    
    ##  EMM_F5   0.8200 0.288  2  -0.4194     2.06  a    
    ##  EMM_F65  0.7880 0.352  2  -0.7276     2.30  a    
    ##  EMM_F34  0.7630 0.352  2  -0.7526     2.28  a    
    ##  EMM_F70  0.7333 0.288  2  -0.5060     1.97  a    
    ##  EMM_F3   0.6533 0.288  2  -0.5860     1.89  a    
    ##  EMM_F66  0.5700 0.288  2  -0.6694     1.81  a    
    ##  Control  0.4333 0.288  2  -0.8060     1.67  a    
    ##  EMM_F49  0.2090 0.352  2  -1.3065     1.72  a    
    ## 
    ## distance_to_yeast = 48:
    ##  Yeast    emmean    SE df lower.CL upper.CL .group
    ##  EMM_F65  1.3830 0.352  2  -0.1326     2.90  a    
    ##  ZAN_F3   1.3780 0.352  2  -0.1376     2.89  a    
    ##  EMM_F64  1.0667 0.288  2  -0.1727     2.31  a    
    ##  EMM_F63  1.0200 0.288  2  -0.2194     2.26  a    
    ##  EMM_F70  1.0133 0.288  2  -0.2260     2.25  a    
    ##  SP_F14   0.9700 0.288  2  -0.2694     2.21  a    
    ##  ZAN_F4   0.9667 0.288  2  -0.2727     2.21  a    
    ##  EMM_F5   0.9633 0.288  2  -0.2760     2.20  a    
    ##  EMM_F66  0.9433 0.288  2  -0.2960     2.18  a    
    ##  EMM_F89  0.9430 0.352  2  -0.5726     2.46  a    
    ##  EMM_F34  0.9280 0.352  2  -0.5876     2.44  a    
    ##  EMM_F48  0.8680 0.352  2  -0.6476     2.38  a    
    ##  EMM_F7   0.7900 0.288  2  -0.4494     2.03  a    
    ##  EMM_F3   0.7200 0.288  2  -0.5194     1.96  a    
    ##  EMM_F49  0.4440 0.352  2  -1.0715     1.96  a    
    ##  Control  0.3100 0.288  2  -0.9294     1.55  a    
    ## 
    ## distance_to_yeast = 55:
    ##  Yeast    emmean    SE df lower.CL upper.CL .group
    ##  EMM_F64  1.5367 0.288  2   0.2973     2.78  a    
    ##  EMM_F48  1.4480 0.352  2  -0.0676     2.96  ab   
    ##  SP_F14   1.2200 0.288  2  -0.0194     2.46  ab   
    ##  ZAN_F4   1.1500 0.288  2  -0.0894     2.39  ab   
    ##  EMM_F5   1.0900 0.288  2  -0.1494     2.33  ab   
    ##  EMM_F66  1.0000 0.288  2  -0.2394     2.24  ab   
    ##  EMM_F65  0.9580 0.352  2  -0.5576     2.47  ab   
    ##  EMM_F7   0.9133 0.288  2  -0.3260     2.15  ab   
    ##  EMM_F63  0.8933 0.288  2  -0.3460     2.13  ab   
    ##  EMM_F3   0.7767 0.288  2  -0.4627     2.02  ab   
    ##  EMM_F70  0.7567 0.288  2  -0.4827     2.00  ab   
    ##  EMM_F34  0.6580 0.352  2  -0.8576     2.17  ab   
    ##  EMM_F89  0.6530 0.352  2  -0.8626     2.17  ab   
    ##  Control  0.3000 0.288  2  -0.9394     1.54  ab   
    ##  EMM_F49  0.2940 0.352  2  -1.2215     1.81  ab   
    ##  ZAN_F3  -0.0970 0.352  2  -1.6126     1.42   b   
    ## 
    ## Degrees-of-freedom method: containment 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F65 - SP_F14  0.134695 0.452 180   0.298  1.0000
    ##  ZAN_F4 - SP_F14   0.310000 0.404 180   0.767  1.0000
    ##  ZAN_F4 - EMM_F65  0.175305 0.452 180   0.388  1.0000
    ##  EMM_F64 - SP_F14  0.410000 0.404 180   1.014  0.9997
    ##  EMM_F64 - EMM_F65 0.275305 0.452 180   0.609  1.0000
    ##  EMM_F64 - ZAN_F4  0.100000 0.404 180   0.247  1.0000
    ##  EMM_F70 - SP_F14  0.416667 0.404 180   1.031  0.9997
    ##  EMM_F70 - EMM_F65 0.281972 0.452 180   0.624  1.0000
    ##  EMM_F70 - ZAN_F4  0.106667 0.404 180   0.264  1.0000
    ##  EMM_F70 - EMM_F64 0.006667 0.404 180   0.016  1.0000
    ##  EMM_F89 - SP_F14  0.459695 0.452 180   1.017  0.9997
    ##  EMM_F89 - EMM_F65 0.325000 0.495 180   0.657  1.0000
    ##  EMM_F89 - ZAN_F4  0.149695 0.452 180   0.331  1.0000
    ##  EMM_F89 - EMM_F64 0.049695 0.452 180   0.110  1.0000
    ##  EMM_F89 - EMM_F70 0.043028 0.452 180   0.095  1.0000
    ##  EMM_F34 - SP_F14  0.514623 0.452 180   1.138  0.9990
    ##  EMM_F34 - EMM_F65 0.379928 0.496 180   0.766  1.0000
    ##  EMM_F34 - ZAN_F4  0.204623 0.452 180   0.453  1.0000
    ##  EMM_F34 - EMM_F64 0.104623 0.452 180   0.231  1.0000
    ##  EMM_F34 - EMM_F70 0.097956 0.452 180   0.217  1.0000
    ##  EMM_F34 - EMM_F89 0.054928 0.496 180   0.111  1.0000
    ##  Control - SP_F14  0.546667 0.404 180   1.353  0.9932
    ##  Control - EMM_F65 0.411972 0.452 180   0.911  0.9999
    ##  Control - ZAN_F4  0.236667 0.404 180   0.586  1.0000
    ##  Control - EMM_F64 0.136667 0.404 180   0.338  1.0000
    ##  Control - EMM_F70 0.130000 0.404 180   0.322  1.0000
    ##  Control - EMM_F89 0.086972 0.452 180   0.192  1.0000
    ##  Control - EMM_F34 0.032044 0.452 180   0.071  1.0000
    ##  EMM_F63 - SP_F14  0.550000 0.404 180   1.361  0.9928
    ##  EMM_F63 - EMM_F65 0.415305 0.452 180   0.918  0.9999
    ##  EMM_F63 - ZAN_F4  0.240000 0.404 180   0.594  1.0000
    ##  EMM_F63 - EMM_F64 0.140000 0.404 180   0.346  1.0000
    ##  EMM_F63 - EMM_F70 0.133333 0.404 180   0.330  1.0000
    ##  EMM_F63 - EMM_F89 0.090305 0.452 180   0.200  1.0000
    ##  EMM_F63 - EMM_F34 0.035377 0.452 180   0.078  1.0000
    ##  EMM_F63 - Control 0.003333 0.404 180   0.008  1.0000
    ##  EMM_F66 - SP_F14  0.646667 0.404 180   1.600  0.9665
    ##  EMM_F66 - EMM_F65 0.511972 0.452 180   1.132  0.9990
    ##  EMM_F66 - ZAN_F4  0.336667 0.404 180   0.833  1.0000
    ##  EMM_F66 - EMM_F64 0.236667 0.404 180   0.586  1.0000
    ##  EMM_F66 - EMM_F70 0.230000 0.404 180   0.569  1.0000
    ##  EMM_F66 - EMM_F89 0.186972 0.452 180   0.413  1.0000
    ##  EMM_F66 - EMM_F34 0.132044 0.452 180   0.292  1.0000
    ##  EMM_F66 - Control 0.100000 0.404 180   0.247  1.0000
    ##  EMM_F66 - EMM_F63 0.096667 0.404 180   0.239  1.0000
    ##  EMM_F49 - SP_F14  0.655683 0.452 180   1.450  0.9864
    ##  EMM_F49 - EMM_F65 0.520988 0.496 180   1.051  0.9996
    ##  EMM_F49 - ZAN_F4  0.345683 0.452 180   0.764  1.0000
    ##  EMM_F49 - EMM_F64 0.245683 0.452 180   0.543  1.0000
    ##  EMM_F49 - EMM_F70 0.239016 0.452 180   0.529  1.0000
    ##  EMM_F49 - EMM_F89 0.195988 0.496 180   0.395  1.0000
    ##  EMM_F49 - EMM_F34 0.141060 0.496 180   0.285  1.0000
    ##  EMM_F49 - Control 0.109016 0.452 180   0.241  1.0000
    ##  EMM_F49 - EMM_F63 0.105683 0.452 180   0.234  1.0000
    ##  EMM_F49 - EMM_F66 0.009016 0.452 180   0.020  1.0000
    ##  EMM_F3 - SP_F14   0.713333 0.404 180   1.765  0.9253
    ##  EMM_F3 - EMM_F65  0.578639 0.452 180   1.280  0.9962
    ##  EMM_F3 - ZAN_F4   0.403333 0.404 180   0.998  0.9998
    ##  EMM_F3 - EMM_F64  0.303333 0.404 180   0.750  1.0000
    ##  EMM_F3 - EMM_F70  0.296667 0.404 180   0.734  1.0000
    ##  EMM_F3 - EMM_F89  0.253639 0.452 180   0.561  1.0000
    ##  EMM_F3 - EMM_F34  0.198711 0.452 180   0.439  1.0000
    ##  EMM_F3 - Control  0.166667 0.404 180   0.412  1.0000
    ##  EMM_F3 - EMM_F63  0.163333 0.404 180   0.404  1.0000
    ##  EMM_F3 - EMM_F66  0.066667 0.404 180   0.165  1.0000
    ##  EMM_F3 - EMM_F49  0.057651 0.452 180   0.127  1.0000
    ##  EMM_F7 - SP_F14   0.840000 0.404 180   2.078  0.7763
    ##  EMM_F7 - EMM_F65  0.705305 0.452 180   1.560  0.9732
    ##  EMM_F7 - ZAN_F4   0.530000 0.404 180   1.311  0.9951
    ##  EMM_F7 - EMM_F64  0.430000 0.404 180   1.064  0.9995
    ##  EMM_F7 - EMM_F70  0.423333 0.404 180   1.047  0.9996
    ##  EMM_F7 - EMM_F89  0.380305 0.452 180   0.841  1.0000
    ##  EMM_F7 - EMM_F34  0.325377 0.452 180   0.720  1.0000
    ##  EMM_F7 - Control  0.293333 0.404 180   0.726  1.0000
    ##  EMM_F7 - EMM_F63  0.290000 0.404 180   0.717  1.0000
    ##  EMM_F7 - EMM_F66  0.193333 0.404 180   0.478  1.0000
    ##  EMM_F7 - EMM_F49  0.184317 0.452 180   0.408  1.0000
    ##  EMM_F7 - EMM_F3   0.126667 0.404 180   0.313  1.0000
    ##  EMM_F48 - SP_F14  1.059695 0.452 180   2.343  0.5922
    ##  EMM_F48 - EMM_F65 0.925000 0.495 180   1.869  0.8865
    ##  EMM_F48 - ZAN_F4  0.749695 0.452 180   1.658  0.9547
    ##  EMM_F48 - EMM_F64 0.649695 0.452 180   1.437  0.9876
    ##  EMM_F48 - EMM_F70 0.643028 0.452 180   1.422  0.9888
    ##  EMM_F48 - EMM_F89 0.600000 0.495 180   1.212  0.9979
    ##  EMM_F48 - EMM_F34 0.545072 0.496 180   1.099  0.9993
    ##  EMM_F48 - Control 0.513028 0.452 180   1.135  0.9990
    ##  EMM_F48 - EMM_F63 0.509695 0.452 180   1.127  0.9991
    ##  EMM_F48 - EMM_F66 0.413028 0.452 180   0.913  0.9999
    ##  EMM_F48 - EMM_F49 0.404012 0.496 180   0.815  1.0000
    ##  EMM_F48 - EMM_F3  0.346361 0.452 180   0.766  1.0000
    ##  EMM_F48 - EMM_F7  0.219695 0.452 180   0.486  1.0000
    ##  EMM_F5 - SP_F14   1.163333 0.404 180   2.878  0.2366
    ##  EMM_F5 - EMM_F65  1.028639 0.452 180   2.275  0.6427
    ##  EMM_F5 - ZAN_F4   0.853333 0.404 180   2.111  0.7555
    ##  EMM_F5 - EMM_F64  0.753333 0.404 180   1.864  0.8885
    ##  EMM_F5 - EMM_F70  0.746667 0.404 180   1.847  0.8953
    ##  EMM_F5 - EMM_F89  0.703639 0.452 180   1.556  0.9738
    ##  EMM_F5 - EMM_F34  0.648711 0.452 180   1.435  0.9877
    ##  EMM_F5 - Control  0.616667 0.404 180   1.526  0.9781
    ##  EMM_F5 - EMM_F63  0.613333 0.404 180   1.517  0.9792
    ##  EMM_F5 - EMM_F66  0.516667 0.404 180   1.278  0.9962
    ##  EMM_F5 - EMM_F49  0.507651 0.452 180   1.123  0.9991
    ##  EMM_F5 - EMM_F3   0.450000 0.404 180   1.113  0.9992
    ##  EMM_F5 - EMM_F7   0.323333 0.404 180   0.800  1.0000
    ##  EMM_F5 - EMM_F48  0.103639 0.452 180   0.229  1.0000
    ##  ZAN_F3 - SP_F14   1.324623 0.452 180   2.929  0.2114
    ##  ZAN_F3 - EMM_F65  1.189928 0.496 180   2.400  0.5503
    ##  ZAN_F3 - ZAN_F4   1.014623 0.452 180   2.244  0.6650
    ##  ZAN_F3 - EMM_F64  0.914623 0.452 180   2.023  0.8094
    ##  ZAN_F3 - EMM_F70  0.907956 0.452 180   2.008  0.8177
    ##  ZAN_F3 - EMM_F89  0.864928 0.496 180   1.744  0.9318
    ##  ZAN_F3 - EMM_F34  0.810000 0.495 180   1.636  0.9594
    ##  ZAN_F3 - Control  0.777956 0.452 180   1.720  0.9388
    ##  ZAN_F3 - EMM_F63  0.774623 0.452 180   1.713  0.9408
    ##  ZAN_F3 - EMM_F66  0.677956 0.452 180   1.499  0.9814
    ##  ZAN_F3 - EMM_F49  0.668940 0.496 180   1.349  0.9934
    ##  ZAN_F3 - EMM_F3   0.611289 0.452 180   1.352  0.9932
    ##  ZAN_F3 - EMM_F7   0.484623 0.452 180   1.072  0.9995
    ##  ZAN_F3 - EMM_F48  0.264928 0.496 180   0.534  1.0000
    ##  ZAN_F3 - EMM_F5   0.161289 0.452 180   0.357  1.0000
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F3  0.065377 0.452 180   0.145  1.0000
    ##  EMM_F34 - ZAN_F3  0.130000 0.495 180   0.263  1.0000
    ##  EMM_F34 - Control 0.064623 0.452 180   0.143  1.0000
    ##  EMM_F3 - ZAN_F3   0.138711 0.452 180   0.307  1.0000
    ##  EMM_F3 - Control  0.073333 0.404 180   0.181  1.0000
    ##  EMM_F3 - EMM_F34  0.008711 0.452 180   0.019  1.0000
    ##  EMM_F66 - ZAN_F3  0.182044 0.452 180   0.403  1.0000
    ##  EMM_F66 - Control 0.116667 0.404 180   0.289  1.0000
    ##  EMM_F66 - EMM_F34 0.052044 0.452 180   0.115  1.0000
    ##  EMM_F66 - EMM_F3  0.043333 0.404 180   0.107  1.0000
    ##  EMM_F70 - ZAN_F3  0.195377 0.452 180   0.432  1.0000
    ##  EMM_F70 - Control 0.130000 0.404 180   0.322  1.0000
    ##  EMM_F70 - EMM_F34 0.065377 0.452 180   0.145  1.0000
    ##  EMM_F70 - EMM_F3  0.056667 0.404 180   0.140  1.0000
    ##  EMM_F70 - EMM_F66 0.013333 0.404 180   0.033  1.0000
    ##  ZAN_F4 - ZAN_F3   0.215377 0.452 180   0.476  1.0000
    ##  ZAN_F4 - Control  0.150000 0.404 180   0.371  1.0000
    ##  ZAN_F4 - EMM_F34  0.085377 0.452 180   0.189  1.0000
    ##  ZAN_F4 - EMM_F3   0.076667 0.404 180   0.190  1.0000
    ##  ZAN_F4 - EMM_F66  0.033333 0.404 180   0.082  1.0000
    ##  ZAN_F4 - EMM_F70  0.020000 0.404 180   0.049  1.0000
    ##  EMM_F89 - ZAN_F3  0.225072 0.496 180   0.454  1.0000
    ##  EMM_F89 - Control 0.159695 0.452 180   0.353  1.0000
    ##  EMM_F89 - EMM_F34 0.095072 0.496 180   0.192  1.0000
    ##  EMM_F89 - EMM_F3  0.086361 0.452 180   0.191  1.0000
    ##  EMM_F89 - EMM_F66 0.043028 0.452 180   0.095  1.0000
    ##  EMM_F89 - EMM_F70 0.029695 0.452 180   0.066  1.0000
    ##  EMM_F89 - ZAN_F4  0.009695 0.452 180   0.021  1.0000
    ##  EMM_F65 - ZAN_F3  0.330072 0.496 180   0.666  1.0000
    ##  EMM_F65 - Control 0.264695 0.452 180   0.585  1.0000
    ##  EMM_F65 - EMM_F34 0.200072 0.496 180   0.403  1.0000
    ##  EMM_F65 - EMM_F3  0.191361 0.452 180   0.423  1.0000
    ##  EMM_F65 - EMM_F66 0.148028 0.452 180   0.327  1.0000
    ##  EMM_F65 - EMM_F70 0.134695 0.452 180   0.298  1.0000
    ##  EMM_F65 - ZAN_F4  0.114695 0.452 180   0.254  1.0000
    ##  EMM_F65 - EMM_F89 0.105000 0.495 180   0.212  1.0000
    ##  EMM_F64 - ZAN_F3  0.348711 0.452 180   0.771  1.0000
    ##  EMM_F64 - Control 0.283333 0.404 180   0.701  1.0000
    ##  EMM_F64 - EMM_F34 0.218711 0.452 180   0.484  1.0000
    ##  EMM_F64 - EMM_F3  0.210000 0.404 180   0.520  1.0000
    ##  EMM_F64 - EMM_F66 0.166667 0.404 180   0.412  1.0000
    ##  EMM_F64 - EMM_F70 0.153333 0.404 180   0.379  1.0000
    ##  EMM_F64 - ZAN_F4  0.133333 0.404 180   0.330  1.0000
    ##  EMM_F64 - EMM_F89 0.123639 0.452 180   0.273  1.0000
    ##  EMM_F64 - EMM_F65 0.018639 0.452 180   0.041  1.0000
    ##  EMM_F7 - ZAN_F3   0.355377 0.452 180   0.786  1.0000
    ##  EMM_F7 - Control  0.290000 0.404 180   0.717  1.0000
    ##  EMM_F7 - EMM_F34  0.225377 0.452 180   0.498  1.0000
    ##  EMM_F7 - EMM_F3   0.216667 0.404 180   0.536  1.0000
    ##  EMM_F7 - EMM_F66  0.173333 0.404 180   0.429  1.0000
    ##  EMM_F7 - EMM_F70  0.160000 0.404 180   0.396  1.0000
    ##  EMM_F7 - ZAN_F4   0.140000 0.404 180   0.346  1.0000
    ##  EMM_F7 - EMM_F89  0.130305 0.452 180   0.288  1.0000
    ##  EMM_F7 - EMM_F65  0.025305 0.452 180   0.056  1.0000
    ##  EMM_F7 - EMM_F64  0.006667 0.404 180   0.016  1.0000
    ##  EMM_F63 - ZAN_F3  0.372044 0.452 180   0.823  1.0000
    ##  EMM_F63 - Control 0.306667 0.404 180   0.759  1.0000
    ##  EMM_F63 - EMM_F34 0.242044 0.452 180   0.535  1.0000
    ##  EMM_F63 - EMM_F3  0.233333 0.404 180   0.577  1.0000
    ##  EMM_F63 - EMM_F66 0.190000 0.404 180   0.470  1.0000
    ##  EMM_F63 - EMM_F70 0.176667 0.404 180   0.437  1.0000
    ##  EMM_F63 - ZAN_F4  0.156667 0.404 180   0.388  1.0000
    ##  EMM_F63 - EMM_F89 0.146972 0.452 180   0.325  1.0000
    ##  EMM_F63 - EMM_F65 0.041972 0.452 180   0.093  1.0000
    ##  EMM_F63 - EMM_F64 0.023333 0.404 180   0.058  1.0000
    ##  EMM_F63 - EMM_F7  0.016667 0.404 180   0.041  1.0000
    ##  EMM_F49 - ZAN_F3  0.471060 0.496 180   0.950  0.9999
    ##  EMM_F49 - Control 0.405683 0.452 180   0.897  0.9999
    ##  EMM_F49 - EMM_F34 0.341060 0.496 180   0.688  1.0000
    ##  EMM_F49 - EMM_F3  0.332349 0.452 180   0.735  1.0000
    ##  EMM_F49 - EMM_F66 0.289016 0.452 180   0.639  1.0000
    ##  EMM_F49 - EMM_F70 0.275683 0.452 180   0.610  1.0000
    ##  EMM_F49 - ZAN_F4  0.255683 0.452 180   0.565  1.0000
    ##  EMM_F49 - EMM_F89 0.245988 0.496 180   0.496  1.0000
    ##  EMM_F49 - EMM_F65 0.140988 0.496 180   0.284  1.0000
    ##  EMM_F49 - EMM_F64 0.122349 0.452 180   0.271  1.0000
    ##  EMM_F49 - EMM_F7  0.115683 0.452 180   0.256  1.0000
    ##  EMM_F49 - EMM_F63 0.099016 0.452 180   0.219  1.0000
    ##  SP_F14 - ZAN_F3   0.482044 0.452 180   1.066  0.9995
    ##  SP_F14 - Control  0.416667 0.404 180   1.031  0.9997
    ##  SP_F14 - EMM_F34  0.352044 0.452 180   0.779  1.0000
    ##  SP_F14 - EMM_F3   0.343333 0.404 180   0.849  1.0000
    ##  SP_F14 - EMM_F66  0.300000 0.404 180   0.742  1.0000
    ##  SP_F14 - EMM_F70  0.286667 0.404 180   0.709  1.0000
    ##  SP_F14 - ZAN_F4   0.266667 0.404 180   0.660  1.0000
    ##  SP_F14 - EMM_F89  0.256972 0.452 180   0.568  1.0000
    ##  SP_F14 - EMM_F65  0.151972 0.452 180   0.336  1.0000
    ##  SP_F14 - EMM_F64  0.133333 0.404 180   0.330  1.0000
    ##  SP_F14 - EMM_F7   0.126667 0.404 180   0.313  1.0000
    ##  SP_F14 - EMM_F63  0.110000 0.404 180   0.272  1.0000
    ##  SP_F14 - EMM_F49  0.010984 0.452 180   0.024  1.0000
    ##  EMM_F5 - ZAN_F3   0.662044 0.452 180   1.464  0.9851
    ##  EMM_F5 - Control  0.596667 0.404 180   1.476  0.9839
    ##  EMM_F5 - EMM_F34  0.532044 0.452 180   1.177  0.9985
    ##  EMM_F5 - EMM_F3   0.523333 0.404 180   1.295  0.9957
    ##  EMM_F5 - EMM_F66  0.480000 0.404 180   1.188  0.9983
    ##  EMM_F5 - EMM_F70  0.466667 0.404 180   1.155  0.9988
    ##  EMM_F5 - ZAN_F4   0.446667 0.404 180   1.105  0.9993
    ##  EMM_F5 - EMM_F89  0.436972 0.452 180   0.966  0.9999
    ##  EMM_F5 - EMM_F65  0.331972 0.452 180   0.734  1.0000
    ##  EMM_F5 - EMM_F64  0.313333 0.404 180   0.775  1.0000
    ##  EMM_F5 - EMM_F7   0.306667 0.404 180   0.759  1.0000
    ##  EMM_F5 - EMM_F63  0.290000 0.404 180   0.717  1.0000
    ##  EMM_F5 - EMM_F49  0.190984 0.452 180   0.422  1.0000
    ##  EMM_F5 - SP_F14   0.180000 0.404 180   0.445  1.0000
    ##  EMM_F48 - ZAN_F3  0.930072 0.496 180   1.876  0.8834
    ##  EMM_F48 - Control 0.864695 0.452 180   1.912  0.8670
    ##  EMM_F48 - EMM_F34 0.800072 0.496 180   1.614  0.9640
    ##  EMM_F48 - EMM_F3  0.791361 0.452 180   1.750  0.9300
    ##  EMM_F48 - EMM_F66 0.748028 0.452 180   1.654  0.9555
    ##  EMM_F48 - EMM_F70 0.734695 0.452 180   1.625  0.9618
    ##  EMM_F48 - ZAN_F4  0.714695 0.452 180   1.581  0.9699
    ##  EMM_F48 - EMM_F89 0.705000 0.495 180   1.424  0.9886
    ##  EMM_F48 - EMM_F65 0.600000 0.495 180   1.212  0.9979
    ##  EMM_F48 - EMM_F64 0.581361 0.452 180   1.286  0.9960
    ##  EMM_F48 - EMM_F7  0.574695 0.452 180   1.271  0.9964
    ##  EMM_F48 - EMM_F63 0.558028 0.452 180   1.234  0.9974
    ##  EMM_F48 - EMM_F49 0.459012 0.496 180   0.926  0.9999
    ##  EMM_F48 - SP_F14  0.448028 0.452 180   0.991  0.9998
    ##  EMM_F48 - EMM_F5  0.268028 0.452 180   0.593  1.0000
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F66 - ZAN_F3  0.160377 0.452 180   0.355  1.0000
    ##  Control - ZAN_F3  0.160377 0.452 180   0.355  1.0000
    ##  Control - EMM_F66 0.000000 0.404 180   0.000  1.0000
    ##  EMM_F34 - ZAN_F3  0.165000 0.495 180   0.333  1.0000
    ##  EMM_F34 - EMM_F66 0.004623 0.452 180   0.010  1.0000
    ##  EMM_F34 - Control 0.004623 0.452 180   0.010  1.0000
    ##  EMM_F89 - ZAN_F3  0.170072 0.496 180   0.343  1.0000
    ##  EMM_F89 - EMM_F66 0.009695 0.452 180   0.021  1.0000
    ##  EMM_F89 - Control 0.009695 0.452 180   0.021  1.0000
    ##  EMM_F89 - EMM_F34 0.005072 0.496 180   0.010  1.0000
    ##  EMM_F70 - ZAN_F3  0.207044 0.452 180   0.458  1.0000
    ##  EMM_F70 - EMM_F66 0.046667 0.404 180   0.115  1.0000
    ##  EMM_F70 - Control 0.046667 0.404 180   0.115  1.0000
    ##  EMM_F70 - EMM_F34 0.042044 0.452 180   0.093  1.0000
    ##  EMM_F70 - EMM_F89 0.036972 0.452 180   0.082  1.0000
    ##  EMM_F3 - ZAN_F3   0.230377 0.452 180   0.509  1.0000
    ##  EMM_F3 - EMM_F66  0.070000 0.404 180   0.173  1.0000
    ##  EMM_F3 - Control  0.070000 0.404 180   0.173  1.0000
    ##  EMM_F3 - EMM_F34  0.065377 0.452 180   0.145  1.0000
    ##  EMM_F3 - EMM_F89  0.060305 0.452 180   0.133  1.0000
    ##  EMM_F3 - EMM_F70  0.023333 0.404 180   0.058  1.0000
    ##  EMM_F49 - ZAN_F3  0.246060 0.496 180   0.496  1.0000
    ##  EMM_F49 - EMM_F66 0.085683 0.452 180   0.189  1.0000
    ##  EMM_F49 - Control 0.085683 0.452 180   0.189  1.0000
    ##  EMM_F49 - EMM_F34 0.081060 0.496 180   0.163  1.0000
    ##  EMM_F49 - EMM_F89 0.075988 0.496 180   0.153  1.0000
    ##  EMM_F49 - EMM_F70 0.039016 0.452 180   0.086  1.0000
    ##  EMM_F49 - EMM_F3  0.015683 0.452 180   0.035  1.0000
    ##  ZAN_F4 - ZAN_F3   0.410377 0.452 180   0.908  0.9999
    ##  ZAN_F4 - EMM_F66  0.250000 0.404 180   0.619  1.0000
    ##  ZAN_F4 - Control  0.250000 0.404 180   0.619  1.0000
    ##  ZAN_F4 - EMM_F34  0.245377 0.452 180   0.543  1.0000
    ##  ZAN_F4 - EMM_F89  0.240305 0.452 180   0.531  1.0000
    ##  ZAN_F4 - EMM_F70  0.203333 0.404 180   0.503  1.0000
    ##  ZAN_F4 - EMM_F3   0.180000 0.404 180   0.445  1.0000
    ##  ZAN_F4 - EMM_F49  0.164317 0.452 180   0.363  1.0000
    ##  EMM_F48 - ZAN_F3  0.505072 0.496 180   1.019  0.9997
    ##  EMM_F48 - EMM_F66 0.344695 0.452 180   0.762  1.0000
    ##  EMM_F48 - Control 0.344695 0.452 180   0.762  1.0000
    ##  EMM_F48 - EMM_F34 0.340072 0.496 180   0.686  1.0000
    ##  EMM_F48 - EMM_F89 0.335000 0.495 180   0.677  1.0000
    ##  EMM_F48 - EMM_F70 0.298028 0.452 180   0.659  1.0000
    ##  EMM_F48 - EMM_F3  0.274695 0.452 180   0.607  1.0000
    ##  EMM_F48 - EMM_F49 0.259012 0.496 180   0.522  1.0000
    ##  EMM_F48 - ZAN_F4  0.094695 0.452 180   0.209  1.0000
    ##  EMM_F65 - ZAN_F3  0.615072 0.496 180   1.240  0.9973
    ##  EMM_F65 - EMM_F66 0.454695 0.452 180   1.006  0.9998
    ##  EMM_F65 - Control 0.454695 0.452 180   1.006  0.9998
    ##  EMM_F65 - EMM_F34 0.450072 0.496 180   0.908  0.9999
    ##  EMM_F65 - EMM_F89 0.445000 0.495 180   0.899  0.9999
    ##  EMM_F65 - EMM_F70 0.408028 0.452 180   0.902  0.9999
    ##  EMM_F65 - EMM_F3  0.384695 0.452 180   0.851  1.0000
    ##  EMM_F65 - EMM_F49 0.369012 0.496 180   0.744  1.0000
    ##  EMM_F65 - ZAN_F4  0.204695 0.452 180   0.453  1.0000
    ##  EMM_F65 - EMM_F48 0.110000 0.495 180   0.222  1.0000
    ##  SP_F14 - ZAN_F3   0.617044 0.452 180   1.365  0.9925
    ##  SP_F14 - EMM_F66  0.456667 0.404 180   1.130  0.9990
    ##  SP_F14 - Control  0.456667 0.404 180   1.130  0.9990
    ##  SP_F14 - EMM_F34  0.452044 0.452 180   1.000  0.9998
    ##  SP_F14 - EMM_F89  0.446972 0.452 180   0.988  0.9998
    ##  SP_F14 - EMM_F70  0.410000 0.404 180   1.014  0.9997
    ##  SP_F14 - EMM_F3   0.386667 0.404 180   0.957  0.9999
    ##  SP_F14 - EMM_F49  0.370984 0.452 180   0.820  1.0000
    ##  SP_F14 - ZAN_F4   0.206667 0.404 180   0.511  1.0000
    ##  SP_F14 - EMM_F48  0.111972 0.452 180   0.248  1.0000
    ##  SP_F14 - EMM_F65  0.001972 0.452 180   0.004  1.0000
    ##  EMM_F63 - ZAN_F3  0.663711 0.452 180   1.468  0.9847
    ##  EMM_F63 - EMM_F66 0.503333 0.404 180   1.245  0.9971
    ##  EMM_F63 - Control 0.503333 0.404 180   1.245  0.9971
    ##  EMM_F63 - EMM_F34 0.498711 0.452 180   1.103  0.9993
    ##  EMM_F63 - EMM_F89 0.493639 0.452 180   1.092  0.9994
    ##  EMM_F63 - EMM_F70 0.456667 0.404 180   1.130  0.9990
    ##  EMM_F63 - EMM_F3  0.433333 0.404 180   1.072  0.9995
    ##  EMM_F63 - EMM_F49 0.417651 0.452 180   0.924  0.9999
    ##  EMM_F63 - ZAN_F4  0.253333 0.404 180   0.627  1.0000
    ##  EMM_F63 - EMM_F48 0.158639 0.452 180   0.351  1.0000
    ##  EMM_F63 - EMM_F65 0.048639 0.452 180   0.108  1.0000
    ##  EMM_F63 - SP_F14  0.046667 0.404 180   0.115  1.0000
    ##  EMM_F5 - ZAN_F3   0.707044 0.452 180   1.564  0.9727
    ##  EMM_F5 - EMM_F66  0.546667 0.404 180   1.353  0.9932
    ##  EMM_F5 - Control  0.546667 0.404 180   1.353  0.9932
    ##  EMM_F5 - EMM_F34  0.542044 0.452 180   1.199  0.9981
    ##  EMM_F5 - EMM_F89  0.536972 0.452 180   1.187  0.9983
    ##  EMM_F5 - EMM_F70  0.500000 0.404 180   1.237  0.9973
    ##  EMM_F5 - EMM_F3   0.476667 0.404 180   1.179  0.9984
    ##  EMM_F5 - EMM_F49  0.460984 0.452 180   1.019  0.9997
    ##  EMM_F5 - ZAN_F4   0.296667 0.404 180   0.734  1.0000
    ##  EMM_F5 - EMM_F48  0.201972 0.452 180   0.447  1.0000
    ##  EMM_F5 - EMM_F65  0.091972 0.452 180   0.203  1.0000
    ##  EMM_F5 - SP_F14   0.090000 0.404 180   0.223  1.0000
    ##  EMM_F5 - EMM_F63  0.043333 0.404 180   0.107  1.0000
    ##  EMM_F64 - ZAN_F3  0.740377 0.452 180   1.637  0.9592
    ##  EMM_F64 - EMM_F66 0.580000 0.404 180   1.435  0.9877
    ##  EMM_F64 - Control 0.580000 0.404 180   1.435  0.9877
    ##  EMM_F64 - EMM_F34 0.575377 0.452 180   1.272  0.9964
    ##  EMM_F64 - EMM_F89 0.570305 0.452 180   1.261  0.9967
    ##  EMM_F64 - EMM_F70 0.533333 0.404 180   1.320  0.9947
    ##  EMM_F64 - EMM_F3  0.510000 0.404 180   1.262  0.9967
    ##  EMM_F64 - EMM_F49 0.494317 0.452 180   1.093  0.9993
    ##  EMM_F64 - ZAN_F4  0.330000 0.404 180   0.816  1.0000
    ##  EMM_F64 - EMM_F48 0.235305 0.452 180   0.520  1.0000
    ##  EMM_F64 - EMM_F65 0.125305 0.452 180   0.277  1.0000
    ##  EMM_F64 - SP_F14  0.123333 0.404 180   0.305  1.0000
    ##  EMM_F64 - EMM_F63 0.076667 0.404 180   0.190  1.0000
    ##  EMM_F64 - EMM_F5  0.033333 0.404 180   0.082  1.0000
    ##  EMM_F7 - ZAN_F3   0.827044 0.452 180   1.829  0.9026
    ##  EMM_F7 - EMM_F66  0.666667 0.404 180   1.649  0.9566
    ##  EMM_F7 - Control  0.666667 0.404 180   1.649  0.9566
    ##  EMM_F7 - EMM_F34  0.662044 0.452 180   1.464  0.9851
    ##  EMM_F7 - EMM_F89  0.656972 0.452 180   1.453  0.9862
    ##  EMM_F7 - EMM_F70  0.620000 0.404 180   1.534  0.9770
    ##  EMM_F7 - EMM_F3   0.596667 0.404 180   1.476  0.9839
    ##  EMM_F7 - EMM_F49  0.580984 0.452 180   1.285  0.9960
    ##  EMM_F7 - ZAN_F4   0.416667 0.404 180   1.031  0.9997
    ##  EMM_F7 - EMM_F48  0.321972 0.452 180   0.712  1.0000
    ##  EMM_F7 - EMM_F65  0.211972 0.452 180   0.469  1.0000
    ##  EMM_F7 - SP_F14   0.210000 0.404 180   0.520  1.0000
    ##  EMM_F7 - EMM_F63  0.163333 0.404 180   0.404  1.0000
    ##  EMM_F7 - EMM_F5   0.120000 0.404 180   0.297  1.0000
    ##  EMM_F7 - EMM_F64  0.086667 0.404 180   0.214  1.0000
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F49 - EMM_F34 0.001060 0.496 180   0.002  1.0000
    ##  Control - EMM_F34 0.050377 0.452 180   0.111  1.0000
    ##  Control - EMM_F49 0.049317 0.452 180   0.109  1.0000
    ##  EMM_F3 - EMM_F34  0.080377 0.452 180   0.178  1.0000
    ##  EMM_F3 - EMM_F49  0.079317 0.452 180   0.175  1.0000
    ##  EMM_F3 - Control  0.030000 0.404 180   0.074  1.0000
    ##  EMM_F65 - EMM_F34 0.150072 0.496 180   0.303  1.0000
    ##  EMM_F65 - EMM_F49 0.149012 0.496 180   0.301  1.0000
    ##  EMM_F65 - Control 0.099695 0.452 180   0.220  1.0000
    ##  EMM_F65 - EMM_F3  0.069695 0.452 180   0.154  1.0000
    ##  EMM_F70 - EMM_F34 0.207044 0.452 180   0.458  1.0000
    ##  EMM_F70 - EMM_F49 0.205984 0.452 180   0.456  1.0000
    ##  EMM_F70 - Control 0.156667 0.404 180   0.388  1.0000
    ##  EMM_F70 - EMM_F3  0.126667 0.404 180   0.313  1.0000
    ##  EMM_F70 - EMM_F65 0.056972 0.452 180   0.126  1.0000
    ##  EMM_F66 - EMM_F34 0.287044 0.452 180   0.635  1.0000
    ##  EMM_F66 - EMM_F49 0.285984 0.452 180   0.632  1.0000
    ##  EMM_F66 - Control 0.236667 0.404 180   0.586  1.0000
    ##  EMM_F66 - EMM_F3  0.206667 0.404 180   0.511  1.0000
    ##  EMM_F66 - EMM_F65 0.136972 0.452 180   0.303  1.0000
    ##  EMM_F66 - EMM_F70 0.080000 0.404 180   0.198  1.0000
    ##  ZAN_F3 - EMM_F34  0.295000 0.495 180   0.596  1.0000
    ##  ZAN_F3 - EMM_F49  0.293940 0.496 180   0.593  1.0000
    ##  ZAN_F3 - Control  0.244623 0.452 180   0.541  1.0000
    ##  ZAN_F3 - EMM_F3   0.214623 0.452 180   0.475  1.0000
    ##  ZAN_F3 - EMM_F65  0.144928 0.496 180   0.292  1.0000
    ##  ZAN_F3 - EMM_F70  0.087956 0.452 180   0.195  1.0000
    ##  ZAN_F3 - EMM_F66  0.007956 0.452 180   0.018  1.0000
    ##  SP_F14 - EMM_F34  0.313711 0.452 180   0.694  1.0000
    ##  SP_F14 - EMM_F49  0.312651 0.452 180   0.691  1.0000
    ##  SP_F14 - Control  0.263333 0.404 180   0.652  1.0000
    ##  SP_F14 - EMM_F3   0.233333 0.404 180   0.577  1.0000
    ##  SP_F14 - EMM_F65  0.163639 0.452 180   0.362  1.0000
    ##  SP_F14 - EMM_F70  0.106667 0.404 180   0.264  1.0000
    ##  SP_F14 - EMM_F66  0.026667 0.404 180   0.066  1.0000
    ##  SP_F14 - ZAN_F3   0.018711 0.452 180   0.041  1.0000
    ##  EMM_F5 - EMM_F34  0.393711 0.452 180   0.871  1.0000
    ##  EMM_F5 - EMM_F49  0.392651 0.452 180   0.868  1.0000
    ##  EMM_F5 - Control  0.343333 0.404 180   0.849  1.0000
    ##  EMM_F5 - EMM_F3   0.313333 0.404 180   0.775  1.0000
    ##  EMM_F5 - EMM_F65  0.243639 0.452 180   0.539  1.0000
    ##  EMM_F5 - EMM_F70  0.186667 0.404 180   0.462  1.0000
    ##  EMM_F5 - EMM_F66  0.106667 0.404 180   0.264  1.0000
    ##  EMM_F5 - ZAN_F3   0.098711 0.452 180   0.218  1.0000
    ##  EMM_F5 - SP_F14   0.080000 0.404 180   0.198  1.0000
    ##  EMM_F48 - EMM_F34 0.510072 0.496 180   1.029  0.9997
    ##  EMM_F48 - EMM_F49 0.509012 0.496 180   1.027  0.9997
    ##  EMM_F48 - Control 0.459695 0.452 180   1.017  0.9997
    ##  EMM_F48 - EMM_F3  0.429695 0.452 180   0.950  0.9999
    ##  EMM_F48 - EMM_F65 0.360000 0.495 180   0.727  1.0000
    ##  EMM_F48 - EMM_F70 0.303028 0.452 180   0.670  1.0000
    ##  EMM_F48 - EMM_F66 0.223028 0.452 180   0.493  1.0000
    ##  EMM_F48 - ZAN_F3  0.215072 0.496 180   0.434  1.0000
    ##  EMM_F48 - SP_F14  0.196361 0.452 180   0.434  1.0000
    ##  EMM_F48 - EMM_F5  0.116361 0.452 180   0.257  1.0000
    ##  EMM_F7 - EMM_F34  0.577044 0.452 180   1.276  0.9963
    ##  EMM_F7 - EMM_F49  0.575984 0.452 180   1.274  0.9964
    ##  EMM_F7 - Control  0.526667 0.404 180   1.303  0.9954
    ##  EMM_F7 - EMM_F3   0.496667 0.404 180   1.229  0.9975
    ##  EMM_F7 - EMM_F65  0.426972 0.452 180   0.944  0.9999
    ##  EMM_F7 - EMM_F70  0.370000 0.404 180   0.915  0.9999
    ##  EMM_F7 - EMM_F66  0.290000 0.404 180   0.717  1.0000
    ##  EMM_F7 - ZAN_F3   0.282044 0.452 180   0.624  1.0000
    ##  EMM_F7 - SP_F14   0.263333 0.404 180   0.652  1.0000
    ##  EMM_F7 - EMM_F5   0.183333 0.404 180   0.454  1.0000
    ##  EMM_F7 - EMM_F48  0.066972 0.452 180   0.148  1.0000
    ##  EMM_F89 - EMM_F34 0.595072 0.496 180   1.200  0.9981
    ##  EMM_F89 - EMM_F49 0.594012 0.496 180   1.198  0.9981
    ##  EMM_F89 - Control 0.544695 0.452 180   1.205  0.9980
    ##  EMM_F89 - EMM_F3  0.514695 0.452 180   1.138  0.9990
    ##  EMM_F89 - EMM_F65 0.445000 0.495 180   0.899  0.9999
    ##  EMM_F89 - EMM_F70 0.388028 0.452 180   0.858  1.0000
    ##  EMM_F89 - EMM_F66 0.308028 0.452 180   0.681  1.0000
    ##  EMM_F89 - ZAN_F3  0.300072 0.496 180   0.605  1.0000
    ##  EMM_F89 - SP_F14  0.281361 0.452 180   0.622  1.0000
    ##  EMM_F89 - EMM_F5  0.201361 0.452 180   0.445  1.0000
    ##  EMM_F89 - EMM_F48 0.085000 0.495 180   0.172  1.0000
    ##  EMM_F89 - EMM_F7  0.018028 0.452 180   0.040  1.0000
    ##  ZAN_F4 - EMM_F34  0.623711 0.452 180   1.379  0.9917
    ##  ZAN_F4 - EMM_F49  0.622651 0.452 180   1.377  0.9918
    ##  ZAN_F4 - Control  0.573333 0.404 180   1.418  0.9890
    ##  ZAN_F4 - EMM_F3   0.543333 0.404 180   1.344  0.9936
    ##  ZAN_F4 - EMM_F65  0.473639 0.452 180   1.047  0.9996
    ##  ZAN_F4 - EMM_F70  0.416667 0.404 180   1.031  0.9997
    ##  ZAN_F4 - EMM_F66  0.336667 0.404 180   0.833  1.0000
    ##  ZAN_F4 - ZAN_F3   0.328711 0.452 180   0.727  1.0000
    ##  ZAN_F4 - SP_F14   0.310000 0.404 180   0.767  1.0000
    ##  ZAN_F4 - EMM_F5   0.230000 0.404 180   0.569  1.0000
    ##  ZAN_F4 - EMM_F48  0.113639 0.452 180   0.251  1.0000
    ##  ZAN_F4 - EMM_F7   0.046667 0.404 180   0.115  1.0000
    ##  ZAN_F4 - EMM_F89  0.028639 0.452 180   0.063  1.0000
    ##  EMM_F63 - EMM_F34 0.833711 0.452 180   1.844  0.8968
    ##  EMM_F63 - EMM_F49 0.832651 0.452 180   1.841  0.8977
    ##  EMM_F63 - Control 0.783333 0.404 180   1.938  0.8545
    ##  EMM_F63 - EMM_F3  0.753333 0.404 180   1.864  0.8885
    ##  EMM_F63 - EMM_F65 0.683639 0.452 180   1.512  0.9799
    ##  EMM_F63 - EMM_F70 0.626667 0.404 180   1.550  0.9746
    ##  EMM_F63 - EMM_F66 0.546667 0.404 180   1.353  0.9932
    ##  EMM_F63 - ZAN_F3  0.538711 0.452 180   1.191  0.9982
    ##  EMM_F63 - SP_F14  0.520000 0.404 180   1.287  0.9960
    ##  EMM_F63 - EMM_F5  0.440000 0.404 180   1.089  0.9994
    ##  EMM_F63 - EMM_F48 0.323639 0.452 180   0.716  1.0000
    ##  EMM_F63 - EMM_F7  0.256667 0.404 180   0.635  1.0000
    ##  EMM_F63 - EMM_F89 0.238639 0.452 180   0.528  1.0000
    ##  EMM_F63 - ZAN_F4  0.210000 0.404 180   0.520  1.0000
    ##  EMM_F64 - EMM_F34 0.867044 0.452 180   1.917  0.8645
    ##  EMM_F64 - EMM_F49 0.865984 0.452 180   1.915  0.8656
    ##  EMM_F64 - Control 0.816667 0.404 180   2.021  0.8106
    ##  EMM_F64 - EMM_F3  0.786667 0.404 180   1.946  0.8504
    ##  EMM_F64 - EMM_F65 0.716972 0.452 180   1.586  0.9691
    ##  EMM_F64 - EMM_F70 0.660000 0.404 180   1.633  0.9601
    ##  EMM_F64 - EMM_F66 0.580000 0.404 180   1.435  0.9877
    ##  EMM_F64 - ZAN_F3  0.572044 0.452 180   1.265  0.9966
    ##  EMM_F64 - SP_F14  0.553333 0.404 180   1.369  0.9923
    ##  EMM_F64 - EMM_F5  0.473333 0.404 180   1.171  0.9986
    ##  EMM_F64 - EMM_F48 0.356972 0.452 180   0.789  1.0000
    ##  EMM_F64 - EMM_F7  0.290000 0.404 180   0.717  1.0000
    ##  EMM_F64 - EMM_F89 0.271972 0.452 180   0.601  1.0000
    ##  EMM_F64 - ZAN_F4  0.243333 0.404 180   0.602  1.0000
    ##  EMM_F64 - EMM_F63 0.033333 0.404 180   0.082  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49 0.224317 0.452 180   0.496  1.0000
    ##  EMM_F66 - EMM_F49 0.360984 0.452 180   0.798  1.0000
    ##  EMM_F66 - Control 0.136667 0.404 180   0.338  1.0000
    ##  EMM_F3 - EMM_F49  0.444317 0.452 180   0.983  0.9998
    ##  EMM_F3 - Control  0.220000 0.404 180   0.544  1.0000
    ##  EMM_F3 - EMM_F66  0.083333 0.404 180   0.206  1.0000
    ##  EMM_F70 - EMM_F49 0.524317 0.452 180   1.160  0.9987
    ##  EMM_F70 - Control 0.300000 0.404 180   0.742  1.0000
    ##  EMM_F70 - EMM_F66 0.163333 0.404 180   0.404  1.0000
    ##  EMM_F70 - EMM_F3  0.080000 0.404 180   0.198  1.0000
    ##  EMM_F34 - EMM_F49 0.553940 0.496 180   1.117  0.9992
    ##  EMM_F34 - Control 0.329623 0.452 180   0.729  1.0000
    ##  EMM_F34 - EMM_F66 0.192956 0.452 180   0.427  1.0000
    ##  EMM_F34 - EMM_F3  0.109623 0.452 180   0.242  1.0000
    ##  EMM_F34 - EMM_F70 0.029623 0.452 180   0.066  1.0000
    ##  EMM_F65 - EMM_F49 0.579012 0.496 180   1.168  0.9986
    ##  EMM_F65 - Control 0.354695 0.452 180   0.784  1.0000
    ##  EMM_F65 - EMM_F66 0.218028 0.452 180   0.482  1.0000
    ##  EMM_F65 - EMM_F3  0.134695 0.452 180   0.298  1.0000
    ##  EMM_F65 - EMM_F70 0.054695 0.452 180   0.121  1.0000
    ##  EMM_F65 - EMM_F34 0.025072 0.496 180   0.051  1.0000
    ##  EMM_F5 - EMM_F49  0.610984 0.452 180   1.351  0.9933
    ##  EMM_F5 - Control  0.386667 0.404 180   0.957  0.9999
    ##  EMM_F5 - EMM_F66  0.250000 0.404 180   0.619  1.0000
    ##  EMM_F5 - EMM_F3   0.166667 0.404 180   0.412  1.0000
    ##  EMM_F5 - EMM_F70  0.086667 0.404 180   0.214  1.0000
    ##  EMM_F5 - EMM_F34  0.057044 0.452 180   0.126  1.0000
    ##  EMM_F5 - EMM_F65  0.031972 0.452 180   0.071  1.0000
    ##  ZAN_F3 - EMM_F49  0.658940 0.496 180   1.329  0.9943
    ##  ZAN_F3 - Control  0.434623 0.452 180   0.961  0.9999
    ##  ZAN_F3 - EMM_F66  0.297956 0.452 180   0.659  1.0000
    ##  ZAN_F3 - EMM_F3   0.214623 0.452 180   0.475  1.0000
    ##  ZAN_F3 - EMM_F70  0.134623 0.452 180   0.298  1.0000
    ##  ZAN_F3 - EMM_F34  0.105000 0.495 180   0.212  1.0000
    ##  ZAN_F3 - EMM_F65  0.079928 0.496 180   0.161  1.0000
    ##  ZAN_F3 - EMM_F5   0.047956 0.452 180   0.106  1.0000
    ##  EMM_F89 - EMM_F49 0.679012 0.496 180   1.369  0.9923
    ##  EMM_F89 - Control 0.454695 0.452 180   1.006  0.9998
    ##  EMM_F89 - EMM_F66 0.318028 0.452 180   0.703  1.0000
    ##  EMM_F89 - EMM_F3  0.234695 0.452 180   0.519  1.0000
    ##  EMM_F89 - EMM_F70 0.154695 0.452 180   0.342  1.0000
    ##  EMM_F89 - EMM_F34 0.125072 0.496 180   0.252  1.0000
    ##  EMM_F89 - EMM_F65 0.100000 0.495 180   0.202  1.0000
    ##  EMM_F89 - EMM_F5  0.068028 0.452 180   0.150  1.0000
    ##  EMM_F89 - ZAN_F3  0.020072 0.496 180   0.040  1.0000
    ##  EMM_F48 - EMM_F49 0.734012 0.496 180   1.480  0.9834
    ##  EMM_F48 - Control 0.509695 0.452 180   1.127  0.9991
    ##  EMM_F48 - EMM_F66 0.373028 0.452 180   0.825  1.0000
    ##  EMM_F48 - EMM_F3  0.289695 0.452 180   0.641  1.0000
    ##  EMM_F48 - EMM_F70 0.209695 0.452 180   0.464  1.0000
    ##  EMM_F48 - EMM_F34 0.180072 0.496 180   0.363  1.0000
    ##  EMM_F48 - EMM_F65 0.155000 0.495 180   0.313  1.0000
    ##  EMM_F48 - EMM_F5  0.123028 0.452 180   0.272  1.0000
    ##  EMM_F48 - ZAN_F3  0.075072 0.496 180   0.151  1.0000
    ##  EMM_F48 - EMM_F89 0.055000 0.495 180   0.111  1.0000
    ##  SP_F14 - EMM_F49  0.750984 0.452 180   1.661  0.9540
    ##  SP_F14 - Control  0.526667 0.404 180   1.303  0.9954
    ##  SP_F14 - EMM_F66  0.390000 0.404 180   0.965  0.9999
    ##  SP_F14 - EMM_F3   0.306667 0.404 180   0.759  1.0000
    ##  SP_F14 - EMM_F70  0.226667 0.404 180   0.561  1.0000
    ##  SP_F14 - EMM_F34  0.197044 0.452 180   0.436  1.0000
    ##  SP_F14 - EMM_F65  0.171972 0.452 180   0.380  1.0000
    ##  SP_F14 - EMM_F5   0.140000 0.404 180   0.346  1.0000
    ##  SP_F14 - ZAN_F3   0.092044 0.452 180   0.204  1.0000
    ##  SP_F14 - EMM_F89  0.071972 0.452 180   0.159  1.0000
    ##  SP_F14 - EMM_F48  0.016972 0.452 180   0.038  1.0000
    ##  ZAN_F4 - EMM_F49  0.780984 0.452 180   1.727  0.9369
    ##  ZAN_F4 - Control  0.556667 0.404 180   1.377  0.9918
    ##  ZAN_F4 - EMM_F66  0.420000 0.404 180   1.039  0.9996
    ##  ZAN_F4 - EMM_F3   0.336667 0.404 180   0.833  1.0000
    ##  ZAN_F4 - EMM_F70  0.256667 0.404 180   0.635  1.0000
    ##  ZAN_F4 - EMM_F34  0.227044 0.452 180   0.502  1.0000
    ##  ZAN_F4 - EMM_F65  0.201972 0.452 180   0.447  1.0000
    ##  ZAN_F4 - EMM_F5   0.170000 0.404 180   0.421  1.0000
    ##  ZAN_F4 - ZAN_F3   0.122044 0.452 180   0.270  1.0000
    ##  ZAN_F4 - EMM_F89  0.101972 0.452 180   0.226  1.0000
    ##  ZAN_F4 - EMM_F48  0.046972 0.452 180   0.104  1.0000
    ##  ZAN_F4 - SP_F14   0.030000 0.404 180   0.074  1.0000
    ##  EMM_F7 - EMM_F49  0.850984 0.452 180   1.882  0.8807
    ##  EMM_F7 - Control  0.626667 0.404 180   1.550  0.9746
    ##  EMM_F7 - EMM_F66  0.490000 0.404 180   1.212  0.9979
    ##  EMM_F7 - EMM_F3   0.406667 0.404 180   1.006  0.9998
    ##  EMM_F7 - EMM_F70  0.326667 0.404 180   0.808  1.0000
    ##  EMM_F7 - EMM_F34  0.297044 0.452 180   0.657  1.0000
    ##  EMM_F7 - EMM_F65  0.271972 0.452 180   0.601  1.0000
    ##  EMM_F7 - EMM_F5   0.240000 0.404 180   0.594  1.0000
    ##  EMM_F7 - ZAN_F3   0.192044 0.452 180   0.425  1.0000
    ##  EMM_F7 - EMM_F89  0.171972 0.452 180   0.380  1.0000
    ##  EMM_F7 - EMM_F48  0.116972 0.452 180   0.259  1.0000
    ##  EMM_F7 - SP_F14   0.100000 0.404 180   0.247  1.0000
    ##  EMM_F7 - ZAN_F4   0.070000 0.404 180   0.173  1.0000
    ##  EMM_F64 - EMM_F49 0.914317 0.452 180   2.022  0.8097
    ##  EMM_F64 - Control 0.690000 0.404 180   1.707  0.9424
    ##  EMM_F64 - EMM_F66 0.553333 0.404 180   1.369  0.9923
    ##  EMM_F64 - EMM_F3  0.470000 0.404 180   1.163  0.9987
    ##  EMM_F64 - EMM_F70 0.390000 0.404 180   0.965  0.9999
    ##  EMM_F64 - EMM_F34 0.360377 0.452 180   0.797  1.0000
    ##  EMM_F64 - EMM_F65 0.335305 0.452 180   0.742  1.0000
    ##  EMM_F64 - EMM_F5  0.303333 0.404 180   0.750  1.0000
    ##  EMM_F64 - ZAN_F3  0.255377 0.452 180   0.565  1.0000
    ##  EMM_F64 - EMM_F89 0.235305 0.452 180   0.520  1.0000
    ##  EMM_F64 - EMM_F48 0.180305 0.452 180   0.399  1.0000
    ##  EMM_F64 - SP_F14  0.163333 0.404 180   0.404  1.0000
    ##  EMM_F64 - ZAN_F4  0.133333 0.404 180   0.330  1.0000
    ##  EMM_F64 - EMM_F7  0.063333 0.404 180   0.157  1.0000
    ##  EMM_F63 - EMM_F49 0.977651 0.452 180   2.162  0.7219
    ##  EMM_F63 - Control 0.753333 0.404 180   1.864  0.8885
    ##  EMM_F63 - EMM_F66 0.616667 0.404 180   1.526  0.9781
    ##  EMM_F63 - EMM_F3  0.533333 0.404 180   1.320  0.9947
    ##  EMM_F63 - EMM_F70 0.453333 0.404 180   1.122  0.9991
    ##  EMM_F63 - EMM_F34 0.423711 0.452 180   0.937  0.9999
    ##  EMM_F63 - EMM_F65 0.398639 0.452 180   0.882  1.0000
    ##  EMM_F63 - EMM_F5  0.366667 0.404 180   0.907  0.9999
    ##  EMM_F63 - ZAN_F3  0.318711 0.452 180   0.705  1.0000
    ##  EMM_F63 - EMM_F89 0.298639 0.452 180   0.660  1.0000
    ##  EMM_F63 - EMM_F48 0.243639 0.452 180   0.539  1.0000
    ##  EMM_F63 - SP_F14  0.226667 0.404 180   0.561  1.0000
    ##  EMM_F63 - ZAN_F4  0.196667 0.404 180   0.487  1.0000
    ##  EMM_F63 - EMM_F7  0.126667 0.404 180   0.313  1.0000
    ##  EMM_F63 - EMM_F64 0.063333 0.404 180   0.157  1.0000
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F49 - Control 0.134016 0.452 180   0.296  1.0000
    ##  EMM_F3 - Control  0.410000 0.404 180   1.014  0.9997
    ##  EMM_F3 - EMM_F49  0.275984 0.452 180   0.610  1.0000
    ##  EMM_F7 - Control  0.480000 0.404 180   1.188  0.9983
    ##  EMM_F7 - EMM_F49  0.345984 0.452 180   0.765  1.0000
    ##  EMM_F7 - EMM_F3   0.070000 0.404 180   0.173  1.0000
    ##  EMM_F48 - Control 0.558028 0.452 180   1.234  0.9974
    ##  EMM_F48 - EMM_F49 0.424012 0.496 180   0.855  1.0000
    ##  EMM_F48 - EMM_F3  0.148028 0.452 180   0.327  1.0000
    ##  EMM_F48 - EMM_F7  0.078028 0.452 180   0.173  1.0000
    ##  EMM_F34 - Control 0.617956 0.452 180   1.367  0.9924
    ##  EMM_F34 - EMM_F49 0.483940 0.496 180   0.976  0.9998
    ##  EMM_F34 - EMM_F3  0.207956 0.452 180   0.460  1.0000
    ##  EMM_F34 - EMM_F7  0.137956 0.452 180   0.305  1.0000
    ##  EMM_F34 - EMM_F48 0.059928 0.496 180   0.121  1.0000
    ##  EMM_F89 - Control 0.633028 0.452 180   1.400  0.9904
    ##  EMM_F89 - EMM_F49 0.499012 0.496 180   1.006  0.9998
    ##  EMM_F89 - EMM_F3  0.223028 0.452 180   0.493  1.0000
    ##  EMM_F89 - EMM_F7  0.153028 0.452 180   0.338  1.0000
    ##  EMM_F89 - EMM_F48 0.075000 0.495 180   0.152  1.0000
    ##  EMM_F89 - EMM_F34 0.015072 0.496 180   0.030  1.0000
    ##  EMM_F66 - Control 0.633333 0.404 180   1.567  0.9721
    ##  EMM_F66 - EMM_F49 0.499317 0.452 180   1.104  0.9993
    ##  EMM_F66 - EMM_F3  0.223333 0.404 180   0.553  1.0000
    ##  EMM_F66 - EMM_F7  0.153333 0.404 180   0.379  1.0000
    ##  EMM_F66 - EMM_F48 0.075305 0.452 180   0.167  1.0000
    ##  EMM_F66 - EMM_F34 0.015377 0.452 180   0.034  1.0000
    ##  EMM_F66 - EMM_F89 0.000305 0.452 180   0.001  1.0000
    ##  EMM_F5 - Control  0.653333 0.404 180   1.616  0.9634
    ##  EMM_F5 - EMM_F49  0.519317 0.452 180   1.148  0.9988
    ##  EMM_F5 - EMM_F3   0.243333 0.404 180   0.602  1.0000
    ##  EMM_F5 - EMM_F7   0.173333 0.404 180   0.429  1.0000
    ##  EMM_F5 - EMM_F48  0.095305 0.452 180   0.211  1.0000
    ##  EMM_F5 - EMM_F34  0.035377 0.452 180   0.078  1.0000
    ##  EMM_F5 - EMM_F89  0.020305 0.452 180   0.045  1.0000
    ##  EMM_F5 - EMM_F66  0.020000 0.404 180   0.049  1.0000
    ##  ZAN_F4 - Control  0.656667 0.404 180   1.625  0.9618
    ##  ZAN_F4 - EMM_F49  0.522651 0.452 180   1.156  0.9988
    ##  ZAN_F4 - EMM_F3   0.246667 0.404 180   0.610  1.0000
    ##  ZAN_F4 - EMM_F7   0.176667 0.404 180   0.437  1.0000
    ##  ZAN_F4 - EMM_F48  0.098639 0.452 180   0.218  1.0000
    ##  ZAN_F4 - EMM_F34  0.038711 0.452 180   0.086  1.0000
    ##  ZAN_F4 - EMM_F89  0.023639 0.452 180   0.052  1.0000
    ##  ZAN_F4 - EMM_F66  0.023333 0.404 180   0.058  1.0000
    ##  ZAN_F4 - EMM_F5   0.003333 0.404 180   0.008  1.0000
    ##  SP_F14 - Control  0.660000 0.404 180   1.633  0.9601
    ##  SP_F14 - EMM_F49  0.525984 0.452 180   1.163  0.9987
    ##  SP_F14 - EMM_F3   0.250000 0.404 180   0.619  1.0000
    ##  SP_F14 - EMM_F7   0.180000 0.404 180   0.445  1.0000
    ##  SP_F14 - EMM_F48  0.101972 0.452 180   0.226  1.0000
    ##  SP_F14 - EMM_F34  0.042044 0.452 180   0.093  1.0000
    ##  SP_F14 - EMM_F89  0.026972 0.452 180   0.060  1.0000
    ##  SP_F14 - EMM_F66  0.026667 0.404 180   0.066  1.0000
    ##  SP_F14 - EMM_F5   0.006667 0.404 180   0.016  1.0000
    ##  SP_F14 - ZAN_F4   0.003333 0.404 180   0.008  1.0000
    ##  EMM_F70 - Control 0.703333 0.404 180   1.740  0.9330
    ##  EMM_F70 - EMM_F49 0.569317 0.452 180   1.259  0.9968
    ##  EMM_F70 - EMM_F3  0.293333 0.404 180   0.726  1.0000
    ##  EMM_F70 - EMM_F7  0.223333 0.404 180   0.553  1.0000
    ##  EMM_F70 - EMM_F48 0.145305 0.452 180   0.321  1.0000
    ##  EMM_F70 - EMM_F34 0.085377 0.452 180   0.189  1.0000
    ##  EMM_F70 - EMM_F89 0.070305 0.452 180   0.155  1.0000
    ##  EMM_F70 - EMM_F66 0.070000 0.404 180   0.173  1.0000
    ##  EMM_F70 - EMM_F5  0.050000 0.404 180   0.124  1.0000
    ##  EMM_F70 - ZAN_F4  0.046667 0.404 180   0.115  1.0000
    ##  EMM_F70 - SP_F14  0.043333 0.404 180   0.107  1.0000
    ##  EMM_F63 - Control 0.710000 0.404 180   1.757  0.9280
    ##  EMM_F63 - EMM_F49 0.575984 0.452 180   1.274  0.9964
    ##  EMM_F63 - EMM_F3  0.300000 0.404 180   0.742  1.0000
    ##  EMM_F63 - EMM_F7  0.230000 0.404 180   0.569  1.0000
    ##  EMM_F63 - EMM_F48 0.151972 0.452 180   0.336  1.0000
    ##  EMM_F63 - EMM_F34 0.092044 0.452 180   0.204  1.0000
    ##  EMM_F63 - EMM_F89 0.076972 0.452 180   0.170  1.0000
    ##  EMM_F63 - EMM_F66 0.076667 0.404 180   0.190  1.0000
    ##  EMM_F63 - EMM_F5  0.056667 0.404 180   0.140  1.0000
    ##  EMM_F63 - ZAN_F4  0.053333 0.404 180   0.132  1.0000
    ##  EMM_F63 - SP_F14  0.050000 0.404 180   0.124  1.0000
    ##  EMM_F63 - EMM_F70 0.006667 0.404 180   0.016  1.0000
    ##  EMM_F64 - Control 0.756667 0.404 180   1.872  0.8850
    ##  EMM_F64 - EMM_F49 0.622651 0.452 180   1.377  0.9918
    ##  EMM_F64 - EMM_F3  0.346667 0.404 180   0.858  1.0000
    ##  EMM_F64 - EMM_F7  0.276667 0.404 180   0.685  1.0000
    ##  EMM_F64 - EMM_F48 0.198639 0.452 180   0.439  1.0000
    ##  EMM_F64 - EMM_F34 0.138711 0.452 180   0.307  1.0000
    ##  EMM_F64 - EMM_F89 0.123639 0.452 180   0.273  1.0000
    ##  EMM_F64 - EMM_F66 0.123333 0.404 180   0.305  1.0000
    ##  EMM_F64 - EMM_F5  0.103333 0.404 180   0.256  1.0000
    ##  EMM_F64 - ZAN_F4  0.100000 0.404 180   0.247  1.0000
    ##  EMM_F64 - SP_F14  0.096667 0.404 180   0.239  1.0000
    ##  EMM_F64 - EMM_F70 0.053333 0.404 180   0.132  1.0000
    ##  EMM_F64 - EMM_F63 0.046667 0.404 180   0.115  1.0000
    ##  ZAN_F3 - Control  1.067956 0.452 180   2.362  0.5786
    ##  ZAN_F3 - EMM_F49  0.933940 0.496 180   1.884  0.8799
    ##  ZAN_F3 - EMM_F3   0.657956 0.452 180   1.455  0.9859
    ##  ZAN_F3 - EMM_F7   0.587956 0.452 180   1.300  0.9955
    ##  ZAN_F3 - EMM_F48  0.509928 0.496 180   1.028  0.9997
    ##  ZAN_F3 - EMM_F34  0.450000 0.495 180   0.909  0.9999
    ##  ZAN_F3 - EMM_F89  0.434928 0.496 180   0.877  1.0000
    ##  ZAN_F3 - EMM_F66  0.434623 0.452 180   0.961  0.9999
    ##  ZAN_F3 - EMM_F5   0.414623 0.452 180   0.917  0.9999
    ##  ZAN_F3 - ZAN_F4   0.411289 0.452 180   0.910  0.9999
    ##  ZAN_F3 - SP_F14   0.407956 0.452 180   0.902  0.9999
    ##  ZAN_F3 - EMM_F70  0.364623 0.452 180   0.806  1.0000
    ##  ZAN_F3 - EMM_F63  0.357956 0.452 180   0.792  1.0000
    ##  ZAN_F3 - EMM_F64  0.311289 0.452 180   0.688  1.0000
    ##  EMM_F65 - Control 1.073028 0.452 180   2.373  0.5703
    ##  EMM_F65 - EMM_F49 0.939012 0.496 180   1.894  0.8754
    ##  EMM_F65 - EMM_F3  0.663028 0.452 180   1.466  0.9849
    ##  EMM_F65 - EMM_F7  0.593028 0.452 180   1.311  0.9950
    ##  EMM_F65 - EMM_F48 0.515000 0.495 180   1.040  0.9996
    ##  EMM_F65 - EMM_F34 0.455072 0.496 180   0.918  0.9999
    ##  EMM_F65 - EMM_F89 0.440000 0.495 180   0.889  0.9999
    ##  EMM_F65 - EMM_F66 0.439695 0.452 180   0.972  0.9998
    ##  EMM_F65 - EMM_F5  0.419695 0.452 180   0.928  0.9999
    ##  EMM_F65 - ZAN_F4  0.416361 0.452 180   0.921  0.9999
    ##  EMM_F65 - SP_F14  0.413028 0.452 180   0.913  0.9999
    ##  EMM_F65 - EMM_F70 0.369695 0.452 180   0.818  1.0000
    ##  EMM_F65 - EMM_F63 0.363028 0.452 180   0.803  1.0000
    ##  EMM_F65 - EMM_F64 0.316361 0.452 180   0.700  1.0000
    ##  EMM_F65 - ZAN_F3  0.005072 0.496 180   0.010  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F49 - ZAN_F3  0.391060 0.496 180   0.789  1.0000
    ##  Control - ZAN_F3  0.397044 0.452 180   0.878  1.0000
    ##  Control - EMM_F49 0.005984 0.452 180   0.013  1.0000
    ##  EMM_F89 - ZAN_F3  0.750072 0.496 180   1.513  0.9797
    ##  EMM_F89 - EMM_F49 0.359012 0.496 180   0.724  1.0000
    ##  EMM_F89 - Control 0.353028 0.452 180   0.781  1.0000
    ##  EMM_F34 - ZAN_F3  0.755000 0.495 180   1.525  0.9782
    ##  EMM_F34 - EMM_F49 0.363940 0.496 180   0.734  1.0000
    ##  EMM_F34 - Control 0.357956 0.452 180   0.792  1.0000
    ##  EMM_F34 - EMM_F89 0.004928 0.496 180   0.010  1.0000
    ##  EMM_F70 - ZAN_F3  0.853711 0.452 180   1.888  0.8780
    ##  EMM_F70 - EMM_F49 0.462651 0.452 180   1.023  0.9997
    ##  EMM_F70 - Control 0.456667 0.404 180   1.130  0.9990
    ##  EMM_F70 - EMM_F89 0.103639 0.452 180   0.229  1.0000
    ##  EMM_F70 - EMM_F34 0.098711 0.452 180   0.218  1.0000
    ##  EMM_F3 - ZAN_F3   0.873711 0.452 180   1.932  0.8574
    ##  EMM_F3 - EMM_F49  0.482651 0.452 180   1.067  0.9995
    ##  EMM_F3 - Control  0.476667 0.404 180   1.179  0.9984
    ##  EMM_F3 - EMM_F89  0.123639 0.452 180   0.273  1.0000
    ##  EMM_F3 - EMM_F34  0.118711 0.452 180   0.263  1.0000
    ##  EMM_F3 - EMM_F70  0.020000 0.404 180   0.049  1.0000
    ##  EMM_F63 - ZAN_F3  0.990377 0.452 180   2.190  0.7027
    ##  EMM_F63 - EMM_F49 0.599317 0.452 180   1.325  0.9945
    ##  EMM_F63 - Control 0.593333 0.404 180   1.468  0.9847
    ##  EMM_F63 - EMM_F89 0.240305 0.452 180   0.531  1.0000
    ##  EMM_F63 - EMM_F34 0.235377 0.452 180   0.521  1.0000
    ##  EMM_F63 - EMM_F70 0.136667 0.404 180   0.338  1.0000
    ##  EMM_F63 - EMM_F3  0.116667 0.404 180   0.289  1.0000
    ##  EMM_F7 - ZAN_F3   1.010377 0.452 180   2.234  0.6717
    ##  EMM_F7 - EMM_F49  0.619317 0.452 180   1.370  0.9923
    ##  EMM_F7 - Control  0.613333 0.404 180   1.517  0.9792
    ##  EMM_F7 - EMM_F89  0.260305 0.452 180   0.576  1.0000
    ##  EMM_F7 - EMM_F34  0.255377 0.452 180   0.565  1.0000
    ##  EMM_F7 - EMM_F70  0.156667 0.404 180   0.388  1.0000
    ##  EMM_F7 - EMM_F3   0.136667 0.404 180   0.338  1.0000
    ##  EMM_F7 - EMM_F63  0.020000 0.404 180   0.049  1.0000
    ##  EMM_F65 - ZAN_F3  1.055072 0.496 180   2.128  0.7448
    ##  EMM_F65 - EMM_F49 0.664012 0.496 180   1.339  0.9938
    ##  EMM_F65 - Control 0.658028 0.452 180   1.455  0.9859
    ##  EMM_F65 - EMM_F89 0.305000 0.495 180   0.616  1.0000
    ##  EMM_F65 - EMM_F34 0.300072 0.496 180   0.605  1.0000
    ##  EMM_F65 - EMM_F70 0.201361 0.452 180   0.445  1.0000
    ##  EMM_F65 - EMM_F3  0.181361 0.452 180   0.401  1.0000
    ##  EMM_F65 - EMM_F63 0.064695 0.452 180   0.143  1.0000
    ##  EMM_F65 - EMM_F7  0.044695 0.452 180   0.099  1.0000
    ##  EMM_F66 - ZAN_F3  1.097044 0.452 180   2.426  0.5308
    ##  EMM_F66 - EMM_F49 0.705984 0.452 180   1.561  0.9730
    ##  EMM_F66 - Control 0.700000 0.404 180   1.732  0.9355
    ##  EMM_F66 - EMM_F89 0.346972 0.452 180   0.767  1.0000
    ##  EMM_F66 - EMM_F34 0.342044 0.452 180   0.756  1.0000
    ##  EMM_F66 - EMM_F70 0.243333 0.404 180   0.602  1.0000
    ##  EMM_F66 - EMM_F3  0.223333 0.404 180   0.553  1.0000
    ##  EMM_F66 - EMM_F63 0.106667 0.404 180   0.264  1.0000
    ##  EMM_F66 - EMM_F7  0.086667 0.404 180   0.214  1.0000
    ##  EMM_F66 - EMM_F65 0.041972 0.452 180   0.093  1.0000
    ##  EMM_F5 - ZAN_F3   1.187044 0.452 180   2.625  0.3882
    ##  EMM_F5 - EMM_F49  0.795984 0.452 180   1.760  0.9268
    ##  EMM_F5 - Control  0.790000 0.404 180   1.955  0.8462
    ##  EMM_F5 - EMM_F89  0.436972 0.452 180   0.966  0.9999
    ##  EMM_F5 - EMM_F34  0.432044 0.452 180   0.955  0.9999
    ##  EMM_F5 - EMM_F70  0.333333 0.404 180   0.825  1.0000
    ##  EMM_F5 - EMM_F3   0.313333 0.404 180   0.775  1.0000
    ##  EMM_F5 - EMM_F63  0.196667 0.404 180   0.487  1.0000
    ##  EMM_F5 - EMM_F7   0.176667 0.404 180   0.437  1.0000
    ##  EMM_F5 - EMM_F65  0.131972 0.452 180   0.292  1.0000
    ##  EMM_F5 - EMM_F66  0.090000 0.404 180   0.223  1.0000
    ##  ZAN_F4 - ZAN_F3   1.247044 0.452 180   2.758  0.3034
    ##  ZAN_F4 - EMM_F49  0.855984 0.452 180   1.893  0.8758
    ##  ZAN_F4 - Control  0.850000 0.404 180   2.103  0.7608
    ##  ZAN_F4 - EMM_F89  0.496972 0.452 180   1.099  0.9993
    ##  ZAN_F4 - EMM_F34  0.492044 0.452 180   1.088  0.9994
    ##  ZAN_F4 - EMM_F70  0.393333 0.404 180   0.973  0.9998
    ##  ZAN_F4 - EMM_F3   0.373333 0.404 180   0.924  0.9999
    ##  ZAN_F4 - EMM_F63  0.256667 0.404 180   0.635  1.0000
    ##  ZAN_F4 - EMM_F7   0.236667 0.404 180   0.586  1.0000
    ##  ZAN_F4 - EMM_F65  0.191972 0.452 180   0.425  1.0000
    ##  ZAN_F4 - EMM_F66  0.150000 0.404 180   0.371  1.0000
    ##  ZAN_F4 - EMM_F5   0.060000 0.404 180   0.148  1.0000
    ##  SP_F14 - ZAN_F3   1.317044 0.452 180   2.913  0.2194
    ##  SP_F14 - EMM_F49  0.925984 0.452 180   2.048  0.7947
    ##  SP_F14 - Control  0.920000 0.404 180   2.276  0.6416
    ##  SP_F14 - EMM_F89  0.566972 0.452 180   1.254  0.9969
    ##  SP_F14 - EMM_F34  0.562044 0.452 180   1.243  0.9972
    ##  SP_F14 - EMM_F70  0.463333 0.404 180   1.146  0.9989
    ##  SP_F14 - EMM_F3   0.443333 0.404 180   1.097  0.9993
    ##  SP_F14 - EMM_F63  0.326667 0.404 180   0.808  1.0000
    ##  SP_F14 - EMM_F7   0.306667 0.404 180   0.759  1.0000
    ##  SP_F14 - EMM_F65  0.261972 0.452 180   0.579  1.0000
    ##  SP_F14 - EMM_F66  0.220000 0.404 180   0.544  1.0000
    ##  SP_F14 - EMM_F5   0.130000 0.404 180   0.322  1.0000
    ##  SP_F14 - ZAN_F4   0.070000 0.404 180   0.173  1.0000
    ##  EMM_F48 - ZAN_F3  1.545072 0.496 180   3.116  0.1354
    ##  EMM_F48 - EMM_F49 1.154012 0.496 180   2.327  0.6040
    ##  EMM_F48 - Control 1.148028 0.452 180   2.539  0.4484
    ##  EMM_F48 - EMM_F89 0.795000 0.495 180   1.606  0.9654
    ##  EMM_F48 - EMM_F34 0.790072 0.496 180   1.593  0.9677
    ##  EMM_F48 - EMM_F70 0.691361 0.452 180   1.529  0.9777
    ##  EMM_F48 - EMM_F3  0.671361 0.452 180   1.485  0.9830
    ##  EMM_F48 - EMM_F63 0.554695 0.452 180   1.227  0.9976
    ##  EMM_F48 - EMM_F7  0.534695 0.452 180   1.182  0.9984
    ##  EMM_F48 - EMM_F65 0.490000 0.495 180   0.990  0.9998
    ##  EMM_F48 - EMM_F66 0.448028 0.452 180   0.991  0.9998
    ##  EMM_F48 - EMM_F5  0.358028 0.452 180   0.792  1.0000
    ##  EMM_F48 - ZAN_F4  0.298028 0.452 180   0.659  1.0000
    ##  EMM_F48 - SP_F14  0.228028 0.452 180   0.504  1.0000
    ##  EMM_F64 - ZAN_F3  1.633711 0.452 180   3.613  0.0327
    ##  EMM_F64 - EMM_F49 1.242651 0.452 180   2.748  0.3092
    ##  EMM_F64 - Control 1.236667 0.404 180   3.060  0.1558
    ##  EMM_F64 - EMM_F89 0.883639 0.452 180   1.954  0.8465
    ##  EMM_F64 - EMM_F34 0.878711 0.452 180   1.943  0.8519
    ##  EMM_F64 - EMM_F70 0.780000 0.404 180   1.930  0.8586
    ##  EMM_F64 - EMM_F3  0.760000 0.404 180   1.880  0.8814
    ##  EMM_F64 - EMM_F63 0.643333 0.404 180   1.592  0.9680
    ##  EMM_F64 - EMM_F7  0.623333 0.404 180   1.542  0.9758
    ##  EMM_F64 - EMM_F65 0.578639 0.452 180   1.280  0.9962
    ##  EMM_F64 - EMM_F66 0.536667 0.404 180   1.328  0.9944
    ##  EMM_F64 - EMM_F5  0.446667 0.404 180   1.105  0.9993
    ##  EMM_F64 - ZAN_F4  0.386667 0.404 180   0.957  0.9999
    ##  EMM_F64 - SP_F14  0.316667 0.404 180   0.783  1.0000
    ##  EMM_F64 - EMM_F48 0.088639 0.452 180   0.196  1.0000
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.54667 0.404 180   1.353  0.9932
    ##  Control - EMM_F65  0.41197 0.452 180   0.911  0.9999
    ##  Control - ZAN_F4   0.23667 0.404 180   0.586  1.0000
    ##  Control - EMM_F64  0.13667 0.404 180   0.338  1.0000
    ##  Control - EMM_F70  0.13000 0.404 180   0.322  1.0000
    ##  Control - EMM_F89  0.08697 0.452 180   0.192  1.0000
    ##  Control - EMM_F34  0.03204 0.452 180   0.071  1.0000
    ##  EMM_F63 - Control  0.00333 0.404 180   0.008  1.0000
    ##  EMM_F66 - Control  0.10000 0.404 180   0.247  1.0000
    ##  EMM_F49 - Control  0.10902 0.452 180   0.241  1.0000
    ##  EMM_F3 - Control   0.16667 0.404 180   0.412  1.0000
    ##  EMM_F7 - Control   0.29333 0.404 180   0.726  1.0000
    ##  EMM_F48 - Control  0.51303 0.452 180   1.135  0.9990
    ##  EMM_F5 - Control   0.61667 0.404 180   1.526  0.9781
    ##  ZAN_F3 - Control   0.77796 0.452 180   1.720  0.9388
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F3   0.06538 0.452 180   0.145  1.0000
    ##  EMM_F34 - Control  0.06462 0.452 180   0.143  1.0000
    ##  EMM_F3 - Control   0.07333 0.404 180   0.181  1.0000
    ##  EMM_F66 - Control  0.11667 0.404 180   0.289  1.0000
    ##  EMM_F70 - Control  0.13000 0.404 180   0.322  1.0000
    ##  ZAN_F4 - Control   0.15000 0.404 180   0.371  1.0000
    ##  EMM_F89 - Control  0.15969 0.452 180   0.353  1.0000
    ##  EMM_F65 - Control  0.26469 0.452 180   0.585  1.0000
    ##  EMM_F64 - Control  0.28333 0.404 180   0.701  1.0000
    ##  EMM_F7 - Control   0.29000 0.404 180   0.717  1.0000
    ##  EMM_F63 - Control  0.30667 0.404 180   0.759  1.0000
    ##  EMM_F49 - Control  0.40568 0.452 180   0.897  0.9999
    ##  SP_F14 - Control   0.41667 0.404 180   1.031  0.9997
    ##  EMM_F5 - Control   0.59667 0.404 180   1.476  0.9839
    ##  EMM_F48 - Control  0.86469 0.452 180   1.912  0.8670
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F3   0.16038 0.452 180   0.355  1.0000
    ##  Control - EMM_F66  0.00000 0.404 180   0.000  1.0000
    ##  EMM_F34 - Control  0.00462 0.452 180   0.010  1.0000
    ##  EMM_F89 - Control  0.00969 0.452 180   0.021  1.0000
    ##  EMM_F70 - Control  0.04667 0.404 180   0.115  1.0000
    ##  EMM_F3 - Control   0.07000 0.404 180   0.173  1.0000
    ##  EMM_F49 - Control  0.08568 0.452 180   0.189  1.0000
    ##  ZAN_F4 - Control   0.25000 0.404 180   0.619  1.0000
    ##  EMM_F48 - Control  0.34469 0.452 180   0.762  1.0000
    ##  EMM_F65 - Control  0.45469 0.452 180   1.006  0.9998
    ##  SP_F14 - Control   0.45667 0.404 180   1.130  0.9990
    ##  EMM_F63 - Control  0.50333 0.404 180   1.245  0.9971
    ##  EMM_F5 - Control   0.54667 0.404 180   1.353  0.9932
    ##  EMM_F64 - Control  0.58000 0.404 180   1.435  0.9877
    ##  EMM_F7 - Control   0.66667 0.404 180   1.649  0.9566
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34  0.05038 0.452 180   0.111  1.0000
    ##  Control - EMM_F49  0.04932 0.452 180   0.109  1.0000
    ##  EMM_F3 - Control   0.03000 0.404 180   0.074  1.0000
    ##  EMM_F65 - Control  0.09969 0.452 180   0.220  1.0000
    ##  EMM_F70 - Control  0.15667 0.404 180   0.388  1.0000
    ##  EMM_F66 - Control  0.23667 0.404 180   0.586  1.0000
    ##  ZAN_F3 - Control   0.24462 0.452 180   0.541  1.0000
    ##  SP_F14 - Control   0.26333 0.404 180   0.652  1.0000
    ##  EMM_F5 - Control   0.34333 0.404 180   0.849  1.0000
    ##  EMM_F48 - Control  0.45969 0.452 180   1.017  0.9997
    ##  EMM_F7 - Control   0.52667 0.404 180   1.303  0.9954
    ##  EMM_F89 - Control  0.54469 0.452 180   1.205  0.9980
    ##  ZAN_F4 - Control   0.57333 0.404 180   1.418  0.9890
    ##  EMM_F63 - Control  0.78333 0.404 180   1.938  0.8545
    ##  EMM_F64 - Control  0.81667 0.404 180   2.021  0.8106
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F49  0.22432 0.452 180   0.496  1.0000
    ##  EMM_F66 - Control  0.13667 0.404 180   0.338  1.0000
    ##  EMM_F3 - Control   0.22000 0.404 180   0.544  1.0000
    ##  EMM_F70 - Control  0.30000 0.404 180   0.742  1.0000
    ##  EMM_F34 - Control  0.32962 0.452 180   0.729  1.0000
    ##  EMM_F65 - Control  0.35469 0.452 180   0.784  1.0000
    ##  EMM_F5 - Control   0.38667 0.404 180   0.957  0.9999
    ##  ZAN_F3 - Control   0.43462 0.452 180   0.961  0.9999
    ##  EMM_F89 - Control  0.45469 0.452 180   1.006  0.9998
    ##  EMM_F48 - Control  0.50969 0.452 180   1.127  0.9991
    ##  SP_F14 - Control   0.52667 0.404 180   1.303  0.9954
    ##  ZAN_F4 - Control   0.55667 0.404 180   1.377  0.9918
    ##  EMM_F7 - Control   0.62667 0.404 180   1.550  0.9746
    ##  EMM_F64 - Control  0.69000 0.404 180   1.707  0.9424
    ##  EMM_F63 - Control  0.75333 0.404 180   1.864  0.8885
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F49 - Control  0.13402 0.452 180   0.296  1.0000
    ##  EMM_F3 - Control   0.41000 0.404 180   1.014  0.9997
    ##  EMM_F7 - Control   0.48000 0.404 180   1.188  0.9983
    ##  EMM_F48 - Control  0.55803 0.452 180   1.234  0.9974
    ##  EMM_F34 - Control  0.61796 0.452 180   1.367  0.9924
    ##  EMM_F89 - Control  0.63303 0.452 180   1.400  0.9904
    ##  EMM_F66 - Control  0.63333 0.404 180   1.567  0.9721
    ##  EMM_F5 - Control   0.65333 0.404 180   1.616  0.9634
    ##  ZAN_F4 - Control   0.65667 0.404 180   1.625  0.9618
    ##  SP_F14 - Control   0.66000 0.404 180   1.633  0.9601
    ##  EMM_F70 - Control  0.70333 0.404 180   1.740  0.9330
    ##  EMM_F63 - Control  0.71000 0.404 180   1.757  0.9280
    ##  EMM_F64 - Control  0.75667 0.404 180   1.872  0.8850
    ##  ZAN_F3 - Control   1.06796 0.452 180   2.362  0.5786
    ##  EMM_F65 - Control  1.07303 0.452 180   2.373  0.5703
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F3   0.39704 0.452 180   0.878  1.0000
    ##  Control - EMM_F49  0.00598 0.452 180   0.013  1.0000
    ##  EMM_F89 - Control  0.35303 0.452 180   0.781  1.0000
    ##  EMM_F34 - Control  0.35796 0.452 180   0.792  1.0000
    ##  EMM_F70 - Control  0.45667 0.404 180   1.130  0.9990
    ##  EMM_F3 - Control   0.47667 0.404 180   1.179  0.9984
    ##  EMM_F63 - Control  0.59333 0.404 180   1.468  0.9847
    ##  EMM_F7 - Control   0.61333 0.404 180   1.517  0.9792
    ##  EMM_F65 - Control  0.65803 0.452 180   1.455  0.9859
    ##  EMM_F66 - Control  0.70000 0.404 180   1.732  0.9355
    ##  EMM_F5 - Control   0.79000 0.404 180   1.955  0.8462
    ##  ZAN_F4 - Control   0.85000 0.404 180   2.103  0.7608
    ##  SP_F14 - Control   0.92000 0.404 180   2.276  0.6416
    ##  EMM_F48 - Control  1.14803 0.452 180   2.539  0.4484
    ##  EMM_F64 - Control  1.23667 0.404 180   3.060  0.1558
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates

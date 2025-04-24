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
library(gvlma)
library(agricolae)
library(emmeans)
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
library(multcompView)
library(dplyr)
library(stringr)
library(tidyverse)
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

## ***Pseudomonas massiliensis* EMM_B5**

Plot is generated using loop around the 4 different classes of yeast,
coming up with 4 plots as an output which will be combined in one plot.

``` r
#read data
B5 <- read.csv("CoCultureAssay/CoCultureAssayData/2024-08-09_PeaceAssay_B5.csv", na.strings = "na") 


#load cbb color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
str(B5)
```

    ## 'data.frame':    1536 obs. of  8 variables:
    ##  $ Bacteria         : chr  "B5" "B5" "B5" "B5" ...
    ##  $ Yeast            : chr  "EMM_F3" "EMM_F3" "EMM_F3" "EMM_F3" ...
    ##  $ Class            : chr  "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" "Dothideomycetes" ...
    ##  $ Replication      : int  1 1 1 1 1 1 1 1 2 2 ...
    ##  $ DAI              : int  2 2 2 2 2 2 2 2 2 2 ...
    ##  $ distance_to_yeast: num  0 11.4 17.4 24.9 31.8 ...
    ##  $ colony_diameter  : num  7.27 7.42 7.51 7.67 7.31 7.48 7.48 7.28 7.78 7.78 ...
    ##  $ increase         : num  0 0 0 0 0 0 0 0 0 0 ...

``` r
#round off of the distance to yeast so that it is visibly a categorical variable
B5$distance_to_yeast <- round(B5$distance_to_yeast, 0)
#change into categorical dataset
B5$Replication=as.factor(B5$Replication)
B5$Yeast=as.factor(B5$Yeast)
B5$DAI=as.factor(B5$DAI)
B5$distance_to_yeast=as.factor(B5$distance_to_yeast)

#remove the data for the distance zero
B5.no.contact <- B5[B5$distance_to_yeast != 0, ]
head(B5.no.contact)
```

    ##   Bacteria  Yeast           Class Replication DAI distance_to_yeast
    ## 2       B5 EMM_F3 Dothideomycetes           1   2                11
    ## 3       B5 EMM_F3 Dothideomycetes           1   2                17
    ## 4       B5 EMM_F3 Dothideomycetes           1   2                25
    ## 5       B5 EMM_F3 Dothideomycetes           1   2                32
    ## 6       B5 EMM_F3 Dothideomycetes           1   2                41
    ## 7       B5 EMM_F3 Dothideomycetes           1   2                48
    ##   colony_diameter increase
    ## 2            7.42        0
    ## 3            7.51        0
    ## 4            7.67        0
    ## 5            7.31        0
    ## 6            7.48        0
    ## 7            7.48        0

``` r
# remove data first data since it is a reference for increase in colony size and donot follow the assumption of linear model i.e. normal distributuion
B5.Day6 <- B5.no.contact[B5.no.contact$DAI == "6",]


plot_listB5 <- list()

# Get unique classes
unique_class <- unique(B5.Day6$Class[B5.Day6$Class != "Control"])

# Loop through each class
for (i in seq_along(unique_class)) {
  
  # Filter for Control and the current class
  data <- B5.Day6 %>%
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
  plot_listB5[[i]] <- p
}

combined_plotB5 <- ggarrange(plotlist = plot_listB5, ncol = 2, nrow = 2, common.legend = FALSE)
```

    ## Warning: Removed 1 row containing non-finite outside the scale range (`stat_summary()`).
    ## Removed 1 row containing non-finite outside the scale range (`stat_summary()`).

``` r
# Annotation for the title
final_plotB5 <- annotate_figure(combined_plotB5,
                              top = text_grob(
    expression("Impact on growth of"~italic("Pseudomonas massiliensis")~"EMM_B5 by Yeast"), color = "Blue2", face = "bold", size = 14, hjust = 0.5))

print(final_plotB5)
```

![](EMM_B5_files/figure-gfm/Plot%20for%20EMM_B5-1.png)<!-- -->

### Stats *Pseudomonas massiliensis* EMM_B5

We are using linear mixed model. Our dependent variable or y is increase
(increase in colony diameter from 1st data) and independent variables
are different Yeast isolates, days after inoculation (DAI), and distance
to yeast which is the distance between the yeast and bacterial colony in
the plate. Each plate is replicated 3 times.

``` r
#filter data to remove 1st day data since the first data is taken as a base to measure the increase in colony size to rule out the variability that is caused by the drop inoculation. So, initially the increase in the colony diameter for 1st data for all colony is "0" that violates the assumption of normality, thus we remove that from analysis. This would be similar for all the bacterial isolates.
library(nlme)
B5.no.1st.data <- B5.no.contact[B5.no.contact$DAI != "2",]
B5try <- lme(increase~distance_to_yeast*Yeast*DAI, data = B5.no.1st.data, random = ~1|Replication, na.action = na.omit)
anova(B5try)
```

    ##                             numDF denDF   F-value p-value
    ## (Intercept)                     1   662 137.37413  <.0001
    ## distance_to_yeast               6   662   6.86297  <.0001
    ## Yeast                          15   662  11.45146  <.0001
    ## DAI                             2   662  17.31295  <.0001
    ## distance_to_yeast:Yeast        90   662   4.70772  <.0001
    ## distance_to_yeast:DAI          12   662   1.31193  0.2065
    ## Yeast:DAI                      30   662   1.83086  0.0048
    ## distance_to_yeast:Yeast:DAI   180   662   1.42602  0.0010

``` r
resultsB5=lme(increase~distance_to_yeast*Yeast+Yeast*DAI, data = B5.no.1st.data, random = ~1|Replication, na.action = na.omit)
summary(resultsB5)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: B5.no.1st.data 
    ##        AIC      BIC    logLik
    ##   1653.995 2347.826 -680.9974
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept)  Residual
    ## StdDev:  0.07482676 0.4423014
    ## 
    ## Fixed effects:  increase ~ distance_to_yeast * Yeast + Yeast * DAI 
    ##                                       Value Std.Error  DF   t-value p-value
    ## (Intercept)                       0.3134921 0.1726661 854  1.815597  0.0698
    ## distance_to_yeast17              -0.2233333 0.2085029 854 -1.071128  0.2844
    ## distance_to_yeast25              -0.0533333 0.2085029 854 -0.255792  0.7982
    ## distance_to_yeast32               0.1033333 0.2085029 854  0.495597  0.6203
    ## distance_to_yeast41              -0.0166667 0.2085029 854 -0.079935  0.9363
    ## distance_to_yeast48               0.1000000 0.2085029 854  0.479610  0.6316
    ## distance_to_yeast55              -0.1111111 0.2085029 854 -0.532900  0.5942
    ## YeastEMM_F3                      -0.0119048 0.2364201 854 -0.050354  0.9599
    ## YeastEMM_F34                      0.0419048 0.2364201 854  0.177247  0.8594
    ## YeastEMM_F47                      0.1971429 0.2364201 854  0.833867  0.4046
    ## YeastEMM_F48                      0.4533333 0.2364201 854  1.917491  0.0555
    ## YeastEMM_F49                      0.1271429 0.2364201 854  0.537784  0.5909
    ## YeastEMM_F63                     -0.0733333 0.2364201 854 -0.310182  0.7565
    ## YeastEMM_F64                      0.0904762 0.2364201 854  0.382693  0.7020
    ## YeastEMM_F65                     -0.8681824 0.2467511 854 -3.518454  0.0005
    ## YeastEMM_F66                      0.2984127 0.2364201 854  1.262214  0.2072
    ## YeastEMM_F7                       0.4361905 0.2364201 854  1.844981  0.0654
    ## YeastEMM_F70                      0.0639683 0.2364201 854  0.270570  0.7868
    ## YeastEMM_F89                      0.1052381 0.2364201 854  0.445132  0.6563
    ## YeastSP_F14                      -2.3152381 0.2364201 854 -9.792900  0.0000
    ## YeastZAN_F3                       0.1487302 0.2364201 854  0.629093  0.5295
    ## YeastZAN_F4                       0.0376190 0.2364201 854  0.159120  0.8736
    ## DAI6                              0.1457143 0.1364972 854  1.067526  0.2860
    ## DAI8                              0.2138095 0.1364972 854  1.566402  0.1176
    ## distance_to_yeast17:YeastEMM_F3   0.4611111 0.2948676 854  1.563790  0.1182
    ## distance_to_yeast25:YeastEMM_F3   0.3266667 0.2948676 854  1.107842  0.2682
    ## distance_to_yeast32:YeastEMM_F3   0.0311111 0.2948676 854  0.105509  0.9160
    ## distance_to_yeast41:YeastEMM_F3   0.3855556 0.2948676 854  1.307555  0.1914
    ## distance_to_yeast48:YeastEMM_F3  -0.0500000 0.2948676 854 -0.169568  0.8654
    ## distance_to_yeast55:YeastEMM_F3   0.1922222 0.2948676 854  0.651893  0.5146
    ## distance_to_yeast17:YeastEMM_F34  0.5544444 0.2948676 854  1.880316  0.0604
    ## distance_to_yeast25:YeastEMM_F34  0.0733333 0.2948676 854  0.248699  0.8037
    ## distance_to_yeast32:YeastEMM_F34 -0.0388889 0.2948676 854 -0.131886  0.8951
    ## distance_to_yeast41:YeastEMM_F34 -0.1922222 0.2948676 854 -0.651893  0.5146
    ## distance_to_yeast48:YeastEMM_F34 -0.1655556 0.2948676 854 -0.561457  0.5746
    ## distance_to_yeast55:YeastEMM_F34  0.1222222 0.2948676 854  0.414499  0.6786
    ## distance_to_yeast17:YeastEMM_F47  0.2366667 0.2948676 854  0.802620  0.4224
    ## distance_to_yeast25:YeastEMM_F47  0.0633333 0.2948676 854  0.214786  0.8300
    ## distance_to_yeast32:YeastEMM_F47 -0.1388889 0.2948676 854 -0.471021  0.6377
    ## distance_to_yeast41:YeastEMM_F47  0.0477778 0.2948676 854  0.162031  0.8713
    ## distance_to_yeast48:YeastEMM_F47 -0.2333333 0.2948676 854 -0.791316  0.4290
    ## distance_to_yeast55:YeastEMM_F47 -0.1588889 0.2948676 854 -0.538848  0.5901
    ## distance_to_yeast17:YeastEMM_F48  0.4066667 0.2948676 854  1.379150  0.1682
    ## distance_to_yeast25:YeastEMM_F48 -0.0155556 0.2948676 854 -0.052754  0.9579
    ## distance_to_yeast32:YeastEMM_F48 -0.1100000 0.2948676 854 -0.373049  0.7092
    ## distance_to_yeast41:YeastEMM_F48 -0.2255556 0.2948676 854 -0.764938  0.4445
    ## distance_to_yeast48:YeastEMM_F48 -0.2788889 0.2948676 854 -0.945810  0.3445
    ## distance_to_yeast55:YeastEMM_F48  0.0533333 0.2948676 854  0.180872  0.8565
    ## distance_to_yeast17:YeastEMM_F49  0.0533333 0.2948676 854  0.180872  0.8565
    ## distance_to_yeast25:YeastEMM_F49 -0.1400000 0.2948676 854 -0.474789  0.6351
    ## distance_to_yeast32:YeastEMM_F49 -0.2233333 0.2948676 854 -0.757402  0.4490
    ## distance_to_yeast41:YeastEMM_F49  0.0422222 0.2948676 854  0.143190  0.8862
    ## distance_to_yeast48:YeastEMM_F49 -0.0655556 0.2948676 854 -0.222322  0.8241
    ## distance_to_yeast55:YeastEMM_F49  0.2333333 0.2948676 854  0.791316  0.4290
    ## distance_to_yeast17:YeastEMM_F63 -0.1055556 0.2948676 854 -0.357976  0.7204
    ## distance_to_yeast25:YeastEMM_F63  0.0744444 0.2948676 854  0.252467  0.8007
    ## distance_to_yeast32:YeastEMM_F63  0.0777778 0.2948676 854  0.263772  0.7920
    ## distance_to_yeast41:YeastEMM_F63 -0.0022222 0.2948676 854 -0.007536  0.9940
    ## distance_to_yeast48:YeastEMM_F63  0.1144444 0.2948676 854  0.388121  0.6980
    ## distance_to_yeast55:YeastEMM_F63  0.2311111 0.2948676 854  0.783779  0.4334
    ## distance_to_yeast17:YeastEMM_F64  0.1011111 0.2948676 854  0.342903  0.7318
    ## distance_to_yeast25:YeastEMM_F64  0.0666667 0.2948676 854  0.226090  0.8212
    ## distance_to_yeast32:YeastEMM_F64 -0.0222222 0.2948676 854 -0.075363  0.9399
    ## distance_to_yeast41:YeastEMM_F64  0.2700000 0.2948676 854  0.915665  0.3601
    ## distance_to_yeast48:YeastEMM_F64  0.0377778 0.2948676 854  0.128118  0.8981
    ## distance_to_yeast55:YeastEMM_F64  0.3233333 0.2948676 854  1.096537  0.2732
    ## distance_to_yeast17:YeastEMM_F65  1.1576943 0.3098473 854  3.736338  0.0002
    ## distance_to_yeast25:YeastEMM_F65  1.4314443 0.3098473 854  4.619838  0.0000
    ## distance_to_yeast32:YeastEMM_F65  1.2047776 0.3098473 854  3.888294  0.0001
    ## distance_to_yeast41:YeastEMM_F65  1.6697776 0.3098473 854  5.389034  0.0000
    ## distance_to_yeast48:YeastEMM_F65  1.1643610 0.3098473 854  3.757854  0.0002
    ## distance_to_yeast55:YeastEMM_F65  1.6692221 0.3098473 854  5.387241  0.0000
    ## distance_to_yeast17:YeastEMM_F66  0.1366667 0.2948676 854  0.463485  0.6431
    ## distance_to_yeast25:YeastEMM_F66 -0.1866667 0.2948676 854 -0.633052  0.5269
    ## distance_to_yeast32:YeastEMM_F66 -0.4400000 0.2948676 854 -1.492195  0.1360
    ## distance_to_yeast41:YeastEMM_F66  0.2800000 0.2948676 854  0.949579  0.3426
    ## distance_to_yeast48:YeastEMM_F66 -0.4266667 0.2948676 854 -1.446977  0.1483
    ## distance_to_yeast55:YeastEMM_F66  0.0044444 0.2948676 854  0.015073  0.9880
    ## distance_to_yeast17:YeastEMM_F7  -0.0022222 0.2948676 854 -0.007536  0.9940
    ## distance_to_yeast25:YeastEMM_F7  -0.1088889 0.2948676 854 -0.369281  0.7120
    ## distance_to_yeast32:YeastEMM_F7  -0.3666667 0.2948676 854 -1.243496  0.2140
    ## distance_to_yeast41:YeastEMM_F7  -0.2800000 0.2948676 854 -0.949579  0.3426
    ## distance_to_yeast48:YeastEMM_F7  -0.4055556 0.2948676 854 -1.375382  0.1694
    ## distance_to_yeast55:YeastEMM_F7  -0.2466667 0.2948676 854 -0.836534  0.4031
    ## distance_to_yeast17:YeastEMM_F70  0.3555556 0.2948676 854  1.205814  0.2282
    ## distance_to_yeast25:YeastEMM_F70  0.2366667 0.2948676 854  0.802620  0.4224
    ## distance_to_yeast32:YeastEMM_F70 -0.1555556 0.2948676 854 -0.527544  0.5980
    ## distance_to_yeast41:YeastEMM_F70  0.3400000 0.2948676 854  1.153060  0.2492
    ## distance_to_yeast48:YeastEMM_F70 -0.2700000 0.2948676 854 -0.915665  0.3601
    ## distance_to_yeast55:YeastEMM_F70 -0.0777778 0.2948676 854 -0.263772  0.7920
    ## distance_to_yeast17:YeastEMM_F89  0.1533333 0.2948676 854  0.520007  0.6032
    ## distance_to_yeast25:YeastEMM_F89 -0.0511111 0.2948676 854 -0.173336  0.8624
    ## distance_to_yeast32:YeastEMM_F89 -0.1300000 0.2948676 854 -0.440876  0.6594
    ## distance_to_yeast41:YeastEMM_F89 -0.0333333 0.2948676 854 -0.113045  0.9100
    ## distance_to_yeast48:YeastEMM_F89 -0.1944444 0.2948676 854 -0.659430  0.5098
    ## distance_to_yeast55:YeastEMM_F89  0.4322222 0.2948676 854  1.465818  0.1431
    ## distance_to_yeast17:YeastSP_F14   2.0744444 0.2948676 854  7.035172  0.0000
    ## distance_to_yeast25:YeastSP_F14   2.6955556 0.2948676 854  9.141578  0.0000
    ## distance_to_yeast32:YeastSP_F14   2.3055556 0.2948676 854  7.818951  0.0000
    ## distance_to_yeast41:YeastSP_F14   2.7922222 0.2948676 854  9.469409  0.0000
    ## distance_to_yeast48:YeastSP_F14   2.4766667 0.2948676 854  8.399249  0.0000
    ## distance_to_yeast55:YeastSP_F14   2.8155556 0.2948676 854  9.548541  0.0000
    ## distance_to_yeast17:YeastZAN_F3   0.2266667 0.2948676 854  0.768707  0.4423
    ## distance_to_yeast25:YeastZAN_F3   0.2644444 0.2948676 854  0.896824  0.3701
    ## distance_to_yeast32:YeastZAN_F3   0.1544444 0.2948676 854  0.523776  0.6006
    ## distance_to_yeast41:YeastZAN_F3  -0.1444444 0.2948676 854 -0.489862  0.6244
    ## distance_to_yeast48:YeastZAN_F3  -0.1233333 0.2948676 854 -0.418267  0.6759
    ## distance_to_yeast55:YeastZAN_F3   0.1111111 0.2948676 854  0.376817  0.7064
    ## distance_to_yeast17:YeastZAN_F4   0.3133333 0.2948676 854  1.062624  0.2883
    ## distance_to_yeast25:YeastZAN_F4   0.0255556 0.2948676 854  0.086668  0.9310
    ## distance_to_yeast32:YeastZAN_F4   0.0577778 0.2948676 854  0.195945  0.8447
    ## distance_to_yeast41:YeastZAN_F4   0.1300000 0.2948676 854  0.440876  0.6594
    ## distance_to_yeast48:YeastZAN_F4   0.0744444 0.2948676 854  0.252467  0.8007
    ## distance_to_yeast55:YeastZAN_F4   0.0622222 0.2948676 854  0.211017  0.8329
    ## YeastEMM_F3:DAI6                  0.0547619 0.1930362 854  0.283687  0.7767
    ## YeastEMM_F34:DAI6                 0.0266667 0.1930362 854  0.138143  0.8902
    ## YeastEMM_F47:DAI6                 0.0600000 0.1930362 854  0.310823  0.7560
    ## YeastEMM_F48:DAI6                 0.0295238 0.1930362 854  0.152944  0.8785
    ## YeastEMM_F49:DAI6                -0.0152381 0.1930362 854 -0.078939  0.9371
    ## YeastEMM_F63:DAI6                 0.0452381 0.1930362 854  0.234350  0.8148
    ## YeastEMM_F64:DAI6                 0.0023810 0.1930362 854  0.012334  0.9902
    ## YeastEMM_F65:DAI6                -0.1521375 0.1943938 854 -0.782625  0.4341
    ## YeastEMM_F66:DAI6                -0.0047619 0.1930362 854 -0.024668  0.9803
    ## YeastEMM_F7:DAI6                  0.0871429 0.1930362 854  0.451433  0.6518
    ## YeastEMM_F70:DAI6                -0.0380952 0.1930362 854 -0.197348  0.8436
    ## YeastEMM_F89:DAI6                -0.0242857 0.1930362 854 -0.125809  0.8999
    ## YeastSP_F14:DAI6                 -0.4057143 0.1930362 854 -2.101753  0.0359
    ## YeastZAN_F3:DAI6                  0.1633333 0.1930362 854  0.846128  0.3977
    ## YeastZAN_F4:DAI6                 -0.0771429 0.1930362 854 -0.399629  0.6895
    ## YeastEMM_F3:DAI8                  0.1609524 0.1930362 854  0.833794  0.4046
    ## YeastEMM_F34:DAI8                 0.0576190 0.1930362 854  0.298488  0.7654
    ## YeastEMM_F47:DAI8                 0.0885714 0.1930362 854  0.458833  0.6465
    ## YeastEMM_F48:DAI8                 0.0638095 0.1930362 854  0.330557  0.7411
    ## YeastEMM_F49:DAI8                 0.1004762 0.1930362 854  0.520504  0.6028
    ## YeastEMM_F63:DAI8                 0.0880952 0.1930362 854  0.456366  0.6482
    ## YeastEMM_F64:DAI8                 0.0361905 0.1930362 854  0.187480  0.8513
    ## YeastEMM_F65:DAI8                -0.8497443 0.2049649 854 -4.145803  0.0000
    ## YeastEMM_F66:DAI8                -0.0138095 0.1930362 854 -0.071539  0.9430
    ## YeastEMM_F7:DAI8                 -0.0523810 0.1930362 854 -0.271353  0.7862
    ## YeastEMM_F70:DAI8                 0.0728571 0.1930362 854  0.377427  0.7059
    ## YeastEMM_F89:DAI8                -0.1347619 0.1930362 854 -0.698117  0.4853
    ## YeastSP_F14:DAI8                 -0.2985714 0.1930362 854 -1.546712  0.1223
    ## YeastZAN_F3:DAI8                  0.0804762 0.1930362 854  0.416897  0.6769
    ## YeastZAN_F4:DAI8                 -0.0190476 0.1930362 854 -0.098674  0.9214
    ##  Correlation: 
    ##                                  (Intr) ds__17 ds__25 ds__32 ds__41 ds__48
    ## distance_to_yeast17              -0.604                                   
    ## distance_to_yeast25              -0.604  0.500                            
    ## distance_to_yeast32              -0.604  0.500  0.500                     
    ## distance_to_yeast41              -0.604  0.500  0.500  0.500              
    ## distance_to_yeast48              -0.604  0.500  0.500  0.500  0.500       
    ## distance_to_yeast55              -0.604  0.500  0.500  0.500  0.500  0.500
    ## YeastEMM_F3                      -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F34                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F47                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F48                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F49                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F63                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F64                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F65                     -0.656  0.422  0.422  0.422  0.422  0.422
    ## YeastEMM_F66                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F7                      -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F70                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastEMM_F89                     -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastSP_F14                      -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastZAN_F3                      -0.685  0.441  0.441  0.441  0.441  0.441
    ## YeastZAN_F4                      -0.685  0.441  0.441  0.441  0.441  0.441
    ## DAI6                             -0.395  0.000  0.000  0.000  0.000  0.000
    ## DAI8                             -0.395  0.000  0.000  0.000  0.000  0.000
    ## distance_to_yeast17:YeastEMM_F3   0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F3   0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F3   0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F3   0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F3   0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F3   0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F34  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F34  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F34  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F34  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F34  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F34  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F47  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F47  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F47  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F47  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F47  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F47  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F48  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F48  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F48  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F48  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F48  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F48  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F49  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F49  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F49  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F49  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F49  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F49  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F63  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F63  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F63  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F63  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F63  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F63  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F64  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F64  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F64  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F64  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F64  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F64  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F65  0.406 -0.673 -0.336 -0.336 -0.336 -0.336
    ## distance_to_yeast25:YeastEMM_F65  0.406 -0.336 -0.673 -0.336 -0.336 -0.336
    ## distance_to_yeast32:YeastEMM_F65  0.406 -0.336 -0.336 -0.673 -0.336 -0.336
    ## distance_to_yeast41:YeastEMM_F65  0.406 -0.336 -0.336 -0.336 -0.673 -0.336
    ## distance_to_yeast48:YeastEMM_F65  0.406 -0.336 -0.336 -0.336 -0.336 -0.673
    ## distance_to_yeast55:YeastEMM_F65  0.406 -0.336 -0.336 -0.336 -0.336 -0.336
    ## distance_to_yeast17:YeastEMM_F66  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F66  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F66  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F66  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F66  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F66  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F7   0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F7   0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F7   0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F7   0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F7   0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F7   0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F70  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F70  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F70  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F70  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F70  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F70  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastEMM_F89  0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastEMM_F89  0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastEMM_F89  0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastEMM_F89  0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastEMM_F89  0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastEMM_F89  0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastSP_F14   0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastSP_F14   0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastSP_F14   0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastSP_F14   0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastSP_F14   0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastSP_F14   0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastZAN_F3   0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastZAN_F3   0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastZAN_F3   0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastZAN_F3   0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastZAN_F3   0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastZAN_F3   0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast17:YeastZAN_F4   0.427 -0.707 -0.354 -0.354 -0.354 -0.354
    ## distance_to_yeast25:YeastZAN_F4   0.427 -0.354 -0.707 -0.354 -0.354 -0.354
    ## distance_to_yeast32:YeastZAN_F4   0.427 -0.354 -0.354 -0.707 -0.354 -0.354
    ## distance_to_yeast41:YeastZAN_F4   0.427 -0.354 -0.354 -0.354 -0.707 -0.354
    ## distance_to_yeast48:YeastZAN_F4   0.427 -0.354 -0.354 -0.354 -0.354 -0.707
    ## distance_to_yeast55:YeastZAN_F4   0.427 -0.354 -0.354 -0.354 -0.354 -0.354
    ## YeastEMM_F3:DAI6                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F34:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F47:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F48:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F49:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F63:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F64:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F65:DAI6                 0.278  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F66:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F7:DAI6                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F70:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F89:DAI6                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastSP_F14:DAI6                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastZAN_F3:DAI6                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastZAN_F4:DAI6                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F3:DAI8                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F34:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F47:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F48:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F49:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F63:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F64:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F65:DAI8                 0.263  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F66:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F7:DAI8                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F70:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastEMM_F89:DAI8                 0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastSP_F14:DAI8                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastZAN_F3:DAI8                  0.279  0.000  0.000  0.000  0.000  0.000
    ## YeastZAN_F4:DAI8                  0.279  0.000  0.000  0.000  0.000  0.000
    ##                                  ds__55 YsEMM_F3 YsEMM_F34 YsEMM_F47 YsEMM_F48
    ## distance_to_yeast17                                                           
    ## distance_to_yeast25                                                           
    ## distance_to_yeast32                                                           
    ## distance_to_yeast41                                                           
    ## distance_to_yeast48                                                           
    ## distance_to_yeast55                                                           
    ## YeastEMM_F3                       0.441                                       
    ## YeastEMM_F34                      0.441  0.500                                
    ## YeastEMM_F47                      0.441  0.500    0.500                       
    ## YeastEMM_F48                      0.441  0.500    0.500     0.500             
    ## YeastEMM_F49                      0.441  0.500    0.500     0.500     0.500   
    ## YeastEMM_F63                      0.441  0.500    0.500     0.500     0.500   
    ## YeastEMM_F64                      0.441  0.500    0.500     0.500     0.500   
    ## YeastEMM_F65                      0.422  0.479    0.479     0.479     0.479   
    ## YeastEMM_F66                      0.441  0.500    0.500     0.500     0.500   
    ## YeastEMM_F7                       0.441  0.500    0.500     0.500     0.500   
    ## YeastEMM_F70                      0.441  0.500    0.500     0.500     0.500   
    ## YeastEMM_F89                      0.441  0.500    0.500     0.500     0.500   
    ## YeastSP_F14                       0.441  0.500    0.500     0.500     0.500   
    ## YeastZAN_F3                       0.441  0.500    0.500     0.500     0.500   
    ## YeastZAN_F4                       0.441  0.500    0.500     0.500     0.500   
    ## DAI6                              0.000  0.289    0.289     0.289     0.289   
    ## DAI8                              0.000  0.289    0.289     0.289     0.289   
    ## distance_to_yeast17:YeastEMM_F3  -0.354 -0.624   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F3  -0.354 -0.624   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F3  -0.354 -0.624   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F3  -0.354 -0.624   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F3  -0.354 -0.624   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F3  -0.707 -0.624   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F34 -0.354 -0.312   -0.624    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F34 -0.354 -0.312   -0.624    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F34 -0.354 -0.312   -0.624    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F34 -0.354 -0.312   -0.624    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F34 -0.354 -0.312   -0.624    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F34 -0.707 -0.312   -0.624    -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F47 -0.354 -0.312   -0.312    -0.624    -0.312   
    ## distance_to_yeast25:YeastEMM_F47 -0.354 -0.312   -0.312    -0.624    -0.312   
    ## distance_to_yeast32:YeastEMM_F47 -0.354 -0.312   -0.312    -0.624    -0.312   
    ## distance_to_yeast41:YeastEMM_F47 -0.354 -0.312   -0.312    -0.624    -0.312   
    ## distance_to_yeast48:YeastEMM_F47 -0.354 -0.312   -0.312    -0.624    -0.312   
    ## distance_to_yeast55:YeastEMM_F47 -0.707 -0.312   -0.312    -0.624    -0.312   
    ## distance_to_yeast17:YeastEMM_F48 -0.354 -0.312   -0.312    -0.312    -0.624   
    ## distance_to_yeast25:YeastEMM_F48 -0.354 -0.312   -0.312    -0.312    -0.624   
    ## distance_to_yeast32:YeastEMM_F48 -0.354 -0.312   -0.312    -0.312    -0.624   
    ## distance_to_yeast41:YeastEMM_F48 -0.354 -0.312   -0.312    -0.312    -0.624   
    ## distance_to_yeast48:YeastEMM_F48 -0.354 -0.312   -0.312    -0.312    -0.624   
    ## distance_to_yeast55:YeastEMM_F48 -0.707 -0.312   -0.312    -0.312    -0.624   
    ## distance_to_yeast17:YeastEMM_F49 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F49 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F49 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F49 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F49 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F49 -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F63 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F63 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F63 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F63 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F63 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F63 -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F64 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F64 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F64 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F64 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F64 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F64 -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F65 -0.336 -0.297   -0.297    -0.297    -0.297   
    ## distance_to_yeast25:YeastEMM_F65 -0.336 -0.297   -0.297    -0.297    -0.297   
    ## distance_to_yeast32:YeastEMM_F65 -0.336 -0.297   -0.297    -0.297    -0.297   
    ## distance_to_yeast41:YeastEMM_F65 -0.336 -0.297   -0.297    -0.297    -0.297   
    ## distance_to_yeast48:YeastEMM_F65 -0.336 -0.297   -0.297    -0.297    -0.297   
    ## distance_to_yeast55:YeastEMM_F65 -0.673 -0.297   -0.297    -0.297    -0.297   
    ## distance_to_yeast17:YeastEMM_F66 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F66 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F66 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F66 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F66 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F66 -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F7  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F7  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F7  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F7  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F7  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F7  -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F70 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F70 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F70 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F70 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F70 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F70 -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F89 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F89 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F89 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F89 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F89 -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F89 -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastSP_F14  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastSP_F14  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastSP_F14  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastSP_F14  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastSP_F14  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastSP_F14  -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastZAN_F3  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastZAN_F3  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastZAN_F3  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastZAN_F3  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastZAN_F3  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastZAN_F3  -0.707 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast17:YeastZAN_F4  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast25:YeastZAN_F4  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast32:YeastZAN_F4  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast41:YeastZAN_F4  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast48:YeastZAN_F4  -0.354 -0.312   -0.312    -0.312    -0.312   
    ## distance_to_yeast55:YeastZAN_F4  -0.707 -0.312   -0.312    -0.312    -0.312   
    ## YeastEMM_F3:DAI6                  0.000 -0.408   -0.204    -0.204    -0.204   
    ## YeastEMM_F34:DAI6                 0.000 -0.204   -0.408    -0.204    -0.204   
    ## YeastEMM_F47:DAI6                 0.000 -0.204   -0.204    -0.408    -0.204   
    ## YeastEMM_F48:DAI6                 0.000 -0.204   -0.204    -0.204    -0.408   
    ## YeastEMM_F49:DAI6                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F63:DAI6                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F64:DAI6                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F65:DAI6                 0.000 -0.203   -0.203    -0.203    -0.203   
    ## YeastEMM_F66:DAI6                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F7:DAI6                  0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F70:DAI6                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F89:DAI6                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastSP_F14:DAI6                  0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastZAN_F3:DAI6                  0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastZAN_F4:DAI6                  0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F3:DAI8                  0.000 -0.408   -0.204    -0.204    -0.204   
    ## YeastEMM_F34:DAI8                 0.000 -0.204   -0.408    -0.204    -0.204   
    ## YeastEMM_F47:DAI8                 0.000 -0.204   -0.204    -0.408    -0.204   
    ## YeastEMM_F48:DAI8                 0.000 -0.204   -0.204    -0.204    -0.408   
    ## YeastEMM_F49:DAI8                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F63:DAI8                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F64:DAI8                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F65:DAI8                 0.000 -0.192   -0.192    -0.192    -0.192   
    ## YeastEMM_F66:DAI8                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F7:DAI8                  0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F70:DAI8                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastEMM_F89:DAI8                 0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastSP_F14:DAI8                  0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastZAN_F3:DAI8                  0.000 -0.204   -0.204    -0.204    -0.204   
    ## YeastZAN_F4:DAI8                  0.000 -0.204   -0.204    -0.204    -0.204   
    ##                                  YsEMM_F49 YsEMM_F63 YsEMM_F64 YsEMM_F65
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
    ## YeastEMM_F63                      0.500                                 
    ## YeastEMM_F64                      0.500     0.500                       
    ## YeastEMM_F65                      0.479     0.479     0.479             
    ## YeastEMM_F66                      0.500     0.500     0.500     0.479   
    ## YeastEMM_F7                       0.500     0.500     0.500     0.479   
    ## YeastEMM_F70                      0.500     0.500     0.500     0.479   
    ## YeastEMM_F89                      0.500     0.500     0.500     0.479   
    ## YeastSP_F14                       0.500     0.500     0.500     0.479   
    ## YeastZAN_F3                       0.500     0.500     0.500     0.479   
    ## YeastZAN_F4                       0.500     0.500     0.500     0.479   
    ## DAI6                              0.289     0.289     0.289     0.277   
    ## DAI8                              0.289     0.289     0.289     0.277   
    ## distance_to_yeast17:YeastEMM_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F34 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F34 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F34 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F34 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F34 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F34 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F47 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F47 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F47 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F47 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F47 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F47 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F48 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F48 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F48 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F48 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F48 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F48 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F49 -0.624    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F49 -0.624    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F49 -0.624    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F49 -0.624    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F49 -0.624    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F49 -0.624    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F63 -0.312    -0.624    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F63 -0.312    -0.624    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F63 -0.312    -0.624    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F63 -0.312    -0.624    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F63 -0.312    -0.624    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F63 -0.312    -0.624    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F64 -0.312    -0.312    -0.624    -0.299   
    ## distance_to_yeast25:YeastEMM_F64 -0.312    -0.312    -0.624    -0.299   
    ## distance_to_yeast32:YeastEMM_F64 -0.312    -0.312    -0.624    -0.299   
    ## distance_to_yeast41:YeastEMM_F64 -0.312    -0.312    -0.624    -0.299   
    ## distance_to_yeast48:YeastEMM_F64 -0.312    -0.312    -0.624    -0.299   
    ## distance_to_yeast55:YeastEMM_F64 -0.312    -0.312    -0.624    -0.299   
    ## distance_to_yeast17:YeastEMM_F65 -0.297    -0.297    -0.297    -0.645   
    ## distance_to_yeast25:YeastEMM_F65 -0.297    -0.297    -0.297    -0.645   
    ## distance_to_yeast32:YeastEMM_F65 -0.297    -0.297    -0.297    -0.645   
    ## distance_to_yeast41:YeastEMM_F65 -0.297    -0.297    -0.297    -0.645   
    ## distance_to_yeast48:YeastEMM_F65 -0.297    -0.297    -0.297    -0.645   
    ## distance_to_yeast55:YeastEMM_F65 -0.297    -0.297    -0.297    -0.645   
    ## distance_to_yeast17:YeastEMM_F66 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F66 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F66 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F66 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F66 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F66 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F7  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F7  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F7  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F7  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F7  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F7  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F70 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F70 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F70 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F70 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F70 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F70 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastEMM_F89 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastEMM_F89 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastEMM_F89 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastEMM_F89 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastEMM_F89 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastEMM_F89 -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastSP_F14  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastSP_F14  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastSP_F14  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastSP_F14  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastSP_F14  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastSP_F14  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastZAN_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastZAN_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastZAN_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastZAN_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastZAN_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastZAN_F3  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast17:YeastZAN_F4  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast25:YeastZAN_F4  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast32:YeastZAN_F4  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast41:YeastZAN_F4  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast48:YeastZAN_F4  -0.312    -0.312    -0.312    -0.299   
    ## distance_to_yeast55:YeastZAN_F4  -0.312    -0.312    -0.312    -0.299   
    ## YeastEMM_F3:DAI6                 -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F34:DAI6                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F47:DAI6                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F48:DAI6                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F49:DAI6                -0.408    -0.204    -0.204    -0.196   
    ## YeastEMM_F63:DAI6                -0.204    -0.408    -0.204    -0.196   
    ## YeastEMM_F64:DAI6                -0.204    -0.204    -0.408    -0.196   
    ## YeastEMM_F65:DAI6                -0.203    -0.203    -0.203    -0.364   
    ## YeastEMM_F66:DAI6                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F7:DAI6                 -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F70:DAI6                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F89:DAI6                -0.204    -0.204    -0.204    -0.196   
    ## YeastSP_F14:DAI6                 -0.204    -0.204    -0.204    -0.196   
    ## YeastZAN_F3:DAI6                 -0.204    -0.204    -0.204    -0.196   
    ## YeastZAN_F4:DAI6                 -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F3:DAI8                 -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F34:DAI8                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F47:DAI8                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F48:DAI8                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F49:DAI8                -0.408    -0.204    -0.204    -0.196   
    ## YeastEMM_F63:DAI8                -0.204    -0.408    -0.204    -0.196   
    ## YeastEMM_F64:DAI8                -0.204    -0.204    -0.408    -0.196   
    ## YeastEMM_F65:DAI8                -0.192    -0.192    -0.192    -0.368   
    ## YeastEMM_F66:DAI8                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F7:DAI8                 -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F70:DAI8                -0.204    -0.204    -0.204    -0.196   
    ## YeastEMM_F89:DAI8                -0.204    -0.204    -0.204    -0.196   
    ## YeastSP_F14:DAI8                 -0.204    -0.204    -0.204    -0.196   
    ## YeastZAN_F3:DAI8                 -0.204    -0.204    -0.204    -0.196   
    ## YeastZAN_F4:DAI8                 -0.204    -0.204    -0.204    -0.196   
    ##                                  YsEMM_F66 YsEMM_F7 YsEMM_F70 YsEMM_F89
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
    ## YeastEMM_F7                       0.500                                
    ## YeastEMM_F70                      0.500     0.500                      
    ## YeastEMM_F89                      0.500     0.500    0.500             
    ## YeastSP_F14                       0.500     0.500    0.500     0.500   
    ## YeastZAN_F3                       0.500     0.500    0.500     0.500   
    ## YeastZAN_F4                       0.500     0.500    0.500     0.500   
    ## DAI6                              0.289     0.289    0.289     0.289   
    ## DAI8                              0.289     0.289    0.289     0.289   
    ## distance_to_yeast17:YeastEMM_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F34 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F34 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F34 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F34 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F34 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F34 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F47 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F47 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F47 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F47 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F47 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F47 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F48 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F48 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F48 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F48 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F48 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F48 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F49 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F49 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F49 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F49 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F49 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F49 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F63 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F63 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F63 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F63 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F63 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F63 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F64 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F64 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F64 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F64 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F64 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F64 -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F65 -0.297    -0.297   -0.297    -0.297   
    ## distance_to_yeast25:YeastEMM_F65 -0.297    -0.297   -0.297    -0.297   
    ## distance_to_yeast32:YeastEMM_F65 -0.297    -0.297   -0.297    -0.297   
    ## distance_to_yeast41:YeastEMM_F65 -0.297    -0.297   -0.297    -0.297   
    ## distance_to_yeast48:YeastEMM_F65 -0.297    -0.297   -0.297    -0.297   
    ## distance_to_yeast55:YeastEMM_F65 -0.297    -0.297   -0.297    -0.297   
    ## distance_to_yeast17:YeastEMM_F66 -0.624    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F66 -0.624    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F66 -0.624    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F66 -0.624    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F66 -0.624    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F66 -0.624    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F7  -0.312    -0.624   -0.312    -0.312   
    ## distance_to_yeast25:YeastEMM_F7  -0.312    -0.624   -0.312    -0.312   
    ## distance_to_yeast32:YeastEMM_F7  -0.312    -0.624   -0.312    -0.312   
    ## distance_to_yeast41:YeastEMM_F7  -0.312    -0.624   -0.312    -0.312   
    ## distance_to_yeast48:YeastEMM_F7  -0.312    -0.624   -0.312    -0.312   
    ## distance_to_yeast55:YeastEMM_F7  -0.312    -0.624   -0.312    -0.312   
    ## distance_to_yeast17:YeastEMM_F70 -0.312    -0.312   -0.624    -0.312   
    ## distance_to_yeast25:YeastEMM_F70 -0.312    -0.312   -0.624    -0.312   
    ## distance_to_yeast32:YeastEMM_F70 -0.312    -0.312   -0.624    -0.312   
    ## distance_to_yeast41:YeastEMM_F70 -0.312    -0.312   -0.624    -0.312   
    ## distance_to_yeast48:YeastEMM_F70 -0.312    -0.312   -0.624    -0.312   
    ## distance_to_yeast55:YeastEMM_F70 -0.312    -0.312   -0.624    -0.312   
    ## distance_to_yeast17:YeastEMM_F89 -0.312    -0.312   -0.312    -0.624   
    ## distance_to_yeast25:YeastEMM_F89 -0.312    -0.312   -0.312    -0.624   
    ## distance_to_yeast32:YeastEMM_F89 -0.312    -0.312   -0.312    -0.624   
    ## distance_to_yeast41:YeastEMM_F89 -0.312    -0.312   -0.312    -0.624   
    ## distance_to_yeast48:YeastEMM_F89 -0.312    -0.312   -0.312    -0.624   
    ## distance_to_yeast55:YeastEMM_F89 -0.312    -0.312   -0.312    -0.624   
    ## distance_to_yeast17:YeastSP_F14  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastSP_F14  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastSP_F14  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastSP_F14  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastSP_F14  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastSP_F14  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastZAN_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastZAN_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastZAN_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastZAN_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastZAN_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastZAN_F3  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast17:YeastZAN_F4  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast25:YeastZAN_F4  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast32:YeastZAN_F4  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast41:YeastZAN_F4  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast48:YeastZAN_F4  -0.312    -0.312   -0.312    -0.312   
    ## distance_to_yeast55:YeastZAN_F4  -0.312    -0.312   -0.312    -0.312   
    ## YeastEMM_F3:DAI6                 -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F34:DAI6                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F47:DAI6                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F48:DAI6                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F49:DAI6                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F63:DAI6                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F64:DAI6                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F65:DAI6                -0.203    -0.203   -0.203    -0.203   
    ## YeastEMM_F66:DAI6                -0.408    -0.204   -0.204    -0.204   
    ## YeastEMM_F7:DAI6                 -0.204    -0.408   -0.204    -0.204   
    ## YeastEMM_F70:DAI6                -0.204    -0.204   -0.408    -0.204   
    ## YeastEMM_F89:DAI6                -0.204    -0.204   -0.204    -0.408   
    ## YeastSP_F14:DAI6                 -0.204    -0.204   -0.204    -0.204   
    ## YeastZAN_F3:DAI6                 -0.204    -0.204   -0.204    -0.204   
    ## YeastZAN_F4:DAI6                 -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F3:DAI8                 -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F34:DAI8                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F47:DAI8                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F48:DAI8                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F49:DAI8                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F63:DAI8                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F64:DAI8                -0.204    -0.204   -0.204    -0.204   
    ## YeastEMM_F65:DAI8                -0.192    -0.192   -0.192    -0.192   
    ## YeastEMM_F66:DAI8                -0.408    -0.204   -0.204    -0.204   
    ## YeastEMM_F7:DAI8                 -0.204    -0.408   -0.204    -0.204   
    ## YeastEMM_F70:DAI8                -0.204    -0.204   -0.408    -0.204   
    ## YeastEMM_F89:DAI8                -0.204    -0.204   -0.204    -0.408   
    ## YeastSP_F14:DAI8                 -0.204    -0.204   -0.204    -0.204   
    ## YeastZAN_F3:DAI8                 -0.204    -0.204   -0.204    -0.204   
    ## YeastZAN_F4:DAI8                 -0.204    -0.204   -0.204    -0.204   
    ##                                  YsSP_F14 YsZAN_F3 YsZAN_F4 DAI6   DAI8  
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
    ## YeastZAN_F4                       0.500    0.500                         
    ## DAI6                              0.289    0.289    0.289                
    ## DAI8                              0.289    0.289    0.289    0.500       
    ## distance_to_yeast17:YeastEMM_F3  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F3  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F3  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F3  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F3  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F3  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F34 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F34 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F34 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F34 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F34 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F34 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F47 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F47 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F47 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F47 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F47 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F47 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F48 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F48 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F48 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F48 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F48 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F48 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F49 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F49 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F49 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F49 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F49 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F49 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F63 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F63 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F63 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F63 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F63 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F63 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F64 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F64 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F64 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F64 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F64 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F64 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F65 -0.297   -0.297   -0.297    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F65 -0.297   -0.297   -0.297    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F65 -0.297   -0.297   -0.297    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F65 -0.297   -0.297   -0.297    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F65 -0.297   -0.297   -0.297    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F65 -0.297   -0.297   -0.297    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F66 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F66 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F66 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F66 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F66 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F66 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F7  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F7  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F7  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F7  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F7  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F7  -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F70 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F70 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F70 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F70 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F70 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F70 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastEMM_F89 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastEMM_F89 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastEMM_F89 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastEMM_F89 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastEMM_F89 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastEMM_F89 -0.312   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastSP_F14  -0.624   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastSP_F14  -0.624   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastSP_F14  -0.624   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastSP_F14  -0.624   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastSP_F14  -0.624   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastSP_F14  -0.624   -0.312   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastZAN_F3  -0.312   -0.624   -0.312    0.000  0.000
    ## distance_to_yeast25:YeastZAN_F3  -0.312   -0.624   -0.312    0.000  0.000
    ## distance_to_yeast32:YeastZAN_F3  -0.312   -0.624   -0.312    0.000  0.000
    ## distance_to_yeast41:YeastZAN_F3  -0.312   -0.624   -0.312    0.000  0.000
    ## distance_to_yeast48:YeastZAN_F3  -0.312   -0.624   -0.312    0.000  0.000
    ## distance_to_yeast55:YeastZAN_F3  -0.312   -0.624   -0.312    0.000  0.000
    ## distance_to_yeast17:YeastZAN_F4  -0.312   -0.312   -0.624    0.000  0.000
    ## distance_to_yeast25:YeastZAN_F4  -0.312   -0.312   -0.624    0.000  0.000
    ## distance_to_yeast32:YeastZAN_F4  -0.312   -0.312   -0.624    0.000  0.000
    ## distance_to_yeast41:YeastZAN_F4  -0.312   -0.312   -0.624    0.000  0.000
    ## distance_to_yeast48:YeastZAN_F4  -0.312   -0.312   -0.624    0.000  0.000
    ## distance_to_yeast55:YeastZAN_F4  -0.312   -0.312   -0.624    0.000  0.000
    ## YeastEMM_F3:DAI6                 -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F34:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F47:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F48:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F49:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F63:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F64:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F65:DAI6                -0.203   -0.203   -0.203   -0.702 -0.351
    ## YeastEMM_F66:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F7:DAI6                 -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F70:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastEMM_F89:DAI6                -0.204   -0.204   -0.204   -0.707 -0.354
    ## YeastSP_F14:DAI6                 -0.408   -0.204   -0.204   -0.707 -0.354
    ## YeastZAN_F3:DAI6                 -0.204   -0.408   -0.204   -0.707 -0.354
    ## YeastZAN_F4:DAI6                 -0.204   -0.204   -0.408   -0.707 -0.354
    ## YeastEMM_F3:DAI8                 -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F34:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F47:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F48:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F49:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F63:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F64:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F65:DAI8                -0.192   -0.192   -0.192   -0.333 -0.666
    ## YeastEMM_F66:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F7:DAI8                 -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F70:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastEMM_F89:DAI8                -0.204   -0.204   -0.204   -0.354 -0.707
    ## YeastSP_F14:DAI8                 -0.408   -0.204   -0.204   -0.354 -0.707
    ## YeastZAN_F3:DAI8                 -0.204   -0.408   -0.204   -0.354 -0.707
    ## YeastZAN_F4:DAI8                 -0.204   -0.204   -0.408   -0.354 -0.707
    ##                                  ds__17:YEMM_F3 ds__25:YEMM_F3 ds__32:YEMM_F3
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
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3   0.500                                      
    ## distance_to_yeast32:YeastEMM_F3   0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F3   0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F3   0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F3   0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F34  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F34  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F34  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F47  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F47  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  ds__41:YEMM_F3 ds__48:YEMM_F3 ds__55:YEMM_F3
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
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3   0.500                                      
    ## distance_to_yeast55:YeastEMM_F3   0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F34  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F34  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F34  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F34  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F47  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F47  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F34 d__25:YEMM_F34 d__32:YEMM_F34
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
    ## DAI6                                                                         
    ## DAI8                                                                         
    ## distance_to_yeast17:YeastEMM_F3                                              
    ## distance_to_yeast25:YeastEMM_F3                                              
    ## distance_to_yeast32:YeastEMM_F3                                              
    ## distance_to_yeast41:YeastEMM_F3                                              
    ## distance_to_yeast48:YeastEMM_F3                                              
    ## distance_to_yeast55:YeastEMM_F3                                              
    ## distance_to_yeast17:YeastEMM_F34                                             
    ## distance_to_yeast25:YeastEMM_F34  0.500                                      
    ## distance_to_yeast32:YeastEMM_F34  0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F34  0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F34  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F34  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F47  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F47  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F34 d__48:YEMM_F34 d__55:YEMM_F34
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F34  0.500                                      
    ## distance_to_yeast55:YeastEMM_F34  0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F47  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F47  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F47  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F47  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F47 d__25:YEMM_F47 d__32:YEMM_F47
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F47  0.500                                      
    ## distance_to_yeast32:YeastEMM_F47  0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F47  0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F47  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F47  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F47 d__48:YEMM_F47 d__55:YEMM_F47
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F47  0.500                                      
    ## distance_to_yeast55:YeastEMM_F47  0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F48  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F48  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F48  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F48  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F48 d__25:YEMM_F48 d__32:YEMM_F48
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F48  0.500                                      
    ## distance_to_yeast32:YeastEMM_F48  0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F48  0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F48  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F48  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F48 d__48:YEMM_F48 d__55:YEMM_F48
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F48  0.500                                      
    ## distance_to_yeast55:YeastEMM_F48  0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F49  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F49  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F49  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F49  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F49 d__25:YEMM_F49 d__32:YEMM_F49
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F49  0.500                                      
    ## distance_to_yeast32:YeastEMM_F49  0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F49  0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F49  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F49  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F49 d__48:YEMM_F49 d__55:YEMM_F49
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F49  0.500                                      
    ## distance_to_yeast55:YeastEMM_F49  0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F63  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F63  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F63  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F63  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F63 d__25:YEMM_F63 d__32:YEMM_F63
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F63  0.500                                      
    ## distance_to_yeast32:YeastEMM_F63  0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F63  0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F63  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F63  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F63 d__48:YEMM_F63 d__55:YEMM_F63
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F63  0.500                                      
    ## distance_to_yeast55:YeastEMM_F63  0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F64  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F64  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F64  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F64  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F64 d__25:YEMM_F64 d__32:YEMM_F64
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F64  0.500                                      
    ## distance_to_yeast32:YeastEMM_F64  0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F64  0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F64  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F64  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F64 d__48:YEMM_F64 d__55:YEMM_F64
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F64  0.500                                      
    ## distance_to_yeast55:YeastEMM_F64  0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F65  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F65  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F65  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F65  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F66  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F66  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F66  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F65 d__25:YEMM_F65 d__32:YEMM_F65
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F65  0.519                                      
    ## distance_to_yeast32:YeastEMM_F65  0.519          0.519                       
    ## distance_to_yeast41:YeastEMM_F65  0.519          0.519          0.519        
    ## distance_to_yeast48:YeastEMM_F65  0.519          0.519          0.519        
    ## distance_to_yeast55:YeastEMM_F65  0.519          0.519          0.519        
    ## distance_to_yeast17:YeastEMM_F66  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F66  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F66  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F66  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F66  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F66  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F7   0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F7   0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F7   0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F7   0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F7   0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F7   0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F70  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F70  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F70  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F70  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F70  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F70  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastEMM_F89  0.476          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F89  0.238          0.476          0.238        
    ## distance_to_yeast32:YeastEMM_F89  0.238          0.238          0.476        
    ## distance_to_yeast41:YeastEMM_F89  0.238          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F89  0.238          0.238          0.238        
    ## distance_to_yeast55:YeastEMM_F89  0.238          0.238          0.238        
    ## distance_to_yeast17:YeastSP_F14   0.476          0.238          0.238        
    ## distance_to_yeast25:YeastSP_F14   0.238          0.476          0.238        
    ## distance_to_yeast32:YeastSP_F14   0.238          0.238          0.476        
    ## distance_to_yeast41:YeastSP_F14   0.238          0.238          0.238        
    ## distance_to_yeast48:YeastSP_F14   0.238          0.238          0.238        
    ## distance_to_yeast55:YeastSP_F14   0.238          0.238          0.238        
    ## distance_to_yeast17:YeastZAN_F3   0.476          0.238          0.238        
    ## distance_to_yeast25:YeastZAN_F3   0.238          0.476          0.238        
    ## distance_to_yeast32:YeastZAN_F3   0.238          0.238          0.476        
    ## distance_to_yeast41:YeastZAN_F3   0.238          0.238          0.238        
    ## distance_to_yeast48:YeastZAN_F3   0.238          0.238          0.238        
    ## distance_to_yeast55:YeastZAN_F3   0.238          0.238          0.238        
    ## distance_to_yeast17:YeastZAN_F4   0.476          0.238          0.238        
    ## distance_to_yeast25:YeastZAN_F4   0.238          0.476          0.238        
    ## distance_to_yeast32:YeastZAN_F4   0.238          0.238          0.476        
    ## distance_to_yeast41:YeastZAN_F4   0.238          0.238          0.238        
    ## distance_to_yeast48:YeastZAN_F4   0.238          0.238          0.238        
    ## distance_to_yeast55:YeastZAN_F4   0.238          0.238          0.238        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                -0.023         -0.023         -0.023        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F65 d__48:YEMM_F65 d__55:YEMM_F65
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F65  0.519                                      
    ## distance_to_yeast55:YeastEMM_F65  0.519          0.519                       
    ## distance_to_yeast17:YeastEMM_F66  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F66  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F66  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F66  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F66  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F66  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F7   0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F7   0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F7   0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F7   0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F7   0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F7   0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F70  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F70  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F70  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F70  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F70  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F70  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastEMM_F89  0.238          0.238          0.238        
    ## distance_to_yeast25:YeastEMM_F89  0.238          0.238          0.238        
    ## distance_to_yeast32:YeastEMM_F89  0.238          0.238          0.238        
    ## distance_to_yeast41:YeastEMM_F89  0.476          0.238          0.238        
    ## distance_to_yeast48:YeastEMM_F89  0.238          0.476          0.238        
    ## distance_to_yeast55:YeastEMM_F89  0.238          0.238          0.476        
    ## distance_to_yeast17:YeastSP_F14   0.238          0.238          0.238        
    ## distance_to_yeast25:YeastSP_F14   0.238          0.238          0.238        
    ## distance_to_yeast32:YeastSP_F14   0.238          0.238          0.238        
    ## distance_to_yeast41:YeastSP_F14   0.476          0.238          0.238        
    ## distance_to_yeast48:YeastSP_F14   0.238          0.476          0.238        
    ## distance_to_yeast55:YeastSP_F14   0.238          0.238          0.476        
    ## distance_to_yeast17:YeastZAN_F3   0.238          0.238          0.238        
    ## distance_to_yeast25:YeastZAN_F3   0.238          0.238          0.238        
    ## distance_to_yeast32:YeastZAN_F3   0.238          0.238          0.238        
    ## distance_to_yeast41:YeastZAN_F3   0.476          0.238          0.238        
    ## distance_to_yeast48:YeastZAN_F3   0.238          0.476          0.238        
    ## distance_to_yeast55:YeastZAN_F3   0.238          0.238          0.476        
    ## distance_to_yeast17:YeastZAN_F4   0.238          0.238          0.238        
    ## distance_to_yeast25:YeastZAN_F4   0.238          0.238          0.238        
    ## distance_to_yeast32:YeastZAN_F4   0.238          0.238          0.238        
    ## distance_to_yeast41:YeastZAN_F4   0.476          0.238          0.238        
    ## distance_to_yeast48:YeastZAN_F4   0.238          0.476          0.238        
    ## distance_to_yeast55:YeastZAN_F4   0.238          0.238          0.476        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                -0.023         -0.023         -0.023        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F66 d__25:YEMM_F66 d__32:YEMM_F66
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F66  0.500                                      
    ## distance_to_yeast32:YeastEMM_F66  0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F66  0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F66  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F66  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F66 d__48:YEMM_F66 d__55:YEMM_F66
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F66  0.500                                      
    ## distance_to_yeast55:YeastEMM_F66  0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F7   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F7   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F7   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  ds__17:YEMM_F7 ds__25:YEMM_F7 ds__32:YEMM_F7
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F7   0.500                                      
    ## distance_to_yeast32:YeastEMM_F7   0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F7   0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F7   0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F7   0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  ds__41:YEMM_F7 ds__48:YEMM_F7 ds__55:YEMM_F7
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F7   0.500                                      
    ## distance_to_yeast55:YeastEMM_F7   0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F70  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F70  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F70  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F70 d__25:YEMM_F70 d__32:YEMM_F70
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastEMM_F70  0.500                                      
    ## distance_to_yeast32:YeastEMM_F70  0.500          0.500                       
    ## distance_to_yeast41:YeastEMM_F70  0.500          0.500          0.500        
    ## distance_to_yeast48:YeastEMM_F70  0.500          0.500          0.500        
    ## distance_to_yeast55:YeastEMM_F70  0.500          0.500          0.500        
    ## distance_to_yeast17:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast41:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast17:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast17:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.500        
    ## distance_to_yeast41:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.250        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__41:YEMM_F70 d__48:YEMM_F70 d__55:YEMM_F70
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast48:YeastEMM_F70  0.500                                      
    ## distance_to_yeast55:YeastEMM_F70  0.500          0.500                       
    ## distance_to_yeast17:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast25:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast32:YeastEMM_F89  0.250          0.250          0.250        
    ## distance_to_yeast41:YeastEMM_F89  0.500          0.250          0.250        
    ## distance_to_yeast48:YeastEMM_F89  0.250          0.500          0.250        
    ## distance_to_yeast55:YeastEMM_F89  0.250          0.250          0.500        
    ## distance_to_yeast17:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastSP_F14   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastSP_F14   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastSP_F14   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastSP_F14   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F3   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F3   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F3   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F3   0.250          0.250          0.500        
    ## distance_to_yeast17:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast25:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast32:YeastZAN_F4   0.250          0.250          0.250        
    ## distance_to_yeast41:YeastZAN_F4   0.500          0.250          0.250        
    ## distance_to_yeast48:YeastZAN_F4   0.250          0.500          0.250        
    ## distance_to_yeast55:YeastZAN_F4   0.250          0.250          0.500        
    ## YeastEMM_F3:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI6                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI6                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI6                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI6                  0.000          0.000          0.000        
    ## YeastEMM_F3:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F34:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F47:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F48:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F49:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F63:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F64:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F65:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F66:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F7:DAI8                  0.000          0.000          0.000        
    ## YeastEMM_F70:DAI8                 0.000          0.000          0.000        
    ## YeastEMM_F89:DAI8                 0.000          0.000          0.000        
    ## YeastSP_F14:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F3:DAI8                  0.000          0.000          0.000        
    ## YeastZAN_F4:DAI8                  0.000          0.000          0.000        
    ##                                  d__17:YEMM_F8 d__25:YEMM_F8 d__32:YEMM_F8
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast25:YeastEMM_F89  0.500                                   
    ## distance_to_yeast32:YeastEMM_F89  0.500         0.500                     
    ## distance_to_yeast41:YeastEMM_F89  0.500         0.500         0.500       
    ## distance_to_yeast48:YeastEMM_F89  0.500         0.500         0.500       
    ## distance_to_yeast55:YeastEMM_F89  0.500         0.500         0.500       
    ## distance_to_yeast17:YeastSP_F14   0.500         0.250         0.250       
    ## distance_to_yeast25:YeastSP_F14   0.250         0.500         0.250       
    ## distance_to_yeast32:YeastSP_F14   0.250         0.250         0.500       
    ## distance_to_yeast41:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast48:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast17:YeastZAN_F3   0.500         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F3   0.250         0.500         0.250       
    ## distance_to_yeast32:YeastZAN_F3   0.250         0.250         0.500       
    ## distance_to_yeast41:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast48:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast17:YeastZAN_F4   0.500         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F4   0.250         0.500         0.250       
    ## distance_to_yeast32:YeastZAN_F4   0.250         0.250         0.500       
    ## distance_to_yeast41:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast48:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast55:YeastZAN_F4   0.250         0.250         0.250       
    ## YeastEMM_F3:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI6                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F3:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI8                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI8                  0.000         0.000         0.000       
    ##                                  d__41:YEMM_F8 d__48:YEMM_F8 d__55:YEMM_F8
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast48:YeastEMM_F89  0.500                                   
    ## distance_to_yeast55:YeastEMM_F89  0.500         0.500                     
    ## distance_to_yeast17:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast32:YeastSP_F14   0.250         0.250         0.250       
    ## distance_to_yeast41:YeastSP_F14   0.500         0.250         0.250       
    ## distance_to_yeast48:YeastSP_F14   0.250         0.500         0.250       
    ## distance_to_yeast55:YeastSP_F14   0.250         0.250         0.500       
    ## distance_to_yeast17:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast32:YeastZAN_F3   0.250         0.250         0.250       
    ## distance_to_yeast41:YeastZAN_F3   0.500         0.250         0.250       
    ## distance_to_yeast48:YeastZAN_F3   0.250         0.500         0.250       
    ## distance_to_yeast55:YeastZAN_F3   0.250         0.250         0.500       
    ## distance_to_yeast17:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast32:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast41:YeastZAN_F4   0.500         0.250         0.250       
    ## distance_to_yeast48:YeastZAN_F4   0.250         0.500         0.250       
    ## distance_to_yeast55:YeastZAN_F4   0.250         0.250         0.500       
    ## YeastEMM_F3:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI6                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F3:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI8                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI8                  0.000         0.000         0.000       
    ##                                  d__17:YS d__25:YS d__32:YS d__41:YS d__48:YS
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
    ## DAI6                                                                         
    ## DAI8                                                                         
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
    ## distance_to_yeast25:YeastSP_F14   0.500                                      
    ## distance_to_yeast32:YeastSP_F14   0.500    0.500                             
    ## distance_to_yeast41:YeastSP_F14   0.500    0.500    0.500                    
    ## distance_to_yeast48:YeastSP_F14   0.500    0.500    0.500    0.500           
    ## distance_to_yeast55:YeastSP_F14   0.500    0.500    0.500    0.500    0.500  
    ## distance_to_yeast17:YeastZAN_F3   0.500    0.250    0.250    0.250    0.250  
    ## distance_to_yeast25:YeastZAN_F3   0.250    0.500    0.250    0.250    0.250  
    ## distance_to_yeast32:YeastZAN_F3   0.250    0.250    0.500    0.250    0.250  
    ## distance_to_yeast41:YeastZAN_F3   0.250    0.250    0.250    0.500    0.250  
    ## distance_to_yeast48:YeastZAN_F3   0.250    0.250    0.250    0.250    0.500  
    ## distance_to_yeast55:YeastZAN_F3   0.250    0.250    0.250    0.250    0.250  
    ## distance_to_yeast17:YeastZAN_F4   0.500    0.250    0.250    0.250    0.250  
    ## distance_to_yeast25:YeastZAN_F4   0.250    0.500    0.250    0.250    0.250  
    ## distance_to_yeast32:YeastZAN_F4   0.250    0.250    0.500    0.250    0.250  
    ## distance_to_yeast41:YeastZAN_F4   0.250    0.250    0.250    0.500    0.250  
    ## distance_to_yeast48:YeastZAN_F4   0.250    0.250    0.250    0.250    0.500  
    ## distance_to_yeast55:YeastZAN_F4   0.250    0.250    0.250    0.250    0.250  
    ## YeastEMM_F3:DAI6                  0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F34:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F47:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F48:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F49:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F63:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F64:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F65:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F66:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F7:DAI6                  0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F70:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F89:DAI6                 0.000    0.000    0.000    0.000    0.000  
    ## YeastSP_F14:DAI6                  0.000    0.000    0.000    0.000    0.000  
    ## YeastZAN_F3:DAI6                  0.000    0.000    0.000    0.000    0.000  
    ## YeastZAN_F4:DAI6                  0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F3:DAI8                  0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F34:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F47:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F48:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F49:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F63:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F64:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F65:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F66:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F7:DAI8                  0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F70:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastEMM_F89:DAI8                 0.000    0.000    0.000    0.000    0.000  
    ## YeastSP_F14:DAI8                  0.000    0.000    0.000    0.000    0.000  
    ## YeastZAN_F3:DAI8                  0.000    0.000    0.000    0.000    0.000  
    ## YeastZAN_F4:DAI8                  0.000    0.000    0.000    0.000    0.000  
    ##                                  d__55:YS d__17:YZAN_F3 d__25:YZAN_F3
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
    ## DAI6                                                                 
    ## DAI8                                                                 
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
    ## distance_to_yeast17:YeastZAN_F3   0.250                              
    ## distance_to_yeast25:YeastZAN_F3   0.250    0.500                     
    ## distance_to_yeast32:YeastZAN_F3   0.250    0.500         0.500       
    ## distance_to_yeast41:YeastZAN_F3   0.250    0.500         0.500       
    ## distance_to_yeast48:YeastZAN_F3   0.250    0.500         0.500       
    ## distance_to_yeast55:YeastZAN_F3   0.500    0.500         0.500       
    ## distance_to_yeast17:YeastZAN_F4   0.250    0.500         0.250       
    ## distance_to_yeast25:YeastZAN_F4   0.250    0.250         0.500       
    ## distance_to_yeast32:YeastZAN_F4   0.250    0.250         0.250       
    ## distance_to_yeast41:YeastZAN_F4   0.250    0.250         0.250       
    ## distance_to_yeast48:YeastZAN_F4   0.250    0.250         0.250       
    ## distance_to_yeast55:YeastZAN_F4   0.500    0.250         0.250       
    ## YeastEMM_F3:DAI6                  0.000    0.000         0.000       
    ## YeastEMM_F34:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F47:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F48:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F49:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F63:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F64:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F65:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F66:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F7:DAI6                  0.000    0.000         0.000       
    ## YeastEMM_F70:DAI6                 0.000    0.000         0.000       
    ## YeastEMM_F89:DAI6                 0.000    0.000         0.000       
    ## YeastSP_F14:DAI6                  0.000    0.000         0.000       
    ## YeastZAN_F3:DAI6                  0.000    0.000         0.000       
    ## YeastZAN_F4:DAI6                  0.000    0.000         0.000       
    ## YeastEMM_F3:DAI8                  0.000    0.000         0.000       
    ## YeastEMM_F34:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F47:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F48:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F49:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F63:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F64:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F65:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F66:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F7:DAI8                  0.000    0.000         0.000       
    ## YeastEMM_F70:DAI8                 0.000    0.000         0.000       
    ## YeastEMM_F89:DAI8                 0.000    0.000         0.000       
    ## YeastSP_F14:DAI8                  0.000    0.000         0.000       
    ## YeastZAN_F3:DAI8                  0.000    0.000         0.000       
    ## YeastZAN_F4:DAI8                  0.000    0.000         0.000       
    ##                                  d__32:YZAN_F3 d__41:YZAN_F3 d__48:YZAN_F3
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast41:YeastZAN_F3   0.500                                   
    ## distance_to_yeast48:YeastZAN_F3   0.500         0.500                     
    ## distance_to_yeast55:YeastZAN_F3   0.500         0.500         0.500       
    ## distance_to_yeast17:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast25:YeastZAN_F4   0.250         0.250         0.250       
    ## distance_to_yeast32:YeastZAN_F4   0.500         0.250         0.250       
    ## distance_to_yeast41:YeastZAN_F4   0.250         0.500         0.250       
    ## distance_to_yeast48:YeastZAN_F4   0.250         0.250         0.500       
    ## distance_to_yeast55:YeastZAN_F4   0.250         0.250         0.250       
    ## YeastEMM_F3:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI6                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F3:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI8                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI8                  0.000         0.000         0.000       
    ##                                  d__55:YZAN_F3 d__17:YZAN_F4 d__25:YZAN_F4
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast17:YeastZAN_F4   0.250                                   
    ## distance_to_yeast25:YeastZAN_F4   0.250         0.500                     
    ## distance_to_yeast32:YeastZAN_F4   0.250         0.500         0.500       
    ## distance_to_yeast41:YeastZAN_F4   0.250         0.500         0.500       
    ## distance_to_yeast48:YeastZAN_F4   0.250         0.500         0.500       
    ## distance_to_yeast55:YeastZAN_F4   0.500         0.500         0.500       
    ## YeastEMM_F3:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI6                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F3:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI8                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI8                  0.000         0.000         0.000       
    ##                                  d__32:YZAN_F4 d__41:YZAN_F4 d__48:YZAN_F4
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast41:YeastZAN_F4   0.500                                   
    ## distance_to_yeast48:YeastZAN_F4   0.500         0.500                     
    ## distance_to_yeast55:YeastZAN_F4   0.500         0.500         0.500       
    ## YeastEMM_F3:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI6                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI6                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI6                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI6                  0.000         0.000         0.000       
    ## YeastEMM_F3:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F34:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F47:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F48:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F49:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F63:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F64:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F65:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F66:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F7:DAI8                  0.000         0.000         0.000       
    ## YeastEMM_F70:DAI8                 0.000         0.000         0.000       
    ## YeastEMM_F89:DAI8                 0.000         0.000         0.000       
    ## YeastSP_F14:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F3:DAI8                  0.000         0.000         0.000       
    ## YeastZAN_F4:DAI8                  0.000         0.000         0.000       
    ##                                  d__55:YZAN_F4 YEMM_F3:DAI6 YEMM_F34:DAI6
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
    ## DAI6                                                                     
    ## DAI8                                                                     
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
    ## distance_to_yeast48:YeastZAN_F4                                          
    ## distance_to_yeast55:YeastZAN_F4                                          
    ## YeastEMM_F3:DAI6                  0.000                                  
    ## YeastEMM_F34:DAI6                 0.000         0.500                    
    ## YeastEMM_F47:DAI6                 0.000         0.500        0.500       
    ## YeastEMM_F48:DAI6                 0.000         0.500        0.500       
    ## YeastEMM_F49:DAI6                 0.000         0.500        0.500       
    ## YeastEMM_F63:DAI6                 0.000         0.500        0.500       
    ## YeastEMM_F64:DAI6                 0.000         0.500        0.500       
    ## YeastEMM_F65:DAI6                 0.000         0.497        0.497       
    ## YeastEMM_F66:DAI6                 0.000         0.500        0.500       
    ## YeastEMM_F7:DAI6                  0.000         0.500        0.500       
    ## YeastEMM_F70:DAI6                 0.000         0.500        0.500       
    ## YeastEMM_F89:DAI6                 0.000         0.500        0.500       
    ## YeastSP_F14:DAI6                  0.000         0.500        0.500       
    ## YeastZAN_F3:DAI6                  0.000         0.500        0.500       
    ## YeastZAN_F4:DAI6                  0.000         0.500        0.500       
    ## YeastEMM_F3:DAI8                  0.000         0.500        0.250       
    ## YeastEMM_F34:DAI8                 0.000         0.250        0.500       
    ## YeastEMM_F47:DAI8                 0.000         0.250        0.250       
    ## YeastEMM_F48:DAI8                 0.000         0.250        0.250       
    ## YeastEMM_F49:DAI8                 0.000         0.250        0.250       
    ## YeastEMM_F63:DAI8                 0.000         0.250        0.250       
    ## YeastEMM_F64:DAI8                 0.000         0.250        0.250       
    ## YeastEMM_F65:DAI8                 0.000         0.235        0.235       
    ## YeastEMM_F66:DAI8                 0.000         0.250        0.250       
    ## YeastEMM_F7:DAI8                  0.000         0.250        0.250       
    ## YeastEMM_F70:DAI8                 0.000         0.250        0.250       
    ## YeastEMM_F89:DAI8                 0.000         0.250        0.250       
    ## YeastSP_F14:DAI8                  0.000         0.250        0.250       
    ## YeastZAN_F3:DAI8                  0.000         0.250        0.250       
    ## YeastZAN_F4:DAI8                  0.000         0.250        0.250       
    ##                                  YEMM_F47:DAI6 YEMM_F48:DAI6 YEMM_F49:DAI6
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## YeastEMM_F3:DAI6                                                          
    ## YeastEMM_F34:DAI6                                                         
    ## YeastEMM_F47:DAI6                                                         
    ## YeastEMM_F48:DAI6                 0.500                                   
    ## YeastEMM_F49:DAI6                 0.500         0.500                     
    ## YeastEMM_F63:DAI6                 0.500         0.500         0.500       
    ## YeastEMM_F64:DAI6                 0.500         0.500         0.500       
    ## YeastEMM_F65:DAI6                 0.497         0.497         0.497       
    ## YeastEMM_F66:DAI6                 0.500         0.500         0.500       
    ## YeastEMM_F7:DAI6                  0.500         0.500         0.500       
    ## YeastEMM_F70:DAI6                 0.500         0.500         0.500       
    ## YeastEMM_F89:DAI6                 0.500         0.500         0.500       
    ## YeastSP_F14:DAI6                  0.500         0.500         0.500       
    ## YeastZAN_F3:DAI6                  0.500         0.500         0.500       
    ## YeastZAN_F4:DAI6                  0.500         0.500         0.500       
    ## YeastEMM_F3:DAI8                  0.250         0.250         0.250       
    ## YeastEMM_F34:DAI8                 0.250         0.250         0.250       
    ## YeastEMM_F47:DAI8                 0.500         0.250         0.250       
    ## YeastEMM_F48:DAI8                 0.250         0.500         0.250       
    ## YeastEMM_F49:DAI8                 0.250         0.250         0.500       
    ## YeastEMM_F63:DAI8                 0.250         0.250         0.250       
    ## YeastEMM_F64:DAI8                 0.250         0.250         0.250       
    ## YeastEMM_F65:DAI8                 0.235         0.235         0.235       
    ## YeastEMM_F66:DAI8                 0.250         0.250         0.250       
    ## YeastEMM_F7:DAI8                  0.250         0.250         0.250       
    ## YeastEMM_F70:DAI8                 0.250         0.250         0.250       
    ## YeastEMM_F89:DAI8                 0.250         0.250         0.250       
    ## YeastSP_F14:DAI8                  0.250         0.250         0.250       
    ## YeastZAN_F3:DAI8                  0.250         0.250         0.250       
    ## YeastZAN_F4:DAI8                  0.250         0.250         0.250       
    ##                                  YEMM_F63:DAI6 YEMM_F64:DAI6 YEMM_F65:DAI6
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## YeastEMM_F3:DAI6                                                          
    ## YeastEMM_F34:DAI6                                                         
    ## YeastEMM_F47:DAI6                                                         
    ## YeastEMM_F48:DAI6                                                         
    ## YeastEMM_F49:DAI6                                                         
    ## YeastEMM_F63:DAI6                                                         
    ## YeastEMM_F64:DAI6                 0.500                                   
    ## YeastEMM_F65:DAI6                 0.497         0.497                     
    ## YeastEMM_F66:DAI6                 0.500         0.500         0.497       
    ## YeastEMM_F7:DAI6                  0.500         0.500         0.497       
    ## YeastEMM_F70:DAI6                 0.500         0.500         0.497       
    ## YeastEMM_F89:DAI6                 0.500         0.500         0.497       
    ## YeastSP_F14:DAI6                  0.500         0.500         0.497       
    ## YeastZAN_F3:DAI6                  0.500         0.500         0.497       
    ## YeastZAN_F4:DAI6                  0.500         0.500         0.497       
    ## YeastEMM_F3:DAI8                  0.250         0.250         0.248       
    ## YeastEMM_F34:DAI8                 0.250         0.250         0.248       
    ## YeastEMM_F47:DAI8                 0.250         0.250         0.248       
    ## YeastEMM_F48:DAI8                 0.250         0.250         0.248       
    ## YeastEMM_F49:DAI8                 0.250         0.250         0.248       
    ## YeastEMM_F63:DAI8                 0.500         0.250         0.248       
    ## YeastEMM_F64:DAI8                 0.250         0.500         0.248       
    ## YeastEMM_F65:DAI8                 0.235         0.235         0.468       
    ## YeastEMM_F66:DAI8                 0.250         0.250         0.248       
    ## YeastEMM_F7:DAI8                  0.250         0.250         0.248       
    ## YeastEMM_F70:DAI8                 0.250         0.250         0.248       
    ## YeastEMM_F89:DAI8                 0.250         0.250         0.248       
    ## YeastSP_F14:DAI8                  0.250         0.250         0.248       
    ## YeastZAN_F3:DAI8                  0.250         0.250         0.248       
    ## YeastZAN_F4:DAI8                  0.250         0.250         0.248       
    ##                                  YEMM_F66:DAI6 YEMM_F7:DAI6 YEMM_F70:DAI6
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
    ## DAI6                                                                     
    ## DAI8                                                                     
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
    ## distance_to_yeast48:YeastZAN_F4                                          
    ## distance_to_yeast55:YeastZAN_F4                                          
    ## YeastEMM_F3:DAI6                                                         
    ## YeastEMM_F34:DAI6                                                        
    ## YeastEMM_F47:DAI6                                                        
    ## YeastEMM_F48:DAI6                                                        
    ## YeastEMM_F49:DAI6                                                        
    ## YeastEMM_F63:DAI6                                                        
    ## YeastEMM_F64:DAI6                                                        
    ## YeastEMM_F65:DAI6                                                        
    ## YeastEMM_F66:DAI6                                                        
    ## YeastEMM_F7:DAI6                  0.500                                  
    ## YeastEMM_F70:DAI6                 0.500         0.500                    
    ## YeastEMM_F89:DAI6                 0.500         0.500        0.500       
    ## YeastSP_F14:DAI6                  0.500         0.500        0.500       
    ## YeastZAN_F3:DAI6                  0.500         0.500        0.500       
    ## YeastZAN_F4:DAI6                  0.500         0.500        0.500       
    ## YeastEMM_F3:DAI8                  0.250         0.250        0.250       
    ## YeastEMM_F34:DAI8                 0.250         0.250        0.250       
    ## YeastEMM_F47:DAI8                 0.250         0.250        0.250       
    ## YeastEMM_F48:DAI8                 0.250         0.250        0.250       
    ## YeastEMM_F49:DAI8                 0.250         0.250        0.250       
    ## YeastEMM_F63:DAI8                 0.250         0.250        0.250       
    ## YeastEMM_F64:DAI8                 0.250         0.250        0.250       
    ## YeastEMM_F65:DAI8                 0.235         0.235        0.235       
    ## YeastEMM_F66:DAI8                 0.500         0.250        0.250       
    ## YeastEMM_F7:DAI8                  0.250         0.500        0.250       
    ## YeastEMM_F70:DAI8                 0.250         0.250        0.500       
    ## YeastEMM_F89:DAI8                 0.250         0.250        0.250       
    ## YeastSP_F14:DAI8                  0.250         0.250        0.250       
    ## YeastZAN_F3:DAI8                  0.250         0.250        0.250       
    ## YeastZAN_F4:DAI8                  0.250         0.250        0.250       
    ##                                  YEMM_F89:DAI6 YSP_F14:DAI6 YZAN_F3:DAI6
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
    ## DAI6                                                                    
    ## DAI8                                                                    
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
    ## distance_to_yeast48:YeastZAN_F4                                         
    ## distance_to_yeast55:YeastZAN_F4                                         
    ## YeastEMM_F3:DAI6                                                        
    ## YeastEMM_F34:DAI6                                                       
    ## YeastEMM_F47:DAI6                                                       
    ## YeastEMM_F48:DAI6                                                       
    ## YeastEMM_F49:DAI6                                                       
    ## YeastEMM_F63:DAI6                                                       
    ## YeastEMM_F64:DAI6                                                       
    ## YeastEMM_F65:DAI6                                                       
    ## YeastEMM_F66:DAI6                                                       
    ## YeastEMM_F7:DAI6                                                        
    ## YeastEMM_F70:DAI6                                                       
    ## YeastEMM_F89:DAI6                                                       
    ## YeastSP_F14:DAI6                  0.500                                 
    ## YeastZAN_F3:DAI6                  0.500         0.500                   
    ## YeastZAN_F4:DAI6                  0.500         0.500        0.500      
    ## YeastEMM_F3:DAI8                  0.250         0.250        0.250      
    ## YeastEMM_F34:DAI8                 0.250         0.250        0.250      
    ## YeastEMM_F47:DAI8                 0.250         0.250        0.250      
    ## YeastEMM_F48:DAI8                 0.250         0.250        0.250      
    ## YeastEMM_F49:DAI8                 0.250         0.250        0.250      
    ## YeastEMM_F63:DAI8                 0.250         0.250        0.250      
    ## YeastEMM_F64:DAI8                 0.250         0.250        0.250      
    ## YeastEMM_F65:DAI8                 0.235         0.235        0.235      
    ## YeastEMM_F66:DAI8                 0.250         0.250        0.250      
    ## YeastEMM_F7:DAI8                  0.250         0.250        0.250      
    ## YeastEMM_F70:DAI8                 0.250         0.250        0.250      
    ## YeastEMM_F89:DAI8                 0.500         0.250        0.250      
    ## YeastSP_F14:DAI8                  0.250         0.500        0.250      
    ## YeastZAN_F3:DAI8                  0.250         0.250        0.500      
    ## YeastZAN_F4:DAI8                  0.250         0.250        0.250      
    ##                                  YZAN_F4:DAI6 YEMM_F3:DAI8 YEMM_F34:DAI8
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
    ## DAI6                                                                    
    ## DAI8                                                                    
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
    ## distance_to_yeast48:YeastZAN_F4                                         
    ## distance_to_yeast55:YeastZAN_F4                                         
    ## YeastEMM_F3:DAI6                                                        
    ## YeastEMM_F34:DAI6                                                       
    ## YeastEMM_F47:DAI6                                                       
    ## YeastEMM_F48:DAI6                                                       
    ## YeastEMM_F49:DAI6                                                       
    ## YeastEMM_F63:DAI6                                                       
    ## YeastEMM_F64:DAI6                                                       
    ## YeastEMM_F65:DAI6                                                       
    ## YeastEMM_F66:DAI6                                                       
    ## YeastEMM_F7:DAI6                                                        
    ## YeastEMM_F70:DAI6                                                       
    ## YeastEMM_F89:DAI6                                                       
    ## YeastSP_F14:DAI6                                                        
    ## YeastZAN_F3:DAI6                                                        
    ## YeastZAN_F4:DAI6                                                        
    ## YeastEMM_F3:DAI8                  0.250                                 
    ## YeastEMM_F34:DAI8                 0.250        0.500                    
    ## YeastEMM_F47:DAI8                 0.250        0.500        0.500       
    ## YeastEMM_F48:DAI8                 0.250        0.500        0.500       
    ## YeastEMM_F49:DAI8                 0.250        0.500        0.500       
    ## YeastEMM_F63:DAI8                 0.250        0.500        0.500       
    ## YeastEMM_F64:DAI8                 0.250        0.500        0.500       
    ## YeastEMM_F65:DAI8                 0.235        0.471        0.471       
    ## YeastEMM_F66:DAI8                 0.250        0.500        0.500       
    ## YeastEMM_F7:DAI8                  0.250        0.500        0.500       
    ## YeastEMM_F70:DAI8                 0.250        0.500        0.500       
    ## YeastEMM_F89:DAI8                 0.250        0.500        0.500       
    ## YeastSP_F14:DAI8                  0.250        0.500        0.500       
    ## YeastZAN_F3:DAI8                  0.250        0.500        0.500       
    ## YeastZAN_F4:DAI8                  0.500        0.500        0.500       
    ##                                  YEMM_F47:DAI8 YEMM_F48:DAI8 YEMM_F49:DAI8
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## YeastEMM_F3:DAI6                                                          
    ## YeastEMM_F34:DAI6                                                         
    ## YeastEMM_F47:DAI6                                                         
    ## YeastEMM_F48:DAI6                                                         
    ## YeastEMM_F49:DAI6                                                         
    ## YeastEMM_F63:DAI6                                                         
    ## YeastEMM_F64:DAI6                                                         
    ## YeastEMM_F65:DAI6                                                         
    ## YeastEMM_F66:DAI6                                                         
    ## YeastEMM_F7:DAI6                                                          
    ## YeastEMM_F70:DAI6                                                         
    ## YeastEMM_F89:DAI6                                                         
    ## YeastSP_F14:DAI6                                                          
    ## YeastZAN_F3:DAI6                                                          
    ## YeastZAN_F4:DAI6                                                          
    ## YeastEMM_F3:DAI8                                                          
    ## YeastEMM_F34:DAI8                                                         
    ## YeastEMM_F47:DAI8                                                         
    ## YeastEMM_F48:DAI8                 0.500                                   
    ## YeastEMM_F49:DAI8                 0.500         0.500                     
    ## YeastEMM_F63:DAI8                 0.500         0.500         0.500       
    ## YeastEMM_F64:DAI8                 0.500         0.500         0.500       
    ## YeastEMM_F65:DAI8                 0.471         0.471         0.471       
    ## YeastEMM_F66:DAI8                 0.500         0.500         0.500       
    ## YeastEMM_F7:DAI8                  0.500         0.500         0.500       
    ## YeastEMM_F70:DAI8                 0.500         0.500         0.500       
    ## YeastEMM_F89:DAI8                 0.500         0.500         0.500       
    ## YeastSP_F14:DAI8                  0.500         0.500         0.500       
    ## YeastZAN_F3:DAI8                  0.500         0.500         0.500       
    ## YeastZAN_F4:DAI8                  0.500         0.500         0.500       
    ##                                  YEMM_F63:DAI8 YEMM_F64:DAI8 YEMM_F65:DAI8
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
    ## DAI6                                                                      
    ## DAI8                                                                      
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
    ## distance_to_yeast48:YeastZAN_F4                                           
    ## distance_to_yeast55:YeastZAN_F4                                           
    ## YeastEMM_F3:DAI6                                                          
    ## YeastEMM_F34:DAI6                                                         
    ## YeastEMM_F47:DAI6                                                         
    ## YeastEMM_F48:DAI6                                                         
    ## YeastEMM_F49:DAI6                                                         
    ## YeastEMM_F63:DAI6                                                         
    ## YeastEMM_F64:DAI6                                                         
    ## YeastEMM_F65:DAI6                                                         
    ## YeastEMM_F66:DAI6                                                         
    ## YeastEMM_F7:DAI6                                                          
    ## YeastEMM_F70:DAI6                                                         
    ## YeastEMM_F89:DAI6                                                         
    ## YeastSP_F14:DAI6                                                          
    ## YeastZAN_F3:DAI6                                                          
    ## YeastZAN_F4:DAI6                                                          
    ## YeastEMM_F3:DAI8                                                          
    ## YeastEMM_F34:DAI8                                                         
    ## YeastEMM_F47:DAI8                                                         
    ## YeastEMM_F48:DAI8                                                         
    ## YeastEMM_F49:DAI8                                                         
    ## YeastEMM_F63:DAI8                                                         
    ## YeastEMM_F64:DAI8                 0.500                                   
    ## YeastEMM_F65:DAI8                 0.471         0.471                     
    ## YeastEMM_F66:DAI8                 0.500         0.500         0.471       
    ## YeastEMM_F7:DAI8                  0.500         0.500         0.471       
    ## YeastEMM_F70:DAI8                 0.500         0.500         0.471       
    ## YeastEMM_F89:DAI8                 0.500         0.500         0.471       
    ## YeastSP_F14:DAI8                  0.500         0.500         0.471       
    ## YeastZAN_F3:DAI8                  0.500         0.500         0.471       
    ## YeastZAN_F4:DAI8                  0.500         0.500         0.471       
    ##                                  YEMM_F66:DAI8 YEMM_F7:DAI8 YEMM_F70:DAI8
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
    ## DAI6                                                                     
    ## DAI8                                                                     
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
    ## distance_to_yeast48:YeastZAN_F4                                          
    ## distance_to_yeast55:YeastZAN_F4                                          
    ## YeastEMM_F3:DAI6                                                         
    ## YeastEMM_F34:DAI6                                                        
    ## YeastEMM_F47:DAI6                                                        
    ## YeastEMM_F48:DAI6                                                        
    ## YeastEMM_F49:DAI6                                                        
    ## YeastEMM_F63:DAI6                                                        
    ## YeastEMM_F64:DAI6                                                        
    ## YeastEMM_F65:DAI6                                                        
    ## YeastEMM_F66:DAI6                                                        
    ## YeastEMM_F7:DAI6                                                         
    ## YeastEMM_F70:DAI6                                                        
    ## YeastEMM_F89:DAI6                                                        
    ## YeastSP_F14:DAI6                                                         
    ## YeastZAN_F3:DAI6                                                         
    ## YeastZAN_F4:DAI6                                                         
    ## YeastEMM_F3:DAI8                                                         
    ## YeastEMM_F34:DAI8                                                        
    ## YeastEMM_F47:DAI8                                                        
    ## YeastEMM_F48:DAI8                                                        
    ## YeastEMM_F49:DAI8                                                        
    ## YeastEMM_F63:DAI8                                                        
    ## YeastEMM_F64:DAI8                                                        
    ## YeastEMM_F65:DAI8                                                        
    ## YeastEMM_F66:DAI8                                                        
    ## YeastEMM_F7:DAI8                  0.500                                  
    ## YeastEMM_F70:DAI8                 0.500         0.500                    
    ## YeastEMM_F89:DAI8                 0.500         0.500        0.500       
    ## YeastSP_F14:DAI8                  0.500         0.500        0.500       
    ## YeastZAN_F3:DAI8                  0.500         0.500        0.500       
    ## YeastZAN_F4:DAI8                  0.500         0.500        0.500       
    ##                                  YEMM_F89:DAI8 YSP_F14:DAI8 YZAN_F3:DAI8
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
    ## DAI6                                                                    
    ## DAI8                                                                    
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
    ## distance_to_yeast48:YeastZAN_F4                                         
    ## distance_to_yeast55:YeastZAN_F4                                         
    ## YeastEMM_F3:DAI6                                                        
    ## YeastEMM_F34:DAI6                                                       
    ## YeastEMM_F47:DAI6                                                       
    ## YeastEMM_F48:DAI6                                                       
    ## YeastEMM_F49:DAI6                                                       
    ## YeastEMM_F63:DAI6                                                       
    ## YeastEMM_F64:DAI6                                                       
    ## YeastEMM_F65:DAI6                                                       
    ## YeastEMM_F66:DAI6                                                       
    ## YeastEMM_F7:DAI6                                                        
    ## YeastEMM_F70:DAI6                                                       
    ## YeastEMM_F89:DAI6                                                       
    ## YeastSP_F14:DAI6                                                        
    ## YeastZAN_F3:DAI6                                                        
    ## YeastZAN_F4:DAI6                                                        
    ## YeastEMM_F3:DAI8                                                        
    ## YeastEMM_F34:DAI8                                                       
    ## YeastEMM_F47:DAI8                                                       
    ## YeastEMM_F48:DAI8                                                       
    ## YeastEMM_F49:DAI8                                                       
    ## YeastEMM_F63:DAI8                                                       
    ## YeastEMM_F64:DAI8                                                       
    ## YeastEMM_F65:DAI8                                                       
    ## YeastEMM_F66:DAI8                                                       
    ## YeastEMM_F7:DAI8                                                        
    ## YeastEMM_F70:DAI8                                                       
    ## YeastEMM_F89:DAI8                                                       
    ## YeastSP_F14:DAI8                  0.500                                 
    ## YeastZAN_F3:DAI8                  0.500         0.500                   
    ## YeastZAN_F4:DAI8                  0.500         0.500        0.500      
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -14.463917335  -0.410039417  -0.006705415   0.383877940   4.218861262 
    ## 
    ## Number of Observations: 1000
    ## Number of Groups: 3

### Loop for running analysis for each day separately EMM_B5

Days after inoculation (DAI) as a factor is always significantly
impacting the growth. Also, our plot will represent data for Day 6 thus
we want relevant stats and comparison on Day 6 to present in the plot.
So, loop was made for each Day data and removing DAI from the model and
keeping rest of it present.

``` r
# Unique days
day <- unique(B5.no.1st.data$DAI)
# Initialize empty lists to store results
B5_Day <- list()
meansB5_Day <- list()
resultsB5_Day <- list()
model_summary_B5 <- list()
anova_table_B5 <- list()
Results_lsmeansB5 <- list()
filtered_comparisonsB5 <- list()

for (i in seq_along(day)) {
    # Filter data for the current day
  B5_Day[[i]] <- filter(B5.no.contact, DAI == day[i])
    # Fit the mixed effects model
  resultsB5_Day[[i]] <- lme(increase ~ distance_to_yeast*Yeast,
                            data = B5_Day[[i]],
                            random = ~1 | Replication,
                            na.action = na.omit)
    model_summary_B5[[i]] <- summary(resultsB5_Day[[i]])
  anova_table_B5[[i]] <- anova(resultsB5_Day[[i]])
   # Estimated marginal means
  meansB5_Day[[i]] <- emmeans(resultsB5_Day[[i]], ~ Yeast|distance_to_yeast)
    # Compact letter display with Tukey adjustment
  Results_lsmeansB5[[i]] <- multcomp::cld(meansB5_Day[[i]], alpha = 0.05,
                                          Letters = letters, reversed = TRUE,
                                          details = TRUE)
   # Print lsmeans results
 # print(Results_lsmeansB5[[i]])
  
  # Filter comparisons involving "Control"
  filtered_comparisonsB5[[i]] <- Results_lsmeansB5[[i]]$comparisons %>%
    filter(grepl("Control", contrast))
   # Print filtered contrasts
  print(filtered_comparisonsB5[[i]])
}
```

    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.80667 0.258 222   3.130  0.1286
    ##  Control - EMM_F3   0.03333 0.258 222   0.129  1.0000
    ##  ZAN_F3 - Control   0.06000 0.258 222   0.233  1.0000
    ##  ZAN_F4 - Control   0.06333 0.258 222   0.246  1.0000
    ##  EMM_F70 - Control  0.09000 0.258 222   0.349  1.0000
    ##  EMM_F63 - Control  0.09333 0.258 222   0.362  1.0000
    ##  EMM_F34 - Control  0.15333 0.258 222   0.595  1.0000
    ##  EMM_F64 - Control  0.16667 0.258 222   0.647  1.0000
    ##  EMM_F49 - Control  0.17333 0.258 222   0.673  1.0000
    ##  EMM_F65 - Control  0.18333 0.258 222   0.711  1.0000
    ##  EMM_F89 - Control  0.21333 0.258 222   0.828  1.0000
    ##  EMM_F47 - Control  0.23333 0.258 222   0.905  0.9999
    ##  EMM_F66 - Control  0.32000 0.258 222   1.242  0.9973
    ##  EMM_F7 - Control   0.39333 0.258 222   1.526  0.9783
    ##  EMM_F48 - Control  0.51000 0.258 222   1.979  0.8341
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F63  0.29000 0.258 222   1.125  0.9991
    ##  Control - SP_F14   0.12667 0.258 222   0.492  1.0000
    ##  Control - EMM_F49  0.01333 0.258 222   0.052  1.0000
    ##  EMM_F64 - Control  0.12333 0.258 222   0.479  1.0000
    ##  EMM_F89 - Control  0.14667 0.258 222   0.569  1.0000
    ##  ZAN_F3 - Control   0.19667 0.258 222   0.763  1.0000
    ##  EMM_F47 - Control  0.24667 0.258 222   0.957  0.9999
    ##  ZAN_F4 - Control   0.26333 0.258 222   1.022  0.9997
    ##  EMM_F7 - Control   0.27000 0.258 222   1.048  0.9996
    ##  EMM_F70 - Control  0.29667 0.258 222   1.151  0.9988
    ##  EMM_F34 - Control  0.37667 0.258 222   1.462  0.9856
    ##  EMM_F66 - Control  0.38667 0.258 222   1.500  0.9815
    ##  EMM_F3 - Control   0.43000 0.258 222   1.669  0.9527
    ##  EMM_F65 - Control  0.59333 0.258 222   2.302  0.6223
    ##  EMM_F48 - Control  0.68333 0.258 222   2.652  0.3689
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.16667 0.258 222   0.647  1.0000
    ##  Control - ZAN_F4   0.07000 0.258 222   0.272  1.0000
    ##  Control - EMM_F63  0.05333 0.258 222   0.207  1.0000
    ##  Control - EMM_F34  0.01667 0.258 222   0.065  1.0000
    ##  Control - EMM_F89  0.00333 0.258 222   0.013  1.0000
    ##  EMM_F49 - Control  0.00667 0.258 222   0.026  1.0000
    ##  EMM_F66 - Control  0.01333 0.258 222   0.052  1.0000
    ##  EMM_F65 - Control  0.05333 0.258 222   0.207  1.0000
    ##  EMM_F70 - Control  0.18667 0.258 222   0.724  1.0000
    ##  EMM_F47 - Control  0.19667 0.258 222   0.763  1.0000
    ##  EMM_F7 - Control   0.26667 0.258 222   1.035  0.9997
    ##  EMM_F64 - Control  0.28000 0.258 222   1.087  0.9994
    ##  EMM_F3 - Control   0.28000 0.258 222   1.087  0.9994
    ##  ZAN_F3 - Control   0.28667 0.258 222   1.112  0.9992
    ##  EMM_F48 - Control  0.48333 0.258 222   1.876  0.8842
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.33000 0.258 222   1.281  0.9962
    ##  Control - EMM_F66  0.06667 0.258 222   0.259  1.0000
    ##  Control - EMM_F63  0.06000 0.258 222   0.233  1.0000
    ##  Control - EMM_F49  0.00000 0.258 222   0.000  1.0000
    ##  EMM_F70 - Control  0.04667 0.258 222   0.181  1.0000
    ##  EMM_F89 - Control  0.05667 0.258 222   0.220  1.0000
    ##  EMM_F3 - Control   0.07333 0.258 222   0.285  1.0000
    ##  EMM_F34 - Control  0.08000 0.258 222   0.310  1.0000
    ##  EMM_F7 - Control   0.14333 0.258 222   0.556  1.0000
    ##  EMM_F47 - Control  0.15000 0.258 222   0.582  1.0000
    ##  EMM_F65 - Control  0.18333 0.258 222   0.711  1.0000
    ##  ZAN_F4 - Control   0.19000 0.258 222   0.737  1.0000
    ##  EMM_F64 - Control  0.25333 0.258 222   0.983  0.9998
    ##  ZAN_F3 - Control   0.34333 0.258 222   1.332  0.9943
    ##  EMM_F48 - Control  0.42000 0.258 222   1.630  0.9612
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F63  0.07333 0.258 222   0.285  1.0000
    ##  EMM_F34 - Control  0.02333 0.258 222   0.091  1.0000
    ##  EMM_F89 - Control  0.09333 0.258 222   0.362  1.0000
    ##  ZAN_F3 - Control   0.10333 0.258 222   0.401  1.0000
    ##  ZAN_F4 - Control   0.12000 0.258 222   0.466  1.0000
    ##  EMM_F47 - Control  0.21667 0.258 222   0.841  1.0000
    ##  SP_F14 - Control   0.22000 0.258 222   0.854  1.0000
    ##  EMM_F64 - Control  0.30000 0.258 222   1.164  0.9987
    ##  EMM_F7 - Control   0.31000 0.258 222   1.203  0.9981
    ##  EMM_F48 - Control  0.31667 0.258 222   1.229  0.9976
    ##  EMM_F49 - Control  0.32000 0.258 222   1.242  0.9973
    ##  EMM_F3 - Control   0.34000 0.258 222   1.319  0.9948
    ##  EMM_F70 - Control  0.44667 0.258 222   1.733  0.9357
    ##  EMM_F65 - Control  0.56333 0.258 222   2.186  0.7059
    ##  EMM_F66 - Control  0.59333 0.258 222   2.302  0.6223
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F70  0.21667 0.258 222   0.841  1.0000
    ##  Control - EMM_F89  0.17333 0.258 222   0.673  1.0000
    ##  Control - EMM_F34  0.14667 0.258 222   0.569  1.0000
    ##  Control - EMM_F3   0.10333 0.258 222   0.401  1.0000
    ##  Control - EMM_F66  0.08667 0.258 222   0.336  1.0000
    ##  Control - SP_F14   0.05333 0.258 222   0.207  1.0000
    ##  Control - EMM_F47  0.03667 0.258 222   0.142  1.0000
    ##  Control - EMM_F49  0.03000 0.258 222   0.116  1.0000
    ##  Control - EMM_F65  0.02333 0.258 222   0.091  1.0000
    ##  EMM_F7 - Control   0.00333 0.258 222   0.013  1.0000
    ##  EMM_F64 - Control  0.00667 0.258 222   0.026  1.0000
    ##  ZAN_F3 - Control   0.03667 0.258 222   0.142  1.0000
    ##  EMM_F63 - Control  0.11333 0.258 222   0.440  1.0000
    ##  ZAN_F4 - Control   0.13000 0.258 222   0.504  1.0000
    ##  EMM_F48 - Control  0.13333 0.258 222   0.517  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F70 - Control  0.02667 0.258 222   0.103  1.0000
    ##  EMM_F63 - Control  0.14667 0.258 222   0.569  1.0000
    ##  EMM_F34 - Control  0.17667 0.258 222   0.686  1.0000
    ##  EMM_F47 - Control  0.19000 0.258 222   0.737  1.0000
    ##  SP_F14 - Control   0.21667 0.258 222   0.841  1.0000
    ##  ZAN_F4 - Control   0.23000 0.258 222   0.893  0.9999
    ##  EMM_F7 - Control   0.25667 0.258 222   0.996  0.9998
    ##  EMM_F3 - Control   0.27667 0.258 222   1.074  0.9995
    ##  EMM_F64 - Control  0.28000 0.258 222   1.087  0.9994
    ##  EMM_F66 - Control  0.29667 0.258 222   1.151  0.9988
    ##  EMM_F49 - Control  0.33333 0.258 222   1.294  0.9958
    ##  EMM_F48 - Control  0.45667 0.258 222   1.772  0.9236
    ##  ZAN_F3 - Control   0.50333 0.258 222   1.953  0.8476
    ##  EMM_F89 - Control  0.58000 0.258 222   2.251  0.6601
    ##  EMM_F65 - Control  0.66667 0.258 222   2.587  0.4131
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   4.09333 0.322 221  12.700  <.0001
    ##  Control - EMM_F63  0.16000 0.322 221   0.496  1.0000
    ##  Control - EMM_F34  0.06667 0.322 221   0.207  1.0000
    ##  Control - EMM_F49  0.05000 0.322 221   0.155  1.0000
    ##  Control - ZAN_F4   0.03667 0.322 221   0.114  1.0000
    ##  EMM_F64 - Control  0.01333 0.322 221   0.041  1.0000
    ##  EMM_F89 - Control  0.02000 0.322 221   0.062  1.0000
    ##  EMM_F70 - Control  0.02333 0.322 221   0.072  1.0000
    ##  EMM_F3 - Control   0.09000 0.322 221   0.279  1.0000
    ##  EMM_F47 - Control  0.15333 0.322 221   0.476  1.0000
    ##  EMM_F65 - Control  0.18997 0.361 221   0.527  1.0000
    ##  EMM_F66 - Control  0.23000 0.322 221   0.714  1.0000
    ##  ZAN_F3 - Control   0.34000 0.322 221   1.055  0.9996
    ##  EMM_F48 - Control  0.40667 0.322 221   1.262  0.9968
    ##  EMM_F7 - Control   0.46333 0.322 221   1.437  0.9877
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.33000 0.322 221   1.024  0.9997
    ##  Control - EMM_F63  0.06667 0.322 221   0.207  1.0000
    ##  EMM_F65 - Control  0.06000 0.322 221   0.186  1.0000
    ##  EMM_F64 - Control  0.20667 0.322 221   0.641  1.0000
    ##  EMM_F49 - Control  0.21667 0.322 221   0.672  1.0000
    ##  ZAN_F4 - Control   0.28000 0.322 221   0.869  1.0000
    ##  EMM_F89 - Control  0.29000 0.322 221   0.900  0.9999
    ##  EMM_F70 - Control  0.40667 0.322 221   1.262  0.9968
    ##  EMM_F66 - Control  0.44000 0.322 221   1.365  0.9927
    ##  EMM_F47 - Control  0.51333 0.322 221   1.593  0.9682
    ##  EMM_F3 - Control   0.54333 0.322 221   1.686  0.9486
    ##  EMM_F7 - Control   0.57333 0.322 221   1.779  0.9214
    ##  ZAN_F3 - Control   0.57667 0.322 221   1.789  0.9178
    ##  EMM_F34 - Control  0.79000 0.322 221   2.451  0.5117
    ##  EMM_F48 - Control  0.94000 0.322 221   2.916  0.2156
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F89  0.06000 0.322 221   0.186  1.0000
    ##  Control - EMM_F63  0.02000 0.322 221   0.062  1.0000
    ##  Control - EMM_F49  0.01333 0.322 221   0.041  1.0000
    ##  ZAN_F4 - Control   0.00333 0.322 221   0.010  1.0000
    ##  EMM_F66 - Control  0.03667 0.322 221   0.114  1.0000
    ##  EMM_F64 - Control  0.09667 0.322 221   0.300  1.0000
    ##  SP_F14 - Control   0.15333 0.322 221   0.476  1.0000
    ##  EMM_F34 - Control  0.22333 0.322 221   0.693  1.0000
    ##  EMM_F65 - Control  0.22667 0.322 221   0.703  1.0000
    ##  EMM_F70 - Control  0.25667 0.322 221   0.796  1.0000
    ##  EMM_F47 - Control  0.29333 0.322 221   0.910  0.9999
    ##  EMM_F7 - Control   0.33667 0.322 221   1.045  0.9996
    ##  EMM_F3 - Control   0.35000 0.322 221   1.086  0.9994
    ##  EMM_F48 - Control  0.42333 0.322 221   1.313  0.9951
    ##  ZAN_F3 - Control   0.58667 0.322 221   1.820  0.9066
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F66  0.24333 0.322 221   0.755  1.0000
    ##  Control - EMM_F70  0.19333 0.322 221   0.600  1.0000
    ##  Control - SP_F14   0.18667 0.322 221   0.579  1.0000
    ##  Control - EMM_F49  0.14000 0.322 221   0.434  1.0000
    ##  Control - EMM_F89  0.11667 0.322 221   0.362  1.0000
    ##  Control - EMM_F64  0.06667 0.322 221   0.207  1.0000
    ##  Control - ZAN_F4   0.04333 0.322 221   0.134  1.0000
    ##  Control - EMM_F34  0.01667 0.322 221   0.052  1.0000
    ##  EMM_F65 - Control  0.03333 0.322 221   0.103  1.0000
    ##  EMM_F7 - Control   0.04667 0.322 221   0.145  1.0000
    ##  EMM_F3 - Control   0.10000 0.322 221   0.310  1.0000
    ##  EMM_F63 - Control  0.11000 0.322 221   0.341  1.0000
    ##  EMM_F47 - Control  0.13333 0.322 221   0.414  1.0000
    ##  EMM_F48 - Control  0.34667 0.322 221   1.076  0.9995
    ##  ZAN_F3 - Control   0.45000 0.322 221   1.396  0.9908
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34  0.25667 0.322 221   0.796  1.0000
    ##  Control - EMM_F63  0.01000 0.322 221   0.031  1.0000
    ##  EMM_F89 - Control  0.01667 0.322 221   0.052  1.0000
    ##  ZAN_F3 - Control   0.08000 0.322 221   0.248  1.0000
    ##  ZAN_F4 - Control   0.12333 0.322 221   0.383  1.0000
    ##  EMM_F48 - Control  0.12667 0.322 221   0.393  1.0000
    ##  EMM_F49 - Control  0.16333 0.322 221   0.507  1.0000
    ##  EMM_F47 - Control  0.22000 0.322 221   0.683  1.0000
    ##  EMM_F7 - Control   0.25333 0.322 221   0.786  1.0000
    ##  SP_F14 - Control   0.28000 0.322 221   0.869  1.0000
    ##  EMM_F70 - Control  0.28667 0.322 221   0.889  0.9999
    ##  EMM_F3 - Control   0.34667 0.322 221   1.076  0.9995
    ##  EMM_F64 - Control  0.44000 0.322 221   1.365  0.9927
    ##  EMM_F65 - Control  0.46333 0.322 221   1.437  0.9877
    ##  EMM_F66 - Control  0.57667 0.322 221   1.789  0.9178
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F70  0.20667 0.322 221   0.641  1.0000
    ##  Control - EMM_F34  0.11000 0.322 221   0.341  1.0000
    ##  Control - EMM_F66  0.10000 0.322 221   0.310  1.0000
    ##  EMM_F89 - Control  0.03333 0.322 221   0.103  1.0000
    ##  EMM_F3 - Control   0.03667 0.322 221   0.114  1.0000
    ##  SP_F14 - Control   0.04333 0.322 221   0.134  1.0000
    ##  EMM_F65 - Control  0.05667 0.322 221   0.176  1.0000
    ##  EMM_F63 - Control  0.06000 0.322 221   0.186  1.0000
    ##  EMM_F47 - Control  0.07667 0.322 221   0.238  1.0000
    ##  ZAN_F4 - Control   0.08000 0.322 221   0.248  1.0000
    ##  EMM_F49 - Control  0.11000 0.322 221   0.341  1.0000
    ##  ZAN_F3 - Control   0.13333 0.322 221   0.414  1.0000
    ##  EMM_F7 - Control   0.15000 0.322 221   0.465  1.0000
    ##  EMM_F48 - Control  0.24000 0.322 221   0.745  1.0000
    ##  EMM_F64 - Control  0.24000 0.322 221   0.745  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - ZAN_F4   0.02000 0.322 221   0.062  1.0000
    ##  EMM_F70 - Control  0.03667 0.322 221   0.114  1.0000
    ##  EMM_F3 - Control   0.18000 0.322 221   0.558  1.0000
    ##  EMM_F47 - Control  0.22667 0.322 221   0.703  1.0000
    ##  SP_F14 - Control   0.24667 0.322 221   0.765  1.0000
    ##  EMM_F34 - Control  0.27000 0.322 221   0.838  1.0000
    ##  EMM_F63 - Control  0.28000 0.322 221   0.869  1.0000
    ##  EMM_F49 - Control  0.39667 0.322 221   1.231  0.9975
    ##  EMM_F7 - Control   0.43000 0.322 221   1.334  0.9942
    ##  EMM_F66 - Control  0.48333 0.322 221   1.500  0.9816
    ##  EMM_F64 - Control  0.49667 0.322 221   1.541  0.9763
    ##  ZAN_F3 - Control   0.50667 0.322 221   1.572  0.9717
    ##  EMM_F65 - Control  0.53667 0.322 221   1.665  0.9535
    ##  EMM_F89 - Control  0.56000 0.322 221   1.737  0.9344
    ##  EMM_F48 - Control  0.72667 0.322 221   2.254  0.6574
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates 
    ## distance_to_yeast = 11:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F65  4.47824 0.484 215   9.258  <.0001
    ##  Control - SP_F14   2.75000 0.432 215   6.362  <.0001
    ##  Control - EMM_F89  0.07667 0.432 215   0.177  1.0000
    ##  Control - EMM_F63  0.02000 0.432 215   0.046  1.0000
    ##  Control - ZAN_F4   0.01000 0.432 215   0.023  1.0000
    ##  EMM_F70 - Control  0.11333 0.432 215   0.262  1.0000
    ##  EMM_F3 - Control   0.12333 0.432 215   0.285  1.0000
    ##  EMM_F34 - Control  0.12333 0.432 215   0.285  1.0000
    ##  EMM_F64 - Control  0.13000 0.432 215   0.301  1.0000
    ##  ZAN_F3 - Control   0.29000 0.432 215   0.671  1.0000
    ##  EMM_F66 - Control  0.32667 0.432 215   0.756  1.0000
    ##  EMM_F49 - Control  0.34333 0.432 215   0.794  1.0000
    ##  EMM_F47 - Control  0.35333 0.432 215   0.817  1.0000
    ##  EMM_F7 - Control   0.48667 0.432 215   1.126  0.9991
    ##  EMM_F48 - Control  0.53667 0.432 215   1.241  0.9973
    ## 
    ## distance_to_yeast = 17:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.97000 0.432 215   2.244  0.6650
    ##  Control - EMM_F65  0.91658 0.484 215   1.895  0.8756
    ##  Control - EMM_F63  0.04667 0.432 215   0.108  1.0000
    ##  EMM_F89 - Control  0.18000 0.432 215   0.416  1.0000
    ##  EMM_F64 - Control  0.28333 0.432 215   0.655  1.0000
    ##  ZAN_F4 - Control   0.41333 0.432 215   0.956  0.9999
    ##  EMM_F49 - Control  0.42333 0.432 215   0.979  0.9998
    ##  EMM_F66 - Control  0.46000 0.432 215   1.064  0.9995
    ##  EMM_F7 - Control   0.49333 0.432 215   1.141  0.9989
    ##  EMM_F3 - Control   0.59000 0.432 215   1.365  0.9927
    ##  EMM_F70 - Control  0.59000 0.432 215   1.365  0.9927
    ##  ZAN_F3 - Control   0.59667 0.432 215   1.380  0.9918
    ##  EMM_F47 - Control  0.69000 0.432 215   1.596  0.9676
    ##  EMM_F34 - Control  0.70667 0.432 215   1.635  0.9601
    ##  EMM_F48 - Control  1.05000 0.432 215   2.429  0.5282
    ## 
    ## distance_to_yeast = 25:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F49 - Control  0.05333 0.432 215   0.123  1.0000
    ##  EMM_F89 - Control  0.06667 0.432 215   0.154  1.0000
    ##  EMM_F64 - Control  0.13333 0.432 215   0.308  1.0000
    ##  ZAN_F4 - Control   0.16000 0.432 215   0.370  1.0000
    ##  EMM_F63 - Control  0.21000 0.432 215   0.486  1.0000
    ##  EMM_F34 - Control  0.22333 0.432 215   0.517  1.0000
    ##  EMM_F66 - Control  0.26667 0.432 215   0.617  1.0000
    ##  EMM_F7 - Control   0.41333 0.432 215   0.956  0.9999
    ##  EMM_F47 - Control  0.44000 0.432 215   1.018  0.9997
    ##  SP_F14 - Control   0.45000 0.432 215   1.041  0.9996
    ##  EMM_F70 - Control  0.49333 0.432 215   1.141  0.9989
    ##  EMM_F48 - Control  0.50000 0.432 215   1.157  0.9988
    ##  EMM_F3 - Control   0.53000 0.432 215   1.226  0.9976
    ##  ZAN_F3 - Control   0.61000 0.432 215   1.411  0.9897
    ##  EMM_F65 - Control  0.74009 0.484 215   1.530  0.9778
    ## 
    ## distance_to_yeast = 32:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - SP_F14   0.21667 0.432 215   0.501  1.0000
    ##  Control - EMM_F89  0.17333 0.432 215   0.401  1.0000
    ##  Control - EMM_F66  0.13333 0.432 215   0.308  1.0000
    ##  Control - EMM_F70  0.09333 0.432 215   0.216  1.0000
    ##  Control - EMM_F49  0.06333 0.432 215   0.147  1.0000
    ##  Control - EMM_F65  0.03324 0.484 215   0.069  1.0000
    ##  EMM_F34 - Control  0.03000 0.432 215   0.069  1.0000
    ##  EMM_F47 - Control  0.04000 0.432 215   0.093  1.0000
    ##  ZAN_F4 - Control   0.04333 0.432 215   0.100  1.0000
    ##  EMM_F7 - Control   0.05333 0.432 215   0.123  1.0000
    ##  EMM_F64 - Control  0.05667 0.432 215   0.131  1.0000
    ##  EMM_F63 - Control  0.09667 0.432 215   0.224  1.0000
    ##  EMM_F3 - Control   0.10000 0.432 215   0.231  1.0000
    ##  EMM_F48 - Control  0.35667 0.432 215   0.825  1.0000
    ##  ZAN_F3 - Control   0.36000 0.432 215   0.833  1.0000
    ## 
    ## distance_to_yeast = 41:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F34  0.13333 0.432 215   0.308  1.0000
    ##  Control - EMM_F7   0.06000 0.432 215   0.139  1.0000
    ##  Control - EMM_F89  0.05333 0.432 215   0.123  1.0000
    ##  Control - EMM_F63  0.01000 0.432 215   0.023  1.0000
    ##  ZAN_F3 - Control   0.07333 0.432 215   0.170  1.0000
    ##  EMM_F49 - Control  0.11000 0.432 215   0.254  1.0000
    ##  ZAN_F4 - Control   0.16333 0.432 215   0.378  1.0000
    ##  SP_F14 - Control   0.22667 0.432 215   0.524  1.0000
    ##  EMM_F48 - Control  0.33333 0.432 215   0.771  1.0000
    ##  EMM_F64 - Control  0.38000 0.432 215   0.879  1.0000
    ##  EMM_F47 - Control  0.44667 0.432 215   1.033  0.9997
    ##  EMM_F70 - Control  0.51333 0.432 215   1.187  0.9983
    ##  EMM_F66 - Control  0.54667 0.432 215   1.265  0.9967
    ##  EMM_F65 - Control  0.61009 0.484 215   1.261  0.9968
    ##  EMM_F3 - Control   0.65000 0.432 215   1.504  0.9811
    ## 
    ## distance_to_yeast = 48:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F89  0.28667 0.432 215   0.663  1.0000
    ##  Control - EMM_F66  0.21667 0.432 215   0.501  1.0000
    ##  Control - SP_F14   0.21000 0.432 215   0.486  1.0000
    ##  Control - EMM_F70  0.16000 0.432 215   0.370  1.0000
    ##  Control - EMM_F34  0.03000 0.432 215   0.069  1.0000
    ##  Control - EMM_F7   0.02667 0.432 215   0.062  1.0000
    ##  EMM_F47 - Control  0.00000 0.432 215   0.000  1.0000
    ##  ZAN_F4 - Control   0.03000 0.432 215   0.069  1.0000
    ##  EMM_F63 - Control  0.08333 0.432 215   0.193  1.0000
    ##  EMM_F65 - Control  0.09009 0.484 215   0.186  1.0000
    ##  EMM_F3 - Control   0.09667 0.432 215   0.224  1.0000
    ##  ZAN_F3 - Control   0.15000 0.432 215   0.347  1.0000
    ##  EMM_F64 - Control  0.17667 0.432 215   0.409  1.0000
    ##  EMM_F49 - Control  0.19000 0.432 215   0.440  1.0000
    ##  EMM_F48 - Control  0.24333 0.432 215   0.563  1.0000
    ## 
    ## distance_to_yeast = 55:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  Control - EMM_F47  0.15333 0.432 215   0.355  1.0000
    ##  Control - EMM_F7   0.08333 0.432 215   0.193  1.0000
    ##  Control - EMM_F70  0.07000 0.432 215   0.162  1.0000
    ##  Control - ZAN_F4   0.00667 0.432 215   0.015  1.0000
    ##  ZAN_F3 - Control   0.01333 0.432 215   0.031  1.0000
    ##  EMM_F66 - Control  0.11000 0.432 215   0.254  1.0000
    ##  EMM_F34 - Control  0.13000 0.432 215   0.301  1.0000
    ##  EMM_F63 - Control  0.18000 0.432 215   0.416  1.0000
    ##  EMM_F3 - Control   0.30000 0.432 215   0.694  1.0000
    ##  EMM_F89 - Control  0.31333 0.432 215   0.725  1.0000
    ##  SP_F14 - Control   0.33333 0.432 215   0.771  1.0000
    ##  EMM_F65 - Control  0.39009 0.484 215   0.806  1.0000
    ##  EMM_F48 - Control  0.43000 0.432 215   0.995  0.9998
    ##  EMM_F49 - Control  0.43667 0.432 215   1.010  0.9998
    ##  EMM_F64 - Control  0.50333 0.432 215   1.164  0.9987
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 16 estimates

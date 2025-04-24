## Things to know before analysis

***This is applied to both replicate***

- Independent variable is colony diameter, colony weight and colony
  forming Unit
- Dependent variable is DAI (Days after Inoculation) in case of
  diameter, yeast isolates for all parameters

``` r
#load the library
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

## Biological Replicate 1

``` r
#read data
BioRep1 <- read.csv("Splitplate/SplitPlateData/2024-09-16_SplitPlate_B52.csv", na.strings = "na")
#load the color palettes
cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#CC79A7")
#check the structure and see the data
str(BioRep1)
```

    ## 'data.frame':    384 obs. of  10 variables:
    ##  $ Bacteria       : chr  "B52" "B52" "B52" "B52" ...
    ##  $ Yeast          : chr  "EMM_F49" "EMM_F49" "EMM_F49" "EMM_F49" ...
    ##  $ Media          : chr  "YePDA" "YePDA" "YePDA" "YePDA" ...
    ##  $ Class          : chr  "Microbotryomycetes" "Microbotryomycetes" "Microbotryomycetes" "Microbotryomycetes" ...
    ##  $ Replication    : int  1 1 2 2 3 3 4 4 5 5 ...
    ##  $ DAI            : int  2 2 2 2 2 2 2 2 2 2 ...
    ##  $ colony_diameter: num  6.85 6.85 7.68 7.69 8.74 8.01 7.88 7.25 7.81 7.8 ...
    ##  $ increase       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ colony_weight  : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ CFU            : num  NA NA NA NA NA NA NA NA NA NA ...

``` r
head(BioRep1)
```

    ##   Bacteria   Yeast Media              Class Replication DAI colony_diameter
    ## 1      B52 EMM_F49 YePDA Microbotryomycetes           1   2            6.85
    ## 2      B52 EMM_F49 YePDA Microbotryomycetes           1   2            6.85
    ## 3      B52 EMM_F49 YePDA Microbotryomycetes           2   2            7.68
    ## 4      B52 EMM_F49 YePDA Microbotryomycetes           2   2            7.69
    ## 5      B52 EMM_F49 YePDA Microbotryomycetes           3   2            8.74
    ## 6      B52 EMM_F49 YePDA Microbotryomycetes           3   2            8.01
    ##   increase colony_weight CFU
    ## 1        0            NA  NA
    ## 2        0            NA  NA
    ## 3        0            NA  NA
    ## 4        0            NA  NA
    ## 5        0            NA  NA
    ## 6        0            NA  NA

``` r
#changing the variables into categorical data
BioRep1$Replication=as.factor(BioRep1$Replication)
BioRep1$Media=as.factor(BioRep1$Media)
BioRep1$Yeast=as.factor(BioRep1$Yeast)
BioRep1$DAI=as.factor(BioRep1$DAI)
```

### Impact of Yeast on Colony Size

``` r
#filter 1st data date and removing the 1st data since there was no growth associated with it. The data is considered to be the baseline for growth so since it violates the assumption of normality, we remove it from the analysis

filterdata1 <-BioRep1[BioRep1$Media == "YePDA"& BioRep1$DAI != 2,]
#load library
library(nlme)
```

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
#stats used for size
size1=lme(increase~Yeast*DAI, data = filterdata1,random=~1|Replication) #mixed effect model
size2=lm(increase ~ Yeast*DAI, data = filterdata1) #use simpler model
summary(size1)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: filterdata1 
    ##       AIC      BIC    logLik
    ##   523.579 563.9382 -247.7895
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept) Residual
    ## StdDev:   0.5200207 1.373088
    ## 
    ## Fixed effects:  increase ~ Yeast * DAI 
    ##                       Value Std.Error  DF   t-value p-value
    ## (Intercept)        2.201667 0.4496493 127  4.896409  0.0000
    ## YeastEMM_F48      -1.695833 0.5605608 127 -3.025244  0.0030
    ## YeastEMM_F49      -1.825833 0.5605608 127 -3.257155  0.0014
    ## YeastSP_F14       -1.988333 0.5605608 127 -3.547043  0.0005
    ## DAI6               3.968333 0.5605608 127  7.079220  0.0000
    ## DAI8               6.285000 0.5605608 127 11.211986  0.0000
    ## YeastEMM_F48:DAI6 -3.514167 0.7927527 127 -4.432866  0.0000
    ## YeastEMM_F49:DAI6 -3.617500 0.7927527 127 -4.563214  0.0000
    ## YeastSP_F14:DAI6  -3.405833 0.7927527 127 -4.296211  0.0000
    ## YeastEMM_F48:DAI8 -5.545000 0.7927527 127 -6.994615  0.0000
    ## YeastEMM_F49:DAI8 -5.741667 0.7927527 127 -7.242696  0.0000
    ## YeastSP_F14:DAI8  -5.397500 0.7927527 127 -6.808554  0.0000
    ##  Correlation: 
    ##                   (Intr) YsEMM_F48 YsEMM_F49 YsSP_F14 DAI6   DAI8  
    ## YeastEMM_F48      -0.623                                           
    ## YeastEMM_F49      -0.623  0.500                                    
    ## YeastSP_F14       -0.623  0.500     0.500                          
    ## DAI6              -0.623  0.500     0.500     0.500                
    ## DAI8              -0.623  0.500     0.500     0.500    0.500       
    ## YeastEMM_F48:DAI6  0.441 -0.707    -0.354    -0.354   -0.707 -0.354
    ## YeastEMM_F49:DAI6  0.441 -0.354    -0.707    -0.354   -0.707 -0.354
    ## YeastSP_F14:DAI6   0.441 -0.354    -0.354    -0.707   -0.707 -0.354
    ## YeastEMM_F48:DAI8  0.441 -0.707    -0.354    -0.354   -0.354 -0.707
    ## YeastEMM_F49:DAI8  0.441 -0.354    -0.707    -0.354   -0.354 -0.707
    ## YeastSP_F14:DAI8   0.441 -0.354    -0.354    -0.707   -0.354 -0.707
    ##                   YEMM_F48:DAI6 YEMM_F49:DAI6 YSP_F14:DAI6 YEMM_F48:DAI8
    ## YeastEMM_F48                                                            
    ## YeastEMM_F49                                                            
    ## YeastSP_F14                                                             
    ## DAI6                                                                    
    ## DAI8                                                                    
    ## YeastEMM_F48:DAI6                                                       
    ## YeastEMM_F49:DAI6  0.500                                                
    ## YeastSP_F14:DAI6   0.500         0.500                                  
    ## YeastEMM_F48:DAI8  0.500         0.250         0.250                    
    ## YeastEMM_F49:DAI8  0.250         0.500         0.250        0.500       
    ## YeastSP_F14:DAI8   0.250         0.250         0.500        0.500       
    ##                   YEMM_F49:DAI8
    ## YeastEMM_F48                   
    ## YeastEMM_F49                   
    ## YeastSP_F14                    
    ## DAI6                           
    ## DAI8                           
    ## YeastEMM_F48:DAI6              
    ## YeastEMM_F49:DAI6              
    ## YeastSP_F14:DAI6               
    ## YeastEMM_F48:DAI8              
    ## YeastEMM_F49:DAI8              
    ## YeastSP_F14:DAI8   0.500       
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -3.2300287 -0.4369499 -0.0165479  0.3504212  4.0042258 
    ## 
    ## Number of Observations: 144
    ## Number of Groups: 6

``` r
anova(size1, size2) #to check if complex or simple model is better
```

    ##       Model df      AIC      BIC    logLik   Test  L.Ratio p-value
    ## size1     1 14 523.5790 563.9382 -247.7895                        
    ## size2     2 13 530.3016 567.7780 -252.1508 1 vs 2 8.722618  0.0031

``` r
#pairwise comparison
lsmeans.size1 <- emmeans(size1, ~ Yeast|DAI) 
# estimate lsmeans 
comp.size1 <- multcomp::cld(object = lsmeans.size1, alpha = 0.05, Letters = letters, reversed = TRUE, details = TRUE)
# contrast with Tukey ajustment 
comp.size1
```

    ## $emmeans
    ## DAI = 4:
    ##  Yeast   emmean   SE df lower.CL upper.CL .group
    ##  Control  2.202 0.45  5    1.046     3.36  a    
    ##  EMM_F48  0.506 0.45  5   -0.650     1.66   b   
    ##  EMM_F49  0.376 0.45  5   -0.780     1.53   b   
    ##  SP_F14   0.213 0.45  5   -0.943     1.37   b   
    ## 
    ## DAI = 6:
    ##  Yeast   emmean   SE df lower.CL upper.CL .group
    ##  Control  6.170 0.45  5    5.014     7.33  a    
    ##  EMM_F48  0.960 0.45  5   -0.196     2.12   b   
    ##  SP_F14   0.776 0.45  5   -0.380     1.93   b   
    ##  EMM_F49  0.727 0.45  5   -0.429     1.88   b   
    ## 
    ## DAI = 8:
    ##  Yeast   emmean   SE df lower.CL upper.CL .group
    ##  Control  8.487 0.45  5    7.331     9.64  a    
    ##  EMM_F48  1.246 0.45  5    0.090     2.40   b   
    ##  SP_F14   1.101 0.45  5   -0.055     2.26   b   
    ##  EMM_F49  0.919 0.45  5   -0.237     2.08   b   
    ## 
    ## Degrees-of-freedom method: containment 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ## DAI = 4:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F49 - SP_F14    0.1625 0.561 127   0.290  0.9915
    ##  EMM_F48 - SP_F14    0.2925 0.561 127   0.522  0.9537
    ##  EMM_F48 - EMM_F49   0.1300 0.561 127   0.232  0.9956
    ##  Control - SP_F14    1.9883 0.561 127   3.547  0.0030
    ##  Control - EMM_F49   1.8258 0.561 127   3.257  0.0078
    ##  Control - EMM_F48   1.6958 0.561 127   3.025  0.0157
    ## 
    ## DAI = 6:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  SP_F14 - EMM_F49    0.0492 0.561 127   0.088  0.9998
    ##  EMM_F48 - EMM_F49   0.2333 0.561 127   0.416  0.9756
    ##  EMM_F48 - SP_F14    0.1842 0.561 127   0.329  0.9877
    ##  Control - EMM_F49   5.4433 0.561 127   9.711  <.0001
    ##  Control - SP_F14    5.3942 0.561 127   9.623  <.0001
    ##  Control - EMM_F48   5.2100 0.561 127   9.294  <.0001
    ## 
    ## DAI = 8:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  SP_F14 - EMM_F49    0.1817 0.561 127   0.324  0.9882
    ##  EMM_F48 - EMM_F49   0.3267 0.561 127   0.583  0.9371
    ##  EMM_F48 - SP_F14    0.1450 0.561 127   0.259  0.9939
    ##  Control - EMM_F49   7.5675 0.561 127  13.500  <.0001
    ##  Control - SP_F14    7.3858 0.561 127  13.176  <.0001
    ##  Control - EMM_F48   7.2408 0.561 127  12.917  <.0001
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(comp.size1$emmeans$Yeast, 
                               comp.size1$emmeans$DAI,
                               str_trim(comp.size1$emmeans$.group))
colnames(sig.diff.letters) <- c("Yeast", 
                                "DAI",
                                "Letters")

# for plotting with letters from significance test
size1plot <- filterdata1 %>%
  group_by(Yeast, DAI) %>%
  dplyr::summarize(
    colonydiameter = mean(increase, na.rm=TRUE),
    se = sd(increase)/sqrt(3)) %>%
  left_join(sig.diff.letters) 
```

    ## `summarise()` has grouped output by 'Yeast'. You can override using the
    ## `.groups` argument.

    ## Joining with `by = join_by(Yeast, DAI)`

``` r
colonysize1 <- filterdata1 %>%
  ggplot(aes(x = Yeast, y = increase, group = Yeast, fill = Yeast)) +
  geom_boxplot()+
  geom_point(shape = 21, color = "black", position = position_jitterdodge(dodge.width = 0.9))+
   geom_text(data = size1plot, aes(label = Letters, y = colonydiameter+(3*se)), vjust = -0.75, hjust = 2) +
  ylab("Increase in Colony Diameter (mm)") +   
  xlab("Yeast Isolates") +  
    scale_fill_manual(values = cbbPalette) +
  theme(strip.text.x = element_text(size = 18, face= "italic"))+
  theme(axis.title.x = element_text(size =25, face = "bold"), axis.title.y = element_text(size = 25, face = "bold"), axis.text.x=element_text(size=20, angle= 0), axis.text.y=element_text(size=20, angle= 0))+
  theme_classic() +
  theme(legend.position = "none")+
  ggtitle(expression(paste("Impact of Yeast on Growth of ", italic("Methylobacterium platani"), " EMM_B52 (Biological Replicate 1)"))) +
   theme(plot.title = element_text(face = "bold", color = "Blue2", size = 12, hjust = 0.5))+
  facet_wrap(~DAI)
colonysize1
```

![](Splitplate_files/figure-gfm/Increase%20in%20colony%20size%20Plot%20for%20Biological%20replicate%201-1.png)<!-- -->

### Impact of Yeast on Colony Weight

``` r
#filter data to avoid other data except of last day since it was detructive sample taking procedure and colony was taken out from the plate to measure the weight on the last day of data collection.

filterdata1_2 <- filterdata1 %>% 
  filter(DAI == "8" )

#stats used for colony weight
cw1=lme(colony_weight~Yeast,data = filterdata1_2, random = ~1|Replication) ##replication random effect added by Dr. Steury
summary(cw1)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: filterdata1_2 
    ##         AIC       BIC   logLik
    ##   -217.6423 -206.9372 114.8212
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept)   Residual
    ## StdDev: 0.007487334 0.01493443
    ## 
    ## Fixed effects:  colony_weight ~ Yeast 
    ##                    Value   Std.Error DF   t-value p-value
    ## (Intercept)   0.11754000 0.005284865 39  22.24087       0
    ## YeastEMM_F48 -0.08830833 0.006096956 39 -14.48400       0
    ## YeastEMM_F49 -0.09863000 0.006096956 39 -16.17692       0
    ## YeastSP_F14  -0.09421000 0.006096956 39 -15.45197       0
    ##  Correlation: 
    ##              (Intr) YEMM_F48 YEMM_F49
    ## YeastEMM_F48 -0.577                  
    ## YeastEMM_F49 -0.577  0.500           
    ## YeastSP_F14  -0.577  0.500    0.500  
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.30425553 -0.52026994 -0.04957339  0.36278982  2.14854189 
    ## 
    ## Number of Observations: 48
    ## Number of Groups: 6

``` r
Anova(cw1)
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: colony_weight
    ##        Chisq Df Pr(>Chisq)    
    ## Yeast 357.29  3  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#pairwise comparision
lsmeans.cw1 <- emmeans(cw1, ~ Yeast) 

# contrast with Tukey ajustment 
comp.cw1 <- multcomp::cld(object = lsmeans.cw1, alpha = 0.05, Letters = letters, reversed = TRUE, details = TRUE)
comp.cw1
```

    ## $emmeans
    ##  Yeast   emmean      SE df lower.CL upper.CL .group
    ##  Control 0.1175 0.00528  5  0.10395   0.1311  a    
    ##  EMM_F48 0.0292 0.00528  5  0.01565   0.0428   b   
    ##  SP_F14  0.0233 0.00528  5  0.00974   0.0369   b   
    ##  EMM_F49 0.0189 0.00528  5  0.00532   0.0325   b   
    ## 
    ## Degrees-of-freedom method: containment 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ##  contrast          estimate     SE df t.ratio p.value
    ##  SP_F14 - EMM_F49   0.00442 0.0061 39   0.725  0.8865
    ##  EMM_F48 - EMM_F49  0.01032 0.0061 39   1.693  0.3411
    ##  EMM_F48 - SP_F14   0.00590 0.0061 39   0.968  0.7683
    ##  Control - EMM_F49  0.09863 0.0061 39  16.177  <.0001
    ##  Control - SP_F14   0.09421 0.0061 39  15.452  <.0001
    ##  Control - EMM_F48  0.08831 0.0061 39  14.484  <.0001
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(comp.cw1$emmeans$Yeast, 
                               str_trim(comp.cw1$emmeans$.group))

colnames(sig.diff.letters) <- c("Yeast", "Letters")

# for plotting with letters from significance test
weightplot <- filterdata1_2 %>%
  group_by(Yeast) %>%
  dplyr::summarize(
    colonyweight = mean(colony_weight, na.rm=TRUE),
    se = sd(colony_weight)/sqrt(3)) %>%
  left_join(sig.diff.letters) 
```

    ## Joining with `by = join_by(Yeast)`

``` r
 colonyweightplot <- filterdata1_2 %>%
  ggplot(aes(x = Yeast, y = colony_weight, group = Yeast, fill = Yeast)) +
  geom_boxplot()+
  geom_point(shape = 21, color = "black", position = position_jitterdodge(dodge.width = 0.9))+
  geom_text(data = weightplot, aes(label = Letters, y = colonyweight+(2*se)), vjust = -0.5, hjust = 2) +
  ylab("Colony Weight (gram)") +   
  xlab("Yeast Isolates") +  
  scale_fill_manual(values = cbbPalette) +
  theme(strip.text.x = element_text(size = 18, face= "italic"))+
  theme(axis.title.x = element_text(size =25, face = "bold"), axis.title.y = element_text(size = 25, face = "bold"), axis.text.x=element_text(size=20, angle= 0), axis.text.y=element_text(size=20, angle= 0))+
  theme_classic() +
  theme(legend.position = "none")
  #ggtitle(expression(paste("Impact of Yeast on Growth of ", italic("Methylobacterium platani"), " EMM_B52"))) +
   #theme(plot.title = element_text(face = "bold", color = "Blue2", size = 12, hjust = 0.5))
colonyweightplot
```

![](Splitplate_files/figure-gfm/Colony%20Weight%20Plot%20for%20Biological%20Replicate%201-1.png)<!-- -->

### Impact of Yeast on Colony Forming Unit (CFU)

``` r
#stats used for colony forming unit
cfu1=lme(log10(CFU)~Yeast,data = filterdata1_2, random = ~1|Replication, na.action = na.omit) 
summary(cfu1)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: filterdata1_2 
    ##        AIC      BIC    logLik
    ##   30.88478 36.85918 -9.442392
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept)  Residual
    ## StdDev:    0.548129 0.2721204
    ## 
    ## Fixed effects:  log10(CFU) ~ Yeast 
    ##                  Value Std.Error DF   t-value p-value
    ## (Intercept)  13.040907 0.3353953 18  38.88219       0
    ## YeastEMM_F48 -1.883912 0.1571088 18 -11.99113       0
    ## YeastEMM_F49 -2.013901 0.1571088 18 -12.81851       0
    ## YeastSP_F14  -1.875739 0.1571088 18 -11.93911       0
    ##  Correlation: 
    ##              (Intr) YEMM_F48 YEMM_F49
    ## YeastEMM_F48 -0.234                  
    ## YeastEMM_F49 -0.234  0.500           
    ## YeastSP_F14  -0.234  0.500    0.500  
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -1.2149837 -0.6404125 -0.1836397  0.7017032  1.8113339 
    ## 
    ## Number of Observations: 24
    ## Number of Groups: 3

``` r
Anova(cfu1)
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: log10(CFU)
    ##        Chisq Df Pr(>Chisq)    
    ## Yeast 226.05  3  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#pairwise comparision
lsmeans.cfu1 <- emmeans(cfu1, ~ Yeast) ###interaction added by Dr. Steury
#emmeans(cfu,pairwise~Yeast) ##added by Dr. steury
# estimate lsmeans 
comp.cfu1 <- multcomp::cld(object = lsmeans.cfu1, alpha = 0.05, Letters = letters, reversed = TRUE, details = TRUE)
# contrast with Tukey ajustment 
comp.cfu1
```

    ## $emmeans
    ##  Yeast   emmean    SE df lower.CL upper.CL .group
    ##  Control   13.0 0.335  2    11.60     14.5  a    
    ##  SP_F14    11.2 0.335  2     9.72     12.6   b   
    ##  EMM_F48   11.2 0.335  2     9.71     12.6   b   
    ##  EMM_F49   11.0 0.335  2     9.58     12.5   b   
    ## 
    ## Degrees-of-freedom method: containment 
    ## Results are given on the log10 (not the response) scale. 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ##  contrast          estimate    SE df t.ratio p.value
    ##  EMM_F48 - EMM_F49  0.12999 0.157 18   0.827  0.8408
    ##  SP_F14 - EMM_F49   0.13816 0.157 18   0.879  0.8154
    ##  SP_F14 - EMM_F48   0.00817 0.157 18   0.052  0.9999
    ##  Control - EMM_F49  2.01390 0.157 18  12.819  <.0001
    ##  Control - EMM_F48  1.88391 0.157 18  11.991  <.0001
    ##  Control - SP_F14   1.87574 0.157 18  11.939  <.0001
    ## 
    ## Degrees-of-freedom method: containment 
    ## Results are given on the log10 (not the response) scale. 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(comp.cfu1$emmeans$Yeast, 
                               str_trim(comp.cfu1$emmeans$.group))

colnames(sig.diff.letters) <- c("Yeast", "Letters")

# for plotting with letters from significance test
cfu1plot <- filterdata1_2 %>%
  filter(Replication %in% c(1:3)) %>% 
  group_by(Yeast) %>%
  dplyr::summarize(
    cfu1 = mean(log10(CFU), na.rm=TRUE),
    se = sd(log10(CFU))/sqrt(3)) %>%
  left_join(sig.diff.letters) 
```

    ## Joining with `by = join_by(Yeast)`

``` r
CFU1PLOT <- filterdata1 %>%
  subset(Media=="YePDA") %>%
  #subset(Yeast %in% c("Control", "SP_F14","EMM_F3", "EMM_F34")) %>%
  ggplot(aes(x = Yeast, y = log10(CFU), group = Yeast, fill = Yeast)) +
  geom_boxplot()+
  geom_point(shape = 21, color = "black", position = position_jitterdodge(dodge.width = 0.9))+
  geom_text(data = cfu1plot, aes(label = Letters, y = cfu1+(4*se)), vjust = -0.1, hjust = 2) +
  ylab("log of Colony Forming Unit per ml") +   
  xlab("Yeast Isolates") +  
  scale_fill_manual(values = cbbPalette) +
  theme(strip.text.x = element_text(size = 18, face= "italic"))+
  theme(axis.title.x = element_text(size =25, face = "bold"), axis.title.y = element_text(size = 25, face = "bold"), axis.text.x=element_text(size=10, angle= 0), axis.text.y=element_text(size=20, angle= 0))+
  theme_classic() +
  theme(legend.position = "none")
  #ggtitle("Impact of Yeast on Growth of Methylobacterium platani EMM_B52") +
  #theme(plot.title = element_text(hjust = 0.5))
CFU1PLOT
```

    ## Warning: Removed 120 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 120 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Splitplate_files/figure-gfm/Colony%20Forming%20Unit(CFU)%20Plot%20for%20Biological%20Replicate%201-1.png)<!-- -->

``` r
#combining similar dataset plots into one
cw1andcfu1 <- ggarrange(colonyweightplot, CFU1PLOT, nrow = 1, ncol = 2, common.legend = F)
```

    ## Warning: Removed 120 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 120 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
cw1andcfu1 <- annotate_figure(cw1andcfu1, top = text_grob(expression(paste("Impact of Yeast on Growth of ", italic("Methylobacterium platani"), " EMM_B52 (Biological  Replicate 1)")), color = "Blue2", size = 12, face = "bold", hjust = 0.5))
cw1andcfu1
```

![](Splitplate_files/figure-gfm/CFU%20and%20colony%20weight%20Combined%20Plot%20for%20Biological%20Replicate%201-1.png)<!-- -->

## Biological Replicate 2

``` r
#read data
BioRep2 <- read.csv("Splitplate/SplitPlateData/2024-11-11_SplitPlate_B52-1.csv", na.strings = "na")
#load the color palettes
cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#CC79A7")
#check the structure and see the data
str(BioRep2)
```

    ## 'data.frame':    240 obs. of  8 variables:
    ##  $ Bacteria       : chr  "B52" "B52" "B52" "B52" ...
    ##  $ Yeast          : chr  "EMM_F48" "EMM_F48" "EMM_F48" "EMM_F48" ...
    ##  $ Replication    : int  1 1 2 2 3 3 4 4 5 5 ...
    ##  $ DAI            : int  2 2 2 2 2 2 2 2 2 2 ...
    ##  $ Colony_diameter: num  7.46 7.83 7.83 7.83 7.55 7.65 7.65 7.8 7.6 7.73 ...
    ##  $ increase       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ colony_weight  : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ CFU            : num  NA NA NA NA NA NA NA NA NA NA ...

``` r
head(BioRep2)
```

    ##   Bacteria   Yeast Replication DAI Colony_diameter increase colony_weight CFU
    ## 1      B52 EMM_F48           1   2            7.46        0            NA  NA
    ## 2      B52 EMM_F48           1   2            7.83        0            NA  NA
    ## 3      B52 EMM_F48           2   2            7.83        0            NA  NA
    ## 4      B52 EMM_F48           2   2            7.83        0            NA  NA
    ## 5      B52 EMM_F48           3   2            7.55        0            NA  NA
    ## 6      B52 EMM_F48           3   2            7.65        0            NA  NA

``` r
#changing the variable into categorical data
BioRep2$Replication=as.factor(BioRep2$Replication)
BioRep2$Yeast=as.factor(BioRep2$Yeast)
BioRep2$DAI=as.factor(BioRep2$DAI)
```

### Impact of Yeast on Colony Size

``` r
#removing the 1st data since there was no growth associated with it
#the data is considered to be the baseline for growth so since it violates the assumption of normality, we remove it from the analysis. For this Day 10 data was also collected. Thus, for consistency, we remove the data for visualization and analysis.
filterdata2 <-BioRep2 %>% 
filter(DAI != "2", DAI != "10")

#load library
library(nlme)
#stats used for size
size2=lme(increase~Yeast*DAI,data = filterdata2,random=~1|Replication) ##replication random effect added by Dr. Steury
summary(size2)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: filterdata2 
    ##        AIC      BIC    logLik
    ##   311.0307 351.3899 -141.5154
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept)  Residual
    ## StdDev:   0.1111038 0.6246895
    ## 
    ## Fixed effects:  increase ~ Yeast * DAI 
    ##                        Value Std.Error  DF   t-value p-value
    ## (Intercept)        2.2116667 0.1859491 127 11.893933  0.0000
    ## YeastEMM_F48      -1.0716667 0.2550284 127 -4.202146  0.0000
    ## YeastEMM_F49      -1.5658333 0.2550284 127 -6.139839  0.0000
    ## YeastSP_F14       -1.0233333 0.2550284 127 -4.012625  0.0001
    ## DAI6               1.7233333 0.2550284 127  6.757417  0.0000
    ## DAI8               2.9225000 0.2550284 127 11.459508  0.0000
    ## YeastEMM_F48:DAI6 -1.1375000 0.3606646 127 -3.153900  0.0020
    ## YeastEMM_F49:DAI6 -0.9391667 0.3606646 127 -2.603989  0.0103
    ## YeastSP_F14:DAI6  -0.8116667 0.3606646 127 -2.250475  0.0261
    ## YeastEMM_F48:DAI8 -1.8066667 0.3606646 127 -5.009271  0.0000
    ## YeastEMM_F49:DAI8 -1.8750000 0.3606646 127 -5.198736  0.0000
    ## YeastSP_F14:DAI8  -1.7200000 0.3606646 127 -4.768973  0.0000
    ##  Correlation: 
    ##                   (Intr) YsEMM_F48 YsEMM_F49 YsSP_F14 DAI6   DAI8  
    ## YeastEMM_F48      -0.686                                           
    ## YeastEMM_F49      -0.686  0.500                                    
    ## YeastSP_F14       -0.686  0.500     0.500                          
    ## DAI6              -0.686  0.500     0.500     0.500                
    ## DAI8              -0.686  0.500     0.500     0.500    0.500       
    ## YeastEMM_F48:DAI6  0.485 -0.707    -0.354    -0.354   -0.707 -0.354
    ## YeastEMM_F49:DAI6  0.485 -0.354    -0.707    -0.354   -0.707 -0.354
    ## YeastSP_F14:DAI6   0.485 -0.354    -0.354    -0.707   -0.707 -0.354
    ## YeastEMM_F48:DAI8  0.485 -0.707    -0.354    -0.354   -0.354 -0.707
    ## YeastEMM_F49:DAI8  0.485 -0.354    -0.707    -0.354   -0.354 -0.707
    ## YeastSP_F14:DAI8   0.485 -0.354    -0.354    -0.707   -0.354 -0.707
    ##                   YEMM_F48:DAI6 YEMM_F49:DAI6 YSP_F14:DAI6 YEMM_F48:DAI8
    ## YeastEMM_F48                                                            
    ## YeastEMM_F49                                                            
    ## YeastSP_F14                                                             
    ## DAI6                                                                    
    ## DAI8                                                                    
    ## YeastEMM_F48:DAI6                                                       
    ## YeastEMM_F49:DAI6  0.500                                                
    ## YeastSP_F14:DAI6   0.500         0.500                                  
    ## YeastEMM_F48:DAI8  0.500         0.250         0.250                    
    ## YeastEMM_F49:DAI8  0.250         0.500         0.250        0.500       
    ## YeastSP_F14:DAI8   0.250         0.250         0.500        0.500       
    ##                   YEMM_F49:DAI8
    ## YeastEMM_F48                   
    ## YeastEMM_F49                   
    ## YeastSP_F14                    
    ## DAI6                           
    ## DAI8                           
    ## YeastEMM_F48:DAI6              
    ## YeastEMM_F49:DAI6              
    ## YeastSP_F14:DAI6               
    ## YeastEMM_F48:DAI8              
    ## YeastEMM_F49:DAI8              
    ## YeastSP_F14:DAI8   0.500       
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.54722046 -0.50154811  0.04139578  0.58429578  2.37202267 
    ## 
    ## Number of Observations: 144
    ## Number of Groups: 6

``` r
anova(size2)
```

    ##             numDF denDF  F-value p-value
    ## (Intercept)     1   127 973.4467  <.0001
    ## Yeast           3   127 112.3483  <.0001
    ## DAI             2   127  77.8973  <.0001
    ## Yeast:DAI       6   127   6.4271  <.0001

``` r
#pairwise comparision
lsmeans.size2 <- emmeans(size2, ~ Yeast|DAI) ###interaction added by Dr. Steury
#emmeans(size,pairwise~Yeast) ##added by Dr. steury
# estimate lsmeans 
comp.size2 <- multcomp::cld(object = lsmeans.size2, alpha = 0.05, Letters = letters, reversed = TRUE, details = TRUE)
# contrast with Tukey ajustment 
comp.size2
```

    ## $emmeans
    ## DAI = 4:
    ##  Yeast   emmean    SE df lower.CL upper.CL .group
    ##  Control  2.212 0.186  5    1.734     2.69  a    
    ##  SP_F14   1.188 0.186  5    0.710     1.67   b   
    ##  EMM_F48  1.140 0.186  5    0.662     1.62   b   
    ##  EMM_F49  0.646 0.186  5    0.168     1.12   b   
    ## 
    ## DAI = 6:
    ##  Yeast   emmean    SE df lower.CL upper.CL .group
    ##  Control  3.935 0.186  5    3.457     4.41  a    
    ##  SP_F14   2.100 0.186  5    1.622     2.58   b   
    ##  EMM_F48  1.726 0.186  5    1.248     2.20   bc  
    ##  EMM_F49  1.430 0.186  5    0.952     1.91    c  
    ## 
    ## DAI = 8:
    ##  Yeast   emmean    SE df lower.CL upper.CL .group
    ##  Control  5.134 0.186  5    4.656     5.61  a    
    ##  SP_F14   2.391 0.186  5    1.913     2.87   b   
    ##  EMM_F48  2.256 0.186  5    1.778     2.73   bc  
    ##  EMM_F49  1.693 0.186  5    1.215     2.17    c  
    ## 
    ## Degrees-of-freedom method: containment 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ## DAI = 4:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F48 - EMM_F49   0.4942 0.255 127   1.938  0.2176
    ##  SP_F14 - EMM_F49    0.5425 0.255 127   2.127  0.1501
    ##  SP_F14 - EMM_F48    0.0483 0.255 127   0.190  0.9976
    ##  Control - EMM_F49   1.5658 0.255 127   6.140  <.0001
    ##  Control - EMM_F48   1.0717 0.255 127   4.202  0.0003
    ##  Control - SP_F14    1.0233 0.255 127   4.013  0.0006
    ## 
    ## DAI = 6:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F48 - EMM_F49   0.2958 0.255 127   1.160  0.6531
    ##  SP_F14 - EMM_F49    0.6700 0.255 127   2.627  0.0470
    ##  SP_F14 - EMM_F48    0.3742 0.255 127   1.467  0.4603
    ##  Control - EMM_F49   2.5050 0.255 127   9.822  <.0001
    ##  Control - EMM_F48   2.2092 0.255 127   8.662  <.0001
    ##  Control - SP_F14    1.8350 0.255 127   7.195  <.0001
    ## 
    ## DAI = 8:
    ##  contrast          estimate    SE  df t.ratio p.value
    ##  EMM_F48 - EMM_F49   0.5625 0.255 127   2.206  0.1272
    ##  SP_F14 - EMM_F49    0.6975 0.255 127   2.735  0.0355
    ##  SP_F14 - EMM_F48    0.1350 0.255 127   0.529  0.9518
    ##  Control - EMM_F49   3.4408 0.255 127  13.492  <.0001
    ##  Control - EMM_F48   2.8783 0.255 127  11.286  <.0001
    ##  Control - SP_F14    2.7433 0.255 127  10.757  <.0001
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(comp.size2$emmeans$Yeast, 
                               comp.size2$emmeans$DAI,
                               str_trim(comp.size2$emmeans$.group))
colnames(sig.diff.letters) <- c("Yeast", 
                                "DAI",
                                "Letters")

# for plotting with letters from significance test
size2plot <- filterdata2 %>%
  group_by(Yeast, DAI) %>%
  dplyr::summarize(
    colonydiameter = mean(increase, na.rm=TRUE),
    se = sd(increase)/sqrt(3)) %>%
  left_join(sig.diff.letters) 
```

    ## `summarise()` has grouped output by 'Yeast'. You can override using the
    ## `.groups` argument.
    ## Joining with `by = join_by(Yeast, DAI)`

``` r
colonysize2 <- filterdata2 %>%
  ggplot(aes(x = Yeast, y = increase, group = Yeast, fill = Yeast)) +
  geom_boxplot()+
  geom_point(shape = 21, color = "black", position = position_jitterdodge(dodge.width = 0.9))+
   geom_text(data = size2plot, aes(label = Letters, y = colonydiameter+(3*se)), vjust = -0.75, hjust = 2) +
  ylab("Increase in Colony Diameter (mm)") +   
  xlab("Yeast Isolates") +  
    scale_fill_manual(values = cbbPalette) +
  theme(strip.text.x = element_text(size = 18, face= "italic"))+
  theme(axis.title.x = element_text(size =25, face = "bold"), axis.title.y = element_text(size = 25, face = "bold"), axis.text.x=element_text(size=20, angle= 0), axis.text.y=element_text(size=20, angle= 0))+
  theme_classic() +
  theme(legend.position = "none")+
  ggtitle(expression(paste("Impact of Yeast on Growth of ", italic("Methylobacterium platani"), " EMM_B52 (Biological Replicate 2)"))) +
   theme(plot.title = element_text(face = "bold", color = "Blue2", size = 12, hjust = 0.5))+
  facet_wrap(~DAI)
colonysize2
```

![](Splitplate_files/figure-gfm/Increase%20in%20colony%20size%20Plot%20for%20Biological%20Replicate%202-1.png)<!-- -->

### Impact of Yeast on Colony Weight

``` r
#filter data to avoid other data except 
filterdata2_2 <- BioRep2 %>% 
  filter(DAI == "10" )

#stats used for colony weight
cw2=lme(colony_weight~Yeast,data = filterdata2_2, random = ~1|Replication) 
summary(cw2)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: filterdata2_2 
    ##         AIC       BIC   logLik
    ##   -202.1782 -191.4731 107.0891
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##          (Intercept)   Residual
    ## StdDev: 3.821768e-07 0.01895415
    ## 
    ## Fixed effects:  colony_weight ~ Yeast 
    ##                  Value   Std.Error DF   t-value p-value
    ## (Intercept)   0.094950 0.005471593 39 17.353264   0e+00
    ## YeastEMM_F48 -0.035005 0.007738001 39 -4.523778   1e-04
    ## YeastEMM_F49 -0.037715 0.007738001 39 -4.873998   0e+00
    ## YeastSP_F14  -0.030940 0.007738001 39 -3.998449   3e-04
    ##  Correlation: 
    ##              (Intr) YEMM_F48 YEMM_F49
    ## YeastEMM_F48 -0.707                  
    ## YeastEMM_F49 -0.707  0.500           
    ## YeastSP_F14  -0.707  0.500    0.500  
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.1707642 -0.4761489 -0.0944384  0.6141134  1.7017377 
    ## 
    ## Number of Observations: 48
    ## Number of Groups: 6

``` r
Anova(cw2)
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: colony_weight
    ##        Chisq Df Pr(>Chisq)    
    ## Yeast 30.687  3  9.895e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#pairwise comparision
lsmeans.cw2 <- emmeans(cw2, ~ Yeast) 

# contrast with Tukey ajustment 
comp.cw2 <- multcomp::cld(object = lsmeans.cw2, alpha = 0.05, Letters = letters, reversed = TRUE, details = TRUE)

comp.cw2
```

    ## $emmeans
    ##  Yeast   emmean      SE df lower.CL upper.CL .group
    ##  Control 0.0950 0.00547  5   0.0809   0.1090  a    
    ##  SP_F14  0.0640 0.00547  5   0.0499   0.0781   b   
    ##  EMM_F48 0.0599 0.00547  5   0.0459   0.0740   b   
    ##  EMM_F49 0.0572 0.00547  5   0.0432   0.0713   b   
    ## 
    ## Degrees-of-freedom method: containment 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ##  contrast          estimate      SE df t.ratio p.value
    ##  EMM_F48 - EMM_F49  0.00271 0.00774 39   0.350  0.9850
    ##  SP_F14 - EMM_F49   0.00677 0.00774 39   0.876  0.8174
    ##  SP_F14 - EMM_F48   0.00406 0.00774 39   0.525  0.9525
    ##  Control - EMM_F49  0.03771 0.00774 39   4.874  0.0001
    ##  Control - EMM_F48  0.03501 0.00774 39   4.524  0.0003
    ##  Control - SP_F14   0.03094 0.00774 39   3.998  0.0015
    ## 
    ## Degrees-of-freedom method: containment 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(comp.cw2$emmeans$Yeast, 
                               str_trim(comp.cw2$emmeans$.group))

colnames(sig.diff.letters) <- c("Yeast", "Letters")

# for plotting with letters from significance test
weightplot <- filterdata2_2 %>%
  group_by(Yeast) %>%
  dplyr::summarize(
    colonyweight = mean(colony_weight, na.rm=TRUE),
    se = sd(colony_weight)/sqrt(3)) %>%
  left_join(sig.diff.letters) 
```

    ## Joining with `by = join_by(Yeast)`

``` r
 colonyweightplot2 <- filterdata2_2 %>%
  ggplot(aes(x = Yeast, y = colony_weight, group = Yeast, fill = Yeast)) +
  geom_boxplot()+
  geom_point(shape = 21, color = "black", position = position_jitterdodge(dodge.width = 0.9))+
  geom_text(data = weightplot, aes(label = Letters, y = colonyweight+(2*se)), vjust = -0.5, hjust = 2) +
  ylab("Colony Weight (gram)") +   
  xlab("Yeast Isolates") +  
  scale_fill_manual(values = cbbPalette) +
  theme(strip.text.x = element_text(size = 18, face= "italic"))+
  theme(axis.title.x = element_text(size =25, face = "bold"), axis.title.y = element_text(size = 25, face = "bold"), axis.text.x=element_text(size=20, angle= 0), axis.text.y=element_text(size=20, angle= 0))+
  theme_classic() +
  theme(legend.position = "none")
colonyweightplot2
```

![](Splitplate_files/figure-gfm/Colony%20weight%20plot%20for%20Biological%20Replicate%202-1.png)<!-- -->

### Impact of Yeast on Colony Forming Unit (CFU)

``` r
#stats used for colony forming unit
cfu2=lme(log10(CFU)~Yeast,data = filterdata2_2, random = ~1|Replication, na.action = na.omit) 
summary(cfu2)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: filterdata2_2 
    ##        AIC      BIC    logLik
    ##   21.69873 22.17537 -4.849363
    ## 
    ## Random effects:
    ##  Formula: ~1 | Replication
    ##         (Intercept)  Residual
    ## StdDev:   0.2355652 0.2863152
    ## 
    ## Fixed effects:  log10(CFU) ~ Yeast 
    ##                  Value Std.Error DF  t-value p-value
    ## (Intercept)  10.114355 0.2087090  5 48.46153  0.0000
    ## YeastEMM_F48 -1.561742 0.2764507  5 -5.64926  0.0024
    ## YeastEMM_F49 -1.702838 0.2764507  5 -6.15965  0.0016
    ## YeastSP_F14  -0.582231 0.2337754  5 -2.49056  0.0551
    ##  Correlation: 
    ##              (Intr) YEMM_F48 YEMM_F49
    ## YeastEMM_F48 -0.594                  
    ## YeastEMM_F49 -0.594  0.642           
    ## YeastSP_F14  -0.560  0.423    0.423  
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -0.9023658 -0.6637413 -0.1819320  0.5092896  1.3785294 
    ## 
    ## Number of Observations: 12
    ## Number of Groups: 4

``` r
Anova(cfu2)
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: log10(CFU)
    ##        Chisq Df Pr(>Chisq)    
    ## Yeast 43.202  3   2.23e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#pairwise comparision
lsmeans.cfu2 <- emmeans(cfu2, ~ Yeast) ###interaction added by Dr. Steury
#emmeans(cfu,pairwise~Yeast) ##added by Dr. steury
# estimate lsmeans 
comp.cfu2 <- multcomp::cld(object = lsmeans.cfu2, alpha = 0.05, Letters = letters, reversed = TRUE, details = TRUE)
# contrast with Tukey ajustment 
comp.cfu2
```

    ## $emmeans
    ##  Yeast   emmean    SE df lower.CL upper.CL .group
    ##  Control  10.11 0.209  3     9.45    10.78  a    
    ##  SP_F14    9.53 0.209  3     8.87    10.20  ab   
    ##  EMM_F48   8.55 0.227  3     7.83     9.27   bc  
    ##  EMM_F49   8.41 0.227  3     7.69     9.13    c  
    ## 
    ## Degrees-of-freedom method: containment 
    ## Results are given on the log10 (not the response) scale. 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ##  contrast          estimate    SE df t.ratio p.value
    ##  EMM_F48 - EMM_F49    0.141 0.234  5   0.604  0.9267
    ##  SP_F14 - EMM_F49     1.121 0.276  5   4.054  0.0352
    ##  SP_F14 - EMM_F48     0.980 0.276  5   3.543  0.0579
    ##  Control - EMM_F49    1.703 0.276  5   6.160  0.0062
    ##  Control - EMM_F48    1.562 0.276  5   5.649  0.0090
    ##  Control - SP_F14     0.582 0.234  5   2.491  0.1762
    ## 
    ## Degrees-of-freedom method: containment 
    ## Results are given on the log10 (not the response) scale. 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(comp.cfu2$emmeans$Yeast, 
                               str_trim(comp.cfu2$emmeans$.group))

colnames(sig.diff.letters) <- c("Yeast", "Letters")

# for plotting with letters from significance test
cfu2plot <- filterdata2_2 %>%
  na.omit() %>% 
  group_by(Yeast) %>%
  dplyr::summarize(
    cfu2 = mean(log10(CFU), na.rm=TRUE),
    se = sd(log10(CFU))/sqrt(3)) %>%
  left_join(sig.diff.letters) 
```

    ## Joining with `by = join_by(Yeast)`

``` r
CFU2PLOT <- filterdata2_2 %>%
  ggplot(aes(x = Yeast, y = log10(CFU), group = Yeast, fill = Yeast)) +
  geom_boxplot()+
  geom_point(shape = 21, color = "black", position = position_jitterdodge(dodge.width = 0.9))+
  geom_text(data = cfu2plot, aes(label = Letters, y = cfu2+(4*se)), vjust = -0.1, hjust = 2) +
  ylab("log of Colony Forming Unit per ml") +   
  xlab("Yeast Isolates") +  
  scale_fill_manual(values = cbbPalette) +
  theme(strip.text.x = element_text(size = 18, face= "italic"))+
  theme(axis.title.x = element_text(size =25, face = "bold"), axis.title.y = element_text(size = 25, face = "bold"), axis.text.x=element_text(size=10, angle= 0), axis.text.y=element_text(size=20, angle= 0))+
  theme_classic() +
  theme(legend.position = "none")
CFU2PLOT
```

    ## Warning: Removed 36 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Splitplate_files/figure-gfm/Colony%20Forming%20Unit(CFU)%20Plot%20for%20Biological%20Replicate%202-1.png)<!-- -->

``` r
#combining similar dataset plots into one
cw2andcfu2 <- ggarrange(colonyweightplot2, CFU2PLOT, nrow = 1, ncol = 2, common.legend = F)
```

    ## Warning: Removed 36 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
cw2andcfu2 <- annotate_figure(cw2andcfu2, top = text_grob(expression(paste("Impact of Yeast on Growth of ", italic("Methylobacterium platani"), " EMM_B52 (Biological Replicate 2)")), color = "Blue2", size = 12, face = "bold", hjust = 0.5))
cw2andcfu2
```

![](Splitplate_files/figure-gfm/CFU%20and%20Colony%20weight%20Combined%20Plot%20for%20Biological%20Replicate%202-1.png)<!-- -->

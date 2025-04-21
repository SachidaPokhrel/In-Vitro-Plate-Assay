# 🧫 In-Vitro-Plate-Assay
The interaction between 16 different yeast isolates belonging to various classes of basidiomycetes yeasts and six different bacterial isolates, isolated from various hosts, was studied using a specialized in-vitro plate assay. This assay was designed to identify different types of interactions, including contact-dependent, distance-dependent (via agar-diffusible compounds), and distance-independent (via volatile compounds) interactions. To assess that following type of assays were conducted.  

## Co-Culture Plate Assay
Each plate contained eight 3 μl drops of yeast and bacteria arranged in diagonal rows positioned at increasing distances inside a 100 mm petri dish containing YePD agar medium, forming the shape of a “V”. At the apex, the microbe drops were placed in direct contact to determine any contact-dependent interactions. At the same time, the distance between colonies increased as the drops diverged to determine contact-independent interactions, such as diffusible compounds or volatile interactions. 
I have different dataset for each bacteria. Each dataset is in the following data structure.
 
|  Column Name               | Type        | Description                                                            |
|----------------------------|-------------|------------------------------------------------------------------------|
| **Yeast**                  | Categorical | Yeast isolate used with the bacteria                                   |
| **Class**                  | Categorical | Different Class of Yeasts used                                         |
| **Replication**            | Categorical | Replicate of the treatment (3 replicates used and identified by number)|
| **ColonyDiameter**(mm)     | Continuous  | Measured diameter of bacterial colony (mm)                             |
| **DAI**                    | Categorical | Days after inoculation                                                 |
| **distance_to_yeast**(mm)  | Categorical | Distance between yeast and bacterial colony on the plate               |

***Each row = one colony measured*** 

Each replicate has 8 bacterial colonies measured present at different distances from 8 colonies of yeast.

## In-vitro Split Plate Assay
This was done for the confirmation of the volatile interactions. Since the previous co-culture experiment could not rule out the role of volatiles versus agar diffusible compounds, we repeated selected interactions inside 100 mm Petri dishes with a bifurcation, which physically separated the yeast from the bacterium, thus only allowing interactions via volatile metabolites. Two colonies were inoculated in one compartment and other compartment is mass streaked with the yeast isolate or left uninoculated (in case of Control). Since the *Methylobacterium platanii* EMM_B52 had a prominent impact due to yeast volatiles we moved forward with *Methylobacterium platanii* EMM_B52 that faced most prominent impact in the presence of yeast.
Two biological replicate was conducted for the experiment and each replicate has different dataset. Each dataset is in the following data structure. 

| Column Name                 | Type        | Description                                                        |
|-----------------------------|-------------|--------------------------------------------------------------------|
| **Yeast**                   | Categorical | Yeast isolate used in combination with EMM_B52                     |
| **Replication**             | Categorical | Replicate of the treatment (3 replicates used and identified by number) |
| **DAI**                     | Categorical | Days after inoculation                                             |
| **ColonyDiameter**          | Continuous  | Measured diameter of bacterial colony (mm)                                  |
| **ColonyWeight**            | Continuous  | Weight of bacterial colony (gram)                                    |
| **CFU**                     | Continuous  | Colony Forming Units per ml                     |

***Each row = one colony measured***

Each replicate includes 2 bacterial colonies measured inoculated in the Split-Plate


## **📊 Statistical Analysis Workflow**
The analyses were conducted in R, using tidyverse-based data wrangling and statistical modeling tools:

- **Data Cleaning and Processing:** 
    - [dplyr::filter()](https://cran.r-project.org/web/packages/dplyr/index.html) for cleaning and filtering data
    - [tidyverse (version 2.0.0)](https://github.com/tidyverse/tidyverse/releases/tag/v2.0.0) for %>% function

- **Data visualtization:** 
  - [cbbPallete](https://ghurault.github.io/HuraultMisc/reference/cbbPalette.html) was used for Color blind palette
  - [ggplot2(version 3.5.1)](https://cloud.r-project.org/web/packages/ggplot2/index.html) and [ggpubr(version 0.6.0)](https://cran.r-project.org/web/packages/ggpubr/index.html) for generating boxplot and line graph visualization

- **Modeling:** 
  - lm() from baseR for linear model
  - [nlme(version 3.1-168)](https://cran.r-project.org/web/packages/nlme/index.html) for linear-mixed effect model 

- **Statistical Testing:** 
  - [rstatix(version 0.7.2)](https://cran.r-project.org/web/packages/rstatix/index.html)
  - [car(version 3.1-3)](https://cran.r-project.org/web/packages/car/index.html) for running ANOVA on the model. 
  - [emmeans(version 1.10.7)](https://cran.r-project.org/web/packages/emmeans/index.html) for running post-hoc tests
  - [multcomp(version 1.4-28)](https://cran.r-project.org/web/packages/multcomp/index.html) for pairwise comparison and 
  - [multcompView(version 0.1-10)](https://cran.r-project.org/web/packages/multcompView/index.html) generate significant letters


📎 Citation
Please cite this repository if you use or build upon this work. For details on data collection and experimental design, refer to the methods section in the accompanying publication (forthcoming).


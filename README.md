# 🧫 In-Vitro-Plate-Assay
The interaction between 16 different yeast isolates belonging to various classes of basidiomycetes yeasts and six different bacterial isolates, isolated from various hosts, was studied using a specialized in-vitro plate assay. This assay was designed to identify different types of interactions, including contact-dependent, distance-dependent (via agar-diffusible compounds), and distance-independent (via volatile compounds) interactions. To assess that following type of assays were conducted.  

## Co-Culture Plate Assay
 Each plate contained eight 3 μl drops of yeast and bacteria arranged in diagonal rows positioned at increasing distances inside a 100 mm petri dish containing YePD agar medium, forming the shape of a “V”. At the apex, the microbe drops were placed in direct contact to determine any contact-dependent interactions. At the same time, the distance between colonies increased as the drops diverged to determine contact-independent interactions, such as diffusible compounds or volatile interactions. 
 
 | Column Name                 | Type        | Description                                                        |
|----------------------------|-------------|--------------------------------------------------------------------|
| **Yeast**                  | Categorical | Combination of yeast and bacteria tested                          |
| **Class**                  | Categorical | Different Class of Yeasts used
| **Replication**            | Categorical | Biological replicate identifier                                    |
| **ColonyDiameter**(mm)     | Continuous  | Measured diameter of bacterial colony (mm)                         |
| **DAI**                    | Categorical | Days after inoculation in the plate                                |
| **distance_to_yeast**(mm)  | Categorical | Distance between yeast and bacterial colony on the plate          |

**Each row = one colony** 
Each replicate has 8 colonies measured at different distances

## In-vitro Split Plate Assay
This was done for the confirmation of the volatile interactions. Since the previous co-culture experiment could not rule out the role of volatiles versus agar diffusible compounds, we repeated selected interactions inside 100 mm Petri dishes with a bifurcation, which physically separated the yeast from the bacterium, thus only allowing interactions via volatile metabolites. Since the Methylobacterium platani had a prominent impact due to yeast volatiles we moved forward with Methylobacterum with its most interaction with the yeast.

| Column Name                 | Type        | Description                                                        |
|----------------------------|-------------|--------------------------------------------------------------------|
| `Treatment`                | Categorical | Yeast-bacteria treatment combination                               |
| `Replication`              | Categorical | Biological replicate identifier                                    |
| `DAI`                      | Categorical | Days after inoculation                                             |
| `ColonyDiameter`           | Continuous  | Diameter of bacterial colony (mm)                                  |
| `ColonyWeight`             | Continuous  | Weight of bacterial colony (mg)                                    |
| `CFU`                      | Continuous  | Colony Forming Units, bacterial count estimate                     |
| ⬇️ **Each row = one colony** |             | Each replicate includes 2 colonies measured                        |


## 📊 Statistical Analysis Workflow
The analyses were conducted in R, using tidyverse-based data wrangling and statistical modeling tools:



Data Cleaning: dplyr, tidyr

Statistical Testing: rstatix for ANOVA and post hoc tests

Modeling: lm, emmeans for marginal means

Visualization: ggplot2, ggpubr

Post Hoc Interpretation: multcompView for compact letter displays

Each .Rmd file is reproducible and includes:

Data import and cleaning

Summary statistics

ANOVA models and assumptions

Post hoc comparisons

Figures and tables for publication

🔁 Reproducibility
All Excel files are converted to .csv for compatibility and ease of import.

- Code is fully annotated for clarity.

Versioning is maintained via Git.

Raw data, metadata, and analysis scripts will be archived and made public post-publication.

📎 Citation and Acknowledgements
Please cite this repository if you use or build upon this work. For details on data collection and experimental design, refer to the methods section in the accompanying publication (forthcoming).


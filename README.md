# Introduction

> *"We have become, by the power of a glorious evolutionary accident called intelligence, the stewards of life's continuity on earth. We did not ask for this role, but we cannot abjure it. We may not be suited to it, but here we are."*
>
> --- Stephen J. Gould

# Modern ecological problems

**This is a course about community ecology and not so much about population ecology.** Community ecology underpins the vast fields of biodiversity and biogeography, and concerns spatial scales from squares of meters to all of Earth. We can look at historical, contemporary, and future processes that have been implicated in shaping the distribution of life on our planet.

Community ecologists tend to analyse how multiple environmental factors act as drivers that influence the distribution of tens or hundreds of species. These data tend to often be messy (not in the sense of untidy data as per the 'tidyverse' definition of tidy data, but it can be that too!) and statistical considerations need to be understood within the context of the data available to us. This translates to errors of measurement and errors due to extreme values, the presence of a few very rare or very abundant species, autocorrelated residuals (due to repeated sampling, for example), colinearity, etc. These challenges make to application of 'basic' statistical approaches problematic, and a new branch of inferential and exploratory statistical needs to be followed. These approaches involve techniques that allow us to work with all the data at once, and because it can simultaneously analyse all the variables (multiple environmental drivers acting on multiple species at multiple places and across multiple times), this group of statistics is called 'multivariate statistics.' There are two main groups of multivariate statistics: 'classifications' and 'ordinations.' Classification generally concerns placing samples (species or environments) into groups or hierarchies of groups, while ordination is best suited for analyses that involve arranging samples along gradients. Often they complement each other, but we shall see later that each approach has its own strengths. Irrespective of the analysis, the data share a few characteristics.

# Content

## Week 1

### Introduction

-   **Lecture 1** [Topic 1: 01 Ecological data](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/01-ecological_data.ipynb); [Topic 1: 02 Data properties](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/02-data_properties.ipynb); [Topic 1: 03 What do we do with the data?](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/03-doing_data.ipynb); [Topic 1: 04 Exploring the data](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/04-exploring_data.ipynb); 00 Where do ecological data come from?; 00 What do we do with these data?

### Biodiversity

-   **Lecture 2** [Topic 2: 05 Biodiversity](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/05-biodiversity.ipynb)
-   **Lecture 3** 05 Biodiversity -- Continue

## Week 2

-   **Lecture 4** Topic 2: 05 Biodiversity -- Continue
-   **Self study** Topic 3: 06 What structures biodiversity?
-   **Self study** Topic 3: 07 Historical, neutral, and niche theories

### Matrices

-   **Lecture 5** [Topic 4: 08 Environmental distance](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/08-environmental_distance.ipynb)
-   **Lecture 6** [Topic 5: 09 Species dissimilarity](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/09-species_dissimilarity.ipynb)
-   **Self study** [Topic 3: 10 Deep dive into gradients](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/10-deep_dive_into_gradients.ipynb)

## Week 3

-   **Lecture 7** [Topic 6: 11 Correlations and associations](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/11-correlations_and_associations.ipynb)

### Ordinations

-   **Lecture 8** Topic 7: 12 Introduction to ordination
-   **Lecture 9** Topic 8: 13 Principal components analysis (PCA)

## Week 4

-   **Lecture 10** Topic 8: 13 Principal components analysis (PCA) -- Continue
-   **Lecture 11** Topic 9: 14 Correspondence analysis (CA)
-   **Lecture 12** Topic 10: 15 Principal coordinates analysis (PCoA)

## Week 5

-   **Lecture 13** Topic 11: 16 non-Metric multi-dimensional scaling (nMDS)
-   **Lecture 14** Topic 12: 17 Redundancy analysis (RDA)
-   **Lecture 15** Topic 13: 18 Canonical correspondence analysis (CCA)

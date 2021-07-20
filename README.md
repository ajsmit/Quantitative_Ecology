# Introduction

> *"We have become, by the power of a glorious evolutionary accident called intelligence, the stewards of life's continuity on earth. We did not ask for this role, but we cannot abjure it. We may not be suited to it, but here we are."*
>
> --- Stephen J. Gould

# Modern ecological problems

**This is a course about community ecology and not so much about population ecology.** Community ecology underpins the vast fields of biodiversity and biogeography, and concerns spatial scales from squares of meters to all of Earth. We can look at historical, contemporary, and future processes that have been implicated in shaping the distribution of life on our planet.

Community ecologists tend to analyse how multiple environmental factors act as drivers that influence the distribution of tens or hundreds of species. These data tend to often be messy (not in the sense of untidy data as per the 'tidyverse' definition of tidy data, but it can be that too!) and statistical considerations need to be understood within the context of the data available to us. This translates to errors of measurement and errors due to extreme values, the presence of a few very rare or very abundant species, autocorrelated residuals (due to repeated sampling, for example), colinearity, etc. These challenges make to application of 'basic' statistical approaches problematic, and a new branch of inferential and exploratory statistical needs to be followed. These approaches involve techniques that allow us to work with all the data at once, and because it can simultaneously analyse all the variables (multiple environmental drivers acting on multiple species at multiple places and across multiple times), this group of statistics is called 'multivariate statistics.' There are two main groups of multivariate statistics: 'classifications' and 'ordinations.' Classification generally concerns placing samples (species or environments) into groups or hierarchies of groups, while ordination is best suited for analyses that involve arranging samples along gradients. Often they complement each other, but we shall see later that each approach has its own strengths. Irrespective of the analysis, the data share a few characteristics.

These multivariate datasets have far more information in them than can de detected by the human eye and univariate statistics.

![More than meets the eye](Resources/more_than_meets_the_eye.jpeg)

# Supplementary Content

These links point to online resources such as datasets and R scripts in support of the video and PDF lecture material. It is essential that you work through these examples and workflows.

## Week 1

### Introduction

These materials are incomplete...

-   **[Lecture 1]** Topic 1: 00 Introduction

    -   [01 Ecological data](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/01-ecological_data.ipynb)

    -   [02 Data properties](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/02-data_properties.ipynb)

    -   [03 What do we do with the data?](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/03-doing_data.ipynb)

    -   [04 Exploring the data](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/04-exploring_data.ipynb)

### Biodiversity

-   **[Lecture 2]** Topic 2: [05 Biodiversity](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/05-biodiversity.ipynb)
-   **[Lecture 3]** Topic 2: 05 Biodiversity (continue)

## Week 2

-   **[Lecture 4]** Topic 2: 05 Biodiversity (continue)
-   **[Self study]** Topic 3: 06 What structures biodiversity?
-   **[Self study]** Topic 3: 07 Historical, neutral, and niche theories

### Matrices

-   **[Lecture 5]** Topic 4: [08 Environmental distance](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/08-environmental_distance.ipynb)
-   **[Lecture 6]** Topic 5: [09 Species dissimilarity](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/09-species_dissimilarity.ipynb)
-   **[Self study]** Topic 3 revisited: [10 Deep dive into gradients](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/10-deep_dive_into_gradients.ipynb)

## Week 3

-   **[Self study]** Topic 6: [11 Correlations and associations](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/11-correlations_and_associations.ipynb) (lecture 7 slot unused)

### Ordinations

-   **[Lecture 8]** Topic 7: 12 Introduction to ordination (theory only; refer to PDF slides and lecture video)
-   **[Lecture 9]** Topic 8: [13 Principal component analysis (PCA)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/13-PCA.ipynb)

## Week 4

-   **[Self study]** Topic 8: 13 Principal component analysis (PCA) (continue; lecture slot 10 unused)
-   **[Example]** Topic 8: [13 PCA of World Health Organization data on progress towards attaining SDGs](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/13-PCA-SDG-example.ipynb)
-   **[Lecture 11]** Topic 9: [14 Correspondence analysis (CA)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/14-CA.ipynb)
-   **[Lecture 12]** Topic 10: [15 Principal coordinate analysis (PCoA)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/15-PCoA.ipynb)

## Week 5

-   **[Lecture 13]** Topic 11: [16 non-Metric multi-dimensional scaling (nMDS)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/16-nMDS.ipynb)
-   **[Lecture 14]** Topic 12: [17 Redundancy analysis (RDA)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/17-RDA.ipynb)
-   **[Lecture 15]** Topic 13: 18 Canonical correspondence analysis (CCA)

## Week 6

-   **[Lecture 16]** Topic 14: 19 Cluster analysis

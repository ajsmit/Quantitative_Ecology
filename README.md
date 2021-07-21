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

-   **[Lecture 1]** Topic 1: Introduction

    -   [Ecological data](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_1-ecological_data.ipynb)

    -   [Data properties](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_1-data_properties.ipynb)

    -   [What do we do with the data?](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_1-doing_data.ipynb)

    -   [Exploring the data](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_1-exploring_data.ipynb)

### Biodiversity

-   **[Lecture 2]** Topic 2: [Biodiversity](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_2-biodiversity.ipynb)
-   **[Lecture 3]** Topic 2: Biodiversity (continue)

## Week 2

-   **[Lecture 4]** Topic 2: Biodiversity (continue)
-   **[Self study]** Topic 3: What structures biodiversity?
-   **[Self study]** Topic 3: Historical, neutral, and niche theories

### Matrices

-   **[Lecture 5]** Topic 4: [Environmental distance](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_4-environmental_distance.ipynb)
-   **[Lecture 6]** Topic 5: [Species dissimilarity](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_5-species_dissimilarity.ipynb)
-   **[Self study]** Topic 3 revisited: [Deep dive into gradients](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_3-deep_dive_into_gradients.ipynb)

## Week 3

-   **[Self study]** Topic 6: [Correlations and associations](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_6-correlations_and_associations.ipynb) (lecture 7 slot unused)

### Ordinations

-   **[Lecture 8]** Topic 7: Introduction to ordination (theory only; refer to PDF slides and lecture video)
-   **[Lecture 9]** Topic 8: [Principal component analysis (PCA)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_8-PCA.ipynb)

## Week 4

-   **[Self study]** Topic 8: Principal component analysis (PCA) (continue; lecture slot 10 unused)
-   **[Self study]** Topic 8: [PCA of World Health Organization data on progress towards attaining SDGs](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_8-PCA-SDG-example.ipynb)
-   **[Lecture 11]** Topic 9: [Correspondence analysis (CA)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_9-CA.ipynb)
-   **[Lecture 12]** Topic 10: [Principal coordinate analysis (PCoA)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_10-PCoA.ipynb)

## Week 5

-   **[Lecture 13]** Topic 11: [non-Metric multi-dimensional scaling (nMDS)](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_11-nMDS.ipynb)
-   **[Lecture 14]** Topic 12: [Constrained ordination](https://github.com/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Topic_12-Constrained_ordination.ipynb)
-   **[Self study]** All topics: Review and assignments (lecture slot 15 unused)

## Week 6

-   **[Lecture 16]** Topic 13: Cluster analysis

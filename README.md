# Introduction

> *"We have become, by the power of a glorious evolutionary accident called intelligence, the stewards of life's continuity on earth. We did not ask for this role, but we cannot abjure it. We may not be suited to it, but here we are."*
>
> --- Stephen J. Gould

## Modern ecological problems

**This is a course about community ecology and not so much about population ecology.** Community ecology underpins the vast fields of biodiversity and biogeography, and concerns spatial scales from squares of meters to all of Earth. We can look at historical, contemporary, and future processes that have been implicated in shaping the distribution of life on our planet.

Community ecologists tend to analyse how multiple environmental factors act as drivers that influence the distribution of tens or hundreds of species. These data tend to often be messy (not in the sense of untidy data as per the 'tidyverse' definition of tidy data, but it can be that too!) and statistical considerations need to be understood within the context of the data available to us. This translates to errors of measurement and errors due to extreme values, the presence of a few very rare or very abundant species, autocorrelated residuals (due to repeated sampling, for example), colinearity, etc. These challenges make to application of 'basic' statistical approaches problematic, and a new branch of inferential and exploratory statistical needs to be followed. These approaches involve techniques that allow us to work with all the data at once, and because it can simultaneously analyse all the variables (multiple environmental drivers acting on multiple species at multiple places and across multiple times), this group of statistics is called 'multivariate statistics.' There are two main groups of multivariate statistics: 'classifications' and 'ordinations.' Classification generally concerns placing samples (species or environments) into groups or hierarchies of groups, while ordination is best suited for analyses that involve arranging samples along gradients. Often they complement each other, but we shall see later that each approach has its own strengths. Irrespective of the analysis, the data share a few characteristics.

Biodiversity: patterns and processes...

# Content

## Introduction

-   [Ecological data](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/ecological_data.ipynb) -- incomplete

-   [Data properties](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/data_properties.ipynb) -- incomplete

-   [What do we do with the data?](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/doing_data.ipynb) -- incomplete

-   [Exploring the data](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/exploring_data.ipynb) -- incomplete

-   Where do ecological data come from?

-   What do we do with these data?

## Biodiversity

-   [Biodiversity](https://nbviewer.jupyter.org/github/ajsmit/Quantitative_Ecology/blob/main/jupyter_lab/Beta-diversity.ipynb)

-   What structures biodiversity: patterns in space (and time)?

-   Historical, neutral, and niche theories

## Matrices

-   Distance matrices: environmental variables

-   Dissimilarity matrices: species

-   Correlations and associations

## Ordinations

-   Principal components analysis (PCA)

-   Correspondence analysis (CA)

-   Principal coordinates analysis (PCoA)

-   non-Metric multidimension scaling (nMDS)

-   Redundancy analysis (RDA)

-   Canonical correspondence analysis (CCA)

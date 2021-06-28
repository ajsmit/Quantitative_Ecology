# Pairwise association matrices

Although we typically start our forrays into data exploration using sites × species and sites × environment tables, the formal statistical analyses usually require 'pairwise association matrices.' Such matrices are symmetrical (or lower triangle) square matrices (i.e. $n \times n$). These matrices tell us about how related any sample is with any other sample in our pool of samples (i.e. relatedness amongst rows with respect to whatever populates the columns, be they species information of environmental information).

Let us consider various kinds of distance matrices under the following headings.

## Association (dissimilarity and similarity)

Two samples with similar species composition are ecologically similar, while two samples that share few species are ecologically distant. In Figure 4.1, below, the data displayed in Figure 2.1 have been converted into a dissimilarity distance matrix of dimension $58 \times 58$. The are a few things of interest in this matrix:

* The distance matrix is square and therefore symmetrical. In other words, there are as many rows as there are columns, and this number corresponds to the number of samples in our sites × species matrix.

* The cells of the diagonal running from top-left to bottom-right contain zeros, showing rather obviously that there is no difference between a sample and the sample itself. 

* The 'upper triangle' above the diagonal is an inversion of the 'lower triangle' below the diagonal; because they are identical in terms of the pairwise relationships that they encode, distance matrices are sometimes represented simply by the lower triangular matrix.

* These matrices contain ecological information. For example, between samples (here each of 58 × 50-km long coastal sections) that are geographically close together, the dissimilarity will is generally low (i.e. the samples are similar in their species composition), while the further sites are removed from each other, the greater the dissimilarity will be. (*Note: not all samples are not always related to each other as a function of distance --- this is a characteristic of the data used for this particular example analysis, so be aware of the context when interpreting distance matrices.*)

* All information about the particular species present within each sample are now gone since that information has been collapsed to a dissimilarity measure.

## Dependence matrices







# Preliminaries {#prelim}

> "In the beginning, the universe was created. This has made a lot of people very angry and been widely regarded as a bad move."
>
> --- Douglas Adams
  
> "The history of life thus consists of long periods of boredom interrupted occasionally by panic."
> 
> --- Elizabeth Kolbert, The Sixth Extinction



## Venue, date and time

Quantitative Ecology is scheduled to replace Plant Ecophysiology and will run between May 14th andJune 29th.

This workshop will take place between **May 14th** and **June 29th**, from **9:00--16:00** on one day in the week. There will also be a field component, where you will be taught about ecological field sampling in marine and terrestrial environments; this will also offer an opportunity to collect real data, which we will then analyse using the multivariate methods used in this course.

## Course outline

### Week 1 -- In the Beginning {-}
  
### Week 2 -- Show and tell {-}

### Week 3 -- Going deeper {-}

### Week 4 -- The Enlightened Researcher {-}

### Week 5 -- The world is yours {-}

## About this Workshop

The aim of this five-day introductory workshop is to guide you through the basics of using R via RStudio for analysis of environmental and biological data. It is ideal for people new to R or who have limited experience. This workshop is not comprehensive, but is necessarily selective. We are not hardcore statisticians, but rather ecologists who have an interest in statistics, and use R frequently. Our emphasis is thus on the steps required to analyse and visualise data in R, rather than focusing on the statistical theory.

The workshop is laid out so it begins simply and slowly to impart the basics of using R. It then gathers pace, so that by the end we are doing intermediate level analyses. Day 1 is concerned with becoming familiar with getting data into R, doing some simple descriptive statistics, data manipulation and visualisation. Day 2 takes a more in depth look at manipulating and visualising data. Day 3 focuses on creating maps. Day 4 deals with the fundamentals of reproducible research. Day 5 allows one to utilise all of the skills learned throughout the week by creating a final project. The workshop is case-study driven, using data and examples primarily from our background in the marine sciences and real life situations. There is no homework but there are in class assignments.

Don't worry if you feel overwhelmed and do not follow everything at any time during the Workshop; that is totally natural with learning a new and powerful program. Remember that you have the notes and material to go through the exercises later at your own pace; we will also be walking the room during sessions and breaks so that we can answer questions one on one. We hope that this Workshop gives you the confidence to start incorporating R into your daily workflow, and if you are already a user, we hope that it will expose you to some new ways of doing things.

Finally, bear in mind that we are self-taught when it comes to R. Our methods will work, but you will learn as you gain more experience with programming that there are many ways to get the right answer or to accomplish the same task.

## Why do we use R for this?

Please refer to the [Intro R Workshop: Data Manipulation, Analysis and Graphing](https://robwschlegel.github.io/Intro_R_Workshop/) for why we feel strongly that you use R [@R2017] for the analyses that we will perform here. All of the reasons provided there are valid here, but one reason perhaps more so than others --- R and RStudio promote the principles of *reproducible research*, and in fact make it very easy to implement. We will focus on these principles throughout the workshop, and the assignments will in fact require that you submit a fully functional working script, complete with all the notes, memos, examples, data, executable code, and output that will result from completing the course material.

What other oprions are there for analysing the kinds of data that we will envounter in ecological research? Software packages like the ones you may be familiar with, such as Statistica and SPSS, as well as specilised software such as Primer (**P**lymouth **R**outines **I**n **M**arine **E**cological **R**esearch; equally suited to data of terrestrial origin), are often used to perform many of the analyses we will encounter. They are rather limited with regards to the full scope of modern multivariate statistical methods in use by ecologists today. This is why we prefert to use R as the *engine* within which to do our ecological data analysis. R is used by academic statisticians the world over, as well as those academic ecological statisticians who are responsible for developing the methods used for answering ecological questions. That package is called **`vegan`** [@vegan2017], and it is the one we will use here.

### Some negatives of using `vegan`

Although there are many positives of using **`vegan`**, there are some negatives:

1.	It can have a steep learning curve for those whom do not like statistics or data manipulation, and it does require frequent use to remain familiar with it and to develop advanced skills

2.	Error trapping can be confusing and frustrating

3.	Rudimentary debugging, although there are some packages available to enhance the process

4.	Handles large datasets (100 MB), but can have some trouble with massive datasets (GBs)

5.	Some simple tasks can be tricky to do in R

6.  There are multiple ways of doing the same thing

### The challenge: learning to program in R and vegan

The same challenges as we highlighted in [Intro R Workshop: Data Manipulation, Analysis and Graphing](https://robwschlegel.github.io/Intro_R_Workshop/) apply here, but we still maintain that these are far outweighed by the benefits.

### Installing R and RStudio

We assume that you already have R installed on your computer, as all of you will have already completed the the Intro R Workshop. If you need a refresher, please refer to [Intro R Workshop: Data Manipulation, Analysis and Graphing](https://robwschlegel.github.io/Intro_R_Workshop/) for the installation instructions.

--- Cheers, AJ and Robert

## Resources

Below you can find the source code to some books and other links to websites about R. With some of the technical skills you'll learn in this course you'll be able to download the source code, compile the book on your own computer and arrive at the fully formatted (typeset) copy of the books that you can purchase for lots of money:

* [ggplot2. Elegant Graphics for Data Analysis](https://github.com/hadley/ggplot2-book) --- the R graphics bible
* [A Compendium of Clean Graphs in R. Version 2.0](http://shinyapps.org/apps/RGraphCompendium/index.php) --- using R's base graphics
* [R for Data Science](http://r4ds.had.co.nz/workflow-basics.html) --- data analysis using tidy principles
* [R Markdown](http://rmarkdown.rstudio.com) --- reproducible reports in R
* [bookdown: Authoring Books and Technical Documents with R Markdown](https://bookdown.org/yihui/bookdown) --- writing books in R
* [Shiny](https://shiny.rstudio.com) --- interactive website driven by R

## Style and code conventions

Early on, develop the habit of unambiguous and consistent style and formatting when writing your code, or anything else for that matter. Pay attention to detail and be pedantic. This will benefit your scientific writing in general. Although many R commands rely on precisely formatted statements (code blocks), style can nevertheless to *some extent* have a personal flavour to it. The key is *consistency*. In this book we use certain conventions to improve readability. We use a consistent set of conventions to refer to code, and in particular to typed commands and package names.

  * Package names are shown in a bold font over a grey box, *e.g.* __`tidyr`__.
  * Functions are shown in normal font followed by parentheses and also over a grey box , *e.g.* `plot()`, or `summary()`.
  * Other R objects, such as data, function arguments or variable names are again in normal font over a grey box, but without parentheses, *e.g.* `x` and `apples`.
  * Sometimes we might directly specify the package that contains the function by using two colons, *e.g.* `dplyr::filter()`.
  * Commands entered onto the R command line (console) and the output that is returned will be shown in a code block, which is a light grey background with code font. The commands entered start at the beginning of a line and the output it produces is preceded by `R>`, like so:


```r
rnorm(n = 10, mean = 0, sd = 13)
```

```
R>  [1]   9.655172   4.390506  -1.592293  -4.470382  13.144224   8.625643
R>  [7]   6.026807  -5.712196 -14.437423  -5.994454
```

Consult these resources for more about R code style :

  * [Google's R style guide](https://google.github.io/styleguide/Rguide.xml)
  * [The tidyverse style guide](http://style.tidyverse.org)
  * [Hadley Wickham's advanced R style guide](http://adv-r.had.co.nz/Style.html)

We can also insert maths expressions, like this $f(k) = {n \choose k} p^{k} (1-p)^{n-k}$ or this: $$f(k) = {n \choose k} p^{k} (1-p)^{n-k}$$

## About this document

This document was written in **`bookdown`** and transformed into the 'GitBook' you see here by **`knitr`**, **pandoc** and \LaTeX\ (Figure \@ref(fig:rmarkdown)). All the source code and associated data are available at AJ Smit's [GitHub page](https://github.com/ajsmit/Intro_R_Workshop). You can download the source code and compile this document on your own computer. If you can compile the document yourself you are officially a geek -- welcome to the club! Note that you will need to complete the exercises in the chapter, An R workflow, before this will be possible.



You will notice that this repository uses [GitHub](https://github.com), and you are advised to set up your own repository for R scripts and all your data. We will touch on GitHub and the principles of reproducible research later, and GitHub forms a core ingredient of such a workflow.

The R session information when compiling this book is shown below:

```r
sessionInfo()
```

```
R> R version 3.4.3 (2017-11-30)
R> Platform: x86_64-apple-darwin15.6.0 (64-bit)
R> Running under: macOS Sierra 10.12.6
R> 
R> Matrix products: default
R> BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
R> LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
R> 
R> locale:
R> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
R> 
R> attached base packages:
R> [1] stats     graphics  grDevices utils     datasets  base     
R> 
R> loaded via a namespace (and not attached):
R>  [1] compiler_3.4.3  backports_1.1.2 bookdown_0.5    magrittr_1.5   
R>  [5] rprojroot_1.3-2 tools_3.4.3     htmltools_0.3.6 yaml_2.1.16    
R>  [9] Rcpp_0.12.14    stringi_1.1.6   rmarkdown_1.8   knitr_1.18     
R> [13] methods_3.4.3   stringr_1.2.0   digest_0.6.13   evaluate_0.10.1
```

## Exercise: It which shall not be named
Now that you have heard (and perhaps read) our argument about the merits of using R, let's double down and spend the next hour seeing first-hand why we think this. Please open the file 'data/SACTN_data.csv' in MS Excel. Gasp! Yes I know. After all of that and now we are using MS Excel? But trust us, there is method to this madness. Your mission, should you choose to accept it, is to spend the next hour creating monthly climatologies and plotting them as a line graph. The South African Coastal Temperature Network (SACTN, which will be used several times during this workshop) data are three monthly temperature time series, each about 30 years long. To complete this objective you will need to first split up the three different time series, and then figure out how to create a monthly climatology for each. A monthly climatology is the average temperature for a given month at a given place. So in this instance, because we have three time series, we will want 36 total values comprised of January - December monthly means for each site (if a time series is 30 years long, then a climatological December will be the mean temperature of all of the data within the 30 Decembers for which data are available). Once those values have been calculated, it should be a relatively easy task to plot them as a dot and line graph. Please keep an eye on the time, if you are not done within an hour please stop anyway. Less than a quarter of workshop attendees have completed this task in the past.

After an hour has passed we will take a break. When we return we will see how to complete this task via R as part of 'The New Age' demonstration.

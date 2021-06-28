# Preliminaries {-}

> *"You can know the name of that bird in all the languages of the world, but when you're finished, you'll know absolutely nothing whatever about the bird. You'll only know about humans in different places, and what they call the bird. ... I learned very early the difference between knowing the name of something and knowing something."*
>
> --- Richard Feynman



## Venue, date and time {-}

Quantitative Ecology is scheduled to replace Plant Ecophysiology, and will run between May 14th and June 29th, 2018. This workshop will take place from **9:00--16:00** on one day in each week during this period. There will also be a field component, where you will be taught about ecological field sampling in marine and terrestrial environments; this will also offer an opportunity to collect real data using actual ecological field methods (these will also be covered in the course), which we will then analyse using the multivariate methods used in this workshop.

## Course outline {-}

* Week 1 -- In the Beginning
* Week 2 -- Show and tell
* Week 3 -- Going deeper
* Week 4 -- The Enlightened Researcher
* Week 5 -- The world is yours

## About this Workshop {-}

The aim of this five-day introductory workshop is to guide you through...

## This is biology: why more R coding? {-}

Please refer to the [Intro R Workshop: Data Manipulation, Analysis and Graphing](https://robwschlegel.github.io/Intro_R_Workshop/) for why we feel strongly that you use R [@R2017] for the analyses that we will perform here. All of the reasons provided there are valid here too, but one reason perhaps more so than others --- R and RStudio promote the principles of *reproducible research*, and in fact make it very easy to implement. We will focus on some of these principles throughout the workshop, and the assignments will in fact require that you submit a fully functional working script, complete with all the notes, memos, examples, data, executable code, and output that will result from completing the course material. Full, detailed reproducibility as per the British Ecological Society [guidelines](https://www.britishecologicalsociety.org/wp-content/uploads/2017/12/guide-to-reproducible-code.pdf) is difficult and requires advanced computing skills, but don't panic! We will not go in to excessive amounts of detail about reproducible research at the expense of learning about ecological sampling and data analysis. For an example of fully reproducible research, see a recent publication, [Seaweeds in Two Oceans: Beta-diversity](https://github.com/ajsmit/Seaweeds_in_Two_Oceans).

Perhaps even more fundamental for our purpose (and also in aid of reproducibility) are the skills around data management. A large portion of our 'workflow' will concern getting the data into a format that is suitable for analysis. And an earlier step is just as important --- that is data entry. Sometimes we have the good fortune to enter our data into a 'tidy' format, but frequently we are given a data set that was put together by one person or multiple people, often over many years. These data are not so easy to coerce into a format that can easily by analysed. We will teach you the skills needed to produce a properly formatted data set from scratch, and provide tips for how to transform an existing untidy data set into a tidy one. The philosophy of tidy data is being formalised into a collection of R packages contained within Hadley Wickham's [`tidyverse`](https://www.tidyverse.org/), and this collection of packages will inform the principles behind our data transformations in order to arrive at a tidy data set, which can then be subjected to the quantitative ecological methods (an array of multivariate methods). Thanksfully, many of these data processing principles have also been summarised by the British Ecological Society, and it can be downloaded [here](https://www.britishecologicalsociety.org/wp-content/uploads/2017/06/BES-Data-Guide-2017_web.pdf). For the purpose of this workshop, please make sure that you have access to both of these publications, as we will refer to them frequently, espcially during the first few contact sessions.

What other oprions are there for analysing the kinds of data that we will encounter in ecological research? Software packages like the ones you may be familiar with, such as Statistica and SPSS, as well as specilised software such as Primer (**P**lymouth **R**outines **I**n **M**arine **E**cological **R**esearch; equally suited to data of terrestrial origin), are often used to perform many of the analyses we will encounter. They are rather limited with regards to the full scope of modern multivariate statistical methods in use by ecologists today. This is why we prefer to use R as the *engine* within which to do our ecological data analysis. R is used by academic statisticians the world over, as well as those academic ecological statisticians who are responsible for developing the methods used for answering today's ecological questions. That package is called **`vegan`** [@vegan2017], and it is the one we will use here.

### Some negatives of using `vegan` {-}

Although there are many positives of using **`vegan`**, there are some negatives:

1.	It can have a steep learning curve for those whom do not like statistics or data manipulation, and it does require frequent use to remain familiar with it and to develop advanced skills

2.	Error trapping can be confusing and frustrating

3.	Rudimentary debugging, although there are some packages available to enhance the process

4.	Handles large datasets (100 MB), but can have some trouble with massive datasets (GBs)

5.	Some simple tasks can be tricky to do in R

6.  There are multiple ways of doing the same thing

### The challenge: learning to program in R and vegan {-}

The same challenges as we highlighted in [Intro R Workshop: Data Manipulation, Analysis and Graphing](https://robwschlegel.github.io/Intro_R_Workshop/) apply here, but we still maintain that these are far outweighed by the benefits.

### Installing R and RStudio {-}

We assume that you already have R installed on your computer, as all of you will have already completed the the Intro R Workshop. If you need a refresher, please refer to [Intro R Workshop: Data Manipulation, Analysis and Graphing](https://robwschlegel.github.io/Intro_R_Workshop/) for the installation instructions.

--- Cheers, AJ and Robert

## Resources {-}

### Required reading {-}
A significant amount of self study will be required to complete this module. The content of this course follows the content of the excellent book by Borcard, D., Gillet, F., & Legendre, P, (2011) [Numerical ecology with R](http://www.springer.com/gp/book/9781441979759). There is a hardcopy in the Library.

### Resources about multivatiate ecological methods {-}
To supplement the required reading, please refer to these helpful websites:

* [Ordination Methods: an overview](http://ordination.okstate.edu/overview.htm) --- as the title says, and overview of ordination methods
* [Ordination Methods for Ecologists](http://ordination.okstate.edu/) --- primarily a focus on ordination methods used by ecologists
* [GUide to STatistical Analysis in Microbial Ecology (GUSTA ME)](https://sites.google.com/site/mb3gustame/) --- although about micobiological ecology, the techniques are the same ones that will be used for the study of larger multicellular species

### General resources about R {-}
Below you can find the source code to some books and other links to websites about R. With some of the technical skills you'll learn in this course you'll be able to download the source code, compile the book on your own computer and arrive at the fully formatted (typeset) copy of the books that you can purchase for lots of money:

* [ggplot2. Elegant Graphics for Data Analysis](https://github.com/hadley/ggplot2-book) --- the R graphics bible
* [A Compendium of Clean Graphs in R. Version 2.0](http://shinyapps.org/apps/RGraphCompendium/index.php) --- using R's base graphics
* [R for Data Science](http://r4ds.had.co.nz/workflow-basics.html) --- data analysis using tidy principles
* [R Markdown](http://rmarkdown.rstudio.com) --- reproducible reports in R
* [bookdown: Authoring Books and Technical Documents with R Markdown](https://bookdown.org/yihui/bookdown) --- writing books in R
* [Shiny](https://shiny.rstudio.com) --- interactive website driven by R

## Style and code conventions {-}

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
R>  [1]  0.4245640 -0.5257653 -4.5610590 -2.4811375 -4.3244899 17.7846893
R>  [7] -3.7848909 -1.8796644 -9.6668739  0.2252982
```

Consult these resources for more about R code style :

  * [Google's R style guide](https://google.github.io/styleguide/Rguide.xml)
  * [The tidyverse style guide](http://style.tidyverse.org)
  * [Hadley Wickham's advanced R style guide](http://adv-r.had.co.nz/Style.html)

We can also insert maths expressions, like this $f(k) = {n \choose k} p^{k} (1-p)^{n-k}$ or this: $$f(k) = {n \choose k} p^{k} (1-p)^{n-k}$$

--- 
title: "Quantitative Ecology"
subtitle: 'A primer in ecological field and computer techniques for BCB (Hons) 2018'
author: "AJ Smit"
date: "2018-02-06"
site: "bookdown::bookdown_site"
papersize: A4
fontsize: 10pt
geometry: margin=4cm
language: Australian
always_allow_html: yes
github-repo: "ajsmit/quantitative_ecology"
cover-image: figures/NSFW.jpg
description: 'This is a book about quatitative ecology: theory and practice.'
thanks: Replication files are available on the author's GitHub account
tables: yes
lof: yes
lot: yes
bibliography: [LaTeX/bibliography.bib, LaTeX/packages.bib]
biblio-style: apalike
link-citations: yes
citecolor: green
linkcolor: cyan
urlcolor: cyan
---

# Preface {-}
This workshop concerns the most common classification and ordination methods used by modern-day ecologists interested in extracting information about features, properties, patterns, and processes from ecological data. Our journey will enable insightful ecologists to make inferences about associations between complex species and environmental data, which may support hypothesis-driven inquiries. With "complex data" I mean that we often work with data sets comprised of 10s or 100s of species, and 10s of measured environmental variables that we intend to use to explain where species occur, why they occur there, how they interact with each other, and why they vary seasonally or from year-to-year. These data often don't behave as we would expect from normally-distributed data (i.e. the structural properties of the data are non-Gaussian), and the error structures may be correlated or display some other kind of artefact. We cannot use the usual kinds of statistical methodologies we are used to, and alternative approaches are necessary. 

In this module we will be using the modern approaches available to us, coupled with the numerical abilities offered by today's powerful computers, to ask questions about these patterns, processes, and associations across vast swaths of Earth's surface, or in some instances, about the *entire* surface of Earth.

# Prerequisites {-}

A prerequisite for this course is a basic proficiency in using R [@R2017]. The necessary experience will have been gained from completing the [Intro R Workshop: Data Manipulation, Analysis and Graphing](https://robwschlegel.github.io/Intro_R_Workshop/) Workshop that was part of your BCB Core Honours module (i.e. Biostatistics). You will also need a laptop with R and RStudio installed as per the instructions provided in that workshop. If you do not have a personal laptop, most computers in the 5th floor lab will be correctly set up for this purpose.



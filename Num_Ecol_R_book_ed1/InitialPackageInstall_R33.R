# ============================================================
# Script for first-time users of "Numerical Ecology with R"  #
# by Daniel Borcard, Francois Gillet and Pierre Legendre     #
# Adapted to R 3.3.2 on 29.01.2017                           #
# ============================================================

# This script installs or provides guidelines to install all 
# the packages necessary to run the code provided in the book, 
# but that do not belong to the standard R distribution (steps 2-3). 

# Steps 1 to 3 must be run only once when installing or upgrading R.
# Step 4 is not mandatory.


# 1. Update installed packages
#    -------------------------
update.packages(checkBuilt=TRUE, ask=FALSE)


# 2. Install packages from the main CRAN site
#    ----------------------------------------

install.packages(c("ade4", "adegraphics", "adespatial", "ape", "cocorresp",
  "cluster", "ellipse", "FactoMineR", "FD", "gclus", "labdsv", "MASS",
  "RColorBrewer", "spdep", "tripack", "vegan", "vegan3d"), dependencies=TRUE, 
  type="both")

# Install mvpart and MVPARTwrap that are not available from CRAN anymore:
# On Windows machines, Rtools (3.4 and above) must be installed first. Go to:
# https://cran.r-project.org/bin/windows/Rtools/
# After that: 
install.packages("devtools")
library(devtools)
install_github("cran/mvpart", force=TRUE)
install_github("cran/MVPARTwrap", force=TRUE)


# 3. Install packages from R-Forge.R-project
#    ---------------------------------------

# Package {AEM}
# The command below *may* work:
install.packages("AEM", repos="http://R-Forge.R-project.org", type="both")
# Otherwise, Windows: 
# Go to https://r-forge.r-project.org/R/?group_id=195
# Download the .zip file and install from the R console.


# 4. OPTIONAL (for power users): Install all R packages from Environmetrics,
#    a CRAN Task View for the Analysis of Ecological and Environmental Data
#    See http://cran.r-project.org/web/views/Environmetrics.html
#    ----------------------------------------------------------------------

install.packages("ctv")
library(ctv)
install.views("Environmetrics")

# Other potentially useful CRAN Task Views...
install.views("Cluster")
install.views("Multivariate")
install.views("Spatial")
install.views("MachineLearning")

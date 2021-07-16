# Quantitative Ecology (BCB743)
# Structuring Biodiversity
# Author: AJ Smit
# Date: 29 June 2021


# Load the packages -------------------------------------------------------

library("coenocliner")

# MODEL 1: a simple linear pH gradient, for example
set.seed(2)
M    <- 20                                 # number of species
ming <- 3.5                                # gradient minimum...
maxg <- 7                                  # ...and maximum
locs <- seq(ming, maxg, length = 100)      # gradient locations
opt  <- runif(M, min = ming, max = maxg)   # species optima
tol  <- rep(0.25, M)                       # species tolerances
h    <- ceiling(rlnorm(M, meanlog = 3))    # max abundances
pars <- cbind(opt = opt, tol = tol, h = h) # put in a matrix

# counts (occurrence) are generated using
# random deviates from the specified distribution
mu <- coenocline(locs, responseModel = "gaussian", params = pars,
                 expectation = TRUE)

class(mu)
dim(mu)
head(mu[, 1:6])

matplot(locs, mu, lty = "solid", type = "l", xlab = "pH", ylab = "Abundance")

# counts generated from a poisson model
simp <- coenocline(locs, responseModel = "gaussian", params = pars,
                   countModel = "poisson")

dim(simp)
head(simp[, 1:6])

matplot(locs, simp, lty = "solid", type = "p", pch = 1:10, cex = 0.8,
        xlab = "pH", ylab = "Abundance")

# counts generated from a negative binomial model
simnb <- coenocline(locs, responseModel = "gaussian", params = pars,
                    countModel = "negbin", countParams = list(alpha = 0.5))

matplot(locs, simnb, lty = "solid", type = "p", pch = 1:10, cex = 0.8,
        xlab = "pH", ylab = "Abundance")


# MODEL 2:
A0    <- c(5,4,7,5,9,8) * 10               # max abundance
m     <- c(25,85,10,60,45,60)              # location on gradient of modal abundance
r     <- c(3,3,4,4,6,5) * 10               # species range of occurrence on gradient
alpha <- c(0.1,1,2,4,1.5,1)                # shape parameter
gamma <- c(0.1,1,2,4,0.5,4)                # shape parameter
locs  <- 1:100                             # gradient locations
pars  <- list(m = m, r = r, alpha = alpha,
              gamma = gamma, A0 = A0)      # species parameters, in list form

mu <- coenocline(locs, responseModel = "beta", params = pars, expectation = TRUE)

dim(mu)

matplot(locs, mu, lty = "solid", type = "l", xlab = "Gradient", ylab = "Abundance")

set.seed(10)
N <- 30                                           # number of samples
M <- 20                                           # number of species

# first gradient (e.g. east-west)
ming1 <- 3.5                                      # 1st gradient minimum...
maxg1 <- 7                                        # ...and maximum
loc1 <- seq(ming1, maxg1, length = N)             # 1st gradient locations
opt1 <- runif(M, min = ming1, max = maxg1)        # species optima
tol1 <- rep(0.5, M)                               # species tolerances
h    <- ceiling(rlnorm(M, meanlog = 3))           # max abundances
par1 <- cbind(opt = opt1, tol = tol1, h = h)      # put in a matrix

# second gradient (e.g. north-south)
ming2 <- 1                                        # 2nd gradient minimum...
maxg2 <- 100                                      # ...and maximum
loc2 <- seq(ming2, maxg2, length = N)             # 2nd gradient locations
opt2 <- runif(M, min = ming2, max = maxg2)        # species optima
tol2 <- ceiling(runif(M, min = 5, max = 50))      # species tolerances
par2 <- cbind(opt = opt2, tol = tol2)             # put in a matrix

# last steps...
pars <- list(px = par1, py = par2)                # put parameters into a list
locs <- expand.grid(x = loc1, y = loc2)           # put gradient locations together

dim(locs)

mu2d <- coenocline(locs, responseModel = "gaussian",
                   params = pars, extraParams = list(corr = 0.5),
                   expectation = TRUE)

mu2d
dim(mu2d)

layout(matrix(1:4, ncol = 2))
op <- par(mar = rep(1, 4))
for (i in c(2,8,13,19)) {
  persp(loc1, loc2, matrix(mu2d[, i], ncol = length(loc2)),
        ticktype = "detailed", zlab = "Abundance",
        theta = 45, phi = 30)
}
par(op)
layout(1)

dMCMC
=====

Markov chain Monti Carlo diagnostic (`dMCMC`). A set of functions that create multipanel plots to quickly evaluate MCMC output, and easily convert MCMC chain lists (returned by e.g. rjags) into handy-dandy latex tables. The motor behind these functions are the excellent R packages **xtable** and **coda**. Though, the main goal of dMCMC is to aggregate the most useful information in a single attractive graph or rapidly format a tables for you. Which, if you work with numerous different Bayesian models daily, can be a major waste of time. One innovation that dMCMC adds to the mix is prior v.s. posteriors plots, which can be a crucial consideration, if one wan't to scrutenize choice of priors amonst other things.  

## Installation

Currently there isn't a release on [CRAN](http://cran.r-project.org/),
though there may one day be one. You can still  
download the [zip](https://github.com/MarcoDVisser/choosecolor/zipball/master) 
or [tar ball](https://github.com/MarcoDVisser/choosecolor/tarball/master).
Then decompress and run `R CMD INSTALL` on it, 
or use the **devtools** package to install the development version.

```r
## Make sure your current packages are up to date
update.packages()
## devtools is required
library(devtools)
install_github("dMCMC", "MarcoDVisser")
```

## examples

## Screenshots

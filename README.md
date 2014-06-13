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

If you plan to use dMCMCs ability to find the priors directly in your BUGS/JAGS model file you will also need to 
download and install [r2bugs](https://github.com/dlebauer/r2bugs). This package assist in converting distributions between R and BUGS format. It can be installed as `dMCMC` with;

```r
install_github("r2bugs", "dlebauer")
```

Otherwise you will need to specify the priors yourself (in the correct R format).

## examples

## Screenshots

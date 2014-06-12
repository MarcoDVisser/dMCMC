dMCMC
=====

Markov chain Monti Carlo diagnostic (**dMCMC**). A set of function that create multipanel plots to quickly evaluate MCMC output, and easily convert MCMC chain lists (returned by e.g. rjags) into handy-dandy latex tables. The motor behind these functions are the excellent R packages **xtable** and **coda**. Though, the innovation here is that dMCMC will aggregate the most useful information in a single graph or rapidly format a tables for you. Crucially dMCMC adds prior v.s. posteriors plots to the mix, which can be crucial to compare.  

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

## Screenshots

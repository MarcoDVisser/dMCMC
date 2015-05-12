dMCMC
=====

Markov chain Monte Carlo diagnostic (`dMCMC`). A set of functions that create multipanel plots to quickly evaluate MCMC output, and easily convert MCMC chain lists (returned by e.g. rjags) into handy-dandy latex tables. The motor behind these functions are the excellent R packages **xtable** and **coda**. Though, the main goal of dMCMC is to aggregate the most useful information in a single attractive graph or rapidly format a table for you. Which, if you work with numerous different Bayesian models daily, should prove usefull. One innovation that dMCMC adds to the mix is prior v.s. posteriors plots, which can be a crucial consideration, if one wan't to scrutenize choice of priors amongst other things.  

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

## Examples

```r
 require(dMCMC)
 sink("examp.txt")
     cat( "
     model {
             for (i in 1:N) {
                     x[i] ~ dnorm(mu, tau)
             }
             mu ~ dnorm(0, .0001)
             tau <- pow(sigma, -2)
             sigma ~ dunif(0, 100)
         } ",fill=TRUE)
          sink()
     
     jags <- jags.model('examp.txt',
                       data = list('x' = rnorm(100,2,2),
                                   'N'=100),
                       n.chains = 4,
                       n.adapt = 100)
     
     mysamples <- coda.samples(jags, c('mu', 'tau'),100)
     
     dMCMCplot(mysamples,"mu","examp.txt")
  
```

![](http://i.imgur.com/J60k45r.png)


`dMCMCs` has the ability to find the priors directly in your BUGS/JAGS model file and translate these to their R equivalents. Here for the example above:


```r
findRprior("examp.txt","mu")
```

Which returns the accociated density function and translated parameters:

```
[1] "BUGS prior identified as:  mu~dnorm(0,.0001)"
[1] "BUGS prior translated to: dnorm"
[[1]]
[1] "dnorm"

[[2]]
[[2]][[1]]
[1] 0 1
```


** more to follow **



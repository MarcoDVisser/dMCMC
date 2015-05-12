##' dMCMCplot Create the untimate diagnostic plot for a MCMC chain 
##'
##' Creates a multipanel graphic with diagnostics on parameter samples
##" from a MCMC chain.  See the exampe below for more details.
##' @title Create a multiple diagnostic plot
##' @param mcmcobject a mcmc object (list of samples).
##' @param interest a set of parameters of interest. These should have
##' true priors (e.g. can be directly related to a prior distribution).
##' Lower-level priors, that depend on a hyperprior, will fail
##' as dMCMC does not attempt to understand the hierachical
##' structure of the model file.
##' 
##' @param model the path to the BUGs/JAGS model text file with
##' the priors specified. This model must be the one
##' that generated the MCMC list. 
##' @param Rprior the prior function directly
##' specified by the user (as a standard R
##' probabilty density function e.g. prior=dunif).
##' @param \dots additional parameters (unused). 
##' diagnostic plots on. Names must be equal to the names
##' in the jags/bugs model. If length > 1, one plot is
##' produced per parameter of interest, in this case output to a
##' pdf is advised. Otherwise the function forces par(ask=TRUE).
##' @author Marco D. Visser
##' @examples
##' sink("examp.txt")
##' cat( "
##' model {
##'	for (i in 1:N) {
##'		x[i] ~ dnorm(mu, tau)
##'	}
##'	mu ~ dnorm(0, .0001)
##'	tau <- pow(sigma, -2)
##'	sigma ~ dunif(0, 100)
##'     } ",fill=TRUE)
##'      sink()
##'
##' jags <- jags.model('examp.txt',
##'                   data = list('x' = rnorm(100,2,2),
##'                               'N'=100),
##'                   n.chains = 4,
##'                   n.adapt = 100)
##' 
##' mysamples <- coda.samples(jags, c('mu', 'tau'),100)
##'
##' dMCMCplot(mysamples,"mu","examp.txt")
##'
##'
##'
##' @seealso \code{\link{jags.model}}
##' @return A multipanel plot
##' @concept MCMC diagnostics
##'
##' @export
dMCMCplot <- function(mcmcobject=NULL, interest=NULL,
                      model=NULL,Rprior=NULL,...) {

  opar<-par("mar","bg")
  on.exit(par(opar))

  ## run checks
  if(is.null(mcmcobject)|is.null(interest)) {
    stop("One of the inputs is null. Please check the input objects")}

  chainsamplenames <- attr(mcmcobject[[1]],"dimnames")[[2]]

  if(!is.element(interest,chainsamplenames)) {
    stop(paste(interest,"not found in MCMC samples"))}
 
  
  if(is.null(Rprior)&is.null(model)) {
    stop("Either specify the model OR the Rprior")
   }

  
  ## multiple or singe plots?  
  if(length(interest)>1) {
    stop("multiple parameters not supported yet but soon!") 
      par(ask=TRUE)
  } else{
    
  
  if(is.null(Rprior)){
    Rtranslation <- findRprior(model=model,interest=interest
                               ,silent=TRUE)
    ## Build prior function 
    Rprior <- eval(parse(text=Rtranslation[[1]]))
    rRprior <- eval(parse(text=gsub("d","r",Rtranslation[[1]])))
    param <-  as.numeric(Rtranslation[[2]][[1]])
    paramn <- length(param)
    
  }

  ## find posterior
  postcoln <- grep(interest,colnames(mcmcobject[[1]]))
  posterrange <- range(sapply(mcmcobject,function(X) range(X[,postcoln])))
  
  if(paramn==1){
    priorsamples <- rRprior(2000,param)
    priorrange <- range(priorsamples)
    fpr <- function(x){Rprior(x,param)}
  }
  if(paramn==2){
    priorsamples <- rRprior(2000,param[1],param[2])
    priorrange <- range(priorsamples)
    fpr <- function(x){Rprior(x,param[1],param[2])}
  } else {stop("priors with > 2 parameters not supported")}

  ## Start building the plot
  par(mar=c(0,0,0,0),bg="grey70",mgp = c(0, -1.4, 0),
      lab=c(4,2,3),tck=-.05,las=1,col.axis="white",col.lab='white')
  layoutmat <- matrix(
                 c(3,3,3,3,4,4,4,4,1,1,
                   3,3,3,3,4,4,4,4,1,1,
                   3,3,3,3,4,4,4,4,2,2,
                   3,3,3,3,4,4,4,4,2,2,
                   5,5,5,5,6,6,6,6,7,7,
                   5,5,5,5,6,6,6,6,7,7,
                   5,5,5,5,6,6,6,6,7,7), byrow=TRUE,
                 ncol=10)
  layout(layoutmat)
    
  curve(fpr(x),priorrange[1],priorrange[2],lwd=2,col=rgb(0,0,1,alpha=0.4))
  text(0.6*priorrange[2],0.8*max(fpr(priorsamples)),"Prior",cex=2,
       col='grey40')

  posterhist <- hist(mcmcobject[[1]][,postcoln],main="",freq=FALSE,
                     xlim=priorrange)
  text(0.5*priorrange[2],0.8*max(posterhist$density),"Posterior",cex=2,
       col='grey40')

  ## plot chains
  chaincolors <- colorRampPalette(
                   c(rgb(0,0,1,alpha=0.2),
                     rgb(0,1,0,alpha=0.2),
                     rgb(1,0,0,alpha=0.2)))

  nchains <- length(mcmcobject)
  colorz <- chaincolors(nchains)
  mcpars <- attr(mcmcobject[[1]],"mcpar")
  
  plot(0,0,ylim=posterrange,xlim=mcpars[1:2],xlab="Iterations",
       ylab=interest)

  iters <- seq(mcpars[1],mcpars[2],mcpars[3])
  chainlength <- sapply(mcmcobject,function(X) length(X[,postcoln]))
  posterlong <- lapply(mcmcobject, function(X) as.numeric(X[,postcoln]))

  for ( i in 1:nchains){
    lines(iters, posterlong[[i]],lwd=1,col=colorz[i])
  }
   text(mcpars[1]+(0.5*(mcpars[2]-mcpars[1])),posterrange[1]+
       (posterrange[2]-posterrange[1])*0.5
       ,interest,cex=4,
        col='grey20')
  postmean <- mean(sapply(posterlong,mean))
  abline(h=postmean,lty=2,lwd=2)
  text(mcpars[1]+(0.1*(mcpars[2]-mcpars[1])),postmean,
       "Mean",cex=2, col='white')

  boxplot(posterlong,col=colorz)
  abline(h=postmean,lty=2,lwd=2)

  subsetpost <- as.mcmc.list(lapply(mcmcobject,function(X) X[,1]))
  gelman.plot(subsetpost,auto.layout=FALSE)

  acfcoords <- par("usr")
  text(acfcoords[1]+(0.5*(acfcoords[2]-acfcoords[1])),
       acfcoords[3]+
       (acfcoords[4]-acfcoords[3])*0.5
       ,"Gelman diagnostic",cex=2,
       col='grey20')

  Autocor <- acf(posterlong[[1]])
  acfcoords <- par("usr")
  text(acfcoords[1]+(0.5*(acfcoords[2]-acfcoords[1])),
       acfcoords[3]+
       (acfcoords[4]-acfcoords[3])*0.5
       ,"ACF",cex=4,
       col='grey20')

  plot(0, 0, type="n", xlab="", ylab="",xlim=c(0,1),ylim=c(0,1),
       xaxt='n',yaxt='n')
  text(0.4,0.75,paste("Summary for:", interest))
  text(0.4,0.65,paste("Max iter:", mcpars[2]))
  text(0.4,0.55,paste("Thinning:", mcpars[3]))
  text(0.4,0.45,paste("Mean prior est:", round(mean(priorsamples),3)))
  text(0.4,0.35,paste("Mean post. est:", round(postmean,3)))
  text(0.4,0.25,paste("Prop. difference:",round(
                      100*((mean(priorsamples)-postmean)/postmean),3)))
  
  }
}

##' findRprior seives through your model files to find and
##' translate your priors to proper R probabilty density functions.
##'
##' Given a model file and a the name of the prior
##' findRprior will return the corresponding 
##' functions in R and translate the parameters from their
##' BUGS formulation into a R compatable form.
##' @title Finds priors and translates these from BUGS to R
##' @param interest a set of parameters to produce
##' @param model the path to the BUGs/JAGS model text file with
##' the priors specified. This model must be the one
##' that generated the MCMC list. 
##' @author Marco D. Visser
##' @concept MCMC diagnostics
##'
##' @export
findRprior <- function(model=NULL,interest=NULL,silent=FALSE){

  if(is.null(interest)) {
    stop("Parameter of interest missing. Please check the inputs")}
  
  if(is.null(model)|!file.exists(model)) {
    stop("No model file supplied or the file does not exist")
   }

  ModelFile <- readLines(model)
  ## find all relavent lines mentioning the prior
  filter1 <- grep(interest,ModelFile)

  if(length(filter1)==0){stop("Prior not found in the model file")}

  ## remove white space
  ModelFileNoWS <- sapply(1:length(ModelFile),
                         function(X) gsub(" ", "",ModelFile[X]))

  ## find the prior
  filter2 <- grep(paste(interest,"~",sep=""),ModelFileNoWS)

  if(length(filter2)!=1){
    if(length(filter2)==0){
      stop(paste("Prior",interest,
                 "not found in model file (is it a true prior?), check inputs?"
                 ))
    } else {
    cat("Multiple lines were found to correspond to the prior\n")
    cat("Identified lines were: \n ")
    print(ModelFile[filter2])
    cat("Are these the true priors or aliases? \n ") 
    stop("No unique prior found in model file")}}

  PriorLine <- filter2
    
  ## find the prior and remove whitespace
  BUGSprior <- gsub(" ", "",ModelFile[PriorLine])
  if(!silent){
    print(paste("BUGS prior identified as: ", BUGSprior) )
  }
  
  PriorSplit <-lapply(list("~","\\(",","),function(x) strsplit(BUGSprior,x))
  FulldensityBUGS <- PriorSplit[[1]][[1]][2]

  ## extract coefficients - thanks to Lazar!
  ## one regex to rule them all : '(-|)[0-9]+\\.*[0-9]*|*\\d+([e|E][+\\-]*\\d*)'
  numbermatch <- gregexpr('(-|)[0-9]+\\.*[0-9]*|*\\d+([e|E][+\\-]*\\d*)'
                          ,FulldensityBUGS)
  coefsBUGS <- regmatches(FulldensityBUGS,numbermatch)
  dfBUGS <- strsplit(FulldensityBUGS,"\\(")[[1]][1]
  Rtranslation <- BUGS2RTranslator(dfBUGS,coefsBUGS)
  
  ## check for corresponding R function
  if(exists(Rtranslation[[1]])){
     if(!silent){
     print(paste("BUGS prior translated to:", Rtranslation[[1]]) ) }
  } else {
    stop(paste("No equivalent R function found for the BUGS function ",dfBUGS))
  }

  return(Rtranslation)
}

##' BUGS2Rtranslator translates probabilty density functions
##' from BUGS to R.
##'
##' Given a list of functions and optionally a list of parameters,
##' BUGS2Rtranslator will return the corresponding 
##' functions in R and translate the parameters from their
##' BUGS formulation into a R compatible form.
##' @title translates from BUGS to R
##' @param bugsnames a list/vector of model names (as characters)
##' @param bugsparameters a list of model parameters
##' which correspond to the models in bugsnames
##' and are in the correct and standard order (e.g. c(mu, sigma)).
##' @author Marco D. Visser
##' @examples
##' 	bugs2rTranslator("norm",list(c(10,100)))
##' @export
BUGS2RTranslator <- function(bugsname="norm",bugsparameters=NULL) {

  if(!is.null(bugsparameters)){
  if(!is.element(class(bugsparameters),"list")){
    stop("parameters must be supplied as a list")
  } }
  
  ## Remove potential confusion due to "d" in e.g. dnorm
  bugsname_d<- lapply(bugsname,function(X) gsub("d","",X))
  bugsname_d<- as.character(unlist(bugsname_d))

  ## Search patterns with word boundries for exact matches for nested
  ## words
  NestFunc <- c("\\bnbinom\\b","\\bnbin\\b","\\bbinom\\b","\\bbin\\b",
                "\\bnorm\\b","\\blnorm\\b")
  ## Likely no word boundries needed here
  UniqFunc <- c("weib","gamma","chisq","unif","pois","exp")

  Rnames <- c("dnbinom","dnbinom","dbinom","dbinom",
              "dnorm","dlnorm","dweibull","dgamma","dchisq",
              "dunif","dpois","dexp")
    

  DistHits <- data.frame(sapply(c(NestFunc,UniqFunc),
                                function(x) grepl(x,bugsname_d)))

  if(ncol(DistHits)<=1) {DistHits <- t(DistHits)}
    
  if (max(rowSums(DistHits)) >  1) {
    badrow <- rowSums(DistHits) > 1
    stop(paste(unique(bugsname[badrow])),
         " correspond(s) to > 1 distribution")
      }


    if (min(rowSums(DistHits)) == 0) {
    badrow <- rowSums(DistHits) == 0
    stop(paste(unique(bugsname[badrow])),
         " no R eqv. found, remove from list and run again")
    }
  
  ## Match to R names
  Rhits<-apply(DistHits, 1, function(X) Rnames[X])
  names(Rhits) <- NULL

  ## Now translate parameters, if needed
  if(!is.null(bugsparameters)) {
  bugsparameters <- lapply(bugsparameters,as.numeric)
  Rparameters <- bugsparameters
  norm <- grep("norm",Rhits)
  weib <- grep("weib",Rhits)
  bin <- grep("bin",Rhits)
  
  if(length(norm)>0){
    for(i in norm){

    Rparameters[[i]] <- c(bugsparameters[[i]][1],
                             bugsparameters[[i]][2]^-0.5)
   }
  }

  if(length(weib)>0){
    for(i in weib){

     Rparameters[[i]] <- c(bugsparameters[[i]][1],
         bugsparameters[[i]][2]^(-1/bugsparameters[[i]][1]))

   #   Rparameters[[i]]<- c(Rparameters[[i]][2],
   #      Rparameters[[i]][1])
                             
   }
  }


  return(list(Rhits,Rparameters))
 } else {return(Rhits)}  

}


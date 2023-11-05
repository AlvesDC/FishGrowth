#' Fit a Length and Likelihood-based Growth Model (LLGM) 
#'
#' The fit function, fitAgrowth, used to fit the models to the length-frequency data.
#'  
#' @param data A vector of data containing the length of each individual observed.
#' @param model The model name to be estimated, between the "*HmA*", "*HmAG*", "*HtA*" and "*HtAG*" options.
#' @param gpars Vector containing the four growth parameters \eqn{L_{\infty}, K, t_{0}, \sigma_{VBGE}^2}. 
#' The values specified for the parameters \eqn{L_{\infty}} and \eqn{K} will be used by the *HmA* and *HtA* models, 
#' while the *HmAG* and *HtAG* models will not use this information, since they will estimate values for 
#' these two parameters. The parameters \eqn{t_{0}} and \eqn{\sigma_{VBGE}^2} should be specified as 
#' assumptions for all models.
#' @param lbin Length classes amplitude. 
#' @param lmin Lower limit of the first length class. 
#' @param lmax Upper limit of the last length class.  
#' @param amax Maximum age to be considered in the model.  
#' 
#' @return Return a fit object.
#' 
#' @author Diego C Alves, \email{dcalves@@uem.br}
#' @references \url{https://doi.org/10.1016/j.fishres.2018.02.013}
#' 
#' @examples
#' Pkl(pars=c(60, 0.3), gpars = c(0, 3))
#' 
#' @export
fitAgrowth <- function(data, model="HmA", gpars=c(60,.2,0,3), lbin=1, lmin=0, lmax=100,
                       amax=10){
  Lobs  <<- hist(data,breaks=c(-999,seq(lmin,lmax,lbin)[-1]),plot=F)$counts
  COV   <<- matrix(0,nrow=lmax,ncol=lmax)
  gpars <<- gpars
  lbin  <<- lbin
  lmin  <<- lmin
  lmax  <<- lmax
  amax  <<- amax
  
  if(model=="HmA"){
    P    <<- Pkl(gpars[c(1,2)],gpars,lbin=lbin,lmax=lmax,lmin=lmin,amax=amax)
    tmp  <- DEoptim(fn = nllk.HmA, lower = c(0.000001,rep(0,amax+1)), 
                    upper = c(1e6, rep(length(data), amax+1)), 
                    control = DEoptim.control(trace=F,itermax=9999,c=.4,
                                              reltol=.000000000000001,CR=0,
                                              F=.65, steptol=150,NP=400,
                                              parallelType=1,
                                              packages=list("mnormt"),
                                              parVar=list("P","Lobs","COV"))
    )$optim
    tmp2 <- optimHess(par=tmp$bestmem,fn=nllk.HmA,control=list(maxit=0))
    PD   <- is.positive.definite(tmp2)}
  
  if(model=="HmAG"){
    tmp  <- DEoptim(fn=nllk.HmAG, lower=c(0.000001,1,0.00001,rep(0,amax+1)),
                    upper=c(1e6,lmax,2,rep(length(data),amax+1)),
                    control=DEoptim.control(trace=F,itermax=9999,c=.4,
                                            reltol=.000000000000001,CR=0,
                                            F=.65,steptol=150,NP=400,
                                            parallelType=1,
                                            packages=list("mnormt"),
                                            parVar=list("Lobs","gpars","COV","Pkl","lmin","lbin","lmax","amax")))$optim
    tmp2 <- optimHess(par=tmp$bestmem,fn=nllk.HmAG,control=list(maxit=0))
    PD   <- is.positive.definite(tmp2)}
  
  if(model=="HtA"){
    P    <<- Pkl(gpars[c(1,2)],gpars,lbin=lbin,lmax=lmax,lmin=lmin,amax=amax)
    tmp  <- DEoptim(fn=nllk.HtA, lower=c(0.000001,rep(0,amax+1)), upper=c(1e6,rep(length(data),amax+1)),
                    control = DEoptim.control(trace=F,itermax=9999,c=.4,reltol=.000000000000001,CR=0,F=.65,steptol=150,NP=400,parallelType=1,
                                              packages=list("mnormt"),parVar=list("Lobs","COV","P")))$optim
    tmp2 <- optimHess(par=tmp$bestmem,fn=nllk.HmA,control=list(maxit=0))
    PD   <- is.positive.definite(tmp2)}
  
  if(model=="HtAG"){
    tmp  <- DEoptim(fn=nllk.HtAG, lower=c(0.000001,1,0.00001,rep(0,amax+1)), upper=c(1e6,lmax,2,rep(length(data),amax+1)),
                    control = DEoptim.control(trace=F,itermax=9999,c=.4,reltol=.000000000000001,CR=0,F=.65,steptol=150,NP=400,parallelType=1,
                                              packages=list("mnormt"),parVar=list("Lobs","gpars","COV","Pkl","lmin","lbin","lmax","amax")))$optim
    tmp2 <- optimHess(par=tmp$bestmem,fn=nllk.HtAG,control=list(maxit=0))
    PD   <- is.positive.definite(tmp2)}
  
  if(PD==TRUE){
    SE    <- sqrt(diag(chol2inv(chol(tmp2))))
    
    return(list(nllk=tmp$bestval,convergence=paste("converged = ",PD),Estimativas=cbind(Est=tmp$bestmem,SE=SE)))
  }else{return(list(nllk=tmp$bestval,convergence=paste("converged = ",PD),Estimativas=cbind(Est=tmp$bestmem)))}}
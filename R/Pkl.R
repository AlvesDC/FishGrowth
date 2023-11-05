#' Transition (P) matrix
#'
#' This function, Pkl, computes the transition matrix, based on the growth parameters of the von Bertalanffy equation.
#'  This function will be automatically used within the estimation function, \code{\link{fitAgrowth}}.
#' 
#' @param pars A vector of growth parameters (Linf and K) that can be estimates
#' @param gpars A vector of growth parameters (t0 and sd) that can not be estimates
#' @param lmin A value of minimum length of the P matrix
#' @param lmax A value of the maximum length of the P matrix
#' @param lbin A value of the length bins of the P matrix
#' @param amax A value of the maximum age of the P matrix
#' 
#' @return Return a P matrix [rows = lengths, cols = ages].
#' 
#' @author Diego C Alves, \email{dcalves@@uem.br}
#' @seealso \code{\link{fitAgrowth}}
#' @references \url{https://doi.org/10.1016/j.fishres.2018.02.013}
#' 
#' @examples
#' Pkl(pars=c(60, 0.3), gpars = c(0, 3))
#' 
#' @export
Pkl <-
  function(pars,
           gpars,
           lmin = 0,
           lmax = 100,
           lbin = 1,
           amax = 10) {
    
    Lclass <- seq(lmin, lmax - lbin, by = lbin)
    Aclass <- seq(0, amax) + .5
    Lmean <- pars[1] * (1 - exp(-pars[2] * (Aclass - gpars[3])))
    SD <- Lmean * gpars[4]
    
    P <-
      matrix(NA,
             length(Lclass),
             length(Aclass),
             dimnames = list(Lclass, Aclass))
    
    for (i in 1:length(Aclass)) {
      d1 <- pnorm(q = Lclass + lbin,
                  mean = Lmean[i],
                  sd = SD[i])
      d2 <- pnorm(q = Lclass, mean = Lmean[i], sd = SD[i])
      P[, i] <- d1 - d2
    }
    
    return(P)
  }

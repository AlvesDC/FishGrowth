#' Negative logarithm likelihood of the Homoscedastic Age Model (HmA) 
#'
#' This function computes the Negative logarithm likelihood of the Homoscedastic Age Model (HmA). 
#' This functions is used in the \code{\link{fitAgrowth}} estimation function. 
#' This objective functions present a single argument, which is the vector of parameters to be estimated (pars), to enable parallel computing.
#' 
#' @param pars Vector of the parameters to be estimated (sigma2, A), where A is a age composition vector.
#' 
#' @return Return a negative logarithm likelihood of the Homoscedastic Age Model (HmA). 
#' 
#' @author Diego C Alves, \email{dcalves@@uem.br}
#' @seealso \code{\link{fitAgrowth}}
#' @references \url{https://doi.org/10.1016/j.fishres.2018.02.013}
#' 
#' @export
nllk.HmA <- function(pars) {
  Lest      <- P %*% pars[-1]
  diag(COV) <- pars[1]
  Nllk      <- dmnorm(as.vector(Lobs), as.vector(Lest), COV, log = T)
  return(-sum(Nllk))
}
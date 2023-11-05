#' Virtual population simulation
#'
#' It simulates a simple virtual population with specified characteristics
#'
#' @param Amax A value of the age max in the population
#' @param Z A value of total mortality (Z)
#' @param k A value of growth coeficient of the VBGF
#' @param Linf A value of asymptotic length
#' @param t0 A value of t0 parameter of the VBGF
#' @param SD A value of the standard deviation of the length at age
#' @param lmin A value of minimum length of the P matrix
#' @param lmax A value of the maximum length of the P matrix
#' @param lbin A value of the length bins of the P matrix
#' @param Nyears Number of years that will be created
#' @param Rmin A value of the minimum Recruitment
#' @param Rmax A value of the maximum Recruitment
#' @param Namos Number of years that will be sampling
#'
#' @return Return a P matrix [rows = lengths, cols = ages].
#'
#' @examples
#' Pkl(pars=c(60, 0.3), gpars = c(0, 3))
#'
#' @export
PopSim <- function(Amax = 10,
                   Z = .45,
                   k = .103,
                   Linf = 63,
                   t0 = -.25,
                   SD = 3,
                   lmin = 0,
                   lmax = 100,
                   lbin = 0.5,
                   Nyears = 20,
                   Rmin = 1000,
                   Rmax = 1000000,
                   Namos = 4) {
  
  Nyears <- Nyears + Amax
  R <- runif(Nyears, Rmin, Rmax)
  Lclass <- seq(lmin, lmax, by = lbin) + lbin / 2
  Aclass <- seq(0, Amax) + .5
  Lmean <- Linf * (1 - exp(-k * (Aclass - t0)))

  # Matriz de N[i,j], individuos da idade i no ano j
  Nij <- matrix(NA, (Amax + 1), (Nyears + Amax))
  tmp <- matrix(NA, (Amax + 1), 1)
  for (j in 1:Nyears) {
    for (i in 1:(Amax + 1)) {
      if (i == 1) {
        tmp[1] <- R[j]
      }
      else{
        tmp[i] <- exp(log(tmp[i - 1]) - Z)
      }
    }
    diag(Nij[, j:(j + Amax)]) <- tmp
  }
  
  Nij <- Nij[, (Amax + 1):Nyears]
  Nij <- ceiling(Nij)
  # Estrutura etaria/LS Real para os anos amostrados
  # Ntrue
  Aj <- Nij[,-(1:(ncol(Nij) - Namos))]
  Ntrue <- matrix(NA, sum(Aj), 3)
  colnames(Ntrue) <- c("year", "LS", "age")
  l <- 1
  m <- 0
  for (j in 1:ncol(Aj)) {
    for (i in 1:(Amax + 1)) {
      tmp <- c()
      m <- Aj[i, j] + m
      tmp <- truncnorm::rtruncnorm(Aj[i, j], a=0, b=Inf, mean = Lmean[i], sd = SD)
      Ntrue[l:m,] <-
        cbind(rep(j, Aj[i, j]), tmp, rep((i - 0.5), Aj[i, j]))
      l <- 1 + m
    }
  }
  
  return(list(
    Ntrue = Ntrue,
    Recs = R,
    Nij = Nij
  ))
}

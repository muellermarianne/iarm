#' Gamma Coefficient of two Ordinal Variables
#'
#' Calculates the Gamma coefficient as a measure of association between two ordinal variables.
#' @param x A matrix of counts.
#' @return Gamma coefficient, standard error and p value.
#' @author Marianne Mueller; original version Greg Rodd
#' @references  Agresti, A. \emph{Categorical Data Analysis}. John Wiley & Sons, 2013, pp. 57-59.
#' @import stats
#' @export
#' @examples # Association between raw score and sex of patients
#' score <- apply(amts[, 4:13], 1, sum, na.rm = TRUE)
#' gamma_coef(table(score,amts$sex))
#'
#' # Association between raw score and indication group of patients
#' score <- apply(desc2[, 5:14], 1, sum, na.rm = TRUE)
#' gamma_coef(table(score,desc2$group))
gamma_coef<-function(x){
  n <- nrow(x)
  m <- ncol(x)
  pic <- pid <- matrix(0, nrow=n, ncol=m)
  rowx <- row(x)
  colx <- col(x)
  for(i in 1:n){
    for(j in 1:m){
      pic[i, j] <- sum(x[rowx < i & colx < j]) + sum(x[rowx > i & colx > j])
      pid[i, j] <- sum(x[rowx < i & colx > j]) + sum(x[rowx > i & colx < j])
    }
  }
  C <- sum(pic*x)/2
  D <- sum(pid*x)/2
  psi <- 2*(D*pic - C*pid)/(C + D)^2
  sigma2 <- sum(x*psi^2) - sum(x*psi)^2
  gamma <- (C - D)/(C + D)
  pwert <- 2*(1-pnorm(abs(gamma)/sqrt(sigma2)))
  c(gamma = gamma, se = sqrt(sigma2), pvalue = pwert)
}

#' Conditional and Partial Gamma Coefficients
#'
#' Calculates conditional and partial Gamma coefficients for x and y given z with confidence intervals.
#' @param x,y,z  Three numeric vectors or factors.
#' @param conf.level Confidence level for the returned confidence interval.
#' @return matrix with estimates, standard errors and confidence interval limits.
#' @export
#' @author Marianne Mueller
#' @examples # Partial Gamma coefficient between an item and an exogenuous variable, given the total score
#' score <- apply(amts[, 4:13], 1, sum, na.rm = TRUE)
#' fz <- cut(score,unique(quantile(score,0:10/10)))
#' partgam(amts$firstww,amts$sex,fz)
#' @references Davis, J. A. A Partial coefficient for Goodman and Kruskal's Gamma.
#'  \emph{Journal of the American Statistical Association}, 62 (317), 1967, pp. 189-193.
partgam <- function(x, y, z, conf.level = 0.95){
  xxx <- xtabs(~ x + y + fz)
  n <- dim(xxx)[1]
  m <- dim(xxx)[2]
  k <- dim(xxx)[3]
  sigma2.k <- rep(NA, k)
  gamma.k  <- rep(NA, k)
  C.k <- rep(NA, k)
  D.k <- rep(NA, k)
  pi.c.k <- pi.d.k <- array(0, c(n, m, k))
  for (l in 1:k){
    xx <- xxx[, , l]
    row.x <- row(xx)
    col.x <- col(xx)
    for(i in 1:n){
       for(j in 1:m){
         pi.c.k[i, j, l]<-sum(xx[row.x < i & col.x < j]) + sum(xx[row.x > i & col.x > j])
         pi.d.k[i, j, l]<-sum(xx[row.x < i & col.x > j]) + sum(xx[row.x > i & col.x < j])
       }
    }
    C.k[l] <- sum(pi.c.k[, , l]*xx)/2
    D.k[l] <- sum(pi.d.k[, , l]*xx)/2
    psi.k <- 2*(D.k[l]*pi.c.k[, , l]-C.k[l]*pi.d.k[, , l])/(C.k[l]+D.k[l])^2
    sigma2.k[l] <- sum(xx*psi.k^2)
    gamma.k[l] <- (C.k[l] - D.k[l])/(C.k[l] + D.k[l])
  }
  C <- sum(C.k)
  D <- sum(D.k)
  psi <- 2*(D*pi.c.k[, , ]-C*pi.d.k[, , ])/(C + D)^2
  sigma2 <- c(sigma2.k, sum(xxx*psi^2))
  gamma <- c(gamma.k, (C - D)/(C + D))
  pr2 <- 1 - (1 - conf.level)/2
  CIa <- outer(sqrt(sigma2),qnorm(pr2)*c(-1, 1)) + gamma
  mm <- data.frame(gamma, se = sqrt(sigma2), CI1 = ifelse(CIa[, 1] > -1 ,CIa[, 1], -1),
  CI2 = ifelse(CIa[, 2] < 1, CIa[, 2], 1))
  row.names(mm)= c(paste("conditional", 1:k), "partial")
  mm
}





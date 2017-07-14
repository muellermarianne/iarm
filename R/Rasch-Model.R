#' Construction of design matrix for partial credit model
#'
#' Construction of design matrix for partial credit model  with
#' the usual parameter restriction.
#'
#' @param x vector of length equal to the number of items. The elements of x are the maximal response categories for the items.
#' @return design matrix W for PCM.
#' @export
#' @examples
#' # desc2 contains 10 items with response categories from 0 to 4 each.
#' W_thresh(rep(4,10))
W_thresh <- function(x){
  mi <- x
  W1 <- diag(sum(mi)-1)
  W1 <- rbind(W1,rep(0,sum(mi)-1))
  W1[sum(mi),cumsum(mi[-length(mi)])] <- -1
  W1
}

#' Computation of item targets for polytomous models
#'
#' The item target is the value of the person parameter where
#' item information is maximized.
#' @param  obj object of class "eRm" (but not "dRm"), a fitted partial
#' credit model using  the function PCM in package eRm.
#' @return vector with item targets.
#' @author Marianne Mueller
#' @import eRm
#' @export
#' @examples
#'  pc.mod <- PCM(desc2[, 5:14])
#'  item_target(pc.mod)
item_target <- function(obj){
  if(obj$model == "RM") stop("Item targets are computed only for polytomous models!")
  thresh1 <- thresholds(obj)[[3]][[1]][, -1] - mean(thresholds(obj)[[3]][[1]][, 1])
  betasum.l <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
  var.X <- function(x) {
    function(theta) {
      xvec <- 0:length(x)
      pp <- (exp(theta*xvec)*exp(-c(0, x)))/(exp(theta*xvec)%*%exp(-c(0, x)))
      (xvec - pp%*%xvec)^2%*%pp
    }
  }
  varf <- lapply(betasum.l, var.X)
  target <- unlist(sapply(varf, optimize, interval = c(-4, 4), maximum = T)[1, ])
  cat("\n")
  cat("Item targets:","\n\n")
  print(target, 3)
  cat("\n\n")
  invisible(target)
}

#' Conditional likelihood ratio tests (CLR)
#'
#' The conditional likelihood ratio tests compare item parameters in low and high score groups
#' for an overall test of homogeneity, and in groups defined by the levels of exogenous factors
#' for tests of no differential item functioning (DIF).
#' @author Marianne Mueller
#' @references Andersen, E.B. (1973). A goodness of fit test for the Rasch model. \emph{Psychometrika}, 38, 123-140.
#' @param dat.items data frame with the responses to the items.
#' @param dat.exo  data frame consisting of exogenuous factor variables.
#' @param model If model="RM" a Rasch model will be fitted,
#' if model="PCM" a partial credit model for polytomous items is used.
#' @return matrix with test statistics, df and p values.
#' @import stats
#' @importFrom psychotools pcmodel
#' @export
#' @examples #CLR overall test and test of  no DIF for agegrp and sex
#' clr_tests(amts[,4:13],amts[,2:3])
clr_tests <- function(dat.items, dat.exo, model = c("RM","PCM")) {
  ok <- complete.cases(dat.items)
  dat.i <- dat.items[ok, ]
  sgrp <- score_groups(dat.i)
  dat.exo <- data.frame(dat.exo)
  ok1 <- complete.cases(cbind(dat.items, dat.exo))
  model <- match.arg(model)
  if (model == "RM") {
    clr0 <- -2*(RM(dat.i)$loglik - (RM(dat.i[sgrp==1, ])$loglik + RM(dat.i[sgrp==2, ])$loglik))
    clrhomo <- function(exo, data){
      sum(sapply(split(data, exo, drop = TRUE), function(x){RM(x)$loglik}))
    }
    ll <- RM(dat.items[ok1, ])$loglik
    m <- dim(dat.i)[2]
  } else {
    clr0 <- -2*(pcmodel(dat.i,hessian=F)$loglik - (pcmodel(dat.i[sgrp==1, ],hessian=F)$loglik + pcmodel(dat.i[sgrp==2, ],hessian=F)$loglik))
    clrhomo <- function(exo, data){
      sum(sapply(split(data, exo, drop = TRUE),function(x){pcmodel(x,hessian=F)$loglik}))
    }
    ll <- pcmodel(dat.items[ok1, ],hessian=F)$loglik
    m <- sum(apply(dat.items[ok1, ], 2, max, na.rm = TRUE))
  }
  clr <- c(clr0, -2*(ll - apply(dat.exo[ok1, ,drop=F],2,clrhomo, dat.items[ok1, ])))
  df <- (m-1)* (c(1, sapply(dat.exo[ok1, ,drop=F],function(x){sum(table(x) != 0)}) - 1))
  mm <- cbind(clr, df)
  row.names(mm) <- c("overall", names(dat.exo[ok1, ,drop=F]))
  pvalue <- apply(mm,1,function(x){1-pchisq(x[1],x[2])})
  symp <- symnum(pvalue, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),symbols = c(" ***", " **", " *", " .", " "))
  mm <- noquote(cbind(clr = format(mm[, 1], digits = 3), df = format(mm[, 2], digits = 1), pvalue = format(pvalue, digits = 2), sig = symp))
  cat("\n")
  cat("Conditional Likelihood Ratio Tests:")
  cat("\n")
  mm
}


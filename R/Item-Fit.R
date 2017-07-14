expvar <- function(data, coeff){
  rv <- rowSums(data, na.rm = TRUE)
  rv[rv==0] <- NA
  k <- dim(data)[2]
  n <- dim(data)[1]
  mi <- sapply(coeff, length)
  m <- sum(mi)
  Er <- matrix(rep(0, m*k), ncol=k)
  Varx <- matrix(rep(0, m*k), ncol=k)
  gr <- psychotools::elementary_symmetric_functions(coeff)[[1]]
  grminus <- vector("list",k)
  for (i in 1:k) {
    grminus[[i]] <- psychotools::elementary_symmetric_functions(coeff[-i])[[1]]
  }
  for (r in 1:(m-1)) {
    for (i in 1:k) {
      pscore <- rep(0,mi[i])
      for (x in 1:mi[i]) {
        if (x <= r & m-mi[i]+x >= r) {
          pscore[x] <- exp(-coeff[[i]][x])*(grminus[[i]])[r+1-x]/gr[r+1]
          }
      }
      pscore <- c(1-sum(pscore), pscore)
      Er[r,i] <- (0:mi[i])%*%pscore
      Varx[r,i] <- ((0:mi[i])-Er[r,i])^2%*%pscore
    }
  }
  Er[m,] <- 1
  Ehat <- Er[rv,]
  VarX <- Varx[rv,]
  colnames(Ehat) <- paste("Item", 1:k, sep="")
  colnames(VarX) <- paste("Item", 1:k, sep="")
  list(Ehat=Ehat, VarX=VarX)
}

#' Observed and expected item mean scores
#'
#' Homogeneity of item responses in the low and high score groups is analyzed by looking at observed and expected item mean scores
#' together with standardized residuals. If the Andersen's CLR test has shown some evidence against homogeneity,
#' this comparison can indicate which items might be responsable.
#' @param  object object of class "Rm", a fitted Rasch model oder partial
#' credit model using  the functions RM or PCM in package eRm.
#' @return list with observed and expected mean scores together with standardized residuals for the two score groups.
#' @author Marianne Mueller
#' @importFrom psychotools elementary_symmetric_functions
#' @export
#' @examples
#' rm.mod <- RM(amts[,4:13])
#' item_obsexp(rm.mod)
#'
#' pc.mod <- PCM(desc2[,5:14])
#' item_obsexp(pc.mod)
item_obsexp <- function(object){
  if (!("Rm"%in%class(object))) stop("object must be of class Rm!")
  k <- dim(object$X)[2]
  if (object$model == "RM") {
    coeff <- (-1)*coef(object)
    mi <- rep(1, k)
  } else {
    mi <- apply(object$X, 2, max, na.rm = TRUE)
    thresh1 <- thresholds(object)[[3]][[1]][, -1] - mean(thresholds(object)[[3]][[1]][, 1])
    coeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
  }
  m <- sum(mi)
  score <- apply(object$X, 1, sum, na.rm=T)
  sgrp <- score_groups(object$X)
  resgrp <- function(x){
    EV <- expvar(object$X[sgrp==x & score > 0 & score < m, ], coeff)
    resm <- cbind(apply(object$X[sgrp==x & score > 0 & score < m, ], 2, mean, na.rm=T), colMeans(EV$Ehat, na.rm=T),
                  length(EV$VarX[!is.na(EV$VarX[, 1]), 1]) * (apply(object$X[sgrp==x & score > 0 & score < m, ], 2, mean, na.rm=T)
                                                              - colMeans(EV$Ehat, na.rm=T))/sqrt(colSums(EV$VarX, na.rm = T)))
    colnames(resm)= c("mean obs", "mean exp", "std.res")
    symp <- symnum(resm[, 3], cutpoints=c(-100, -3, -2, 2, 3, 100), symbols=c("--", "-", " ", "+", "++"))
    resm <- noquote(cbind(format(resm, digits = 3), sig = symp))
    resm
  }
  result <- lapply(c(1,2),resgrp)
  cat("Score group 1:","\n")
  print(result[[1]])
  cat("\n")
  cat("Score group 2:","\n")
  print(result[[2]])
  cat("\n")
  invisible(return)
}


#' Item Outfit and Infit Statistics
#'
#' To avoid bias observed item responses are compared to expected responses under
#' the conditional distribution of responses given the total score. This leads to standardized  residuals
#' which can be summarized to outfit and infit statistics in the usual way.
#' @param  object an object of class "Rm", a fitted Rasch model oder partial
#' credit model using  the functions RM or PCM in package eRm.
#' @param se if TRUE the standard errors will be included.
#' @details The fit statistics and their standard errors are calculated as described in Christensen et al.
#' P values are are based on the normal distribution of the standardized fit statistics.
#' @author Marianne Mueller
#' @export
#' @references Christensen, K. B. , Kreiner, S. & Mesbah, M. (Eds.)
#' \emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 86 - 90.
#'
#'  Kreiner, S. & Christensen, K. B. (2011) Exact evaluation of Bias in Rasch model residuals.
#' \emph{Advances in Mathematics Research}, 12, 19-40.
#' @return an object of class outfit containing:
#' \item{outfit}{outfit statistics}
#' \item{outfit.se}{standard errors of outfit statistics}
#' \item{out.pvalue}{p values of outfit statistics}
#' \item{infit}{infit statistics}
#' \item{infit.se}{standard errors of infit statistics}
#' \item{in.pvalue}{p values of infit statistics}
#' @examples
#' rm.mod <- RM(amts[,4:13])
#' out_infit(rm.mod)
out_infit <- function(object,se=TRUE){
  if (!("Rm"%in%class(object))) stop("object must be of class Rm!")
  rv <- rowSums(object$X, na.rm = TRUE)
  k <- dim(object$X)[2]
  if (object$model=="RM") {
    koeff <- (-1)*coef(object)
    mi <- rep(1,k)
  } else {
    mi <- apply(object$X, 2, max, na.rm = TRUE)
    thresh1 <- thresholds(object)[[3]][[1]][, -1] - mean(thresholds(object)[[3]][[1]][, 1])
    koeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
  }
  m <- sum(mi)
  rvstrich <- rv[rv > 0 & rv < m]
  ER <- matrix(rep(0, (m-1)*k), ncol = k)
  VarZ <- matrix(rep(0, (m-1)*k), ncol = k)
  VarX.R <- matrix(rep(0, (m-1)*k), ncol = k)
  gr <- elementary_symmetric_functions(koeff)[[1]]
  grminus <- vector("list",k)
  for (i in 1:k) {
    grminus[[i]] <- elementary_symmetric_functions(koeff[-i])[[1]]
  }
  for (r in 1:(m-1)) {
    for (i in 1:k) {
      pscore <- rep(0, mi[i])
      for (x in 1:mi[i]) {
        if (x <= r  & x >= r-m+mi[i]) {
          pscore[x] <- exp(-koeff[[i]][x]) * (grminus[[i]])[r+1-x] / gr[r+1]
        }
      }
      pscore <- c(1-sum(pscore), pscore)
      ER[r,i] <- as.numeric((0:mi[i]) %*% pscore)
      VarX.R[r,i] <- as.numeric(((0:mi[i])-ER[r,i]) ^ 2 %*% pscore)
      VarZ[r,i] <- as.numeric((((0:mi[i]) - ER[r,i])^2/VarX.R[r, i] - 1)^2 %*% pscore)
    }
  }
  Ehat <- ER[rvstrich,]
  VarX <- VarX.R[rvstrich,]
  Z.vi2 <- ((object$X[rv > 0 & rv < m,] - Ehat) ^ 2) / VarX
  Outfit = colMeans(Z.vi2, na.rm=T)
  X <- as.data.frame(object$X)
  nir <- t(sapply(split(X[rv > 0 & rv < m, ], factor(rvstrich, levels = 1:(m-1))),function(x){apply(x, 2, function(y){length(na.exclude(y))})}))
  if (se==TRUE) {
    Outfit.se <- sqrt(colSums((nir/colSums(nir)^2) * VarZ))
    fit <- cbind(Outfit,Outfit.se)
    pwert <- function(x){ifelse(x[, 1] > 1, 2*(1 - pnorm((x[, 1] - 1)/x[, 2])), 2*(pnorm((x[, 1] - 1)/x[, 2])))}
    out.pwert <- pwert(fit)
  }
  Wri <- VarX.R / matrix(rep(colSums(nir * VarX.R), m-1), ncol = k, byrow = T)
  Wvi <- Wri[rvstrich, ]
  Infit <-  colSums(Wvi * Z.vi2, na.rm = T)
  if (se==TRUE) {
    Infit.se <- sqrt(colSums(nir * Wri^2 * VarZ,na.rm = T))
    names(Infit)=colnames(X)
    fit <- cbind(Infit,Infit.se)
    in.pwert <- pwert(fit)
    result <- list(outfit=Outfit, outfit.se=Outfit.se, out.pvalue=out.pwert, infit=Infit, infit.se=Infit.se, in.pvalue=in.pwert)
  } else {
    result <- list(outfit=Outfit, infit=Infit)
  }
  class(result) <- "outfit"
  result
}

#' Print method for output of out_infit
#' @param x object of class outfit.
#' @param ... arguments passed to other functions.
#' @export
print.outfit <- function(x, ...){
  if(length(x)==2) {
    print(data.frame(Outfit=x[[1]],Infit=x[[2]]),digits=3)
    cat("\n")
  } else {
    symp <- function(x){symnum(x,cutpoints=c(0,0.001,0.01,0.05,0.1,1),symbols=c("***   "," **  "," *   "," .   ","    "))}
    tab2 <- noquote(cbind(Outfit=round(x[[1]],digits=3),se=round(x[[2]],digits=3), pvalue=round(x[[3]],digits=3),sig=symp(x[[3]]),
                                 Infit=round(x[[4]],digits=3),se=round(x[[5]],digits=3),pvalue=round(x[[6]],digits=3) ,sig=symp(x[[6]])))
    cat("\n")
    print(tab2)
    cat("\n")
  }
}




pscore_poly  <- function(i,x,r,coeff){
  pscore <- c()
  m <- sum(sapply(coeff, length))
  mi <- length(coeff[[i]])
  if (x > r | m-mi+x < r) {
    pscore <- 0
  } else {
    pscore <- ifelse(x > 0, exp(-coeff[[i]][x]), 1)*(psychotools::elementary_symmetric_functions(coeff[-i])[[1]])[r+1-x]/psychotools::elementary_symmetric_functions(coeff)[[1]][r+1]
  }
  pscore
}


#' Item restscore association
#'
#' The observed Gamma coefficient between the score of a single item and the total score of the remaining items
#' is compared with the corresponding expected Gamma coefficient under the Rasch model.
#' @param object an object of class "Rm", a fitted Rasch model oder partial
#' credit model using  the functions RM or PCM in package eRm.
#' @export
#' @return a matrix containing:
#' \item{observed}{observed gamma coefficients}
#' \item{expected}{expected gamma coefficients}
#' \item{se}{standard errors}
#' \item{pvalue}{p values (under normal distribution assumption)}
#' \item{sig}{significance stars: 0  " *** "  0.001  " ** "  0.01  " * "  0.05   " . "  0.1  " "  1}
#' @references Kreiner, S. (2011). A note on item-restscore association in Rasch models.
#' \emph{Applied Psychological Measurement}, 35, 557-561.
#' @author Marianne Mueller
#' @examples
#' rm.mod <- RM(amts[,4:13])
#' item_restscore(rm.mod)
item_restscore <- function(object){
  if (!("Rm"%in%class(object))) stop("object must be of class Rm!")
  k <- dim(object$X)[2]
  n <- dim(object$X)[1]
  if (object$model=="RM") {
    coeff <- (-1)*coef(object)
    mi <- rep(1,k)
  } else {
    mi <- apply(object$X, 2, max, na.rm = TRUE)
    thresh1 <- thresholds(object)[[3]][[1]][, -1] - mean(thresholds(object)[[3]][[1]][, 1])
    coeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
  }
  m <- sum(mi)
  score <- apply(object$X,1,sum,na.rm=T)
  restscore <- matrix(rep(score,k),ncol=k)- object$X
  mm <- t(apply(rbind(object$X,restscore),2,function(x){gamma_coef(table(x[1:n],x[(n+1):(2*n)]))[1:2]}))
  #! beob gammas.se sind nicht gleich wie Svend's
  # erwartete gammas
  rvneu <- factor(score,levels = 0:m)
  npmat <- vector("list", k)
  for (i in 1:k) {
    nmat <- matrix(table(rvneu)[matrix(0:(m - mi[i]) + rep(0:mi[i], each = m - mi[i] + 1), ncol = mi[i] + 1) + 1], ncol = mi[i] + 1)
    pmat <- sapply(0:mi[i], function(x){sapply(0:(m - mi[i]),function(r){pscore_poly(i, x, r + x, coeff)})})
    npmat[[i]] <- nmat*pmat
  }
  expected <- sapply(npmat, function(x){gamma_coef(x)[1]})
  pvalue <- round(ifelse((mm[, 1] - expected) > 0, 2*(1 - pnorm((mm[, 1] - expected)/mm[, 2])), 2*(pnorm((mm[, 1] - expected)/mm[, 2]))), digits = 4)
  mm <- cbind(observed = mm[, 1], expected, se = mm[, 2], pvalue)
  symp <- symnum(pvalue, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols=c("*** ", "** ", "* ", "." , " "))
  mm <- noquote(cbind(format(mm, digits = 3), sig = symp))
  cat("\n")
  mm
}


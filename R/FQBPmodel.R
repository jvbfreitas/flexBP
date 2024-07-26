flexBP = function(formula, data, tol = 0.001, maxiter = 25,
                  compm, linkmu = "log", tuning = 0.5){

  nameslink = c("log", "identity")
  if(all(nameslink != linkmu)){
    stop("the link function is not defined")
  }
  formula = as.formula(formula)
  nformula = all.vars(formula)
  fnames = 0
  jaux = 1
  listaux = list(NULL)
  for(i in 1:length(nformula)){
    if(is.factor(data[,nformula[i]])){
      fnames[jaux] = nformula[i]
      listaux[[jaux]] = "contr.sum"
      jaux = jaux+1
    }
  }
  call <- match.call()
  if(jaux>1){
    names(listaux) = fnames
    X = as.matrix(model.matrix(formula, data = data, contrasts = listaux))
  }else{
    X = as.matrix(model.matrix(formula, data = data))
  }
  p = ncol(X)
  y = model.frame(formula, data = data)[,1]
  # fit_cp <- mcglm(linear_pred = list(formula),
  #                 matrix_pred = compm,
  #                 variance = "tweedie",
  #                 data = data,
  #                 control_algorithm = list(verbose = TRUE,
  #                                          tol = tol,
  #                                          tuning = 0.3, max_iter = 1000))
  iniaux = glm(formula, family = Gamma(link = "log"), data = data)
  beta0 = iniaux$coefficients
  tau0 = c(4, rep(1, length(compm)))
  p = tau0[1]
  pp = length(compm)

  crit = 1
  repeat{
    if(linkmu == "identity"){mi = X%*%beta0}
    if(linkmu == "log"){mi = exp(X%*%beta0)}
    u = y-mi
    vmi = (mi*(1+mi))^p
    Vmi = diag(as.vector(vmi))
    Omega = 0
    for(j in 1:pp){
      Omega = Omega + tau0[j+1]*compm[[j]]
    }
    Sigma = sqrt(Vmi)%*%Omega%*%sqrt(Vmi)
    SII = solve(Sigma)
    if(linkmu == "identity"){D = X}
    if(linkmu == "log"){D = diag(as.vector(mi))%*%X}
    scoreb = t(D)%*%SII%*%u
    Sb = -t(D)%*%SII%*%D
    beta1 = beta0-solve(Sb)%*%scoreb
    if(linkmu == "identity"){mi = X%*%beta1}
    if(linkmu == "log"){mi = exp(X%*%beta1)}
    u = y-mi
    vmi = (mi*(1+mi))^p
    Vmi = diag(as.vector(vmi))
    Sigma = sqrt(Vmi)%*%Omega%*%sqrt(Vmi)
    SII = solve(Sigma)
    if(linkmu == "identity"){D = X}
    if(linkmu == "log"){D = diag(as.vector(mi))%*%X}
    dvdp = diag(as.vector(0.5*(log(mi*(1+mi)))*(mi*(1+mi))^(p/2)))
    dsdp = dvdp%*%Omega%*%sqrt(Vmi) + sqrt(Vmi)%*%Omega%*%dvdp
    ds = list(NULL)
    ds[[1]] = dsdp
    for(j in 1:pp){
      ds[[j+1]] = sqrt(Vmi)%*%compm[[j]]%*%sqrt(Vmi)
    }
    wp = list(NULL)
    for(j in 1:(pp+1)){
      wp[[j]] = SII%*%ds[[j]]%*%SII
    }
    scorep = 0
    for(j in 1:(pp+1)){
      scorep[j] = sum(diag(as.matrix(wp[[j]])%*%(u%*%t(u)-as.matrix(Sigma))))
    }
    St = matrix(0,pp+1,pp+1)
    for(i in 1:(pp+1)){
      for(j in 1:(pp+1)){
        St[i,j] = -sum(diag(as.matrix(wp[[i]])%*%as.matrix(Sigma)%*%as.matrix(wp[[j]])%*%as.matrix(Sigma)))
      }
    }
    tau1 = tau0-tuning*solve(St)%*%scorep
    p = tau1[1]
    crit = crit+1
    if(crit==maxiter||abs(sum(beta1-beta0))+abs(sum(tau1-tau0))<=tol){
      break
    }
    beta0 = beta1
    tau0 = tau1
  }

  dvdk = (p/2)*((mi*(1+mi))^(p/2-1))*(1+2*mi)*mi
  scoretaus = wp
  scorebetas = list(NULL)
  n = length(u)
  np = ncol(X)
  for(i in 1:np){
    scorebetas[[i]] = diag(as.vector(dvdk*X[,i]),n)%*%Omega%*%sqrt(Vmi) + sqrt(Vmi)%*%Omega%*%diag(as.vector(dvdk*X[,i]),n)
  }
  Stb = matrix(0,length(tau0), np)
  for(i in 1:np){
    for(j in 1:length(tau0)){
      Stb[j,i] = -sum(diag(as.matrix(scoretaus[[j]]%*%Sigma%*%scorebetas[[i]]%*%Sigma)))
    }
  }
  Stheta = cbind(rbind(Sb, Stb),rbind(matrix(0,nrow(Sb),ncol(St)), St))
  Vb = -Sb
  scoretau = scorep
  Vtau = scoretau%*%t(scoretau)
  Vtb = scoretau%*%t(scoreb)
  Vtheta = cbind(rbind(Vb,Vtb), rbind(t(Vtb),Vtau))
  Js = ginv(as.matrix(Stheta))%*%Vtheta%*%t(ginv(as.matrix(Stheta)))
  varbeta = diag(Js[1:np,1:np])
  if(linkmu == "identity"){mi = X%*%beta0}
  if(linkmu == "log"){mi = exp(X%*%beta0)}
  u = y-mi
  vmi = (mi*(1+mi))^p
  Vmi = diag(as.vector(vmi))
  Omega = 0
  for(j in 1:pp){
    Omega = Omega + tau0[j+1]*compm[[j]]
  }
  Sigma = sqrt(Vmi)%*%Omega%*%sqrt(Vmi)
  fit <- list()
  attr(fit, "class") <- c("FLEXBP")
  fit$title <- "FLEXBP:  FLEXIBLE QUASI-BETA PRIME MODEL"
  fit$model <- list()
  fit$model$link <- linkmu
  fit$model$varfun <- "[mu(1+mu)]^p"
  fit$model$compm <- compm
  fit$call <- call
  fit$formula <- formula
  fit$ndata = n
  fit$iterations <- crit
  fit$coefficients <- beta0
  eta <- as.vector(X %*% fit$coefficients)
  fit$linear.predictors <- eta
  mu <- as.vector(mi)
  fit$fitted.values <- mi
  fit$residuals <- y-mi
  fit$family <- "Quasi-Beta prime"
  fit$y <- as.vector(y)
  fit$tau <- tau0
  fit$robust.variance <- varbeta
  fit$robust.se = sqrt(varbeta)
  fit$comp$X = X
  fit$comp$Sigma = Sigma
  fit$comp$SII = solve(Sigma)
  npa = (ncol(X)+length(tau0))
  GP = -0.5*t(u)%*%SII%*%u - determinant(as.matrix(Sigma))$modulus -n/2*log(2*pi)
  AIC = -2*GP + 2*(ncol(X)+length(tau0));AIC
  BIC = -2*GP + (ncol(X)+length(tau0))*log(n);BIC
  fit$AIC = AIC
  fit$BIC = BIC
  return(fit)
}

print.flexBP = function(x, digits = NULL, quote = FALSE, prefix = "", ...){
  if(is.null(digits)) digits = options()$digits else options(digits =
                                                               digits)
  cat("\n", x$title)
  cat("n"," Model:\n")
  cat(" Link (mean):", x$model$link, "\n")
  cat("\nCall:\n")
  dput(x$call)
  ys = matrix(0, nrow = x$ndata, ncol = 4)
  ys[, 1] = x$y
  ys[, 2] = x$linear.predictors
  ys[, 3] = x$fitted.values
  ys[, 4] = x$residuals
  dimnames(ys) = list(1:length(x$y), c("Y", "LP", "fitted",
                                       "Residual")) #       cat("\nFitted Values:\n")
  cat("\nNumber of observations : ", x$nobs, "\n")
  cat("\n\nCoefficients (mean):\n")
  print(t(x$coefficients), digits = digits)
  cat("\n\nCoefficients (covariance):\n")
  print(t(x$tau), digits = digits)
  cat("\nNumber of Iterations: ", x$iterations)
  invisible(x)
}

residuals.flexBP <- function (x, type = c("standard", "pearson"), ...)
{
  type <- match.arg(type)
  y   <- x$y
  rs   <- x$residuals
  mi = x$fitted.values
  Lamb = diag(as.vector(mi))
  SII = x$comp$SII
  W = as.matrix(Lamb%*%SII%*%t(Lamb))
  sw = sqrtm(W)$B
  X = x$comp$X
  H = sw%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sw
  hij = diag(H)
  r = sw%*%solve(Lamb)%*%rs
  rij = r/sqrt(1-hij)
  rij = as.vector(rij)
  res <- switch(type,
                standard = rs,
                pearson = rij)
  res
}


LI.flexBP = function(model, pert = c("case-weight", "response")){
  y   <- model$y
  rs   <- as.vector(model$residuals)
  mi = model$fitted.values
  Lamb = diag(as.vector(mi))
  SII = model$comp$SII
  W = as.matrix(Lamb%*%SII%*%t(Lamb))
  sw = sqrtm(W)$B
  X = model$comp$X
  H = sw%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sw
  hij = diag(H)
  if(pert == "case-weight"){
    aux = diag(rs)%*%solve(Lamb)%*%t(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%solve(Lamb)%*%diag(rs)
    vaux = eigen(aux)$vectors[,1]
  }
  if(pert == "response"){
    B = sqrt(diag(as.matrix(model$comp$Sigma)))
    aux = diag(B)%*%solve(Lamb)%*%t(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%solve(Lamb)%*%diag(B)
    vaux = eigen(aux)$vectors[,1]
  }
  return(vaux)
}

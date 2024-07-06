# This function evalutates the estimation accuracy
# You can set the correct number of factors r 
# otherwise r will be selected via the eigen-ratio method
analysis = function(mMatArr, Q, r = NULL) {
  if(length(dim(mMatArr)) == 3)
    mMat = apply(mMatArr, c(1,2), sum)
  else
    mMat = mMatArr %*% t(mMatArr)
  p = dim(mMatArr)[1]
  ## threshold on each slice of mMatArr
  mEig = eigen(mMat)
  val = mEig$values
  vec = mEig$vectors
  
  R = p/2
  ratio = val[2:p] / val[1:p-1]
  rr = which.min(ratio[1:R])
  if(!is.null(r)) rhat = r
  else rhat = rr
  hatA = vec[, 1:rhat]
  
  temp = hatA%*%t(hatA)%*%Q%*%t(Q)
  err = sqrt(1 - sum(diag(temp)) / max(ncol(Q), rhat))
  
  c(rr, err)
}
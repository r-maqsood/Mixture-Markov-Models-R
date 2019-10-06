
## calls new code for estimation step ... sequence classification
## uses same alpha (or mixture weights) for all mixture components ... a vector of size K
classifyMixtureMarkov2 <- function (testData, K, mixMarkovModel, states, missingStates){  

 if (K < 1) stop("Wrong number of mixture components K...\n")

  lenStates <- length(states)
  lenMissingStates <- length(missingStates)

  if(lenMissingStates > lenStates) stop("Missing states vector can't be more than original states...\n")

  ## create empty vector, find and store indices of missing states ... this vector will be used to fill x arrays ... Pi resolved automatically
  vec = c()
  for(i in 1: lenMissingStates){
    vec[i] = which(states==missingStates[i])
  }
  print("missing States ....")
  print(paste("Index: ",vec))
  print(paste("lenMissing: ",lenMissingStates))


  xdims <- dim(testData)
  p <- xdims[1] + length(vec)   ## states equal p + no. of missing states
  n <- xdims[3]
  
  MMdims <- dim(mixMarkovModel)
  MMp <- MMdims[1]
  MMK <- MMdims[3]


  ## same weight will be assigned to each mixture in libEM.c ... classify function
  alpha <- rep(0, K)
  

  if (K != MMK) stop("K doesn't match with mixture markov model components\n")


  x1 <- as.vector(testData)  # converts all transition matrices into vector... stores column wise
  id <- rep(0, n)

  MM <- as.vector(mixMarkovModel)  # converts all transition matrices into vector... stores column wise
  

  gamma <- rep(-1, p*p*K)
  z <- rep(0, n*K)

 
  Q <- .C("classifier_MMM2", p1 = as.integer(p), K1 = as.integer(K), n1 = as.integer(n), x1 = as.double(x1), alpha = as.double(alpha), Pi1 = as.double(MM), gamma1 = as.double(z), id = as.integer(id), missingStates = as.integer(vec), lenMissing1 = as.integer(lenMissingStates), PACKAGE = "ClickClust")


  a <- array(Q$Pi1, c(K, p, p))
  b <- array(NA, c(p, p, K))
  for (i in 1:K){
    b[,,i] <- t(a[i,,])
  }
  b[b == -1] <- NA
  
  
  ret <- list(z = t(matrix(Q$gamma1, nrow = K)), id = Q$id + 1, gamma = b) 
  
  class(ret) <- "EM"
  return(ret)
  
}

rm(list=ls())

library('ltmle')

get.data <- function(n, include.counterfactuals=F){
  
  # confounders
  W1 <- rnorm(n, 0, .5)
  W2 <- rnorm(n, 0, .5)
  W3 <- rnorm(n, 0, .5)
  
  # instrumental variable
  I <- rnorm(n, 0, 0.5)
  
  # unmeasured common causes
  U1 <- runif(n)  # M & A
  U2 <- runif(n)  # M & Y
  
  # M variable
  M <- runif(n, -U1, U2)
  
  # propensity score & exposure
  pA <- plogis(.7*(W1+W2+W3) +2.5*I)
  A <- as.numeric( pA < U1)
  
  # outcome
  UY <- runif(n, 0, 1)
  Y1 <- generate.Y(W1, W2, W3, A=1, U2, UY)
  Y0 <- generate.Y(W1, W2, W3, A=0, U2, UY)
  
  if(include.counterfactuals){
    print(summary(pA) )
    hist(pA)
    df <- data.frame(cbind(Y1,Y0))
  } else{
    Y <- rep(NA, n)
    Y[A==1] <- Y1[A==1]
    Y[A==0] <- Y0[A==0]
    df <- data.frame(cbind(W1,W2,W3,I,M,A,Y))
  }
  df
}

generate.Y <- function(W1,W2,W3,A, U2, UY){
  ct <- 1*W1 + 1.5*W2 +  1*W3 + 0.1*UY
  ind <- U2>.5
  ind*plogis(2.5*A+ ct) + (1-ind)*plogis(-.5*A - ct)
}

get.gcomp <- function(X, ATE, Qform=NULL, this){
  
  est <-ltmle(data=X, Anodes='A', Ynodes='Y', abar=list(1,0), 
              Qform=Qform, gcomp=T,
              estimate.time=F, variance.method='ic') 
  est <- summary(est)$effect.measures$ATE
  estimate <- est$estimate
  bias <- estimate - ATE
  cover <- (est$CI[1] <= ATE & ATE <= est$CI[2])
  get.col.names(estimate, bias, cover, this)
}


get.ipw <- function(df, pA, ATE, this){
  wt <- (df$A==1)/pA - (df$A==0)/(1-pA)
  estimate <- mean(wt*df$Y)
  IC <- wt*df$Y - estimate
  se <- sqrt(var(IC)/nrow(df))
  cover <- (estimate-1.96*se <= ATE & ATE <= estimate+1.96*se)
  bias <- estimate - ATE
  get.col.names(estimate, bias, cover, this)
}

get.col.names <- function(estimate, bias, cover, this){
  out <- data.frame(cbind(est=estimate, bias=bias, cov=cover))
  colnames(out) <- paste(this, colnames(out), sep='.')
  out
}

do.estimation <- function(n, ATE){
  df <- get.data(n)
  
  # unadjusted
  unadj <- get.gcomp(X=subset(df, select=c(A,Y)), ATE=ATE, this='unadj' )
  
  # ipw naive 
  pA <- predict(glm(A~W1+W2+W3+I+M, family='binomial', data=df), type='response')
  ipw.naive <- get.ipw(df=df, pA=pA, ATE=ATE, this='ipw.naive')
  
  # ipw informed
  pA <- predict(glm(A~W1+W2+W3, family='binomial', data=df), type='response')
  ipw.roadmap <- get.ipw(df=df, pA=pA,ATE=ATE, this='ipw.roadmap')
  
  # gcomp naive
  gcomp.naive <- get.gcomp(X=df, ATE=ATE, this='gcomp.naive')
  
  # gcomp informed
  gcomp.roadmap <- get.gcomp(X=subset(df, select=-c(M,I)), ATE=ATE, this='gcomp.roadmap')
  
  # output
  yay <- data.frame(cbind(unadj,  gcomp.naive, ipw.naive,  gcomp.roadmap, ipw.roadmap))
  yay
}

# true effect
pop <- get.data(n=500000, include.counterfactuals=T)
ATE <- mean(pop$Y1) - mean(pop$Y0)
ATE

R <- 1000 # iterations
n <- 1000 # sample size
out <- data.frame(matrix(NA, nrow=R, ncol=15))
for(r in 1:R){
  temp <- do.estimation(n=n, ATE=ATE)
  out[r,] <- temp
}
colnames(out) <- colnames(temp)
x <- colMeans(out)
y <- data.frame(matrix(x, ncol=3, byrow=T))
rownames(y) <- c('Unadjusted', 
                 'Gcomp-naive', 'IPW-naive', 'Gcomp-roadmap','IPW-roadmap')
colnames(y) <- c('Ave. Est', 'Bias', 'Coverage')
write.csv(y, file=paste('out',Sys.Date(), 'csv', sep='.'))
save(out, ATE, file=paste('out',Sys.Date(), 'Rdata', sep='.'))
round(ATE*100,1)
round(y*100,1)

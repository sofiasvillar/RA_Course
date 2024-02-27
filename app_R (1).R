#Credits:
#We would like to extend our gratitude to Lukas Pin and Colin Starr from the MRC Biostatistics Unit (MRC BSU) for their contributions to the development of this App. 
#Lukas Pin and Colin Starr's involvement ranged from conceptualizing the analytical frameworks to providing detailed feedback on script functionality. 

#Acknowledgements:
#Lukas Pin: For his expertise in the statistical methods and the implementation of the approaches encoded in these scripts.
#Colin Starr: For detailed feedback and improvement of the script functionality, which enhanced the quality and usability of the scripts.

wald.test.binary <- function(x0,x1, measure="simple mean difference"){
  p0 <- mean(x0); p1 <- mean(x1)
  q0 <- 1-p0 ; q1 <- 1-p1
  sdx0 <- p0*q0; sdx1 <- p1*(1-p1)
  n0 <- length(x0); n1 <- length(x1)

  if(measure=="simple mean difference"){
    if((sdx0==0 && sdx1==0)){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- (p0-p1)/sqrt(sdx0/n0+sdx1/n1)
    }
  } else if(measure=="rr") {
    if(p0==1 || (p0==0 && p1*q1==0)){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    }
    Z <- (q1/q0-1)/sqrt(p0/(n0*q0^3)+p1*q1/n1)
  } else if(measure=="or"){
    if((q0*p1==1) || (p0==1 && p1==1) || (p0==0 && p1==0)){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    }
    else{
      Z <- (p0*q1/(q0*p1)-1)/sqrt(p0/(n0*q0^3)+p1/(n1*q1^3))
    }
  } else if(measure=="log relative risk"){
    if(p0==1 || p1==1 || (p0==0 && p1==0) ){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- log(q1/q0)/ sqrt(p0/(n0*q0)+p1/(n1*q1))
    }
  } else if(measure=="lor"){
    if(p0==0 || q0==0 || p1==0 || q1==0){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- log(q1/q0)/ sqrt(p0/(n0*q0)+p1/(n1*q1))
    }
  }
  return(2*pnorm(-1*abs(Z)))
}

superior <- function(par1=c(0,0), par2=FALSE, dist="bern"){
    if(par1[1]<par1[2]){
      return(1)
    }
    if(par1[1]>par1[2]){
      return(0)
    }
    if(par1[1]==par1[2]){
      return(2)
    }
}


### Equal Randomisation ###
two_arm_ER <- function(N=200, p = c(0,0),  burnin=3, measure="simple mean difference", Z_corrected=FALSE){
  # Generate Samples 
  n1 = min(max(rbinom(1, N, 0.5),burnin), N-burnin)
  n0 = N - n1
  x0 <- rbinom(n0,1, p[1]) 
  x1 <- rbinom(n1, 1, p[2])
  
  # Response
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist=dist,par1 = p)
  if(sup==2){ #2 means that no arm is theoretically superior 
    #   imb <- NA
    per.sup <- NA
  } else {
    if(sup==0){ 
      #    low <- n1 ; lar <- n0
      per.sup <- n0/N
    }
    if(sup==1){ 
      #   low <- n0 ; lar <- n1
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test 
  Z_P <-wald.test.binary(x0,x1, measure = measure)
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, per.sup))
}

### PW ###
PW = function(N=200, p = c(0,0),  burnin=3, measure="simple mean difference", Z_corrected=FALSE){
  
  allocation = response = rep(0,N)
  if(burnin>0){
    allocation[1:burnin] = 0
    response[1:burnin] = rbinom(burnin, 1, p[1])
    allocation[(burnin+1):(2*burnin)] = 1
    response[(burnin+1):(2*burnin)] = rbinom(burnin, 1, p[2])
  }
  
  for (i in (2*burnin+1):N){
    if(response[i-1]==1){
      allocation[i]=allocation[i-1]
    }
    else{ 
      allocation[i]=abs(allocation[i-1]-1)
    }
    response[i] = rbinom(1,1,p[allocation[i]+1])
  }
  
  # Response
  x0 <- response[allocation==0]; x1 <- response[allocation==1]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist="bern",par1 = p)
  if(sup==2){ #2 means that no arm is theoretically superior 
    #   imb <- NA
    per.sup <- NA
  } else {
    if(sup==0){ 
      #    low <- n1 ; lar <- n0
      per.sup <- n0/N
    }
    if(sup==1){ 
      #   low <- n0 ; lar <- n1
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test 
  Z_P <-wald.test.binary(x0,x1, measure = measure)
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, per.sup))
}

### RPW ###
randomiseRPW = function(N=200, p = c(0,0),  burnin=3, measure="simple mean difference", Z_corrected=FALSE){
  
  allocation = response = rep(0,N)
  if(burnin>0){
    allocation[1:burnin] = 0
    response[1:burnin] = rbinom(burnin, 1, p[1])
    allocation[(burnin+1):(2*burnin)] = 1
    response[(burnin+1):(2*burnin)] = rbinom(burnin, 1, p[2])
  }
  #urn = c(0,1, rep(NA,N)) without burnin
  urn = c(rep(0, burnin),rep(1, burnin), rep(NA,(N-2*burnin)))
  
  for (i in (2*burnin+1):N){
    ball = sample(urn[1:(i-1)], 1) # chose one ball
    allocation[i] = ball
    response[i] = rbinom(1,1,p[allocation[i]+1])
    urn[i] = ifelse(response[i], ball, 1-ball)
  }
  
  # Response
  x0 <- response[allocation==0];  x1 <- response[allocation==1]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist="bern",par1 = p)
  if(sup==2){ #2 means that no arm is theoretically superior 
    #   imb <- NA
    per.sup <- NA
  } else {
    if(sup==0){ 
      #    low <- n1 ; lar <- n0
      per.sup <- n0/N
    }
    if(sup==1){ 
      #   low <- n0 ; lar <- n1
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test 
  Z_P <-wald.test.binary(x0,x1, measure = measure)
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, per.sup))
}

### Bayesian RAR ###

# define a function to select the arm with the highest expected value
select_arm <- function(num_arms, success_counts, failure_counts) {
  # generate a beta distribution for each arm
  values <- sapply(1:num_arms, function(i) generate_beta(success_counts[i] + 1, failure_counts[i] + 1))
  # select the arm with the highest expected value
  return(which.max(values))
}

generate_beta <- function(alpha, beta) {
  rbeta(1, alpha, beta)
}

BRAR = function(N=200, p = c(0,0),burnin=3, measure="simple mean difference", Z_corrected=FALSE){
  
  A = rep(NA, N)   # Vector of allocations
  X = rep(NA, N)   # Vector of responses
  
  #Burnin
  if(burnin>0){
    A[1:burnin] = 1
    X[1:burnin] = rbinom(burnin, 1, p[1])
    A[(burnin+1):(2*burnin)] = 2
    X[(burnin+1):(2*burnin)] = rbinom(burnin, 1, p[2])
  }

  for (i in (2*burnin+1):N){
    #Count Successes and Failures
    success_counts <- sapply(1:2, function(i) sum(X[A==i], na.rm=TRUE))
    failure_counts <- sapply(1:2, function(i) length(X[A==i])) - success_counts
    A[i] = select_arm(num_arms=2, success_counts, failure_counts)
    X[i] = rbinom(1,1,p[A[i]])
  }
  
  # Response
  x0 <- X[A==1];  x1 <- X[A==2]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist="bern",par1 = p)
  if(sup==2){ #2 means that no arm is theoretically superior 
    #   imb <- NA
    per.sup <- NA
  } else {
    if(sup==0){ 
      #    low <- n1 ; lar <- n0
      per.sup <- n0/N
    }
    if(sup==1){ 
      #   low <- n0 ; lar <- n1
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test 
  Z_P <-wald.test.binary(x0,x1, measure = measure)
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, per.sup))
}

two_arm_SMLE <- function(N=200, dist="bern",  par1 = c(0,0), measure="simple mean difference", burnin=3, ar="WMW", Z_corrected=FALSE, first=FALSE){
  
  A = rep(NA, N)   # Vector of allocations
  X = rep(NA, N)   # Vector of responses
  est <- rep(NA, N) 
  
  #Burnin
  if(burnin>0){
    A[1:burnin] = 1
    X[1:burnin] = rbinom(burnin, 1, par1[1])
    A[(burnin+1):(2*burnin)] = 2
    X[(burnin+1):(2*burnin)] = rbinom(burnin, 1, par1[2])
  }
  
  
  for (i in (2*burnin+1):N){
    n1 = sum(A==1, na.rm = TRUE)
    n2 = sum(A==2, na.rm = TRUE)
    
    if(ar=="minF"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if(measure=="simple mean difference"){
          if((sqrt(p1.hat) + sqrt(p2.hat)) == 0){
            rho1.hat = 0.5
          } else {
            rho1.hat = sqrt(p1.hat)/(sqrt(p1.hat) + sqrt(p2.hat))
          }
        } else if(measure=="log relative risk"){
          if((sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat)) == 0){
            rho1.hat = 0.5
          } else {
            rho1.hat = 1- sqrt(p2.hat)*(1-p1.hat) / (sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat))
          }
        }
        p = 1-rho1.hat
    }
    if(ar=="Neyman"){
      if(measure=="simple mean difference"){
        sig1_N <- sd(X[A==1], na.rm = TRUE)
        sig2_N <- sd(X[A==2], na.rm = TRUE)
        if((sig1_N+sig2_N)==0){
          p <- 0.5
        } else {
          p <- 1 - sig1_N/(sig1_N+sig2_N)
        }
      } else if(measure=="log relative risk"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if((sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))==0){
          p <- 0.5
        } else {
          p <- sqrt((1-p1.hat)*p2.hat) / (sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))
        }
      }
      est[i] <- p
    }
    if(ar=="AD" && measure =="simple mean difference"){
      s1 = sum(X[A==1], na.rm = TRUE)
      p1.hat = (s1)/(n1)
      s2 = sum(X[A==2], na.rm = TRUE)
      p2.hat = (s2)/(n2)
      if((p1.hat + p2.hat) == 0){
        rho1.hat = 0.5
      } else {
        rho1.hat = p1.hat/(p1.hat + p2.hat)
      }
      p = 1-rho1.hat
    }
    if(p==0){p <- 1/N}
    if(p==1){p <- 1-1/N}
    est[i] <- p
    A[i] = rbinom(1, 1, p) + 1
    X[i]  <- rbinom(1,1,par1[A[i]]) 
  }
  
  # Response
  x0 <- X[A==1];  x1 <- X[A==2]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist=dist,par1 = par1)
  if(sup==2){ #2 means that no arm is theoretically superior 
    #   imb <- NA
    per.sup <- NA
  } else {
    if(sup==0){ 
      #    low <- n1 ; lar <- n0
      per.sup <- n0/N
    }
    if(sup==1){ 
      #   low <- n0 ; lar <- n1
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test 
  Z_P <-wald.test.binary(x0,x1, measure = measure)
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, per.sup))
}

two_arm_DBCD = function(N=200, dist="bern",  par1 = c(0,0), measure="simple mean difference", burnin=3, ar="minF", Z_corrected=FALSE, first=FALSE){
  
  A = rep(NA, N)   # Vector of allocations
  X = rep(NA, N)   # Vector of responses
  est <- rep(NA, N) 
  
  #Burnin
  if(burnin>0){
    A[1:burnin] = 1
    X[1:burnin] = rbinom(burnin, 1, par1[1])
    A[(burnin+1):(2*burnin)] = 2
    X[(burnin+1):(2*burnin)] = rbinom(burnin, 1, par1[2])
    est[1:(2*burnin)] <- 0.5
  }
  
  
  for (i in (2*burnin+1):N){
    n1 = sum(A==1, na.rm = TRUE)
    n2 = sum(A==2, na.rm = TRUE)
    
    if(ar=="minF"){
      s1 = sum(X[A==1], na.rm = TRUE)
      p1.hat = (s1)/(n1)
      s2 = sum(X[A==2], na.rm = TRUE)
      p2.hat = (s2)/(n2)
      if(measure=="simple mean difference"){
        if((sqrt(p1.hat) + sqrt(p2.hat)) == 0){
          rho1.hat = 0.5
        } else {
          rho1.hat = sqrt(p1.hat)/(sqrt(p1.hat) + sqrt(p2.hat))
        }
      } else if(measure=="log relative risk"){
        if((sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat)) == 0){
          rho1.hat = 0.5
        } else {
          rho1.hat = 1- sqrt(p2.hat)*(1-p1.hat) / (sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat))
        }
      }
      p = 1-rho1.hat
    }
    if(ar=="Neyman"){
      if(measure=="simple mean difference"){
        sig1_N <- sd(X[A==1], na.rm = TRUE)
        sig2_N <- sd(X[A==2], na.rm = TRUE)
        if((sig1_N+sig2_N)==0){
          p <- 0.5
        } else {
          p <- 1 - sig1_N/(sig1_N+sig2_N)
        }
      } else if(measure=="log relative risk"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if((sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))==0){
          p <- 0.5
        } else {
          p <- sqrt((1-p1.hat)*p2.hat) / (sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))
        }
      }
      est[i] <- p
    }
    if(ar=="AD" && measure =="simple mean difference"){
      s1 = sum(X[A==1], na.rm = TRUE)
      p1.hat = (s1)/(n1)
      s2 = sum(X[A==2], na.rm = TRUE)
      p2.hat = (s2)/(n2)
      if((p1.hat + p2.hat) == 0){
        rho1.hat = 0.5
      } else {
        rho1.hat = p1.hat/(p1.hat + p2.hat)
      }
      p = 1-rho1.hat
    }
    if(p==0){p <- 1/N}
    if(p==1){p <- 1-1/N}
    est[i] <- p
    rho.hat <- c(1-p,p)
    
    alloc.prop = c(n1/i, n2/i)
    
    phi1 = DBCD(rho.hat, alloc.prop, gamma=2)
    
    
    A[i] = rbinom(1, 1, phi1[2]) + 1
    X[i]  <- rbinom(1,1,par1[A[i]]) 
  }
  
  # Response
  x0 <- X[A==1];  x1 <- X[A==2]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist="bern",par1 = par1)
  if(sup==2){ #2 means that no arm is theoretically superior 
    #   imb <- NA
    per.sup <- NA
  } else {
    if(sup==0){ 
      #    low <- n1 ; lar <- n0
      per.sup <- n0/N
    }
    if(sup==1){ 
      #   low <- n0 ; lar <- n1
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test 
  Z_P <-wald.test.binary(x0,x1, measure = measure)
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, per.sup))
}

DBCD = function(rho.hat, alloc.prop, gamma=2){
  
  num = rho.hat*(rho.hat/alloc.prop)^gamma
  den = sum(num)
  
  return(num/den)
  
}

two_arm_ERADE = function(N=200, dist="bern",  par1 = c(0,0), measure="simple mean difference", burnin=3, ar="minF", Z_corrected=FALSE, first=FALSE){
  
  A = rep(NA, N)   # Vector of allocations
  X = rep(NA, N)   # Vector of responses
  est <- rep(NA, N) 
  
  #Burnin
  if(burnin>0){
    A[1:burnin] = 1
    X[1:burnin] = rbinom(burnin, 1, par1[1])
    A[(burnin+1):(2*burnin)] = 2
    X[(burnin+1):(2*burnin)] = rbinom(burnin, 1, par1[2])
    est[1:(2*burnin)] <- 0.5
  }
  
  for (i in (2*burnin+1):N){
    n1 = sum(A==1, na.rm = TRUE)
    n2 = sum(A==2, na.rm = TRUE)
    
    if(ar=="minF"){
      s1 = sum(X[A==1], na.rm = TRUE)
      p1.hat = (s1)/(n1)
      s2 = sum(X[A==2], na.rm = TRUE)
      p2.hat = (s2)/(n2)
      if(measure=="simple mean difference"){
        if((sqrt(p1.hat) + sqrt(p2.hat)) == 0){
          rho1.hat = 0.5
        } else {
          rho1.hat = sqrt(p1.hat)/(sqrt(p1.hat) + sqrt(p2.hat))
        }
      } else if(measure=="log relative risk"){
        if((sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat)) == 0){
          rho1.hat = 0.5
        } else {
          rho1.hat = 1- sqrt(p2.hat)*(1-p1.hat) / (sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat))
        }
      }
      p = 1-rho1.hat
    }
    if(ar=="Neyman"){
      if(measure=="simple mean difference"){
        sig1_N <- sd(X[A==1], na.rm = TRUE)
        sig2_N <- sd(X[A==2], na.rm = TRUE)
        if((sig1_N+sig2_N)==0){
          p <- 0.5
        } else {
          p <- 1 - sig1_N/(sig1_N+sig2_N)
        }
      } else if(measure=="log relative risk"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if((sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))==0){
          p <- 0.5
        } else {
          p <- sqrt((1-p1.hat)*p2.hat) / (sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))
        }
      }
      est[i] <- p
    }
    if(ar=="AD" && measure =="simple mean difference"){
      s1 = sum(X[A==1], na.rm = TRUE)
      p1.hat = (s1)/(n1)
      s2 = sum(X[A==2], na.rm = TRUE)
      p2.hat = (s2)/(n2)
      if((p1.hat + p2.hat) == 0){
        rho1.hat = 0.5
      } else {
        rho1.hat = p1.hat/(p1.hat + p2.hat)
      }
      p = 1-rho1.hat
    }
    if(p==0){p <- 1/N}
    if(p==1){p <- 1-1/N}
    est[i] <- p
    alloc.prop = n1/i
    
    phi1 = ERADE(1-p, alloc.prop) #USE AUC here instead 
    A[i] = rbinom(1, 1, 1-phi1) + 1
    X[i]  <- rbinom(1,1,par1[A[i]]) 
  }
  
  
  # Response
  x0 <- X[A==1];  x1 <- X[A==2]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist="bern",par1 = par1)
  if(sup==2){ #2 means that no arm is theoretically superior 
    #   imb <- NA
    per.sup <- NA
  } else {
    if(sup==0){ 
      #    low <- n1 ; lar <- n0
      per.sup <- n0/N
    }
    if(sup==1){ 
      #   low <- n0 ; lar <- n1
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test 
  Z_P <-wald.test.binary(x0,x1, measure = measure)
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, per.sup))
}


###############################################################################

ERADE = function(rho.hat, alloc.prop, alpha=0.5){
  
  if(alloc.prop > rho.hat){
    p = alpha*rho.hat
  } else if(alloc.prop < rho.hat){
    p = 1 - alpha*(1 - rho.hat)
  } else{
    p = rho.hat
  }
  
  return(p)
}


### Simulation ###
simulation <- function(N=366, p = c(0.941,0.991), burnin=3, nsim= 10^4, method="ER", measure="simple mean difference", proportion=NULL){
  
  ########## Theoretical Variance Calculation ############################
  fun <- function(x){x*(1-x)}
  variance <-  data.frame(lapply(p,fun))
  
  sim = matrix(nrow = nsim, ncol = 3)
  ########## Select Method ############################
  for(i in 1:nsim){
    if(method=="ER"){
      sim[i,] <- two_arm_ER(N=N,  p =p, measure=measure, burnin = burnin, Z_corrected=FALSE)
    } 
    else if(method=="PW"){
      sim[i,] <- PW(N=N,  p =p,  burnin = burnin, measure=measure, Z_corrected=FALSE)
    }
    else if(method=="RPW"){
      sim[i,] <- randomiseRPW(N=N,  p =p,  burnin = burnin, measure=measure, Z_corrected=FALSE)
    }
    else if(method=="BRAR"){
      sim[i,] <- BRAR(N=N,  p =p,  burnin = burnin, measure=measure, Z_corrected=FALSE)
    }
    else if(method=="SMLE Neyman"){
      sim[i,] <- two_arm_SMLE(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="Neyman")
    }
    else if(method=="SMLE minF"){
      sim[i,] <- two_arm_SMLE(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="minF")
    }
    else if(method=="SMLE AD"){
      if(measure=="simple mean difference"){
        sim[i,] <- two_arm_SMLE(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="AD")
      } else { }
    }
    else if(method=="DBCD Neyman"){
      sim[i,] <- two_arm_DBCD(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="Neyman")
    }
    else if(method=="DBCD minF"){
      sim[i,] <- two_arm_DBCD(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="minF")
    }
    else if(method=="DBCD AD"){
      if(measure=="simple mean difference"){
        sim[i,] <- two_arm_DBCD(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="AD")
      } else { }
    }
    else if(method=="ERADE Neyman"){
      sim[i,] <- two_arm_ERADE(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="Neyman")
    }
    else if(method=="ERADE minF"){
      sim[i,] <- two_arm_ERADE(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="minF")
    }
    else if(method=="ERADE AD"){
      if(measure=="simple mean difference"){
        sim[i,] <- two_arm_ERADE(N=N,  par1 =p,  measure=measure, burnin = burnin, ar="AD")
      } else { }
    }
  }
  
  ### Output
  reject_Z <- sum(sim[,2] < 0.05)
  power <- reject_Z/nsim
  power_err <- sqrt(power*(1-power)/nsim)
  power <- round(power*100,3)
  per.sup <- round(mean(sim[,3])*100,2)
  emr <- round(mean(sim[,1]),4)
  var.per.sup <- round(var(sim[,3])*10000)
  var.emr <- round(var(sim[,1]),4)

  list(power, per.sup, var.per.sup, emr, var.emr, round(power_err*100,2))
}

### Results ###
# Enter Burn-In per ARM!
#simulation(N=120, p = c(0.1,0.3), burnin=10, nsim= 10^4, method="ERADE AD", measure="simple mean difference") 

library(shiny)
library(shinyjs)

#library(Rcpp)
#sourceCpp("code.cpp")

# Necessary to reproduce R behaviour in Rcpp::sample
#RNGkind(sample.kind = "Rounding")

ui <- fluidPage(
  
  useShinyjs(),
  
  titlePanel("RAR App for binary endpoints"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("N", "Sample Size", min = 50, max = 2000, value = 366),
      sliderInput("burnin", "Burn-In per Arm", min = 1, max = 500, value = 30),
      sliderInput("nsim", "Number of Simulations", min = 0, max = 10000, value = 5000, step = 100),
      selectInput("method", "Method", choices = c("ER", "PW", "RPW", "BRAR", 
                                                  "SMLE Neyman", "SMLE minF", "SMLE AD", 
                                                  "DBCD Neyman", "DBCD minF", "DBCD AD",
                                                  "ERADE Neyman", "ERADE minF", "ERADE AD")),
      selectInput("measure", "Measure", choices = c("simple mean difference", "log relative risk")),
      numericInput("p0", "Success Probability Arm 0", value = 0.941, min = 0, max = 1, step = 0.001),
      numericInput("p1", "Success Probability Arm 1", value = 0.991, min = 0, max = 1, step = 0.001),
      actionButton("runSimulation", "Run Simulation")
    ),
    mainPanel(
      uiOutput("resultsText")
    )
  )
)

server <- function(input, output, session) {
  # Define a reactive expression that triggers when the 'runSimulation' button is clicked
  results <- eventReactive(input$runSimulation, {
    
    shinyjs::disable('runSimulation')
    
    # Use p0 and p1 directly
    probabilities <- c(input$p0, input$p1)
    
    # Call the simulation function and capture the output
    
    method = switch(input$method,
                    "ER"=0,
                    "PW"=1,
                    "RPW"=2,
                    "BRAR"=3,
                    "SMLE Neyman"=4,
                    "SMLE minF"=5,
                    "SMLE AD"=6,
                    "DBCD Neyman"=7,
                    "DBCD minF"=8,
                    "DBCD AD"=9,
                    "ERADE Neyman"=10,
                    "ERADE minF"=11,
                    "ERADE AD"=12
                    )
    measure = switch(input$measure,
                     "simple mean difference"=0,
                     "log relative risk"=1
                     )

    out = simulation(N = input$N, p = probabilities, burnin = input$burnin, nsim = input$nsim, method = input$method, measure=input$measure)
    
    #outCpp = simulationCpp(input$N, input$p0, input$p1, input$burnin, input$nsim, method, measure)
        
    shinyjs::enable('runSimulation')
    
    list(round(out[[1]], 3),
         round(out[[2]], 3),
         round(out[[3]]),
         round(out[[4]], 4),
         round(out[[5]], 4),
         round(out[[6]], 2))
    
  }, ignoreNULL = TRUE) # ignore initial NULL value before any action
  
  # Render the results when 'results' reactive expression is triggered
  output$resultsText <- renderUI({
    # Convert the results to HTML format
    if(input$p0==input$p1){
      HTML(paste0("<div>",
                                "<h3>Simulation Results</h3>",
                                "<p><strong>Type-I error:</strong> ", results()[1], "%</p>",
                                "<p><strong>Monte Carlo Error Type-I error:</strong> ", results()[6], "%</p>",
                                #"<p><strong>Percentage of Patients assigned to treatment arm 1:</strong> ", results()[2], "%</p>",
                                #"<p><strong>Variance of percentage of Patients assigned to treatment arm 1:</strong> ", results()[3], "</p>",
                                "<p><strong>Expected Mean Response:</strong> ", results()[4], "</p>",
                                "<p><strong>Variance of Expected Mean Response:</strong> ", results()[5], "</p>",
                                "</div>"))
    } else{
        HTML(paste0("<div>",
                                "<h3>Simulation Results</h3>",
                                "<p><strong>Power:</strong> ", results()[1], "%</p>",
                                "<p><strong>Monte Carlo Error Power:</strong> ", results()[6], "%</p>",
                                "<p><strong>Percentage of Patients assigned to superior treatment arm:</strong> ", results()[2], "%</p>",
                                "<p><strong>Variance of percentage of Patients assigned to superior treatment arm:</strong> ", results()[3], "</p>",
                                "<p><strong>Expected Mean Response:</strong> ", results()[4], "</p>",
                                "<p><strong>Variance of Expected Mean Response:</strong> ", results()[5], "</p>",
                                "</div>"))
    }
  })
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

shinyApp(ui = ui, server = server)


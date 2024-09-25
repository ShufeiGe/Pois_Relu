#Input:
#@ data.train: training dataset.
#@ Rn.idx: number of break points in interval [0,1] (the length of the annealling power series)
#@ K: number of MCMC iterations within each SMC loop.
#@ N: number of hyperplanes (i.e. neurons in the first layer).
#@ seed: random seed.
#@ no.cores: number of CPUs used; the default value is 1.
#@a,b, mu_0, sigma_2_0: hyperparameters in the model.
#@save.immediate.result: output the immediate result or not.
Annealed_SMC <- function(data.train,Rn.idx=FALSE,K=FALSE,N=FALSE,seed=FALSE,no.cores=FALSE,a=FALSE,b=FALSE,mu_0=FALSE,sigma2_0=FALSE,save.immediate.result=FALSE){
  if(is.null(data.train)){
    stop("input data must be providided", call.=FALSE)
  }
  if(!Rn.idx){
    Rn.idx <- 10
  }
  
  if(!K){
    K <- 1000
  }
  
  if(!N){
    N <- 10
  }
  
  if(!seed){
    seed <- 123
  }
  
  if(!no.cores){
    no.cores <- 1
  }
  
  if(!a){
    a <- 100
  }
  
  
  
  if(!b){
    b <- 1/2
  }
  
  if(!mu_0){
    mu_0 <- 0
  }
  
  
  
  if(!sigma2_0){
    sigma2_0 <- 1
  }
  
  
  phir_seq <- seq(0,1, length.out =(Rn.idx+1))
  phi.add <- 1/Rn.idx  
  
  dm <- dim(data.train)
  n <- dm[1]
  np2 <- dm[2]
  np <- np2-1
  y <- data.train[,np2]
  
  set.seed(seed)
  #initialize particles and weights
  temp.all <- vector(mode="list",length = K)
  strt.tm <- Sys.time()
  
  
  if(!save.immediate.result){
    result.temp<-vector(mode="list",length = K)
    for(i in 1:K){
      result.temp[[i]] <-  intlz_function(data.train,y)
    }
    
    for (r in 1:Rn.idx){
      
      temp.r <- result.temp
      temp.r.update <-mclapply(temp.r,function(temp.k.r){
        temp.k.r <-rMH_kernel(data.train,temp.k.r,(r+1),y,phir_seq)
        return(temp.k.r)
      },mc.cores=no.cores)
      result.temp <- temp.r.update
      
      #if(r==1){
      #  rnt.tm <-Sys.time() - strt.tm 
      #  cat("Estimated running times is",rnt.tm/60*Rn.idx, "(minutes) \n")
      #}
      #resample
      if(r<Rn.idx){
        logl.r1 <- sapply(result.temp,"[[","logl.r1")
        log.w <- phi.add*logl.r1
        wgt.temp <- exp(log.w-max(log.w))/sum(exp(log.w-max(log.w)))
        id.x <- sample(1:K,size=K,replace = TRUE,prob=wgt.temp)
        for(k in 1:K){
          temp.all[[k]] <- result.temp[[id.x[k]]]
        }
        result.temp <- temp.all
      }
    }
    return(result.temp) 
  }else{
    result.all<-vector(mode="list",length = (Rn.idx+1))
    for(i in 1:K){
      result.all[[1]][[i]] <-  intlz_function(y)
    }
    for (r in 1:Rn.idx){
      
      #for(k in 1:K){
      #  temp.k.r <-rMH_kernel(result.all[[r]][[k]],(r+1))
      #   result.all[[(r+1)]][[k]] <- temp.k.r
      # }
      temp.r <- result.all[[r]]
      temp.r.update <-mclapply(temp.r,function(temp.k.r){
        temp.k.r <-rMH_kernel(data.train,temp.k.r,(r+1),y,phir_seq)
        return(temp.k.r)
      },mc.cores=no.cores)
      result.all[[r+1]] <- temp.r.update
      
      #if(r==1){
      #  rnt.tm <-Sys.time() - strt.tm 
      #  cat("Estimated running times is",rnt.tm/60*Rn.idx, "(minutes) \n")
      #}
      
      #resample
      if(r<Rn.idx){
        logl.r1 <- sapply(result.all[[r+1]],"[[","logl.r1")
        log.w <- phi.add*logl.r1
        wgt.temp <- exp(log.w-max(log.w))/sum(exp(log.w-max(log.w)))
        id.x <- sample(1:K,size=K,replace = TRUE,prob=wgt.temp)
        for(k in 1:K){
          temp.all[[k]] <- result.all[[r+1]][[id.x[k]]]
        }
        result.all[[r+1]] <- temp.all
      }
      
      #if(r%%10==0){
      #  save.image(output)
      #}
    }
    return(result.all) 
  }
  
  
  
}



#initialize particles.
intlz_function <- function(data.train,y){
  
  norms <- matrix(NA,nrow=np,ncol=N)
  for(i in 1:N){
    norm.temp <- rnorm(np)
    norms[,i] <- norm.temp/sqrt(sum(norm.temp^2))
  }
  
  ells <- runif(N)*sqrt(2)
  
  logl.r1 <- log(1)
  
  Zij <- matrix(NA,nrow=n,ncol=N2)
  Zij[,N2] <- 1 #column N+1 corresponds the constant term
  
  proj <- data.train[,1:np]%*%norms - matrix(ells,nrow=n,ncol=N,byrow = TRUE)
  Zij[,1:N] <-  (proj>0)*proj
  
  
  w <- matrix(rnorm(N2),N2,1)
  
  err.v <- y-Zij%*%w
  err.sq <- t(err.v)%*%err.v
  
  
  sigma2 <- invgamma::rinvgamma(1, shape=a, rate = b)
  
  
  input.ls <- list(Zij=Zij,norms=norms,ells=ells, sigma2=sigma2,w=w,logl.r1=logl.r1,accpet=NA)
  
  return(input.ls)
}

#the MCMC kernel at r step: pi(x)^{phi_r}
rMH_kernel <- function(data.train,input.ls,r,y,phir_seq){
  #set.seed(123)
  Zij <- input.ls$Zij
  norms <- input.ls$norms
  sigma2 <- input.ls$sigma2
  w <- input.ls$w
  ells <- input.ls$ells
  #weight for next step.  
  err.v <- y-Zij%*%w
  err.sq <-t(err.v)%*%err.v 
  logl.r1 <- -err.sq/2/sigma2-n*log(sigma2)/2
  
  
  #update a hyperplane by randomly propose a new one. 
  j <- sample(1:N,1)  # N: no. of hyper-planes 
  
  norm.temp <- rnorm(np)
  norm.temp <- norm.temp/sqrt(sum(norm.temp^2))
  
  #ell <- runif(1)*sqrt(2)
  ell <- runif(1)*sqrt(np)
  proj <- data.train[,1:np]%*%norm.temp - ell
  
  
  Zj <-  (proj>0)*proj
  
  Zj.chg <- Zij[,j]-Zj
  wj <- w[j]
  
  logr <- wj^2*t(Zj.chg)%*%(Zj.chg)+2*wj*t(Zj.chg)%*%y-2*wj%*%t(Zj.chg)%*%Zij%*%w
  logr <- -phir_seq[r]*logr/2/sigma2 #logr/2/sigma2
  mu_r <- log(runif(1))
  accpet <- (mu_r < logr)
  
  if(accpet){
    Zij[,j] <- Zj 
    norms[,j] <- norm.temp
    ells[j] <- ell
    err.v <- y-Zij%*%w
    err.sq <- t(err.v)%*%err.v
  }
  
  a.star <- phir_seq[r]*(a+n/2) #a.star <- a+n/2
  b.star <- phir_seq[r]*(b + err.sq/2) #b.star <-  b + t(err.v)%*%err.v/2
  
  sigma2 <- invgamma::rinvgamma(1, shape=a.star, rate = b.star)
  sigma12<- sigma2*sigma2_0
  
  for(i in 1:N2){
    Zij.i <- Zij[,i]
    t1j <- sigma2+sigma2_0*t(Zij.i)%*%Zij.i 
    t2j <- mu_0*sigma2+sigma2_0*t((y-Zij[,-i]%*%matrix(w[-i,],N,1)))%*%Zij.i
    muj <- t2j/t1j
    sigmaj <-sqrt(sigma12/t1j/phir_seq[r]) #sqrt(sigma12/t1j)
    w[i] <- rnorm(1,mean=muj, sd = sigmaj)
  }
  result <- list(Zij=Zij,norms=norms,ells=ells,sigma2=sigma2,w=w,logl.r1=logl.r1,accept =accpet)
  return(result)
}


#output summarized norms, ells and weights of the annealed SMC.
summary_smc <- function(result,phi.add) {
  
  logl.r1 <- sapply(result,"[[","logl.r1")
  log.w <- phi.add*logl.r1
  
  wgt<- exp(log.w-max(log.w))/sum(exp(log.w-max(log.w)))
  
  ws <- sapply(result,"[[","w")
  sigma2s <- sapply(result,"[[","sigma2")
  norms <- sapply(result,"[[","norms")
  ells <- sapply(result,"[[","ells")
  
  
  
  sigma2_mean_sd <- weighted_mean_sd(sigma2s,wgt) 
  ws_mean_sd <- t(as.matrix(apply(ws,1,weighted_mean_sd,wgt=wgt)))
  norms_mean_sd <- t(as.matrix(apply(norms,1,weighted_mean_sd,wgt=wgt)))
  ells_mean_sd <- t(as.matrix(apply(ells,1,weighted_mean_sd,wgt=wgt)))
  
  ws <- ws_mean_sd[,1]
  norms <- matrix(norms_mean_sd[,1],np,N,byrow = FALSE)
  ells <- ells_mean_sd[,1]
  return(list(norms=norms,ells=ells,ws=ws))
}


#calculate weighted mean and standard deviation
weighted_mean_sd <- function(vec,wgt){
  v1 <- sum(wgt)
  v2 <- sum(wgt^2)
  mean_wgt <- sum(vec*wgt)/v1
  sd_wgt <- sqrt(sum((vec-mean_wgt)^2*wgt)*v1/(v1^2-v2))
  return(c(mean_wgt,sd_wgt))
}

#compute the RMSE of prediction
rmse_smc_f <- function(data.in,norms,ells,ws){
  dm <- dim(data.in)
  n <- dm[1]
  np2 <- dm[2]
  np <- np2 -1
  y <- data.in[,np2]
  
  N2 <- length(ws)
  N <- N2 -1
  
  Zij <- matrix(NA,nrow=n,ncol=N2)
  Zij[,N2] <- 1 #column N+1 corresponds the constant term
  
  proj <- data.in[,1:np]%*%norms - matrix(ells,nrow=n,ncol=N,byrow = TRUE)
  Zij[,1:N] <-  (proj>0)*proj
  
  y.fit <- Zij%*%ws
  rmse <- sqrt(mean((y.fit -y)^2))
  return(rmse)
}





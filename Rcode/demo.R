require(invgamma)
require(parallel)
source("main_functions.R")


Rn.idx <- 100 #length of annealing power series
K <- 100      #number of particles
t.ratio <- 0.75 #training proportion in the train/test split 
N <- 2   #number of hyperplanes (lines/planes) in the simulated data
N2 <- N+1  
seed <- 123
no.cores <- 1 #number of cores used.

n_org <- 500  #number of obserzations
np <- 2       #number of predictors
np2 <- np+1   #number of predictors + 1 constant term

Npp <- 2#number of hyperplanes (lines/planes) in the model
Npp2 <- Npp+1

simulate_data <- function(seed,n,p,K){
  set.seed(seed)
  norms <- matrix(NA,nrow=p,ncol=K)
  ells <- runif(K,0,sqrt(p))
  Zij <- NULL
  for (i in 1:K){
    normal.temp <-rnorm(p)
    norms[,i] <- normal.temp/sqrt(sum(normal.temp^2))
  }
  
  data.x <- matrix(runif((n*p),-1,1),nrow=n,ncol=p)
  
  
  sigma2.true <- 0.1^2
  w <- matrix(rnorm(K),nrow=K,ncol=1)
  w0 <- rnorm(1)*0.1
  err <-sqrt(sigma2.true)*rnorm(n) #0.05*rnorm(n)
  
  proj <- data.x%*%norms - matrix(ells,nrow=n,ncol=K,byrow = TRUE)
  Zij<-  (proj>0)*proj
  y.true <- Zij%*%w+w0
  y  <-  y.true+err
  data <- cbind(data.x,y)
  return(data)
}



set.seed(seed)





#set hyperprameters & initialize particles
a <- 100
b <- 1/2

#a.star <- a+n/2
mu_0 <- 0
sigma2_0 <- 1

save.immediate.result <- FALSE


phi.add <- 1/Rn.idx  
(n <- round(n_org*t.ratio)) #number of observation in the training dataset



data <- simulate_data(seed,n_org,np,N) 
if(max(abs(range(data[,1:np])))>1){
  for(j in 1:np){
    data[,j] <- (data[,j] -min(data[,j]))/(max(data[,j])-min(data[,j]))*2-1
  }
}


id <- sample(1:n_org,n)
data.train <- data[id,]
data.test <- data[-id,]


n <- dim(data.train)[1]

#running time
running_time <-NULL


#result_all
strt.tm <- Sys.time()
result <-  Annealed_SMC(data.train,Rn.idx,K,Npp,seed,no.cores,a,b,mu_0,sigma2_0,save.immediate.result=FALSE)
end.tm <- Sys.time()
tm.tk <- end.tm - strt.tm

tm.tk

smc_o_avg <- summary_smc(result,phi.add) 

(rmse.train <- rmse_smc_f(data.train,smc_o_avg$norms,smc_o_avg$ells,smc_o_avg$ws))
(rmse.test <- rmse_smc_f(data.test,smc_o_avg$norms,smc_o_avg$ells,smc_o_avg$ws))

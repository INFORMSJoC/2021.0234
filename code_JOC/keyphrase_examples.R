library(mvtnorm)
library(quanteda)
library(coda)
library(MASS)
library(invgamma)
library(pROC)
library(ROCR)
#####################################################################################
###This code is to try semi-supervised methods on a known statistical paper##########
#####There is no ground truth for this paper, known keyphrases are observed from titles#####
#######################################################################################
big.graph.generate <- function (article){

  long_text <- tolower(readChar(paste("/Users/wangguanshen/Desktop/keyphrase/SemEval2010/",article,sep=""), nchars=10e6, useBytes = FALSE))
  long_text <- gsub("[/=]", " ", long_text)
  long_out <- fcm(long_text, context = "window", window = 2,tri=F)
  n <- nrow(long_out)
  d <- 0.85
  long_test <- as.matrix(cbind(long_out%*%rep(1,n),seq(n)))
  A <- as.matrix(long_out)
  diag(A) <- 0
  D <- diag(as.vector(long_out%*%rep(1,n)))
  G <- solve(D)%*%as.matrix(A)
  B <- diag(n)-d*t(G)
  u_0 <- solve(B)%*%rep(1-d,n)
  
  minus_list <- which(long_test[,1] < 12 | u_0 < sort(u_0,decreasing = TRUE)[200] )
  minus_list <- c(minus_list,87,145,152,153,171,172,181,232,328,333,334,369,411,425,438,523,563,589,591,824)#further remove some meanless units
  A_minus <- A[-minus_list,-minus_list]
  n_minus <- n-length(minus_list)
  D_minus <- diag(as.vector(A_minus%*%rep(1,n_minus)))
  dictionary_minus = cbind(long_test[-minus_list,2],seq(n_minus))
  print(dictionary_minus)

  return(list(n=n,A=A,D=D,dictionary=long_test,A_minus=A_minus,D_minus=D_minus,n_minus=n_minus,minus_list=minus_list,dictionary_minus=dictionary_minus))
}

graph<- big.graph.generate("EM_Preprocessed.txt")
inv_logit <- function(x){
  return(exp(x)/(1 + exp(x)))
}
#Posterior of theta, given sigma^2
posterior_gibbstheta <- function(Y,alpha,theta,u_0,B,sigma2){
  temp <- (1-alpha)*inv_logit(theta);
  C <- t(B%*% (theta-u_0) ) %*% (B%*% (theta-u_0) )  
  lglk <- sum(Y*log(temp)+(1-Y)*log(1-temp))-C/(2*sigma2)
  return(lglk)
}
posterior_gibbstheta2 <- function(Y,alpha,theta){
  temp <- (1-alpha)*inv_logit(theta);
  lglk <- sum(Y*log(temp)+(1-Y)*log(1-temp))
  return(lglk)
}
#Posterior of theta (intergral sigma^2 out)
log.posterior <- function(n,Y,alpha,theta,u_0,B){
  temp <- (1-alpha)*inv_logit(theta);
  C <- t(B%*%(theta-u_0))%*%(B%*%(theta-u_0))
  lglk <- sum(Y*log(temp)+(1-Y)*log(1-temp))-(n/2+0.01)*log(C/2+0.01)
  return(lglk)
}
#semi-supervised score to initial point
Base_to_start <- function(Base_Line){
  ini.point <- Base_Line
  ini.point[which(ini.point >= 1)] = 0.99
  ini.point <- log(ini.point/(1-ini.point))
  return(ini.point)
}
#select FDR cutoff
FDR_calculate <- function(cutoff,poster_md_adjust){
  set <- poster_md_adjust[poster_md_adjust >= cutoff]
  FDR <- sum(1-set)/length(set)
  return(FDR)
}
vec_FDR_cal <- Vectorize(FDR_calculate, "cutoff")
###find FDR cutoff, no ground truth needed
FDR_cutoff_notruth <- function(poster_pi_md,c,Y){
  poster_md_adjust <- force_obs_to_key2(Y,poster_pi_md)
  cutoffs <- unique(sort(poster_md_adjust,decreasing = TRUE))
  FDRs <- vec_FDR_cal(cutoffs,poster_md_adjust)
  index <- max(which(FDRs<c))
  cutoff <- cutoffs[index]
  return(cutoff)
}
#Gibbs and MH sampling function
#T: length of chain. ini,point: starting points.
#u_0: textrank score. B: graph matrix. Y: observed labels.
Gibbs.MH <- function(Burn_in,T,ini,n_minus,graph,Y_minus,B_minus,u_0_minus,alpha_est,grid){
  d <- 0.85
  theta_store <- matrix(0,nrow = T,ncol=n_minus)
  sigma2_store <- rep(NA,T)
  alpha_store <- rep(NA,T)
  accept <- 0 
  theta <- ini
  for (t in 1:T) {
    #sample from sigma
    C <- t(B_minus%*% (theta-u_0_minus) )%*%(B_minus%*% (theta-u_0_minus) )
    sigma2 <- rinvgamma(1,n_minus/2+0.001,rate=C/2+0.001)
    sigma2_store[t] <- sigma2
    theta_star <- mvrnorm(1, theta, solve(t(B_minus)%*%B_minus)*sigma2*5/n_minus)
    theta_star <- sapply(theta_star, function(y) min(max(y,-700),700))
    log_MH_rate <- posterior_gibbstheta(Y_minus,alpha_est,theta_star,u_0_minus,B_minus,sigma2)-posterior_gibbstheta(Y_minus,alpha_est,theta,u_0_minus,B_minus,sigma2)
    MH_rate <- exp(log_MH_rate)
    u <- runif(1)
    if (u < MH_rate) {
      theta <- theta_star
      accept <- accept + 1
    }
    theta_store[t,] <- theta
    alpha_est <- alpha_find(theta,Y_minus,grid)
    alpha_store[t] <- alpha_est
  }
  print(paste0("Accept rate: ", accept/T))
  prob <- inv_logit(theta_store)
  poster_pi_md <- apply(prob[Burn_in:T,],2,median)
  poster_pi_md[which(poster_pi_md==1)] <- 2 + rnorm(length(which(poster_pi_md==1)),0,0.01)
  poster_pi_mn <- apply(prob[Burn_in:T,],2,mean)
  poster_theta_mn <- apply(theta_store[Burn_in:T,],2,mean)
  poster_pi_mnV2 <- inv_logit(poster_theta_mn)
  alpha_mn <- mean(alpha_store)
  return(list(theta_store=theta_store,sigma2_store=sigma2_store,prob=prob,
    poster_pi_md=poster_pi_md,poster_pi_mn=poster_pi_mn,poster_pi_mnV2=poster_pi_mnV2,alpha_mn=alpha_mn,alpha_store=alpha_store,accept=accept))
}


#call function to run chains
semi.keyphrase2 <- function(graph,obs_label,grid){
  
  n_minus <- graph$n_minus
  d <- 0.85
  G_minus <- solve(graph$D_minus)%*%as.matrix(graph$A_minus)
  B_minus <- diag(n_minus)-d*t(G_minus)
  u_0_minus <- solve(B_minus)%*%rep(1-d,n_minus)
  w_minus <- sqrt(solve(graph$D_minus))
  B_star_minus <-  diag(n_minus)-d*w_minus%*%graph$A_minus%*%w_minus
  Y_minus <- rep(0,n_minus)
  Y_minus[obs_label]=1
  Base_Line_minus <- solve(B_star_minus)%*%Y_minus
  T <- 50000
  Burn_in <- 2000
  ini <- Base_to_start(u_0_minus)
  alpha_est <- alpha_find(u_0_minus,Y_minus,grid)
  test.chain <- Gibbs.MH(Burn_in,T,ini,n_minus,graph,Y_minus,B_minus,u_0_minus,alpha_est,grid)
  return(list(poster_pi_md=test.chain$poster_pi_md,poster_pi_mn=test.chain$poster_pi_mn,poster_pi_mnV2=test.chain$poster_pi_mnV2,
              obs_label=obs_label,dictionary_minus=graph$dictionary_minus,
              u_0_minus=c(u_0_minus),Base_Line_minus=c(Base_Line_minus),Y_minus=Y_minus,
              alpha_est = test.chain$alpha_mn, sigma2_store = test.chain$sigma2_store,
              accept_rate=test.chain$accept/T))
}

#enforece observed labels to be positive keyphrase
force_obs_to_key <- function (Y,poster_pi_mean,k){
  poster_pi_mean[which(Y==1)]=10 + rnorm(k,0.01)
  return(poster_pi_mean)
}
force_obs_to_key2 <- function (Y,poster_pi_mean){
  poster_pi_mean[which(Y==1)]=1 
  return(poster_pi_mean)
}
alpha_find <- function(Base_Line,Y,grid){
  alpha_est <- grid[which.max(vec.alpha_lk(Base_Line,Y,grid))] 
  return(alpha_est)
}
alpha_lk <- function(Base_Line,Y,alpha){
  pi <- inv_logit(c(Base_Line))
  return (sum(log((1-alpha)*pi)*Y + log(1-(1-alpha)*pi)*(1-Y)))
}
vec.alpha_lk <- Vectorize(alpha_lk,vectorize.args = 'alpha')
k <- 4
grid <- (seq(10,42)-5)/seq(10,42)
set.seed(12345)
Non_para_result <- semi.keyphrase2(graph,c(2,18,35,75),grid)
cutoff <- FDR_cutoff_notruth(Non_para_result$poster_pi_mn,0.3,Non_para_result$Y_minus) #cutoff is 0.598,
cutoff
#identified keywords by each method
graph$dictionary_minus[Non_para_result$poster_pi_mn>cutoff,]
num_of_positive <- sum(Non_para_result$poster_pi_mn>cutoff)
graph$dictionary_minus[order(force_obs_to_key(Non_para_result$Y_minus,Non_para_result$u_0_minus,k),decreasing = TRUE)[1:num_of_positive],]
graph$dictionary_minus[order(Non_para_result$Base_Line_minus,decreasing = TRUE)[1:num_of_positive],]



##############################################################
####following is example of Amazon review text from Kaggle####
##############################################################
amazon.graph.generate <- function (article){
  long_text <- tolower(readChar(paste("/Users/wangguanshen/Desktop/keyphrase/SemEval2010/",article,sep=""), nchars=10e6, useBytes = FALSE))
  long_text <- gsub("[/=]", " ", long_text)
  long_out <- fcm(long_text, context = "window", window = 2,tri=F)
  n <- nrow(long_out)
  d <- 0.85
  long_test <- as.matrix(cbind(long_out%*%rep(1,n),seq(n)))
  A <- as.matrix(long_out)
  diag(A) <- 0
  D <- diag(as.vector(long_out%*%rep(1,n)))
  G <- solve(D)%*%as.matrix(A)
  B <- diag(n)-d*t(G)
  u_0 <- solve(B)%*%rep(1-d,n)
  A_minus <- A
  D_minus <- D
  n_minus <- n

  return(list(n=n,A=A,D=D,dictionary=long_test,A_minus=A_minus,D_minus=D_minus,n_minus=n_minus
              #,dictionary_minus=dictionary_minus
              ))
}
amazon.graph<- amazon.graph.generate("amazon_pre.txt")

set.seed(12345)
#two observed keywords
Amazon_result <- semi.keyphrase2(amazon.graph,c(3,83),grid)
cutoff2 <- FDR_cutoff_notruth(Amazon_result$poster_pi_mn,0.25,Amazon_result$Y_minus) #cutoff 0.65
#identified keywords by each method
amazon.graph$dictionary[Amazon_result$poster_pi_mn>=cutoff2,]
num_of_positive2 <- sum(Amazon_result$poster_pi_mn>=cutoff2)
amazon.graph$dictionary[order(force_obs_to_key(Amazon_result$Y_minus,Amazon_result$u_0_minus,2),decreasing = TRUE)[1:num_of_positive2],]
amazon.graph$dictionary[order(Amazon_result$Base_Line_minus,decreasing = TRUE)[1:num_of_positive2],]

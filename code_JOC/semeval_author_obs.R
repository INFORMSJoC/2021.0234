library(RcppArmadillo)
library(microbenchmark)
library(mvtnorm)
library(quanteda)
library(MASS)
library(Rcpp)
library(coda)
library(invgamma)
library(inline)
library(pROC)
library(ROCR)
library(rbenchmark)

args <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(args)

file_list <- list.files("/users/guanshenw/scratch/keyphrase/2010_train/pre_process")
author_key_list <- list.files("/users/guanshenw/scratch/keyphrase/2010_train/pre_process_author_truth")
reader_key_list <- list.files("/users/guanshenw/scratch/keyphrase/2010_train/pre_process_reader_truth")

big.graph.generate <- function(i){
  long_text <- tolower(readChar(paste("/users/wangguanshen/Desktop/keyphrase/SemEval2010/pre_process/",file_list[i],sep=""), nchars=10e6, useBytes = FALSE))
  long_text <- gsub("[[:punct:]]", "", long_text)
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
  minus_list <- which(long_test[,1] < 12)
  n_minus1 <- n-length(minus_list) #stemming and freq based filter
  minus_list <- which(long_test[,1] < 12 | u_0 < sort(u_0,decreasing = TRUE)[150]) 
  A_minus <- A[-minus_list,-minus_list]
  n_minus <- n-length(minus_list)
  D_minus <- diag(as.vector(A_minus%*%rep(1,n_minus)))
  
  key_directory <- paste0("/users/guanshenw/scratch/keyphrase/2010_train/pre_process_truth/",gsub(".txt.final","",file_list[i]))
  long_key <- gsub("[\t\r;]","",tolower(readChar(key_directory, nchars=10e6, useBytes = FALSE)))
  long_key <- gsub(",", " ",long_key)
  longkey_list <- strsplit(long_key, " ")[[1]]
  long_truth <- unique(long_test[longkey_list[longkey_list %in% rownames(long_test)],2])
  dictionary_minus = cbind(long_test[-minus_list,2],seq(n_minus))
  truth_minus <- unique(dictionary_minus[longkey_list[longkey_list %in% rownames(dictionary_minus)],2])
  obs_key_directory <- paste0("/users/guanshenw/scratch/keyphrase/2010_train/pre_process_author_truth/",gsub(".txt.final","",author_key_list[i]))
  obs_key <- gsub("[\t\r;]","",tolower(readChar(obs_key_directory, nchars=10e6, useBytes = FALSE)))
  obs_key <- gsub(",", " ",obs_key)
  obskey_list <- strsplit(obs_key, " ")[[1]]
  obs <- unique(dictionary_minus[obskey_list[obskey_list %in% rownames(dictionary_minus)],2])  
  return(list(file=file_list[i],n=n,A=A,D=D,dictionary=long_test,truth=long_truth,A_minus=A_minus,
              D_minus=D_minus,n_minus=n_minus1,minus_list=minus_list,dictionary_minus=dictionary_minus,
              truth_minus=truth_minus,obs=obs))
}
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
#semi-supervised score to initial point
Base_to_start <- function(Base_Line){
  ini.point <- Base_Line
  ini.point[which(ini.point >= 1)] = 0.99
  ini.point <- log(ini.point/(1-ini.point))
  return(ini.point)
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
    #sample from theta
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
  
  return(list(
    theta_store=theta_store,sigma2_store=sigma2_store,prob=prob,
    poster_pi_md=poster_pi_md,poster_pi_mn=poster_pi_mn,poster_pi_mnV2=poster_pi_mnV2,alpha_mn=alpha_mn,accept=accept))
}


#call function to run chains
semi.keyphrase2 <- function(graph,grid){
  file <- graph$file
  n_minus <- graph$n_minus
  truth_minus <- graph$truth_minus
  d <- 0.85
  G_minus <- solve(graph$D_minus)%*%as.matrix(graph$A_minus)
  B_minus <- diag(n_minus)-d*t(G_minus)
  u_0_minus <- solve(B_minus)%*%rep(1-d,n_minus)
  w_minus <- sqrt(solve(graph$D_minus))
  B_star_minus <-  diag(n_minus)-d*w_minus%*%graph$A_minus%*%w_minus
  
  obs_label <- graph$obs
  k <- length(obs_label)
  True_Y_minus <- rep(0,n_minus)
  True_Y_minus[truth_minus]=1
  Y_minus <- rep(0,n_minus)
  Y_minus[obs_label] = 1 
  Base_Line_minus <- solve(B_star_minus)%*%Y_minus
  T <- 50000
  Burn_in <- 2000
  ini <-Base_Line_minus
  alpha_est <- alpha_find(u_0_minus,Y_minus,grid)
  test.chain <- Gibbs.MH(Burn_in,T,ini,n_minus,graph,Y_minus,B_minus,u_0_minus,alpha_est,grid)
  
  return(list(file = file,poster_pi_md=test.chain$poster_pi_md,poster_pi_mn=test.chain$poster_pi_mn,poster_pi_mnV2=test.chain$poster_pi_mnV2,
              truth_minus=truth_minus,obs_label=obs_label,dictionary_minus=graph$dictionary_minus,
              u_0_minus=c(u_0_minus),Base_Line_minus=c(Base_Line_minus),Y_minus=Y_minus,
              alpha_est = test.chain$alpha_mn,n=graph$n,
              accept_rate=test.chain$accept/T))
}


#Given the FDR, find out where the cutoff is. 
FDR_cutoff <- function(poster_pi_md,c,Y,truth){
  poster_md_adjust <- force_obs_to_key2(Y,poster_pi_md)
  cutoffs <- unique(sort(poster_md_adjust,decreasing = TRUE))
  FDRs <- vec_FDR_cal(cutoffs,poster_md_adjust)
  index <- max(which(FDRs<c))
  cutoff <- cutoffs[index]
  FDR_pos <- sum(poster_md_adjust>=cutoff) 
  FDR_tp <- sum(which(poster_md_adjust>=cutoff) %in% truth)
  Real_FDR <- (FDR_pos-FDR_tp)/FDR_pos
  return(list(FDR_pos=FDR_pos,FDR_tp=FDR_tp,Real_FDR=Real_FDR))
}
vec_FDR_cutoff <- Vectorize(FDR_cutoff,"c")

FDR_cutoff2 <- function(poster_pi_md,c,Y,truth){
  poster_md_adjust <- (poster_pi_md - min(poster_pi_md))/ (max(poster_pi_md) - min(poster_pi_md)) - 0.25
  poster_md_adjust <- force_obs_to_key2(Y,poster_md_adjust)
  cutoffs <- unique(sort(poster_md_adjust,decreasing = TRUE))
  FDRs <- vec_FDR_cal(cutoffs,poster_md_adjust)
  index <- max(which(FDRs<c))
  cutoff <- cutoffs[index]
  FDR_pos <- sum(poster_md_adjust>=cutoff) 
  FDR_tp <- sum(which(poster_md_adjust>=cutoff) %in% truth)
  Real_FDR <- (FDR_pos-FDR_tp)/FDR_pos
  return(list(FDR_pos=FDR_pos,FDR_tp=FDR_tp,Real_FDR=Real_FDR))
}
vec_FDR_cutoff2 <- Vectorize(FDR_cutoff2,"c")

FDR_calculate <- function(cutoff,poster_md_adjust){
  set <- poster_md_adjust[poster_md_adjust >= cutoff]
  FDR <- sum(1-set)/length(set)
  return(FDR)
}
vec_FDR_cal <- Vectorize(FDR_calculate, "cutoff")
#given predictions and truth labels, calculate precision and recall
precision.recall.auc <- function (Y,truth_Y,poster_pi_md,per,k){
  roc_md <- roc(truth_Y,poster_pi_md,direction="<")
  auc <- roc_md$auc
  pred_md <- prediction(poster_pi_md, truth_Y)
  pref_md <- performance(pred_md, "prec", "rec");
  recall_md <- pref_md@x.values[[1]][round(seq(0.1,1,0.1)*length(pref_md@x.values[[1]]))]
  precision_md <- pref_md@y.values[[1]][round(seq(0.1,1,0.1)*length(pref_md@x.values[[1]]))]
  poster_md_adjust <- force_obs_to_key(Y,poster_pi_md,k)
  roc_md2 <- roc(truth_Y,poster_md_adjust,direction="<")
  auc2 <- roc_md2$auc
  pred_md2 <- prediction(poster_md_adjust, truth_Y)
  pref_md2 <- performance(pred_md2, "prec", "rec");
  recall_md2<- pref_md2@x.values[[1]][round(seq(0.1,1,0.1)*length(pref_md2@x.values[[1]]))]
  precision_md2 <- pref_md2@y.values[[1]][round(seq(0.1,1,0.1)*length(pref_md2@x.values[[1]]))]  
  return(list(auc=auc,auc_l = auc2, recall=recall_md,recall_l = recall_md2,precision=precision_md,
              precision_l = precision_md2))
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

alpha_lk <- function(Base_Line,Y,alpha){
  pi <- inv_logit(c(Base_Line))
  return (sum(log((1-alpha)*pi)*Y + log(1-(1-alpha)*pi)*(1-Y)))
}
vec.alpha_lk <- Vectorize(alpha_lk,vectorize.args = 'alpha')

alpha_find <- function(Base_Line,Y,grid){
  alpha_est <- grid[which.max(vec.alpha_lk(Base_Line,Y,grid))] 
  return(alpha_est)
}
set.seed(12345)
alpha_min <- (10-5)/10
alpha_max <- (150-5)/150
grid <- (seq(10,150)-5)/seq(10,150)
start_time <- Sys.time()

group_article <- function(j){
  total_ans <- list()
  for (i in ((j-1)*10 + 1): ((j-1)*10 + 10)){
    if (i<=144){
      graph <- big.graph.generate(i)
      tryCatch({Ans1 <- semi.keyphrase2(graph,grid)
      total_ans <- append(total_ans,list(Ans1))},
      error = function(e){
        message(paste("An error occurred for artcile", i,file_list[i],":\n"), e)}) 
    }
  } 
  return(total_ans)
}
total_ans <- group_article(i)
end_time <- Sys.time()
end_time - start_time
save_path <- paste("Semeval_obs_author_",as.numeric(i),".Rdata",sep="")
save(total_ans,file=save_path)      

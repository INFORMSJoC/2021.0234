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
#Jun 23, 2022 update:
#Remove all unnecessary parts. Only mean and median are kept. Pauc and Correlation discarded.
#Oct 3 update:
#estimated alpha in both MCMCs, add time comparison.

#Sep 24 update:
#Main difference with HPC_Gibbs_MH_sym.R:
#1. Only 2 chains, 3 estimators, no MAP. estimated alpha is included.
#2. update FDR tables, Real FDR now is overall FDR not average
#3. update precision, recall calculation.
#This function is used to generate graph, given text and true labels.
#output of this function includes: dimension n; graph matrix A,D; words dictionary; truth index;

#Apr 27, 2021 update:
#For the purpose of comparison
#Save all articles name, all dictionaries, all truths index
sourceCpp("/users/guanshenw/scratch/keyphrase/Component_MCMC.cpp")
file_list <- list.files("/users/guanshenw/scratch/keyphrase/Test/pre_process")
save(file_list,file = "/users/guanshenw/scratch/keyphrase/file_list.Rdata")
graph.generate <- function (i){
  directory <- paste0("/users/guanshenw/scratch/keyphrase/Test/pre_process/",file_list[i],collapse = '')
  text <- tolower(readChar(directory, nchars=10e6, useBytes = FALSE))
  raw_text <- tolower(readChar(directory, nchars=10e6, useBytes = FALSE))
  raw_text<- gsub("[\t\r;]","",raw_text)
  raw_text<- gsub("\n"," ",raw_text)
  units <- length(strsplit(raw_text, " ")[[1]])
  key_directory <- paste0("/users/guanshenw/scratch/keyphrase/Test/",gsub("abstr","uncontr",file_list[i]))
  key <- gsub("[\t\r;]","",tolower(readChar(key_directory, nchars=10e6, useBytes = FALSE)))
  key <- gsub("\n", " ",key)
  key_list <- strsplit(key, " ")[[1]]
  out <- fcm(raw_text, context = "window", window = 2,tri=F)
  n <- nrow(out)
  test <- as.matrix(cbind(out%*%rep(1,n),seq(n)))
  A <- as.matrix(out)
  D <- diag(as.vector(out%*%rep(1,n)))
  degree <- sum((as.matrix(out)!=0)%*%rep(1,n))
  truth <- unique(test[key_list[key_list %in% rownames(test)],2])
  return(list("n"=n,"A"=A,"D"=D,"test"=test,"truth"=truth,"units"=units,"degree"=degree))
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
#Posterior of theta (intergral sigma^2 out)
log.posterior <- function(n,Y,alpha,theta,u_0,B){
  temp <- (1-alpha)*inv_logit(theta);
  C <- t(B%*%(theta-u_0))%*%(B%*%(theta-u_0))
  lglk <- sum(Y*log(temp)+(1-Y)*log(1-temp))-(n/2+0.001)*log(C/2+0.001)
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
Gibbs.MH <- function(Burn_in,T,ini,n,graph,Y,B,u_0,alpha_est,grid){
  d <- 0.85
  theta_store <- matrix(0,nrow = T,ncol=n)
  sigma2_store <- rep(NA,T)
  alpha_store <- rep(NA,T)
  accept <- 0 
  theta <- ini
  for (t in 1:T) {
    #sample from sigma
    C <- t(B%*% (theta-u_0) )%*%(B%*% (theta-u_0) )
    sigma2 <- rinvgamma(1,n/2+0.001,rate=C/2+0.001)
    sigma2_store[t] <- sigma2
    #sample from theta
    theta_star <- mvrnorm(1, theta, solve(t(B)%*%B)*sigma2*4/n)
    theta_star <- sapply(theta_star, function(y) min(max(y,-700),700))
    log_MH_rate <- posterior_gibbstheta(Y,alpha_est,theta_star,u_0,B,sigma2)-posterior_gibbstheta(Y,alpha_est,theta,u_0,B,sigma2)
    MH_rate <- exp(log_MH_rate)
    u <- runif(1)
    if (u < MH_rate) {
      theta <- theta_star
      accept <- accept + 1
    }
    theta_store[t,] <- theta
    alpha_est <- alpha_find(theta,Y,grid)
    alpha_store[t] <- alpha_est
  }
  print(paste0("Accept rate: ", accept/T))
  prob <- inv_logit(theta_store)
  poster_pi_md <- apply(prob[Burn_in:T,],2,median)
  poster_pi_md[which(poster_pi_md==1)] <- 1 + rnorm(length(which(poster_pi_md==1)),0,0.01)
  poster_pi_mn <- apply(prob[Burn_in:T,],2,mean)

  alpha_mn <- mean(alpha_store)
  alpha_md <- median(alpha_store)
  return(list(poster_pi_md=poster_pi_md,poster_pi_mn=poster_pi_mn,
              sigma2_store=sigma2_store,alpha_mn=alpha_mn,alpha_md=alpha_md))
}


#call function to run chains
semi.keyphrase2 <- function(graph,k,grid){
  alpha <- (length(graph$truth) - k) / length(graph$truth)
  n <- graph$n
  truth <- graph$truth
  d <- 0.85
  G <- solve(graph$D)%*%as.matrix(graph$A)
  B <- diag(n)-d*t(G)
  w <- sqrt(solve(graph$D))
  B_star <-  diag(n)-d*w%*%graph$A%*%w
  #observed label Y
  Y <- rep(0,n)
  Y[sample(truth,k,replace=FALSE)]=1
  obs_label <- which(Y==1)
  u_0 <- solve(B)%*%rep(1-d,n)
  Base_Line <- solve(B_star)%*%Y
  ini <- Base_to_start(Base_Line)
  alpha_est <- alpha_find(u_0,Y,grid)
  
  T <- 50000
  Burn_in <- 2000
  bench_test <- benchmark("Gibbs" = {
    chain1 <- Gibbs.MH(Burn_in,T,ini,n,graph,Y,B,u_0,alpha_est,grid)
    poster_pi_md <- chain1$poster_pi_md
    poster_pi_mn <- chain1$poster_pi_mn
    sigma_chain <- chain1$sigma2_store
  },
  "Componentwise" = {
    chain3 <- Componetwise_MCMC(T,ini,n,grid,alpha_est,u_0,B,Y)
    prob3 <- inv_logit(chain3$theta)
    poster_pi_md3 <- apply(prob3[Burn_in:T,],2,median)
    poster_pi_md3[which(poster_pi_md3==1)] <- 1 + rnorm(length(which(poster_pi_md3==1)),0,0.01)
    poster_pi_mn3 <- apply(prob3[Burn_in:T,],2,mean)
  },replications = 1,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self"))
  #print(bench_test)
  #a very easy look at convergence
  z.score <- geweke.diag(mcmc(sigma_chain[Burn_in:T]))
  if (abs(z.score$z) > 8) {
    print("chain may not converge well")
    return(list(z.score=z.score$z))
  }
  truth_Y <- rep(0,n)
  truth_Y[truth] = 1
  
  u_0_adjust <- force_obs_to_key(Y,u_0,k)
  Base_Line_adjust <- force_obs_to_key(Y,Base_Line,k)
  pre.rec.auc.md <- precision.recall.auc(Y,truth,truth_Y,poster_pi_md,k,c=c(0.05,0.1,0.15,0.2,0.25,0.3))
  pre.rec.auc.mn <- precision.recall.auc(Y,truth,truth_Y,poster_pi_mn,k,c=c(0.05,0.1,0.15,0.2,0.25,0.3))

  pre.rec.auc.md3 <- precision.recall.auc(Y,truth,truth_Y,poster_pi_md3,k,c=c(0.05,0.1,0.15,0.2,0.25,0.3))
  pre.rec.auc.mn3 <- precision.recall.auc(Y,truth,truth_Y,poster_pi_mn3,k,c=c(0.05,0.1,0.15,0.2,0.25,0.3))

  pre.rec.auc.u_0 <- precision.recall.auc(Y,truth,truth_Y,c(u_0),k,c=0.15)
  pre.rec.auc.bl <- precision.recall.auc(Y,truth,truth_Y,c(Base_Line),k,c=0.15)
  
  return(list(truth=truth,obs_label=obs_label,u_0_adjust=u_0_adjust,Base_Line_adjust=Base_Line_adjust,alpha=alpha,alpha_mn1=chain1$alpha_mn,
              alpha_md1=chain1$alpha_md,alpha_mn3=chain3$alpha_mn,alpha_md3=chain3$alpha_md,
              poster_pi_md=poster_pi_md,poster_pi_mn=poster_pi_mn,#poster_pi_mnV2=poster_pi_mnV2,
              poster_pi_md3=poster_pi_md3,poster_pi_mn3=poster_pi_mn3,#poster_pi_mn3V2=poster_pi_mn3V2,
              pos_md=pre.rec.auc.md,pos_mn=pre.rec.auc.mn,#pos_mnV2=pre.rec.auc.mnV2,
              pos_md3=pre.rec.auc.md3,pos_mn3=pre.rec.auc.mn3,#pos_mn3V2=pre.rec.auc.mn3V2,
              textrank=pre.rec.auc.u_0,semi=pre.rec.auc.bl,
              FDR_pos_md=pre.rec.auc.md$FDR_pos,FDR_tp_md=pre.rec.auc.md$FDR_tp,Real_FDR_md=pre.rec.auc.md$Real_FDR,
              FDR_pos_mn=pre.rec.auc.mn$FDR_pos,FDR_tp_mn=pre.rec.auc.mn$FDR_tp,Real_FDR_mn=pre.rec.auc.mn$Real_FDR,

              FDR_pos_md3=pre.rec.auc.md3$FDR_pos,FDR_tp_md3=pre.rec.auc.md3$FDR_tp,Real_FDR_md3=pre.rec.auc.md3$Real_FDR,
              FDR_pos_mn3=pre.rec.auc.mn3$FDR_pos,FDR_tp_mn3=pre.rec.auc.mn3$FDR_tp,Real_FDR_mn3=pre.rec.auc.mn3$Real_FDR,
              z.score=z.score$z,time_elapsed=bench_test$elapsed
  ))
} 

#given predictions and truth labels, calculate precision and recall
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
FDR_calculate <- function(cutoff,poster_md_adjust){
  set <- poster_md_adjust[poster_md_adjust >= cutoff]
  FDR <- sum(1-set)/length(set)
  return(FDR)
}
vec_FDR_cal <- Vectorize(FDR_calculate, "cutoff")
#calculate precision and recall, given an output of method and specified FDR level.
precision.recall.auc <- function (Y,truth,truth_Y,poster_pi_md,k,c=c(0.05,0.1,0.15,0.2,0.25,0.3)){
  roc_md <- roc(truth_Y,poster_pi_md,direction="<")
  auc <- roc_md$auc

  pred_md <- prediction(poster_pi_md, truth_Y)
  pref_md <- performance(pred_md, "prec", "rec")
  recall_md <- pref_md@x.values[[1]][round(seq(0.1,1,0.1)*length(pref_md@x.values[[1]]))]
  precision_md <- pref_md@y.values[[1]][round(seq(0.1,1,0.1)*length(pref_md@x.values[[1]]))]
  
  poster_md_adjust <- force_obs_to_key(Y,poster_pi_md,k)
  roc_md2 <- roc(truth_Y,poster_md_adjust,direction="<")
  #Pauc <- c(auc(roc_md2,partial.auc=c(0,0.1)),  auc(roc_md2,partial.auc=c(0,0.3)),auc(roc_md2,partial.auc=c(0,0.5)), auc(roc_md2,partial.auc=c(0,0.7)),auc(roc_md2,partial.auc=c(0,0.9)))
  auc2 <- roc_md2$auc
  pred_md2 <- prediction(poster_md_adjust, truth_Y)
  pref_md2 <- performance(pred_md2, "prec", "rec")
  
  recall_md2<- pref_md2@x.values[[1]][round(seq(0.1,1,0.1)*length(pref_md2@x.values[[1]]))]
  precision_md2 <- pref_md2@y.values[[1]][round(seq(0.1,1,0.1)*length(pref_md2@x.values[[1]]))] 
  FDR_result <- vec_FDR_cutoff(poster_pi_md,c,Y,truth) 
  FDR_pos <- unlist(FDR_result[1,])
  FDR_tp <- unlist(FDR_result[2,])
  Real_FDR <- unlist(FDR_result[3,])
  return(list(auc=auc,auc_l = auc2, 
              recall=recall_md,recall_l = recall_md2,precision=precision_md,
              precision_l = precision_md2,FDR_pos=FDR_pos,FDR_tp=FDR_tp,Real_FDR=Real_FDR))
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
#Given the number of positives found in FDR control, find how baselines method work.
Baseline_finding <- function(u_0_adjust,Base_Line_adjust,FDR_cut_md,truth,total_word,split_points,split_key) {
  Real_FDR_TR <- vector()
  Real_FDR_BL <- vector()
  check_mn_tr <- vector()
  check_mn_bl <- vector()
  for (i in 1:(length(split_points)-1)){
    part <- round(FDR_cut_md[,1]/total_word*(split_points[i+1]-split_points[i]))
    TP_md_TR0.05 <- sum(order(u_0_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[1]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_BL0.05 <- sum(order(Base_Line_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[1]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_TR0.1 <- sum(order(u_0_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[2]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_BL0.1 <- sum(order(Base_Line_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[2]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_TR0.15 <- sum(order(u_0_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[3]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_BL0.15 <- sum(order(Base_Line_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[3]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_TR0.2 <- sum(order(u_0_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[4]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_BL0.2 <- sum(order(Base_Line_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[4]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_TR0.25 <- sum(order(u_0_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[5]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_BL0.25 <- sum(order(Base_Line_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[5]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_TR0.3 <- sum(order(u_0_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[6]] %in% truth[(split_key[i]+1):split_key[i+1]])
    TP_md_BL0.3 <- sum(order(Base_Line_adjust[(split_points[i]+1):split_points[i+1]],decreasing=T)[1:part[6]] %in% truth[(split_key[i]+1):split_key[i+1]])
    check_mn_tr <- rbind(check_mn_tr,c(TP_md_TR0.05,TP_md_TR0.1,TP_md_TR0.15,TP_md_TR0.2,TP_md_TR0.25,TP_md_TR0.3,part))
    check_mn_bl <- rbind(check_mn_bl,c(TP_md_BL0.05,TP_md_BL0.1,TP_md_BL0.15,TP_md_BL0.2,TP_md_BL0.25,TP_md_BL0.3,part))
    FDR_cut_md[1,3] <- FDR_cut_md[1,3] + TP_md_TR0.05
    FDR_cut_md[1,4] <- FDR_cut_md[1,4] + TP_md_BL0.05
    FDR_cut_md[2,3] <- FDR_cut_md[2,3] + TP_md_TR0.1
    FDR_cut_md[2,4] <- FDR_cut_md[2,4] + TP_md_BL0.1
    FDR_cut_md[3,3] <- FDR_cut_md[3,3] + TP_md_TR0.15
    FDR_cut_md[3,4] <- FDR_cut_md[3,4] + TP_md_BL0.15
    FDR_cut_md[4,3] <- FDR_cut_md[4,3] + TP_md_TR0.2
    FDR_cut_md[4,4] <- FDR_cut_md[4,4] + TP_md_BL0.2
    FDR_cut_md[5,3] <- FDR_cut_md[5,3] + TP_md_TR0.25
    FDR_cut_md[5,4] <- FDR_cut_md[5,4] + TP_md_BL0.25
    FDR_cut_md[6,3] <- FDR_cut_md[6,3] + TP_md_TR0.3
    FDR_cut_md[6,4] <- FDR_cut_md[6,4] + TP_md_BL0.3
    Real_FDR_TR <- rbind(Real_FDR_TR, 1-c(TP_md_TR0.05,TP_md_TR0.1,TP_md_TR0.15,TP_md_TR0.2,TP_md_TR0.25,TP_md_TR0.3)/part)
    Real_FDR_BL <- rbind(Real_FDR_BL, 1-c(TP_md_BL0.05,TP_md_BL0.1,TP_md_BL0.15,TP_md_BL0.2,TP_md_BL0.25,TP_md_BL0.3)/part)
  }
  FDR_cut_md[,6] <- 1 - FDR_cut_md[,3]/FDR_cut_md[,1]
  FDR_cut_md[,7] <- 1 - FDR_cut_md[,4]/FDR_cut_md[,1]
  return(list(FDR_cut_md=FDR_cut_md,Real_FDR_TR=Real_FDR_TR,Real_FDR_BL=Real_FDR_BL,check_mn_tr=check_mn_tr,check_mn_bl=check_mn_bl))
}
#Find alpha at current iteration
alpha_find <- function(Base_Line,Y,grid){
  alpha_est <- grid[which.max(vec.alpha_lk(Base_Line,Y,grid))] 
  return(alpha_est)
}
alpha_lk <- function(Base_Line,Y,alpha){
  pi <- inv_logit(c(Base_Line))
  return (sum(log((1-alpha)*pi)*Y + log(1-(1-alpha)*pi)*(1-Y)))
}
vec.alpha_lk <- Vectorize(alpha_lk,vectorize.args = 'alpha')

main.function2 <- function(k) {
  u_0_adjust <- vector()
  total_word <- 0 
  num_of_key <- vector()
  split_points <- vector()
  split_points <- c(split_points,0)
  total_keys <- 0 
  split_key <- vector()
  obs_label <- list()
  files <- c()
  dict <- list()
  truth_list <- list()
  words <- vector()
  units <- vector()
  degree <- vector()
  alpha <- vector()
  time_elapsed <- vector()
  alpha_mn1 <- vector()
  alpha_md1 <- vector()
  alpha_mn3 <- vector()
  alpha_md3 <- vector()
  md_TP_Pos <- vector()
  mn_TP_Pos <- vector()
  split_key <- c(split_key,0)
  #Base_Line <- vector()
  Base_Line_adjust <- vector()
  poster_pi_md <- vector()
  poster_pi_mn <- vector()

  poster_pi_md3 <- vector()
  poster_pi_mn3 <- vector()
  TP_0.15mn_article <- vector()
  truth <- vector()
  FDR_cut_md <- matrix(0,6,7)
  colnames(FDR_cut_md) <- c("Positve adj","True Postive adj","Textrank TP adj","Semi TP adj","Real FDR","Textrank FDR","Semi FDR")
  rownames(FDR_cut_md) <- c("FDR0.05","FDR0.1","FDR0.15","FDR0.2","FDR0.25","FDR0.3")
  FDR_cut_mn <- FDR_cut_md
  FDR_cut_md3 <- FDR_cut_md
  FDR_cut_mn3 <- FDR_cut_md
  article <- 0
  alpha_min <- (10-5)/10
  alpha_max <- (42-5)/42
  grid <- (seq(10,42)-5)/seq(10,42)
  for (i in 201:500){
    graph <- graph.generate(i)
    if (length(graph$truth) > 10) {
      
      Ans <- semi.keyphrase2(graph,k,grid)
      if (abs(Ans$z.score) < 8){
        files <- c(files,file_list[i])
        dict <- append(dict,list(graph$test))
        truth_list <- append(truth_list,list(graph$truth))
        truth <- c(truth, Ans$truth)
        units <- c(units, graph$units)
        degree <- c(degree,graph$degree)
        alpha <- c(alpha, Ans$alpha)
        time_elapsed <- rbind(time_elapsed,Ans$time_elapsed)
        alpha_mn1 <- c(alpha_mn1, Ans$alpha_mn1)
        alpha_md1 <- c(alpha_md1, Ans$alpha_md1)
        alpha_mn3 <- c(alpha_mn3, Ans$alpha_mn3)
        alpha_md3 <- c(alpha_md3, Ans$alpha_md3)
        obs_label <- append(obs_label,list(Ans$obs_label))
        u_0_adjust <- c(u_0_adjust,Ans$u_0_adjust)
        Base_Line_adjust <- c(Base_Line_adjust,Ans$Base_Line_adjust)
        poster_pi_md <- c(poster_pi_md,Ans$poster_pi_md)
        poster_pi_mn <- c(poster_pi_mn,Ans$poster_pi_mn)

        poster_pi_md3 <- c(poster_pi_md3,Ans$poster_pi_md3)
        poster_pi_mn3 <- c(poster_pi_mn3,Ans$poster_pi_mn3)

        FDR_cut_md[,c(1:2,5)] <- FDR_cut_md[,c(1:2,5)] + cbind(Ans$FDR_pos_md,Ans$FDR_tp_md,Ans$Real_FDR_md)
        FDR_cut_mn[,c(1:2,5)] <- FDR_cut_mn[,c(1:2,5)] + cbind(Ans$FDR_pos_mn,Ans$FDR_tp_mn,Ans$Real_FDR_mn)

        FDR_cut_md3[,c(1:2,5)] <- FDR_cut_md3[,c(1:2,5)] + cbind(Ans$FDR_pos_md3,Ans$FDR_tp_md3,Ans$Real_FDR_md3)
        FDR_cut_mn3[,c(1:2,5)] <- FDR_cut_mn3[,c(1:2,5)] + cbind(Ans$FDR_pos_mn3,Ans$FDR_tp_mn3,Ans$Real_FDR_mn3)
        md_TP_Pos <- rbind(md_TP_Pos,c(Ans$FDR_tp_md,Ans$FDR_pos_md))
        mn_TP_Pos <- rbind(mn_TP_Pos,c(Ans$FDR_tp_mn,Ans$FDR_pos_mn))
        
        TP_0.15mn_article <- rbind(TP_0.15mn_article,c(Ans$FDR_pos_mn[3],Ans$FDR_tp_mn[3]))

        article <- article + 1
        words <- c(words,graph$n)
        num_of_key <-c(num_of_key ,length(graph$truth))
        total_word <- total_word + graph$n
        total_keys <- total_keys + length(graph$truth)
        split_points <- c(split_points,total_word)
        split_key <- c(split_key,total_keys)  
      }
    }
    else {print(paste0(file_list[i], " does not have enough keyphrase"))}
  }
  FDR_cut_md[,5] <- 1 - (FDR_cut_md[,2] / FDR_cut_md[,1])
  FDR_cut_mn[,5] <- 1 - (FDR_cut_mn[,2] / FDR_cut_mn[,1])

  FDR_cut_md3[,5] <- 1 - (FDR_cut_md3[,2] / FDR_cut_md3[,1])
  FDR_cut_mn3[,5] <- 1 - (FDR_cut_mn3[,2] / FDR_cut_mn3[,1])

  md_comparison <- Baseline_finding(u_0_adjust,Base_Line_adjust,FDR_cut_md,truth,total_word,split_points,split_key)
  mn_comparison <- Baseline_finding(u_0_adjust,Base_Line_adjust,FDR_cut_mn,truth,total_word,split_points,split_key)
 
  md3_comparison <- Baseline_finding(u_0_adjust,Base_Line_adjust,FDR_cut_md3,truth,total_word,split_points,split_key)
  mn3_comparison <- Baseline_finding(u_0_adjust,Base_Line_adjust,FDR_cut_mn3,truth,total_word,split_points,split_key)
  #FDR0.05 mean
  mn_0.05 <- apply(cbind(mn_TP_Pos,mn_comparison$check_mn_tr,mn_comparison$check_mn_bl[,1:6])[,c(1,7,13,25,19)],2,sum)
  precision_0.05_mn <- c(mn_0.05[1]/mn_0.05[2],mn_0.05[3]/mn_0.05[5],mn_0.05[4]/mn_0.05[5]) 
  recall_0.05_mn <- c(mn_0.05[1]/total_keys,mn_0.05[3]/total_keys,mn_0.05[4]/total_keys)
  F1_0.05_mn <- 2*(precision_0.05_mn*recall_0.05_mn)/(precision_0.05_mn+recall_0.05_mn)
  #FDR0.1 mean
  mn_0.1 <- apply(cbind(mn_TP_Pos,mn_comparison$check_mn_tr,mn_comparison$check_mn_bl[,1:6])[,c(2,8,14,26,20)],2,sum)
  precision_0.1_mn <- c(mn_0.1[1]/mn_0.1[2],mn_0.1[3]/mn_0.1[5],mn_0.1[4]/mn_0.1[5]) 
  recall_0.1_mn <- c(mn_0.1[1]/total_keys,mn_0.1[3]/total_keys,mn_0.1[4]/total_keys)
  F1_0.1_mn <- 2*(precision_0.1_mn*recall_0.1_mn)/(precision_0.1_mn+recall_0.1_mn)
  #FDR0.15 mean
  mn_0.15 <- apply(cbind(mn_TP_Pos,mn_comparison$check_mn_tr,mn_comparison$check_mn_bl[,1:6])[,c(3,9,15,27,21)],2,sum)
  precision_0.15_mn <- c(mn_0.15[1]/mn_0.15[2],mn_0.15[3]/mn_0.15[5],mn_0.15[4]/mn_0.15[5]) 
  recall_0.15_mn <- c(mn_0.15[1]/total_keys,mn_0.15[3]/total_keys,mn_0.15[4]/total_keys)
  F1_0.15_mn <- 2*(precision_0.15_mn*recall_0.15_mn)/(precision_0.15_mn+recall_0.15_mn)
  #FDR0.2 mean
  mn_0.2 <- apply(cbind(mn_TP_Pos,mn_comparison$check_mn_tr,mn_comparison$check_mn_bl[,1:6])[,c(4,10,16,28,22)],2,sum)
  precision_0.2_mn <- c(mn_0.2[1]/mn_0.2[2],mn_0.2[3]/mn_0.2[5],mn_0.2[4]/mn_0.2[5]) 
  recall_0.2_mn <- c(mn_0.2[1]/total_keys,mn_0.2[3]/total_keys,mn_0.2[4]/total_keys)
  F1_0.2_mn <- 2*(precision_0.2_mn*recall_0.2_mn)/(precision_0.2_mn+recall_0.2_mn)

  precision_mn <- rbind(precision_0.05_mn,precision_0.1_mn,precision_0.15_mn,precision_0.2_mn)
  recall_mn <- rbind(recall_0.05_mn,recall_0.1_mn,recall_0.15_mn,recall_0.2_mn)
  F1_mn <- rbind(F1_0.05_mn,F1_0.1_mn,F1_0.15_mn,F1_0.2_mn)
  check_mn_tr <- mn_comparison$check_mn_tr
  check_mn_bl <- mn_comparison$check_mn_bl
  
  FDR_cut_md <- md_comparison$FDR_cut_md
  FDR_cut_mn <- mn_comparison$FDR_cut_md
  
  FDR_cut_md3 <- md3_comparison$FDR_cut_md
  FDR_cut_mn3 <- mn3_comparison$FDR_cut_md
 
  output <- list(u_0_adjust=u_0_adjust,Base_Line_adjust=Base_Line_adjust,
                 poster_pi_md=poster_pi_md,poster_pi_mn=poster_pi_mn,#poster_pi_mnV2=poster_pi_mnV2,
                 poster_pi_md3=poster_pi_md3,poster_pi_mn3=poster_pi_mn3,#poster_pi_mn3V2=poster_pi_mn3V2,
                 split_points=split_points,truth=truth,split_key=split_key,units=units,degree=degree,obs_label=obs_label,article=article,total_keys=total_keys,total_word=total_word,
                 files=files,dict=dict,truth_list=truth_list,
                 alpha=alpha,alpha_mn1=alpha_mn1,alpha_mn3=alpha_mn3,alpha_md1=alpha_md1,alpha_md3=alpha_md3,
                 num_of_key=num_of_key,words=words,time_elapsed=time_elapsed,
                 check_mn_tr=check_mn_tr,check_mn_bl=check_mn_bl,
                 precision_mn=precision_mn,recall_mn=recall_mn,F1_mn=F1_mn,
                 FDR_cut_md=FDR_cut_md,FDR_cut_mn=FDR_cut_mn,
                 FDR_cut_md3=FDR_cut_md3,FDR_cut_mn3=FDR_cut_mn3,
                 TP_0.15mn_article=TP_0.15mn_article,mn_comparison=mn_comparison,mn3_comparison=mn3_comparison
                 )
  return(output)
}


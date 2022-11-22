library(mvtnorm)
library(quanteda)
library(coda)
library(MASS)
library(invgamma)
library(pROC)
library(ROCR)
big.graph.generate <- function (article){
  long_text <- tolower(readChar(paste("/users/wangguanshen/Desktop/keyphrase/SemEval2010/train_preprocess/",article,sep=""), nchars=10e6, useBytes = FALSE))
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
  minus_list <- which(long_test[,1] < 12 | u_0 < sort(u_0,decreasing = TRUE)[150]) 
  A_minus <- A[-minus_list,-minus_list]
  n_minus <- n-length(minus_list)
  D_minus <- diag(as.vector(A_minus%*%rep(1,n_minus)))
  longkeyphrase <- tolower(c('ensembl kalman filter,data assimil methodolog,hydrocarbon reservoir simul,energi explor,tigr grid comput environ,grid comput,cyberinfrastructur develop project,high perform comput,tigr grid middlewar,strateg applic area,gridway metaschedul,pool licens,grid-en,reservoir model,enkf,tigr'))
  long_key <- gsub("[\t\r]","",longkeyphrase)
  long_key <- gsub(",", " ",long_key)
  longkey_list <- strsplit(long_key, " ")[[1]]
  long_truth <- unique(long_test[longkey_list[longkey_list %in% rownames(long_test)],2])
  dictionary_minus = cbind(long_test[-minus_list,2],seq(n_minus))
  print(dictionary_minus)
  truth_minus <- unique(dictionary_minus[longkey_list[longkey_list %in% rownames(dictionary_minus)],2])
  print(truth_minus)
  return(list(n=n,A=A,D=D,dictionary=long_test,truth=long_truth,A_minus=A_minus,D_minus=D_minus,n_minus=n_minus,minus_list=minus_list,dictionary_minus=dictionary_minus,truth_minus=truth_minus))
}

big.graphH41<- big.graph.generate("C-42(2).txt")
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
  MH_rate_store <- rep(NA,T)
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
    MH_rate_store[t] <- MH_rate
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
    theta_store=theta_store,sigma2_store=sigma2_store,prob=prob,MH_rate_store=MH_rate_store,
    poster_pi_md=poster_pi_md,poster_pi_mn=poster_pi_mn,poster_pi_mnV2=poster_pi_mnV2,alpha_mn=alpha_mn,accept=accept))
}


#call function to run chains
semi.keyphrase2 <- function(graph,obs_label,grid){
  k <- length(obs_label)
  n_minus <- graph$n_minus
  truth_minus <- graph$truth_minus
  d <- 0.85
  G_minus <- solve(graph$D_minus)%*%as.matrix(graph$A_minus)
  B_minus <- diag(n_minus)-d*t(G_minus)
  u_0_minus <- solve(B_minus)%*%rep(1-d,n_minus)
  w_minus <- sqrt(solve(graph$D_minus))
  B_star_minus <-  diag(n_minus)-d*w_minus%*%graph$A_minus%*%w_minus
  True_Y_minus <- rep(0,n_minus)
  True_Y_minus[truth_minus]=1
  Y_minus <- rep(0,n_minus)
  Y_minus[obs_label]=1
  Base_Line_minus <- solve(B_star_minus)%*%Y_minus
  T <- 50000
  Burn_in <- 2000
  ini <- Base_to_start(u_0_minus)
  alpha_est <- alpha_find(u_0_minus,Y_minus,grid)

  test.chain <- Gibbs.MH(Burn_in,T,ini,n_minus,graph,Y_minus,B_minus,u_0_minus,alpha_est,grid)

  per <- seq(0.05,1,0.05)
  pre.rec.auc.md <- precision.recall.auc(Y_minus,True_Y_minus,test.chain$poster_pi_md,per,k)
  pre.rec.auc.mn <- precision.recall.auc(Y_minus,True_Y_minus,test.chain$poster_pi_mn,per,k)
  pre.rec.auc.mnV2 <- precision.recall.auc(Y_minus,True_Y_minus,test.chain$poster_pi_mnV2,per,k)
  pre.rec.auc.u_0 <- precision.recall.auc(Y_minus,True_Y_minus,c(u_0_minus),per,k)
  pre.rec.auc.bl <- precision.recall.auc(Y_minus,True_Y_minus,c(Base_Line_minus),per,k)

  return(list(poster_pi_md=test.chain$poster_pi_md,poster_pi_mn=test.chain$poster_pi_mn,poster_pi_mnV2=test.chain$poster_pi_mnV2,
              truth_minus=truth_minus,obs_label=obs_label,dictionary_minus=graph$dictionary_minus,
              u_0_minus=c(u_0_minus),Base_Line_minus=c(Base_Line_minus),Y_minus=Y_minus,
              alpha_est = test.chain$alpha_mn,sigma2_store=test.chain$sigma2_store,theta_store=test.chain$theta_store,MH_rate_store=test.chain$MH_rate_store,
              pos_md=pre.rec.auc.md,pos_mn=pre.rec.auc.mn,pos_mnV2=pre.rec.auc.mnV2,
              textrank=pre.rec.auc.u_0,semi=pre.rec.auc.bl,
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
alpha_find <- function(Base_Line,Y,grid){
  alpha_est <- grid[which.max(vec.alpha_lk(Base_Line,Y,grid))] 
  return(alpha_est)
}
alpha_lk <- function(Base_Line,Y,alpha){
  pi <- inv_logit(c(Base_Line))
  return (sum(log((1-alpha)*pi)*Y + log(1-(1-alpha)*pi)*(1-Y)))
}
vec.alpha_lk <- Vectorize(alpha_lk,vectorize.args = 'alpha')

grid <- (seq(10,42)-5)/seq(10,42)

set.seed(541)
obs_label1 <- c(8,13,10,11,69,29)
k1 <- length(obs_label1)
Ans1 <- semi.keyphrase2(big.graphH41,obs_label1,grid)
save(Ans1,file = "/users/guanshenw/scratch/keyphrase/biggraphC42(sub_title).Rdata")
positive1 <- vec_FDR_cutoff(Ans1$poster_pi_mn,0.25,Ans1$Y_minus,Ans1$truth_minus)[[1]] # of positive identified by BSS method
Ans1$dictionary_minus[Ans1$obs_label,]#check observed keywords

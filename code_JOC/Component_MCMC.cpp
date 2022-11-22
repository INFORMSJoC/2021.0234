#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::export]]
arma::vec inv_logit1(arma::vec x){
  return exp(x)/(1 + exp(x));
}
//[[Rcpp::export]]
arma::vec log_posterior(int n, arma::vec Y, double alpha, arma::rowvec theta,arma::vec u_0, arma::mat B)
{
  double n1 = n;
  arma::vec temp = inv_logit1(theta.t());
  arma::mat C2 = B * (theta.t()-u_0);
  arma::mat C = C2.t() * C2;
  return(sum(Y.t()*log((1-alpha)*temp)+(1-Y.t())*log(1-(1-alpha)*temp))-(n1/2+0.001)*log(C/2+0.001));
  //return(Y.t()*log((1-alpha)*temp)+(1-Y.t())*log(1-(1-alpha)*temp));
  //return(-(n1/2+0.01)*log(C/2+0.01));
}

//[[Rcpp::export]]
double alpha_lk_cpp(arma::rowvec theta, arma::vec Y, double alpha )
{
  arma::vec pi = inv_logit1(theta.t());
  return(sum(Y.t()*log((1-alpha)*pi)+(1-Y.t())*log(1-(1-alpha)*pi)));
}
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
//mvrnormArma(100, mu, sigma)
//[[Rcpp::export]]
Rcpp::List Componetwise_MCMC(int T, arma::vec ini,int n, NumericVector grid, double alpha_est, arma::vec u_0,arma::mat B,arma::vec Y)
{
  int i,t,j;
  //arma::vec theta = ini;
  arma::rowvec theta_star,thetai_star,theta_nomove;
  NumericVector lg_MH_rate;
  NumericVector MH_rate,V = rep(0.5,n);
  double var1;
  NumericVector u;
  arma::vec sp_var;
  arma::mat theta_store(T, n);
  NumericMatrix var_store(T, n);
  NumericVector alpha_store(T);
  int accept = 0;
  int g = grid.size();
  arma::vec alpha_grid(g);
  theta_store.row(0) = ini.t();
  var_store(0,_) = V;
  alpha_store[0] = alpha_est;
  for (t = 1; t < T; t++) 
  {
    for (i = 0; i < n; i++){
      if (t < 10){
        var1 = 1;
      }
      else {
        sp_var = theta_store.submat(0,i,t-1,i);
        //sp_var.print();
        //var1 is the sd.
        var1 = sqrt(2.4*(var(sp_var) + 0.01));
      }
      //Rprintf("%f\n",var1);
      var_store(t,i) = var1;
      //prop_mean = theta_store(i,t-1);
      thetai_star = as<arma::vec>(wrap(rnorm(1,theta_store(t-1,i),var1)));
      if (i==0){
        //theta_star = join_cols(thetai_star,theta_store.submat(t-1,1,t-1,n-1));
        theta_nomove = theta_store.row(t-1);
      }
      //else if (i==(n-1)){ 
      //theta_star = join_cols(theta_store.submat(0,t,n-2,t ), thetai_star);
      //theta_nomove = join_cols(theta_store.submat(0,t,n-2,t),theta_store.submat(n-1,t-1,n-1,t-1));
      //  }
      else {
        //theta_star = join_cols(join_cols(theta_store.submat(0,t,i-1,t ), thetai_star),theta_store.submat(i+1,t-1,n-1,t-1 ));
        theta_nomove = join_rows(theta_store.submat(t,0,t,i-1),theta_store.submat(t-1,i,t-1,n-1));
      }
      theta_star = theta_nomove;
      theta_star[i] = as_scalar(thetai_star);
      //Rprintf("the iteration is %i,%i\n", t,i+1);
      lg_MH_rate = as<NumericVector>(wrap(log_posterior(n,Y,alpha_est,theta_star,u_0,B)-log_posterior(n,Y,alpha_est,theta_nomove,u_0,B)));
      MH_rate = exp(lg_MH_rate);
      u = runif(1);
      if (MH_rate[0]>u[0]) {
        //theta_store.submat(i,t,i,t) = as_scalar(thetai_star);
        theta_store(t,i) = as_scalar(thetai_star);
        accept += 1;
      }
      else {
        //theta_store.submat(i,t,i,t) = theta_store.submat(i,t-1,i,t-1);
        theta_store(t,i) = theta_store(t-1,i);
      }
    }
    for (j = 0; j < g; j++){
      alpha_grid[j] = alpha_lk_cpp(theta_store.row(t),Y,grid[j]);
    }
    alpha_est = grid[alpha_grid.index_max()];
    alpha_store[t] = alpha_est;
  }
  return Rcpp::List::create(Rcpp::Named("theta") = theta_store, Rcpp::Named("accept") = accept, Rcpp::Named("var") = var_store, Rcpp::Named("alpha_store") = alpha_store);
}
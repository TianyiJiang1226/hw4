#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double importance(int n) {
  NumericVector x = rnorm(n);
  double res = mean(pow(x,2) * exp(-(pow(abs(x),3))/3) / dnorm(x))/2.575798633708138;
  return res;
}



// [[Rcpp::export]]
double rejection(int n) {
  NumericVector x = rnorm(n);
  NumericVector u = runif(n);
  NumericVector dec(n,0.0);
  for (int i = 0; i < n; i++) {
    if(u[i] < ((exp(-(pow(abs(x[i]),3))/3)/dnorm(x)[i])/2.575798633708138)){
      dec[i] = 1.0;
    }
  }
  double res = sum(pow(x,2) * dec)/sum(dec);
  return res;
}


// [[Rcpp::export]]
double SIR(int n) {
  NumericVector x = rnorm(n);
  NumericVector f = exp(-(pow(abs(x),3))/3) / dnorm(x);
  double s = sum(f);
  NumericVector standf = f/s;
  NumericVector y = sample(x, n, true,standf);
  double res = mean(pow(y,2));
  return res;
}

// [[Rcpp::export]]
double Philippe(int n) {
  NumericVector x = rnorm(n);
  NumericVector u = runif(n);
  NumericVector y(n,0.0);
  for (int i = 0; i < n; i++) {
    if(u[i] < ((exp(-(pow(abs(x[i]),3))/3)/dnorm(x)[i])/2.575798633708138)){
      y[i] = x[i];
    }
  }
  y = y[y != 0.0];
  y = y.sort();
  int n2 = y.size() - 1;
  double res = 0;
  for (int i = 0; i < n2; i++) {
    res += (y[i+1] - y[i]) * pow(y[i],2) * exp(-(pow(abs(y[i]),3))/3) /2.575798633708138;
  }
  return res;
}


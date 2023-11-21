#include <Rcpp.h>
#include <cstdio>
#include <math.h>

using namespace Rcpp;

// Helper function to calculate normal density
  
double NormalDensity(double x, double mean, double std_dev) {
    // Sample from normal distribution
    double nx = (x - mean)/std_dev;
    double log_sig = log(std_dev * std_dev);
    double sqr_x = nx *nx;
    double densi = -log_sig - .5 * sqr_x ;
    long double result = exp(densi);
    // printf("Value of densi: %Lf\n", result); 
    return result;
    
}

// [[Rcpp::export]]
double lTargetKernelCpp(int d, double newX, double kdeValue, NumericVector tar, double sdConvertor) {
    double s_2 = sdConvertor * kdeValue;
    long double sumDnorm = 0;
    for(int i = 0; i < tar.length(); i++) {
      double mu = tar[i];
      // printf("Value of sum: %f\n", mu); 
      sumDnorm += NormalDensity(newX, mu, s_2);
    }

    // printf("____________________\nValue of sum: %Lf\n", sumDnorm); 
    double result = log(sumDnorm) -  log( sdConvertor * kdeValue);
    return result;
}

NumericVector getColumn(NumericMatrix tar,int col) {
  int nrow = tar.nrow();
  NumericVector result(nrow);
  
  for(int i = 0; i < nrow; i++){
    result[i] = tar(i, col); // 1 = second column 
  }
  
  return result;
}

// [[Rcpp::export]]
double loglikelihood_uqCpp(NumericVector params, NumericVector proxy_val, NumericVector kde_bw, NumericMatrix tar, double sdConvertor, IntegerVector loca ){
  long double sumll = 0;
  int t_i = tar.rows();
  // printf("Value of length of t: %i\n", t_i);
  for(int i = 0; i < proxy_val.length(); i++) {
    int location = loca[i];
    double newx = proxy_val[i];
    double bw_value = kde_bw[location];
    NumericVector tar_vec = getColumn(tar, location);
    // printf("Value of sum: %i\n", tar_vec.length()); 
    sumll += lTargetKernelCpp(i, newx, bw_value, tar_vec, sdConvertor);
  }
  
  return sumll;
} 








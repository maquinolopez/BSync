#include <Rcpp.h>
#include <cstdio>
#include <math.h>
#include <stdexcept>

using namespace Rcpp;



// [[Rcpp::export]]
NumericVector linearInterpolation(NumericVector X, NumericVector y, NumericVector newx) {
  int xy_size = X.size();
  int newx_size = newx.size();
  NumericVector results(newx_size);
  
  for (int i = 0; i < newx_size; ++i) {
    if (newx[i] <= X[0]) {
      results[i] = y[0];
    } else if (newx[i] >= X[xy_size - 1]) {
      results[i] = y[xy_size - 1];
    } else {
      for (int j = 0; j < xy_size - 1; ++j) {
        if (X[j] <= newx[i] && newx[i] <= X[j + 1]) {
          double slope = (y[j + 1] - y[j]) / (X[j + 1] - X[j]);
          results[i] = y[j] + slope * (newx[i] - X[j]);
          break;
        }
      }
    }
  }
  
  return results;
}


// This is the tau function
// [[Rcpp::export]]
NumericVector tauC(NumericVector x, NumericVector param, double b_length, NumericVector breaks) {
  double tau0 = param[0];
  int n_slopes = param.size() - 2; // Excluding tau0 and w
  NumericVector slopes = param[Rcpp::Range(2, param.size() - 1)]; // Equivalent to param[-c(1,2)] in R
  
  NumericVector tau_is(n_slopes + 1);
  tau_is[0] = tau0;
  double cumsum = tau0;
  
  for(int i = 0; i < n_slopes; ++i) {
    cumsum += slopes[i] * b_length;
    tau_is[i + 1] = cumsum;
  }
  NumericVector taus = linearInterpolation(breaks, tau_is, x);
  return taus;
}


//___________________________________
// This functions are for the normal cases (no UQ)

// [[Rcpp::export]]
double tdistroC(NumericVector  X, NumericVector Mu, double sigma, double a, double b) {
  // if (X.size() != Mu.size()) {
  //   throw std::invalid_argument("X and Mu must be of the same size");
  // }
  
  double sigmaSquared = sigma * sigma; // Precompute sigma squared
  double logSigma = 0.5 * log(sigmaSquared); // Precompute .5 * log(sigma)
  double aTerm = (2 * a + 1) / 2.0; // Precompute (2 * a + 1) / 2.0
  
  double sum = 0.0;
  
  for (size_t i = 0; i < X.size(); ++i) {
    double x = X[i];
    double mu = Mu[i];
    
    // if (std::isnan(x) || std::isnan(mu)) continue; // Skip NaN values for both x and mu
    double diffSquared = (x - mu) * (x - mu); // Precompute (x - mu)^2
    double term = aTerm * log(b + diffSquared / (2 * sigmaSquared)) + logSigma;
    sum -= term; // Subtracting because the original R function multiplies by -1
  }
  
  return sum;
} 

// [[Rcpp::export]]
double loglikelihoodC(NumericVector params, NumericVector tar_prox, 
                      NumericVector inp_X, NumericVector inp_ProxyValue, 
                      NumericVector tar_X, double my_sd, double ta, double tb, 
                      double b_length, NumericVector breaks) {
  
  // Calculate t_times using tauC function
  NumericVector t_times = tauC(inp_X, params, b_length, breaks);
  
  // Perform linear interpolation to get new_target values
  NumericVector new_target = linearInterpolation(tar_X, tar_prox, t_times);
  
  // Calculate the likelihood using tdistroC function
  double likelihood = tdistroC(inp_ProxyValue, new_target, my_sd, ta, tb);
  
  return likelihood;
}



//___________________________________
// This functions are for the UQ


// Helper function to calculate normal density

double NormalDensity(double x, double mean, double std_dev) {
  // Sample from normal distribution
  double nx = (x - mean)/std_dev;
  double log_sig = log(std_dev);
  double sqr_x = nx *nx;
  double densi = -log_sig - .5 * sqr_x ;
  long double result = exp(densi);
  // printf("Value of densi: %Lf\n", result); 
  return result;
  
}

// [[Rcpp::export]]
IntegerVector localizador(NumericVector tar_ages, NumericVector x) {
  IntegerVector loc(x.size());
  for (int i = 0; i < x.size(); ++i) {
    // Find the nearest lower index in tar_ages for each element in x
    auto it = std::lower_bound(tar_ages.begin(), tar_ages.end(), x[i]);
    loc[i] = std::distance(tar_ages.begin(), it);
  }
  return loc;
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
double lTargetKernelCpp(int d, double newX, double kdeValue, NumericVector tar, double sdConvertor) {
  double s_2 = sdConvertor * kdeValue;
  long double sumDnorm = 0;
  for(int i = 0; i < tar.length(); i++) {
    double mu = tar[i];
    sumDnorm += NormalDensity(newX, mu, s_2);
  }
  
  double result = log(sumDnorm) ;//-  log( sdConvertor * kdeValue);
  return result;
}


// Function to print an entire IntegerVector
void printVector(const NumericVector& vec) {
  for(int i = 0; i < vec.size(); ++i) {
    Rcpp::Rcout << vec[i] << " ";
  }
  Rcpp::Rcout << std::endl;
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
    double newvalue = lTargetKernelCpp(i, newx, bw_value, tar_vec, sdConvertor);
    // used for debugging    
    // if (std::isinf(newvalue)) {
    //   newvalue = 1;
    // }
    sumll += newvalue;
  }
  

  
  return sumll;
} 




// [[Rcpp::export]]
double loglikelihood_uqCfull(NumericVector params, NumericVector inp_X, 
                                NumericVector inp_ProxyValue, NumericVector bandwidths, 
                                NumericMatrix tar_mt, double sd_convertor, 
                                double b_length, NumericVector breaks, 
                                NumericVector tar_ages) {
  // Calculate t_times using tauC function
  NumericVector t_times = tauC(inp_X, params, b_length, breaks);
  
  // Get the values of the target at the required times
  IntegerVector loca = localizador(tar_ages, t_times);
  // printIntegerVector(loca);
  // Calculate the log likelihood using loglikelihood_uqCpp function
  double ll = loglikelihood_uqCpp(params, inp_ProxyValue, bandwidths, tar_mt, sd_convertor, loca);
  if (std::isinf(ll)) {
    Rcpp::Rcout << "ll is inf" << std::endl;
  }
  return ll;
}


// Function to calculate grouped means
NumericVector calculateGroupedMeans(NumericVector values, IntegerVector groups) {
  std::map<int, std::pair<double, int>> sum_counts;
  int max_group = IntegerVector::is_na(groups[0]) ? 0 : groups[0];
  
  // Sum values and count occurrences for each group
  for (int i = 0; i < values.size(); ++i) {
    if (!NumericVector::is_na(values[i]) && !IntegerVector::is_na(groups[i])) {
      sum_counts[groups[i]].first += values[i];
      sum_counts[groups[i]].second += 1;
      if (groups[i] > max_group) {
        max_group = groups[i];
      }
    }
  }
  
  // Calculate means
  NumericVector means(max_group + 1, NA_REAL); // Initialize all elements to NA
  for (const auto& pair : sum_counts) {
    means[pair.first] = pair.second.first / pair.second.second;
  }
  
  return means;
}



// [[Rcpp::export]]
NumericVector adjustProxyValues(NumericVector tarX, NumericVector tarProxyValue, NumericVector t_times, 
                                IntegerVector window_labels, NumericVector mean_values, NumericVector inpProxyValue) {
  // Step 1: Perform linear interpolation
  NumericVector new_tar = linearInterpolation(tarX, tarProxyValue, t_times);
  
  // Step 2: Calculate grouped means and calculate shifts
  NumericVector grouped_means = calculateGroupedMeans(new_tar, window_labels);
  NumericVector shifts = grouped_means - mean_values;
  
  // Step 3: Adjust inp$ProxyValue
  NumericVector adjustedProxyValue(inpProxyValue.size());
  for (int i = 0; i < inpProxyValue.size(); ++i) {
    adjustedProxyValue[i] = inpProxyValue[i] + shifts[window_labels[i]];
  }
  
  return adjustedProxyValue;
}



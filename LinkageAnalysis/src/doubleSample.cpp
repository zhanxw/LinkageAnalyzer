#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double doubleSample(NumericVector prob, IntegerVector numG3, int obs, int nTrial) {
  const int n = prob.size();

  int hit = 0;
  for (int trial = 0; trial < nTrial; ++trial) {
    int resampleObs = 0;
   for (int i = 0; i < n; ++i) {
// REprintf("prob[%d] = %f\\n", i, prob[i]);
// REprintf("numG3[%d] = %d\\n", i, numG3[i]);

      if (numG3[i] == 0 || prob[i] == 0.) {
        continue;
      }
      const int r = rbinom(1, numG3[i], prob[i])[0];
      resampleObs += r;
    }
    if (resampleObs <= obs) {
      hit++;
    }
  }
  if (nTrial == 0) {
    REprintf("nTrail = 0\\n");
    return 1.0;
  }
  return 1.0 * hit / nTrial;
}

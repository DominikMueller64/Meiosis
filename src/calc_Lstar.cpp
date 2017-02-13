// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <boost/math/tools/roots.hpp>
#include "calc_Lstar.hpp"

// [[Rcpp::interfaces(r,cpp)]]

typedef double (*func_t)(double, double, int, double);  // define function pointer

double func(double Lstar, double L, int m, double p)
{

  double denom = 0;
  if (m == 0) {
    denom = 1 - std::exp(-Lstar / 50);
  } else {
    double lambda1 = Lstar / 50 * (m + 1) * (1 - p);
    double lambda2 = Lstar / 50 * p;
    for (int i = 0; i < m + 1; i++) {
      denom += R::dpois(i, lambda1, false) * (m + 1 - i) / (m + 1);
    }
    denom = 1 - denom * std::exp(-lambda2);
  }

  return 2 * L - 2 * Lstar / denom;
}

// Functor for root finding.
class functor
{
  public:
  functor(double L, int m, double p, func_t f): L(L), m(m), p(p), f(f) {} // constructor
  double operator()(double Lstar) {
    // return func(Lstar, L, m, p);
    return f(Lstar, L, m, p);
  }

  private:
    double L, p;
    int m;
    func_t f;
};

// Functor for determining the terminal condition.
struct TermCond  {
  double epsilon;
  TermCond(double epsilon): epsilon(epsilon) {}
  bool operator() (double min, double max)  {
    return std::abs(min - max) <= epsilon;
  }
};


// [[Rcpp::export]]
double calc_Lstar(double L, int m, double p, double epsilon)
{

  if (L <= 50) {
    throw std::invalid_argument("L must be larger than 50.");
  }
  if (m < 0) {
    throw std::invalid_argument("m must be a non-negative integer.");
  }
  if (p < 0 || p > 1) {
    throw std::invalid_argument("p must be in [0, 1].");
  }
  if (1 - p < p * std::numeric_limits<double>::epsilon()){
    m = 0;
    p = 0;
  }
  if (epsilon <= 0) {
    throw std::invalid_argument("epsilon must be positive.");
  }
 // functor f(L, m, p);
 functor f(L, m, p, &func);
 std::pair<double, double> r = boost::math::tools::bisect(f, 1e-8, L, TermCond(epsilon));
 return r.first;
}

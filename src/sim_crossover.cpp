// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <random>
#include <chrono>
#include "globals.hpp"
#include "sim_crossover.hpp"

// [[Rcpp::interfaces(r,cpp)]]

// const unsigned random_seed = std::chrono::system_clock::now().time_since_epoch().count();
std::random_device rdev;
const unsigned random_seed = rdev();
std::mt19937 engine{random_seed};


std::vector<double> runif(std::size_t n,
                          double a,
                          double b)
{
  std::vector<double> ret(n);
  std::uniform_real_distribution<double> dist (a, b);
  for(auto& x : ret){
    x = dist(engine);
  }
  return ret;
}

// This function is fast for large L, m and small p.
// [[Rcpp::export]]
std::vector<double> crossover(const double L,
                              const int m,
                              const double p,
                              const bool obligate_chiasma,
                              const double Lstar)
{
  if(m == 0) { // no-interference model is a lot easier
    int n_xo = 0;

    if (obligate_chiasma) {
      // rejection sampling to get at least one chiasma
      std::poisson_distribution<int> poisson{Lstar / 50};
      // while ((n_xo = R::rpois(Lstar / 50)) == 0);
      while ((n_xo = poisson(engine)) == 0);

      std::binomial_distribution<int> binomial{n_xo, 0.5};
      n_xo = binomial(engine);
    }
    else {
      std::poisson_distribution<int> poisson{L / 100};
      n_xo = poisson(engine);
    }

    auto tmp = runif(n_xo, 0, L);
    std::sort(tmp.begin(), tmp.end());
    return tmp;
  }

  int n_points, first, n_nichi, n_ichi;
  const double lambda1 = Lstar / 50 * (double) (m + 1) * (1 - p);
  const double lambda2 = Lstar / 50 * p;

  while(true) {
    // chiasma and intermediate points
    std::poisson_distribution<int> poisson1{lambda1};
    n_points = poisson1(engine);

    // which point is the first chiasma?
    std::uniform_int_distribution<int> unif_int(0, m);
    first = unif_int(engine);
    if (first > n_points) n_ichi = 0;
    else n_ichi = n_points / (m + 1) + (int)(first < (n_points % (m + 1)));

    // no. chiasma from no interference process
    std::poisson_distribution<int> poisson2{lambda2};
    if (p > 0) n_nichi = poisson2(engine);
    else n_nichi = 0;

    if (!obligate_chiasma || n_ichi + n_nichi > 0) break;
  }

  // locations of chiasmata and intermediate points for process w/ interference
  auto point_locations = runif(n_points, 0, L);
  std::sort(point_locations.begin(), point_locations.end());

  // move every (m+1)st point back to front
  int n_chi = 0;
  for (std::size_t j = first; j < n_points; j += (m + 1), ++n_chi)
    point_locations.at(n_chi) = point_locations.at(j);

  // chiasma locations from non-interference process
  const auto nichi_locations = runif(n_nichi, 0, L);

  // combine interference and no interference chiasma locations
  std::vector<double> chi_locations(n_chi + n_nichi);
  std::copy(point_locations.begin(), point_locations.begin() + n_chi, chi_locations.begin());
  std::copy(nichi_locations.begin(), nichi_locations.end(), chi_locations.begin() + n_chi);
  std::sort(chi_locations.begin(), chi_locations.end());

  // thining by 0.5
  int n_xo = 0;
  const auto coins = runif(n_chi, 0, L);
  for(std::size_t i = 0; i < n_chi + n_nichi; ++i) {
    if (coins.at(i) < 0.5) {
      chi_locations.at(n_xo) = chi_locations.at(i);
      ++n_xo;
    }
  }
  std::vector<double> xo_locations(n_xo);
  std::copy(chi_locations.begin(), chi_locations.begin() + n_xo, xo_locations.begin());
  return xo_locations;
}

// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <Rcpp.h>
#include <vector>
#include <random>
#include <chrono>
#include "globals.hpp"

#ifndef crossover__hpp
#define crossover__hpp

namespace crossover_ns
{

    template<typename t_xlocations>
    // inline std::vector<double> runif(const std::size_t n,
    inline t_xlocations runif(const std::size_t n,
                              const double a,
                              const double b,
                              std::mt19937& engine)
    {
      // std::vector<double> ret(n);
      t_xlocations ret(n);
      std::uniform_real_distribution<double> distribution (a, b);
      auto dist = std::bind(distribution, std::ref(engine));
      for(auto& x : ret){
        x = dist();
      }
      return ret;
    }

    // This function is fast for large L, m and small p.
    template<typename t_xlocations>
    // inline std::vector<double> crossover(const double L,
    inline t_xlocations crossover(const double L,
                                  const int m,
                                  const double p,
                                  const bool obligate_chiasma,
                                  const double Lstar,
                                  std::mt19937& engine)
    {
      if (m == 0) { // no-interference model is a lot easier
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

        // auto tmp = runif(n_xo, 0, L, engine);
        auto tmp = runif<t_xlocations>(n_xo, 0, L, engine);
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
      // auto point_locations = runif(n_points, 0.0, L, engine);
      auto point_locations = runif<t_xlocations>(n_points, 0.0, L, engine);
      std::sort(point_locations.begin(), point_locations.end());

      // move every (m+1)st point back to front
      int n_chi = 0;
      for (std::size_t j = first; j < n_points; j += (m + 1), ++n_chi)
        point_locations[n_chi] = point_locations[j]; // no range-checks, faster

      // chiasma locations from non-interference process
      // const auto nichi_locations = runif(n_nichi, 0.0, L, engine);
      const auto nichi_locations = runif<t_xlocations>(n_nichi, 0.0, L, engine);

      // combine interference and no interference chiasma locations
      // std::vector<double> chi_locations(n_chi + n_nichi);
      t_xlocations chi_locations(n_chi + n_nichi);
      std::copy(point_locations.begin(), point_locations.begin() + n_chi, chi_locations.begin());
      std::copy(nichi_locations.begin(), nichi_locations.end(), chi_locations.begin() + n_chi);
      std::sort(chi_locations.begin(), chi_locations.end());

      // for (auto i: point_locations)
      //   std::cout << i << ' ';
      // std::cout << std::endl;

      // thining by 0.5
      int n_xo = 0;
      // const auto coins = runif(n_chi + n_nichi, 0.0, 1.0, engine);
      const auto coins = runif<t_xlocations>(n_chi + n_nichi, 0.0, 1.0, engine);
      for(std::size_t i = 0; i < n_chi + n_nichi; ++i) {
        if (coins.at(i) < 0.5) {
          chi_locations[n_xo] = chi_locations[i];
          ++n_xo;
        }
      }

      // std::vector<double> xo_locations(n_xo);
      t_xlocations xo_locations(n_xo);
      std::copy(chi_locations.begin(), chi_locations.begin() + n_xo, xo_locations.begin());

      return xo_locations;
    }

} // crossover_ns

#endif // crossover__hpp

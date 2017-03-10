// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <random>
#include <chrono>
#include "Meiosis.h"

std::random_device rdev2;
std::mt19937 engine2;

typedef std::vector<double> t_alleles;
typedef std::vector<double> t_locations;

double meiosis_geno_sum(const t_alleles& patalle,
                        const t_alleles& matalle,
                        const t_locations& xlocations,
                        const t_locations& pos,
                        const double patsum,
                        const double matsum,
                        std::mt19937& engine)
{
  // Draw the parent where to start from.
  std::uniform_int_distribution<int> dist(0, 1);
  bool cur_allele = (dist(engine) != 0);

  // If no crossover occurs, early exit.
  if (xlocations.size() == 0) {
    if (cur_allele == 0) {
      // return std::accumulate(matalle.begin(), matalle.end(), 0.0);
      return matsum;
    }
    else {
      // return std::accumulate(patalle.begin(), patalle.end(), 0.0);
      return patsum;
    }
  }

  std::size_t n_loci = pos.size();
  t_locations xloc(xlocations);
  xloc.push_back(*(pos.end() - 1) + 1e-6);
  double sum = 0.0;

  std::size_t ix1 = 0; // first slicing position
  std::size_t ix2; // second slicing position
  for (std::size_t j = 0; j < xloc.size(); ++j){
    const auto it = std::upper_bound(pos.begin() + ix1, pos.end(), xloc[j]);
    ix2 = it - pos.begin();
    if (cur_allele){
      for (std::size_t i = ix1; i != ix2; ++i) {
        sum += patalle[i];
      }
    } else{
      for (std::size_t i = ix1; i != ix2; ++i) {
        sum += matalle[i];
      }
    }
    cur_allele = !cur_allele; // switch to other chromatid
    ix1 = ix2;
  }
  return sum;
}


// [[Rcpp::export]]
Rcpp::List bcgv(std::vector<std::vector<std::vector<double>>> individual,
                std::vector<std::vector<double>> positions,
                std::vector<std::vector<double>> locus_effects,
                const int n_gam,
                const double se_level,
                const int min_rep,
                const int max_rep,
                const int m,
                const double p
                )
{
  auto t1 = std::chrono::high_resolution_clock::now();

  if (min_rep > max_rep)
    Rcpp::stop("'min_rep' must be smaller than or equal to 'max_rep'");


  engine2.seed(rdev2());
  const auto& xof = Meiosis::crossover<std::vector<double>>;
  // const auto& meif = Meiosis::meiosis_geno<std::vector<double>, std::vector<double>>;

  // prepare xodat parameters
  // const int m = 0;
  // const double p = 0.0;
  const bool obligate_chiasma = true;
  const std::size_t n_chr = positions.size();
  std::vector<double> L(n_chr), Lstar(n_chr);
  const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
  for (std::size_t i{0}; i != n_chr; ++i) {
    const auto& pos = positions[i];
    L[i] = std::ceil(pos.back());
    Lstar[i] = Meiosis::calc_Lstar(L[i], m, p, epsilon);
  }

  // transform genotypes by multiplication with locus effects and compute sums
  std::vector<std::vector<double>> sums(individual.size(), std::vector<double>(n_chr));
  for (std::size_t g{0}; g != individual.size(); ++g) {
  // for (auto& gamete: individual) {
    auto& gamete = individual[g];
    for (std::size_t i{0}; i != n_chr; ++i) {
      auto& chromatid = gamete[i];
      auto& leff = locus_effects[i];
      std::transform(chromatid.begin(), chromatid.end(),
                     leff.begin(), chromatid.begin(), std::multiplies<double>());
      sums[g][i] = std::accumulate(chromatid.begin(), chromatid.end(), 0.0);
    }
  }

  std::size_t n = 0;
  double mean = 0.0;
  double M2 = 0.0;
  double se = std::numeric_limits<double>::max();
  double delta;
  double delta2;

  while (n < min_rep || (se >= se_level && n < max_rep)) {
    double max_value = std::numeric_limits<double>::lowest();

    for (std::size_t g{0}; g != n_gam; ++g) {
      double value = 0.0;
      for (std::size_t i{0}; i != n_chr; ++i) {
        const auto& xlocations = xof(L[i], m, p, obligate_chiasma, Lstar[i], engine2);

        // const auto& tmp = meif(individual[0][i], individual[1][i], xlocations,
        //                        positions[i], engine2);
        // value = std::accumulate(tmp.begin(), tmp.end(), value);

        value += meiosis_geno_sum(individual[0][i], individual[1][i], xlocations,
                                 positions[i], sums[0][i], sums[1][i], engine2);
      }
      if (value > max_value) max_value = value;
    }

    ++n;
    delta = max_value - mean;
    mean += delta / n;
    delta2 = max_value - mean;
    M2 += delta * delta2;
    se = std::sqrt(M2 / (n * (n - 1)));
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  return Rcpp::List::create(Rcpp::Named("bcgv") = mean,
                            Rcpp::Named("se") = se,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("time (ms)") = duration / 1000.0);
}


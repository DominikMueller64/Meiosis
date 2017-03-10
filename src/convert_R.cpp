// // [[Rcpp::plugins(cpp11)]]
// #include <cstddef> // std::size_t
// #include <vector>
// #include <unordered_map>
// #include <algorithm>
// #include <Rcpp.h>
// #include "Meiosis.h"
// #include "founders.hpp"
// #include "convert_R.hpp"

// Rcpp::List convert_xodat(const Rcpp::List& individual,
//                          const Rcpp::List& pos,
//                          const Founders& founders)
// {
//   Rcpp::List new_ind(individual.size());

//   for (size_t g = 0; g != individual.size(); ++g) {
//     const Rcpp::List& gam = individual[g];
//     Rcpp::List new_gam(gam.size());
//     for (size_t k = 0; k != gam.size(); ++k) {
//       const Rcpp::NumericVector p = pos[k];
//       const Rcpp::List& chr = gam[k];
//       const Rcpp::IntegerVector& alleles = chr[0];
//       const Rcpp::NumericVector& locations = chr[1];

//       Rcpp::IntegerVector out(p.size());
//       std::size_t ix1 = 0;
//       std::size_t ix2;
//       for (std::size_t j = 0; j < locations.size(); ++j) {
//         auto it = std::upper_bound(p.begin() + ix1, p.end(), locations[j]);
//         ix2 = it - p.begin();
//         for (std::size_t i = ix1; i < ix2; ++i) {
//           out[i] = founders.obj.at(alleles[j])[k][i]; // alleles[j] founder, k-th chr, i-th locus
//         }
//         ix1 = ix2;
//       }
//       new_gam[k] = out;
//     }
//     new_ind[g] = new_gam;
//   }
//   return new_ind;
// }


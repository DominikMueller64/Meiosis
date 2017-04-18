// #include <Rcpp.h>
// #include <numeric> // std::accumulate
// #include <cstddef> // std::size_t

// // [[Rcpp::export("to_matrix")]]
// Rcpp::IntegerMatrix to_matrix(const Rcpp::List& x)
// {

//   const Rcpp::List& gam = x[0];
//   Rcpp::IntegerVector cols(gam.size());
//   for (std::size_t i = 0; i != gam.size(); ++i)
//     cols(i) = Rcpp::as<Rcpp::IntegerVector>(gam[i]).size();

//   std::size_t m = std::accumulate(cols.begin(), cols.end(), 0);
//   Rcpp::IntegerMatrix mat(2, m);

//   for (size_t ib = 0; ib != 2; ++ib) {
//     std::size_t cjb = 1;
//     const Rcpp::List& xib = x[ib];
//     for (size_t jb = 0; jb != cols.size(); ++jb) {
//       const Rcpp::IntegerVector& xijb = xib[jb];
//       for (size_t j = 0; j != cols(jb); ++j) {
//         mat(ib - 1, cjb + j - 1) = xijb(j);
//       } // j
//       cjb += cols(jb);
//     } // jb
//   } // ib

//   return mat;
// }


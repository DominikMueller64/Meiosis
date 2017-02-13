// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <tuple>
#include <random>
#include "globals.hpp"
#include "sim_meiosis.hpp"

// [[Rcpp::interfaces(r,cpp)]]

typedef std::vector<int> ivec;
typedef std::vector<double> vec;

// [[Rcpp::export("meiosis_xodat")]]
Rcpp::List meiosis_xodat_R(const Rcpp::IntegerVector& patalle,
                           const Rcpp::NumericVector& patloc,
                           const Rcpp::IntegerVector& matalle,
                           const Rcpp::NumericVector& matloc,
                           const Rcpp::NumericVector& xlocations)
{

  const auto& ret = meiosis_xodat<Rcpp::IntegerVector, Rcpp::NumericVector>(
                               patalle, patloc, matalle, matloc, xlocations);

  return Rcpp::List::create(Rcpp::Named("alleles") = ret.first,
                            Rcpp::Named("locations") = ret.second);
}

std::pair<ivec, vec > meiosis_xodat_cpp(const ivec& patalle,
                                        const vec& patloc,
                                        const ivec& matalle,
                                        const vec& matloc,
                                        const vec& xlocations)
{
  return meiosis_xodat<ivec, vec >(patalle, patloc, matalle, matloc, xlocations);
}

// //' @export
// // [[Rcpp::export("meiosis")]]
// Rcpp::List sim_meiosis_R(const Rcpp::List& paternal,
//                          const Rcpp::List& maternal,
//                          const std::vector<double>& xlocations)
// {

//   const auto& patalle_ = Rcpp::as<std::vector<size_t> >(paternal[0]);
//   const auto& patloc_ = Rcpp::as<std::vector<double> >(paternal[1]);
//   const auto& matalle_ = Rcpp::as<std::vector<size_t> >(maternal[0]);
//   const auto& matloc_ = Rcpp::as<std::vector<double> >(maternal[1]);
//   // const auto& xlocations_ = Rcpp::as<std::vector<double> >(xlocations);


//   // const Rcpp::IntegerVector& patalle = paternal[0];
//   // const Rcpp::NumericVector& patloc = paternal[1];
//   // const Rcpp::NumericVector& matalle = maternal[0];
//   // const Rcpp::NumericVector& matloc = maternal[1];

//   // const std::vector<std::size_t> patalle_(patalle.begin(), patalle.end());
//   // const std::vector<double> patloc_(patloc.begin(), patloc.end());
//   // const std::vector<std::size_t> matalle_(matalle.begin(), matalle.end());
//   // const std::vector<double> matloc_(matloc.begin(), matloc.end());
//   // const std::vector<double> xlocations_(xlocations.begin(), xlocations.end());

//   const auto& chromatid = sim_meiosis(Chromatid(patalle_, patloc_),
//                                       Chromatid(matalle_, matloc_), xlocations);
//   return Rcpp::List::create(Rcpp::Named("alleles") = chromatid.alleles,
//                             Rcpp::Named("locations") = chromatid.locations);
// }

// template<typename t_alleles, typename t_locations>

// t_chromatid meiosis_xodat(const t_alleles& patalle,
//                           const t_locations& patloc,
//                           const t_alleles& matalle,
//                           const t_locations& patloc,
//                           const t_locations& xlocations);
// {
//   // Draw the parent where to start from.
//   std::uniform_int_distribution<int> dist(0, 1);
//   int cur_allele = dist(engine);

//   // If no crossover occurs, early exit.
//   if (tmp.size() == 0) {
//     if (cur_allele == 0) {
//       return std::make_tuple(matalle, matloc);
//       // return maternal;
//     }
//     else {
//       return std::make_tuple(patalle, patloc);
//       // return paternal;
//     }
//   }

//   // std::vector<double> tmp(xlocations);

//   const std::size_t matloc_size = matloc.size();
//   const std::size_t patloc_size = patloc.size();

//   // const std::size_t product_size = tmp.size() + 1;
//   // std::vector<double> product(product_size);
//   // product[0] = -1.0;
//   // std::copy(tmp.begin(), tmp.end(), product.begin() + 1);

//   xlocations.insert(xlocations.begin(), -1.0);
//   const std::size_t product_size = xlocations.size() + 1;


//   const std::size_t biggest_length = product_size + matloc_size + patloc_size;
//   std::vector<double> loc(biggest_length);
//   std::vector<size_t> alle(biggest_length);

//   int curpos = 0;
//   std::size_t i, j;
//   std::size_t jmat = 0;
//   std::size_t jpat = 0;
//   double prod0, prod1;
//   for(i = 1; i < product_size; ++i) {
//     prod0 = product[i - 1];
//     prod1 = product[i];
//     if(cur_allele == 0) { // mat chr
//       while(matloc[jmat] < prod0) {
//         ++jmat;
//       }
//       while(matloc[jmat] < prod1 && jmat < matloc.size()) {
//         loc[curpos] = matloc[jmat];
//         alle[curpos] = matalle[jmat];
//         ++jmat;
//         ++curpos;
//       }
//       loc[curpos] = prod1;
//       alle[curpos] = matalle[jmat];
//       ++curpos;
//     }
//     else { // pat chr
//       while(patloc[jpat] < prod0) {
//         ++jpat;
//       }
//       while(patloc[jpat] < prod1 && jpat < patloc.size()) {
//         loc[curpos] = patloc[jpat];
//         alle[curpos] = patalle[jpat];
//         ++jpat;
//         ++curpos;
//       }
//       loc[curpos] = prod1;
//       alle[curpos] = patalle[jpat];
//       ++curpos;
//     }
//     cur_allele = 1 - cur_allele;
//   }

//   const double lastxo = product[product_size - 1];

//   if(cur_allele == 0) { // mat chr
//     for(j = 0; j < matloc_size; ++j) {
//       if(matloc[j] > lastxo) {
//         loc[curpos] = matloc[j];
//         alle[curpos] = matalle[j];
//         ++curpos;
//       }
//     }
//   }
//   else{ // pat chr
//     for(j = 0; j < patloc_size; ++j) {
//       if(patloc[j] > lastxo) {
//         loc[curpos] = patloc[j];
//         alle[curpos] = patalle[j];
//         ++curpos;
//       }
//     }
//   }

//   if(curpos > 1) { // clean up repeated alleles
//     std::size_t lastpos = 0;
//     for(i = 1; i < curpos; ++i) {
//       if(alle[lastpos] == alle[i]) {
//         loc[lastpos] = loc[i];
//       }
//       else {
//         ++lastpos;
//         loc[lastpos] = loc[i];
//         alle[lastpos] = alle[i];
//       }
//     }
//     curpos = lastpos + 1;
//   }
//   // copy over to short vectors
//   std::vector<double> loc_result(loc.begin(), loc.begin() + curpos);
//   std::vector<std::size_t> alle_result(alle.begin(), alle.begin() + curpos);
//   return std::make_tuple(alle_result, loc_result);
//   // return Chromatid(alle_result, loc_result, paternal, maternal);
// }

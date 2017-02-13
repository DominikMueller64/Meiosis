// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <tuple>
#include <random>
#include "globals.hpp"

#ifndef SIM_MEIOSIS_H
#define SIM_MEIOSIS_H

template<typename t_alleles, typename t_locations>
std::pair<t_alleles, t_locations> meiosis_xodat(const t_alleles& patalle,
                                                const t_locations& patloc,
                                                const t_alleles& matalle,
                                                const t_locations& matloc,
                                                const t_locations& xlocations)
{
  // Draw the parent where to start from.
  std::uniform_int_distribution<int> dist(0, 1);
  int cur_allele = dist(engine);

  // If no crossover occurs, early exit.
  if (xlocations.size() == 0) {
    if (cur_allele == 0) {
      return std::make_pair(matalle, matloc);
      // return maternal;
    }
    else {
      return std::make_pair(patalle, patloc);
      // return paternal;
    }
  }

  // std::vector<double> tmp(xlocations);
  const std::size_t matloc_size = matloc.size();
  const std::size_t patloc_size = patloc.size();
  const std::size_t product_size = xlocations.size() + 1;
  std::vector<double> product(product_size);
  product[0] = -1.0;
  std::copy(xlocations.begin(), xlocations.end(), product.begin() + 1);

  // This modifies xlocations, which is not desirable
  // xlocations.insert(xlocations.begin(), -1.0);
  // const std::size_t product_size = xlocations.size() + 1;


  const std::size_t biggest_length = product_size + matloc_size + patloc_size;
  std::vector<double> loc(biggest_length);
  std::vector<size_t> alle(biggest_length);

  int curpos = 0;
  std::size_t i, j;
  std::size_t jmat = 0;
  std::size_t jpat = 0;
  double prod0, prod1;
  for(i = 1; i < product_size; ++i) {
    prod0 = product[i - 1];
    prod1 = product[i];
    if(cur_allele == 0) { // mat chr
      while(matloc[jmat] < prod0) {
        ++jmat;
      }
      while(matloc[jmat] < prod1 && jmat < matloc.size()) {
        loc[curpos] = matloc[jmat];
        alle[curpos] = matalle[jmat];
        ++jmat;
        ++curpos;
      }
      loc[curpos] = prod1;
      alle[curpos] = matalle[jmat];
      ++curpos;
    }
    else { // pat chr
      while(patloc[jpat] < prod0) {
        ++jpat;
      }
      while(patloc[jpat] < prod1 && jpat < patloc.size()) {
        loc[curpos] = patloc[jpat];
        alle[curpos] = patalle[jpat];
        ++jpat;
        ++curpos;
      }
      loc[curpos] = prod1;
      alle[curpos] = patalle[jpat];
      ++curpos;
    }
    cur_allele = 1 - cur_allele;
  }

  const double lastxo = product[product_size - 1];

  if(cur_allele == 0) { // mat chr
    for(j = 0; j < matloc_size; ++j) {
      if(matloc[j] > lastxo) {
        loc[curpos] = matloc[j];
        alle[curpos] = matalle[j];
        ++curpos;
      }
    }
  }
  else{ // pat chr
    for(j = 0; j < patloc_size; ++j) {
      if(patloc[j] > lastxo) {
        loc[curpos] = patloc[j];
        alle[curpos] = patalle[j];
        ++curpos;
      }
    }
  }

  if(curpos > 1) { // clean up repeated alleles
    std::size_t lastpos = 0;
    for(i = 1; i < curpos; ++i) {
      if(alle[lastpos] == alle[i]) {
        loc[lastpos] = loc[i];
      }
      else {
        ++lastpos;
        loc[lastpos] = loc[i];
        alle[lastpos] = alle[i];
      }
    }
    curpos = lastpos + 1;
  }
  // copy over to short vectors
  t_alleles alle_result(alle.begin(), alle.begin() + curpos);
  t_locations loc_result(loc.begin(), loc.begin() + curpos);
  return std::make_pair(alle_result, loc_result);
  // return Chromatid(alle_result, loc_result, paternal, maternal);
}


#endif

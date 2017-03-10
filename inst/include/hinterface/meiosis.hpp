// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <tuple>
#include <utility>  // std::pair
#include <random>
#include "globals.hpp"
#include "crossover.hpp"

#ifndef meiosis__hpp
#define meiosis__hpp

namespace meiosis_ns
{
    typedef std::vector<int> ivec;
    typedef std::vector<double> vec;


    // Template function of simulating meiosis using genotypes.
    template<typename t_alleles, typename t_locations>
    t_alleles meiosis_geno(const t_alleles& patalle,
                           const t_alleles& matalle,
                           const t_locations& xlocations,
                           const t_locations& pos,
                           std::mt19937& engine)
    {
      // Draw the parent where to start from.
      std::uniform_int_distribution<int> dist(0, 1);
      bool cur_allele = (dist(engine) != 0);

      // If no crossover occurs, early exit.
      if (xlocations.size() == 0) {
        if (cur_allele == 0) {
          return matalle;
        }
        else {
          return patalle;
        }
      }

      std::size_t n_loci = pos.size();
      t_locations xloc(xlocations);
      xloc.push_back(*(pos.end() - 1) + 1e-6);
      t_alleles alleles(n_loci);

      std::size_t ix1 = 0; // first slicing position
      std::size_t ix2; // second slicing position
      for (std::size_t j = 0; j < xloc.size(); ++j){
        const auto it = std::upper_bound(pos.begin() + ix1, pos.end(), xloc[j]);
        ix2 = it - pos.begin();
        if (cur_allele){
          for (std::size_t i = ix1; i != ix2; ++i) {
            alleles[i] += patalle[i];
          }
        } else{
          for (std::size_t i = ix1; i != ix2; ++i) {
            alleles[i] += matalle[i];
          }
        }
        cur_allele = !cur_allele; // switch to other chromatid
        ix1 = ix2;
      }
      return alleles;
    }



  // Template function of simulating meiosis using crossover data.
    template<typename t_alleles, typename t_locations>
    std::pair<t_alleles, t_locations> meiosis_xodat(const t_alleles& patalle,
                                                    const t_locations& patloc,
                                                    const t_alleles& matalle,
                                                    const t_locations& matloc,
                                                    const t_locations& xlocations,
                                                    std::mt19937& engine)
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
    }


    inline std::vector<int> meiosis_geno_cpp(const std::vector<int>& patalle,
                                             const std::vector<int>& matalle,
                                             const std::vector<double>& xlocations,
                                             const std::vector<double>& pos,
                                             std::mt19937& engine)
    {
      return meiosis_geno<std::vector<int>, std::vector<double> >(patalle, matalle,
                                                                  xlocations, pos, engine);
    }


    inline std::pair<ivec,vec> meiosis_xodat_cpp(const ivec& patalle,
                                                 const vec& patloc,
                                                 const ivec& matalle,
                                                 const vec& matloc,
                                                 const vec& xlocations,
                                                 std::mt19937& engine)
    {
      return meiosis_xodat<ivec,vec>(patalle, patloc, matalle, matloc, xlocations, engine);
    }



} // meiosis_ns

#endif // meiosis__hpp

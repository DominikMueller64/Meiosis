// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <utility> // std::make_pair()
#include "Meiosis.h"
#include "Converter.hpp"

Converter::Converter(t_positions positions): positions(positions){}


std::size_t Converter::size(){return hashmapvec.size();}

void Converter::insert_founder(std::vector<int> keys, t_geno geno)
{
  insert_gamete(keys[0], geno[0]);
  insert_gamete(keys[1], geno[1]);
}


void Converter::insert_gamete(int key, t_gamete gamete)
{

  if (size() == 0) { // still empty
    const auto& n_chr = gamete.size();
    hashmapvec.reserve(n_chr); // reserve memory
    for (const auto& v: gamete) {
      t_hashmap hashmap; // instantiate new hashmap
      hashmap[key] = v; // insert first element
      hashmapvec.push_back(hashmap); // add to vector of hashmaps
    }
  }
  else { // not empty
    if (size() != gamete.size())
      Rcpp::stop("The length of 'value' is not conformable.");
    for (std::size_t i{0}; i != size(); ++i) {
      auto& hashmap = hashmapvec[i];
      auto search = hashmap.find(key);
      if (search != hashmap.end()) // element already exists
        hashmap.erase(key); // erase if it exists and insert new one
      hashmap[key] = gamete[i];
    }
  }
}

Rcpp::List Converter::convert_gamete(const Rcpp::List& gamete)
{
  const auto& f = Meiosis::convert<Rcpp::IntegerVector, Rcpp::NumericVector, t_pos, t_alleles>;
  Rcpp::List new_gam(gamete.size());
  for (size_t k = 0; k != gamete.size(); ++k) {
    const Rcpp::List& chr = gamete[k];
    new_gam[k] = f(chr[0], chr[1], positions[k], hashmapvec[k]);
  }
  return new_gam;
}

Rcpp::List Converter::convert(const Rcpp::List& individual)
{
  const auto& n_gam = individual.size();
  Rcpp::List new_ind(n_gam);
  for (size_t g = 0; g != n_gam; ++g) {
    new_ind[g] = convert_gamete(individual[g]);
  }
  return new_ind;
}
      // Rcpp::IntegerVector out(p.size());
      // std::size_t ix1 = 0;
      // std::size_t ix2;
      // for (std::size_t j = 0; j < locations.size(); ++j) {
      //   auto it = std::upper_bound(p.begin() + ix1, p.end(), locations[j]);
      //   ix2 = it - p.begin();
      //   for (std::size_t i = ix1; i < ix2; ++i) {
      //     // out[i] = founders.obj.at(alleles[j])[k][i]; // alleles[j] founder, k-th chr, i-th locus
      //     out[i] = hashmapvec[k][alleles[j]][i];
      //   }
      //   ix1 = ix2;
      // }
      // new_gam[k] = out;

// t_value Converter::at(int key) {return obj.at(key);}
// void Converter::erase(int key) {obj.erase(key);}
// bool Converter::empty(){return obj.empty();}
// void Converter::clear(){obj.clear();}


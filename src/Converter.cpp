// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <map>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <utility> // std::make_pair()
#include "Meiosis.h"
#include "Converter.hpp"

Converter::Converter(t_positions positions): positions(positions){}


std::size_t Converter::size(){return mapvec.size();}

void Converter::insert_founder(ivec keys, t_geno geno)
{
  insert_gamete(keys[0], geno[0]);
  insert_gamete(keys[1], geno[1]);
}


void Converter::insert_gamete(int key, t_gamete gamete)
{

  if (mapvec.empty()) { // still empty
    const auto& n_chr = gamete.size();
    mapvec.reserve(n_chr); // reserve memory
    for (const auto& v: gamete) {
      // boost::container::map<int, ivec> map;
      std::map<int, ivec> map;
      // t_hashmap map; // instantiate new hashmap
      map[key] = v; // insert first element
      mapvec.push_back(map); // add to vector of hashmaps
    }
  }
  else { // not empty
    if (size() != gamete.size())
      Rcpp::stop("The length of 'value' is not conformable.");
    for (std::size_t i{0}; i != size(); ++i) {
      auto& map = mapvec[i];
      if (map.find(key) != map.end()) // element already exists
        map.erase(key); // erase if it exists and insert new one
      map[key] = gamete[i];
    }
  }
}

Rcpp::List Converter::convert_gamete(const Rcpp::List& gamete)
{
  const auto& f = Meiosis::convert<Rcpp::IntegerVector, Rcpp::NumericVector, vec>;
  Rcpp::List new_gam(gamete.size());
  for (size_t k = 0; k != gamete.size(); ++k) {
    const Rcpp::List& chr = gamete[k];
    new_gam[k] = f(chr[0], chr[1], positions[k], mapvec[k]);
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

// Converter 2

Rcpp::IntegerVector conv(const Rcpp::IntegerVector& alleles,
                         const Rcpp::NumericVector& locations,
                         const std::vector<double>& pos,
                         const std::map<int, std::size_t>& map,
                         const std::vector<std::vector<int>>& mat)
{
  Rcpp::IntegerVector out(pos.size());
  std::size_t ix1 = 0;
  std::size_t ix2;
  for (std::size_t j = 0; j < locations.size(); ++j) {
    auto it = std::upper_bound(pos.begin() + ix1, pos.end(), locations[j]);
    ix2 = it - pos.begin();
    for (std::size_t i = ix1; i < ix2; ++i) {
      out[i] = mat[map.at(alleles[j])][i];
    }
    ix1 = ix2;
  }
  return out;
}


Converter2::Converter2(t_positions positions):
  positions(positions), counter(0){
  const auto n_chr = positions.size();
  // std::vector<int> n_loci;
  // n_loci.resize(n_chr);
  // for (const auto& p: positions) n_loci.push_back(p.size());
  matvec = std::vector<std::vector<std::vector<int>>>(n_chr); // initialize with n_chr chrom.
}

void Converter2::insert_gamete(int key, t_gamete gamete)
{
  if (map.find(key) != map.end())
    Rcpp::stop("Founder allele(s) already registered.");

  map[key] = counter++;
  for (std::size_t i{0}; i != gamete.size(); ++i)
    matvec[i].push_back(gamete[i]);
}

void Converter2::insert_founder(std::vector<int> keys, t_geno geno)
{
  insert_gamete(keys[0], geno[0]);
  insert_gamete(keys[1], geno[1]);
}



Rcpp::List Converter2::convert_gamete(const Rcpp::List& gamete)
{
  Rcpp::List new_gam(gamete.size());
  for (size_t k = 0; k != gamete.size(); ++k) {
    const Rcpp::List& chr = gamete[k];
    new_gam[k] = conv(chr[0], chr[1], positions[k], map, matvec[k]);
  }
  return new_gam;
}

Rcpp::List Converter2::convert(const Rcpp::List& individual)
{
  const auto& n_gam = individual.size();
  Rcpp::List new_ind(n_gam);
  for (size_t g = 0; g != n_gam; ++g) {
    new_ind[g] = convert_gamete(individual[g]);
  }
  return new_ind;
}



// Converter 3

Rcpp::IntegerVector conv3(const Rcpp::IntegerVector& alleles,
                         const Rcpp::NumericVector& locations,
                         const std::vector<double>& pos,
                         const std::map<int, std::size_t>& map,
                         const std::vector<int>& allvec)
{
  Rcpp::IntegerVector out(pos.size());
  std::size_t ix1 = 0;
  std::size_t ix2;
  std::size_t n = pos.size();
  for (std::size_t j = 0; j < locations.size(); ++j) {
    auto it = std::upper_bound(pos.begin() + ix1, pos.end(), locations[j]);
    ix2 = it - pos.begin();
    for (std::size_t i = ix1; i < ix2; ++i) {
      out[i] = allvec[map.at(alleles[j]) * n + i];
    }
    ix1 = ix2;
  }
  return out;
}


Converter3::Converter3(t_positions positions):
  positions(positions), counter(0){
  const auto n_chr = positions.size();
  allvec = std::vector<std::vector<int>>(n_chr); // initialize with n_chr chrom.
}

void Converter3::insert_gamete(int key, t_gamete gamete)
{
  if (map.find(key) != map.end())
    Rcpp::stop("Founder allele(s) already registered.");

  map[key] = counter++;
  for (std::size_t i{0}; i != gamete.size(); ++i)
    allvec[i].insert(allvec[i].end(), gamete[i].begin(), gamete[i].end());
}

void Converter3::insert_founder(std::vector<int> keys, t_geno geno)
{
  insert_gamete(keys[0], geno[0]);
  insert_gamete(keys[1], geno[1]);
}

Rcpp::List Converter3::convert_gamete(const Rcpp::List& gamete)
{
  Rcpp::List new_gam(gamete.size());
  for (size_t k = 0; k != gamete.size(); ++k) {
    const Rcpp::List& chr = gamete[k];
    new_gam[k] = conv3(chr[0], chr[1], positions[k], map, allvec[k]);
  }
  return new_gam;
}

Rcpp::List Converter3::convert(const Rcpp::List& individual)
{
  const auto& n_gam = individual.size();
  Rcpp::List new_ind(n_gam);
  for (size_t g = 0; g != n_gam; ++g) {
    new_ind[g] = convert_gamete(individual[g]);
  }
  return new_ind;
}

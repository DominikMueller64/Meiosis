// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cstdlib>

#ifndef __Converter__hpp
#define __Converter__hpp

// convenience typedef
typedef std::vector<int> t_alleles;
typedef std::vector<double> t_locations;
typedef std::vector<double> t_pos;
typedef std::vector<t_pos> t_positions;
// typedef std::vector<int> t_value;
typedef std::unordered_map<int, t_alleles> t_hashmap;
typedef std::vector<t_hashmap> t_hashmapvec;
typedef std::vector<t_alleles> t_gamete;
typedef std::vector<t_gamete> t_geno;

class Converter
{
public:
  Converter(t_positions positions);

  std::size_t size();
  void insert_gamete(int key, t_gamete gamete);
  void insert_founder(std::vector<int> keys, t_geno geno);

  Rcpp::List convert_gamete(const Rcpp::List& gamete);
  Rcpp::List convert(const Rcpp::List& individual);

  t_hashmapvec hashmapvec;
  t_positions positions;
};

#endif

// bool empty();
// void clear();
// t_value at(int key);
// void erase(int key);

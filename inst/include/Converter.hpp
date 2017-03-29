// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <map>
#include <vector>
#include <cstdlib>

#include <boost/container/flat_map.hpp>
#include <boost/container/map.hpp>

#ifndef __Converter__hpp
#define __Converter__hpp

// convenience typedef
// typedef std::vector<int> t_alleles;
// typedef std::vector<double> t_locations;
typedef std::vector<int> ivec;
typedef std::vector<double> vec;
typedef std::vector<vec> t_positions;

// typedef std::map<int, ivec> t_hashmap;
// typedef std::vector<t_hashmap> t_hashmapvec;

typedef std::vector<ivec> t_gamete;
typedef std::vector<t_gamete> t_geno;




class Converter
{
public:
  Converter(t_positions positions);

  std::size_t size();
  void insert_gamete(int key, t_gamete gamete);
  void insert_founder(ivec keys, t_geno geno);

  bool check_xodat_individual(const Rcpp::List& ind);
  Rcpp::List convert_gamete(const Rcpp::List& gamete);
  Rcpp::List convert(const Rcpp::List& individual);

  std::vector<std::map<int, ivec>> mapvec;
  // std::vector<std::unordered_map<int, ivec>> mapvec;
  t_positions positions;
};


// typedef std::vector<std::vector<std::vector<int>>> t_matvec;
// typedef std::map<int, std::size_t> t_map;

// class Converter2
// {
// public:
//   Converter2(t_positions positions);


//   void insert_gamete(int key, t_gamete gamete);
//   void insert_founder(std::vector<int> keys, t_geno geno);

//   Rcpp::List convert_gamete(const Rcpp::List& gamete);
//   Rcpp::List convert(const Rcpp::List& individual);

//   t_matvec matvec;
//   t_map map;
//   t_positions positions;
//   std::size_t counter;
// };


// class Converter3
// {
// public:
//   Converter3(t_positions positions);


//   void insert_gamete(int key, t_gamete gamete);
//   void insert_founder(std::vector<int> keys, t_geno geno);

//   Rcpp::List convert_gamete(const Rcpp::List& gamete);
//   Rcpp::List convert(const Rcpp::List& individual);

//   std::vector<std::vector<int>> allvec;
//   t_map map;
//   t_positions positions;
//   std::size_t counter;
// };

#endif

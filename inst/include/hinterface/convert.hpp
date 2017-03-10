// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <algorithm>

#ifndef convert__hpp
#define convert__hpp

namespace convert_ns
{
    // typedef std::vector<int> t_alleles;
    // typedef std::vector<int> t_genotypes;
    // typedef std::vector<double> t_locations;
    // typedef std::unordered_map<int, std::vector<int> > t_founder;

    // inline t_founder list2map(const Rcpp::List& founder)
    // {
    //   // Convert founder.
    //   const Rcpp::CharacterVector& keys = founder.attr("names");
    //   t_founder founder_;
    //   for (std::size_t i = 0; i < founder.size(); ++i){
    //     const std::string& key = Rcpp::as<std::string>(keys[i]);
    //     std::stringstream sstream(key);
    //     std::size_t key_;
    //     sstream >> key_;
    //     const auto& value_ = Rcpp::as<t_genotypes >(founder[i]);
    //     founder_.insert(std::make_pair(key_, value_));
    //   }
    //   return founder_;
    // }

    template<typename t_alleles, typename t_locations,
             typename t_positions, typename t_geno>
    inline t_alleles convert(const t_alleles& alleles,
                             const t_locations& locations,
                             const t_positions& pos,
                             const std::unordered_map<int, t_geno>& founder)
    {
      t_alleles out(pos.size());
      std::size_t ix1 = 0;
      std::size_t ix2;
      for (std::size_t j = 0; j < locations.size(); ++j) {
        auto it = std::upper_bound(pos.begin() + ix1, pos.end(), locations[j]);
        ix2 = it - pos.begin();
        for (std::size_t i = ix1; i < ix2; ++i) {
          out[i] = founder.at(alleles[j])[i];
        }
        ix1 = ix2;
      }
      return out;
    }

} // convert_ns
#endif



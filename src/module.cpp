// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "Meiosis.h"
#include "Converter.hpp"

RCPP_EXPOSED_CLASS(Converter);

RCPP_MODULE(Module){
  using namespace Rcpp ;
  class_<Converter>("Converter")
    .constructor<t_positions>("Construct converter (use_names = FALSE).")
    .constructor<t_positions, bool>("Construct converter.")
    .method("convert", &Converter::convert)
    .method("size", &Converter::size)
    .method("insert_founder", &Converter::insert_founder)
    .method("[[<-", &Converter::insert_founder)
    // .method("check_xodat_individual", &Converter::check_xodat_individual)
    ;

  // class_<Converter2>("Converter2")
  //   .constructor<t_positions>("Construct converter2.")
  //   .method("convert", &Converter2::convert)
  //   .method( "insert_founder", &Converter2::insert_founder)
  //   ;


  // class_<Converter3>("Converter3")
  //   .constructor<t_positions>("Construct converter3.")
  //   .method("convert", &Converter3::convert)
  //   .method( "insert_founder", &Converter3::insert_founder)
  //   ;
}



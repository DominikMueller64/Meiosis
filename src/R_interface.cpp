// [[Rcpp::plugins(cpp11)]]
#include <algorithm>
#include <cstddef>
#include <set>
#include <iterator>
#include "Meiosis.h"

#include <array>

std::random_device rdev;
std::mt19937 engine;

// [[Rcpp::export]]
double calc_Lstar(double L, int m, double p, double epsilon)
{
  return Meiosis::calc_Lstar(L, m, p, epsilon);
}

// [[Rcpp::export]]
Rcpp::NumericVector crossover(const double L,
                              const int m,
                              const double p,
                              const bool obligate_chiasma,
                              const double Lstar)
{
  engine.seed(rdev());
  return Meiosis::crossover<Rcpp::NumericVector>(L, m, p, obligate_chiasma, Lstar, engine);
}


// [[Rcpp::export]]
Rcpp::IntegerVector meiosis_geno_(const Rcpp::IntegerVector& patalle,
                                 const Rcpp::IntegerVector& matalle,
                                 const Rcpp::NumericVector& pos,
                                 const double L,
                                 const int m,
                                 const double p,
                                 const bool obligate_chiasma,
                                 const double Lstar
                                 )
{
  engine.seed(rdev()); // extremely costly, massive overhead!
  const auto& xlocations = Meiosis::crossover<std::vector<double>>(L, m, p, obligate_chiasma, Lstar, engine);
  return Meiosis::meiosis_geno<Rcpp::IntegerVector, std::vector<double>, Rcpp::NumericVector>(
        patalle, matalle, xlocations, pos, engine);
}

// [[Rcpp::export]]
Rcpp::List meiosis_geno(const Rcpp::List& individual,
                        const Rcpp::List& positions,
                        const Rcpp::List& xodat_param)
{
  engine.seed(rdev()); // costs about 3 microseconds
  const Rcpp::NumericVector& L = xodat_param[0];
  const int m = xodat_param[1];
  const double p = xodat_param[2];
  const bool obligate_chiasma = xodat_param[3];
  const Rcpp::NumericVector& Lstar = xodat_param[4];
  const auto& xo = Meiosis::crossover<std::vector<double>>;
  const auto& mei = Meiosis::meiosis_geno<Rcpp::IntegerVector,
                                          std::vector<double>,
                                          Rcpp::NumericVector>;

  Rcpp::List gamete(L.size());

  const auto& paternal = Rcpp::as<Rcpp::List>(individual[0]);
  const auto& maternal = Rcpp::as<Rcpp::List>(individual[1]);
  for (std::size_t i = 0; i != L.size(); ++i){
    // gamete[i] = meiosis_geno_(Rcpp::as<Rcpp::List>(individual[0])[i],
    //                           Rcpp::as<Rcpp::List>(individual[1])[i],
    //                           positions[i], L[i], m, p, obligate_chiasma, Lstar[i]);
    gamete[i] = mei(paternal[i], maternal[i],
                    xo(L[i], m, p, obligate_chiasma, Lstar[i], engine),
                    positions[i], engine);

  }
  return gamete;
}




// [[Rcpp::export]]
Rcpp::List meiosis_xodat_(const Rcpp::IntegerVector& patalle,
                         const Rcpp::NumericVector& patloc,
                         const Rcpp::IntegerVector& matalle,
                         const Rcpp::NumericVector& matloc,
                         const double L,
                         const int m,
                         const double p,
                         const bool obligate_chiasma,
                         const double Lstar
                         )
{
  engine.seed(rdev()); // costs about 3 microseconds
  const auto& xlocations = Meiosis::crossover<Rcpp::NumericVector>(L, m, p, obligate_chiasma, Lstar, engine);
  const auto& ret = Meiosis::meiosis_xodat<Rcpp::IntegerVector,
                                           Rcpp::NumericVector,
                                           Rcpp::NumericVector>(
       patalle, patloc, matalle, matloc, xlocations, engine);
  return Rcpp::List::create(Rcpp::Named("alleles") = ret.first,
                            Rcpp::Named("locations") = ret.second);
}



// [[Rcpp::export]]
Rcpp::List meiosis_xodat(const Rcpp::List& individual,
                         const Rcpp::List& xodat_param)
{
  engine.seed(rdev()); // costs about 3 microseconds
  const Rcpp::NumericVector& L = xodat_param[0];
  const int m = xodat_param[1];
  const double p = xodat_param[2];
  const bool obligate_chiasma = xodat_param[3];
  const Rcpp::NumericVector& Lstar = xodat_param[4];

  const Rcpp::List& paternal = individual[0];
  const Rcpp::List& maternal = individual[1];

  const auto& xo = Meiosis::crossover<std::vector<double>>;
  const auto& mei = Meiosis::meiosis_xodat<Rcpp::IntegerVector,
                                           Rcpp::NumericVector,
                                           std::vector<double>>;
  Rcpp::List gamete(L.size());
  for (std::size_t i = 0; i != L.size(); ++i){
    const Rcpp::List& pat = paternal[i];
    const Rcpp::List& mat = maternal[i];
    const auto& xlocations = xo(L[i], m, p, obligate_chiasma, Lstar[i], engine);
    const auto& ret = mei(pat[0], pat[1], mat[0], mat[1], xlocations, engine);
    gamete[i] = Rcpp::List::create(ret.first, ret.second);
  }
  return gamete;
}



// [[Rcpp::export]]
Rcpp::List cross_geno(const Rcpp::List& father,
                      const Rcpp::List& mother,
                      const Rcpp::List& positions,
                      const Rcpp::List& xodat_param)
{
  // const auto& pat = meiosis_geno(father, pos, L, m, p, obligate_chiasma, Lstar);
  // const auto& mat = meiosis_geno(mother, pos, L, m, p, obligate_chiasma, Lstar);
  const auto& pat = meiosis_geno(father, positions, xodat_param);
  const auto& mat = meiosis_geno(mother, positions, xodat_param);
  return Rcpp::List::create(Rcpp::Named("paternal") = pat,
                            Rcpp::Named("maternal") = mat);
}


// [[Rcpp::export]]
Rcpp::List cross_xodat(const Rcpp::List& father,
                       const Rcpp::List& mother,
                       const Rcpp::List& xodat_param)
{
  // const auto& pat = meiosis_xodat(father, L, m, p, obligate_chiasma, Lstar);
  // const auto& mat = meiosis_xodat(mother, L, m, p, obligate_chiasma, Lstar);
  const auto& pat = meiosis_xodat(father, xodat_param);
  const auto& mat = meiosis_xodat(mother, xodat_param);
  return Rcpp::List::create(Rcpp::Named("paternal") = pat,
                            Rcpp::Named("maternal") = mat);
}


// [[Rcpp::export]]
Rcpp::List dh_geno(const Rcpp::List& individual,
                   const Rcpp::List& positions,
                   const Rcpp::List& xodat_param)
{
  // const auto& gam = meiosis_geno(individual, pos, L, m, p, obligate_chiasma, Lstar);
  const auto& gam = meiosis_geno(individual, positions, xodat_param);
  return Rcpp::List::create(Rcpp::Named("paternal") = gam,
                            Rcpp::Named("maternal") = gam);
}


// [[Rcpp::export]]
Rcpp::List dh_xodat(const Rcpp::List& individual,
                    const Rcpp::List& xodat_param)
{
  // const auto& gam = meiosis_xodat(individual, L, m, p, obligate_chiasma, Lstar);
  const auto& gam = meiosis_xodat(individual, xodat_param);
  return Rcpp::List::create(Rcpp::Named("paternal") = gam,
                            Rcpp::Named("maternal") = gam);
}

// [[Rcpp::export]]
double realized_ibd(const Rcpp::List& individual_1,
                    const Rcpp::List& individual_2) {

  const Rcpp::List& pat_1 = individual_1[0];
  const Rcpp::List& mat_1 = individual_1[1];
  const Rcpp::List& pat_2 = individual_2[0];
  const Rcpp::List& mat_2 = individual_2[1];
  auto f = Meiosis::realized_ibd<Rcpp::IntegerVector, Rcpp::NumericVector>;

  std::size_t n_chr = pat_1.size();
  double len = 0;
  double tot = 0;

   for (std::size_t i = 0; i < n_chr; i++) {
     const Rcpp::List& pat_1_i = pat_1[i];
     const Rcpp::List& mat_1_i = mat_1[i];
     const Rcpp::List& pat_2_i = pat_2[i];
     const Rcpp::List& mat_2_i = mat_2[i];

     Rcpp::NumericVector tmp = pat_1_i[1];
     len += tmp[tmp.size() - 1];
     tot += f(pat_1_i[0], pat_1_i[1], pat_2_i[0], pat_2_i[1]) +
            f(pat_1_i[0], pat_1_i[1], mat_2_i[0], mat_2_i[1]) +
            f(mat_1_i[0], mat_1_i[1], pat_2_i[0], pat_2_i[1]) +
            f(mat_1_i[0], mat_1_i[1], mat_2_i[0], mat_2_i[1]);
   }
   return tot / (4.0 * len);
}


// [[Rcpp::export]]
double realized_f(const Rcpp::List& individual) {

  const Rcpp::List& pat = individual[0];
  const Rcpp::List& mat = individual[1];
  auto f = Meiosis::realized_ibd<Rcpp::IntegerVector, Rcpp::NumericVector>;

  std::size_t n_chr = pat.size();
  double len = 0;
  double tot = 0;

   for (std::size_t i = 0; i < n_chr; i++) {
     const Rcpp::List& pat_i = pat[i];
     const Rcpp::List& mat_i = mat[i];

     Rcpp::NumericVector tmp = pat_i[1];
     len += tmp[tmp.size() - 1];
     tot += f(pat_i[0], pat_i[1], mat_i[0], mat_i[1]);
   }
   return 0.5 * (1 + tot / len);
}


// // [[Rcpp::export(".meiosis_geno_std")]]
// std::vector<int> meiosis_geno_std(const std::vector<int>& patalle,
//                                   const std::vector<int>& matalle,
//                                   const std::vector<double>& pos,
//                                   const double L,
//                                   const int m,
//                                   const double p,
//                                   const bool obligate_chiasma,
//                                   const double Lstar
//                                   )
// {
//   engine.seed(rdev());
//   const auto& xlocations = Meiosis::crossover<std::vector<double> >(L, m, p, obligate_chiasma, Lstar, engine);
//   return Meiosis::meiosis_geno<std::vector<int>, std::vector<double>>(
//   patalle, matalle, xlocations, pos, engine);

// }


// // [[Rcpp::export("meiosis_geno_std")]]
// std::vector<std::vector<int> > meiosis_geno_ind_std(
//                                  const std::vector<std::vector<std::vector<int> > >& individual,
//                                  const std::vector<std::vector<double> >& pos,
//                                  const std::vector<double>& L,
//                                  const int m,
//                                  const double p,
//                                  const bool obligate_chiasma,
//                                  const std::vector<double>& Lstar)
// {
//   std::vector<std::vector<int> > gamete(L.size());
//   for (std::size_t i = 0; i != L.size(); ++i){
//     gamete[i] = meiosis_geno_std(individual[0][i], individual[1][i],
//                                  pos[i], L[i], m, p, obligate_chiasma, Lstar[i]);
//   }
//   return gamete;
// }


// // This appears to be (always) slower than the .rcpp version!
// // [[Rcpp::export(".meiosis_xodat_std")]]
// Rcpp::List meiosis_xodat_std(const std::vector<int>& patalle,
//                              const std::vector<double>& patloc,
//                              const std::vector<int>& matalle,
//                              const std::vector<double>& matloc,
//                              const double L,
//                              const int m,
//                              const double p,
//                              const bool obligate_chiasma,
//                              const double Lstar
//                              )
// {
//   engine.seed(rdev());
//   const auto& xlocations = Meiosis::crossover<std::vector<double> >(L, m, p, obligate_chiasma, Lstar, engine);
//   const auto& ret = Meiosis::meiosis_xodat<std::vector<int>, std::vector<double> >(
//                        patalle, patloc, matalle, matloc, xlocations, engine);
//   return Rcpp::List::create(Rcpp::Named("alleles") = ret.first,
//                             Rcpp::Named("locations") = ret.second);
// }

// // [[Rcpp::export]]
// Rcpp::List test(const Rcpp::List& individual,
//                 const Rcpp::List& xodat_param)
// {
//   const Rcpp::NumericVector& L = xodat_param[0];
//   int m = xodat_param[1];
//   double p = xodat_param[2];
//   bool obligate_chiasma = xodat_param[3];
//   const Rcpp::NumericVector& Lstar = xodat_param[4];

//   const Rcpp::List& paternal = individual[0];
//   const Rcpp::List& maternal = individual[1];
//   Rcpp::List gamete(L.size());
//   for (std::size_t i = 0; i != L.size(); ++i){
//     const Rcpp::List& pat = paternal[i];
//     const Rcpp::List& mat = maternal[i];
//     gamete[i] = meiosis_xodat(pat[0], pat[1], mat[0], mat[1], L[i],
//                               m, p, obligate_chiasma, Lstar[i]);
//   }
//   return gamete;
// }

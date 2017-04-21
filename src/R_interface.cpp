// [[Rcpp::plugins(cpp11)]]
#include <algorithm>
#include <cstddef>
#include <set>
#include <iterator>
#include <vector>
#include <tuple>
#include <limits>
#include <cmath>
#include "Meiosis.h"

std::random_device rdev;
std::mt19937 engine;
const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

//' @rdname seed_rng
// // [[Rcpp::export]]
int seed_rng(const Rcpp::Nullable<int>& seed = R_NilValue)
{
  int seed_;
  if (seed.isNotNull()) seed_ = Rcpp::as<int>(seed);
  else seed_ = rdev();
  engine.seed(seed_);
  return seed_;
}

//' @rdname calc_Lstar
// // [[Rcpp::export]]
double calc_Lstar(const double L,
                  const int m,
                  const double p,
                  const Rcpp::Nullable<double> epsilon = R_NilValue)
{
  if (epsilon.isNotNull()) {
    return Meiosis::calc_Lstar(L, m, p, Rcpp::as<double>(epsilon));
  }
  return Meiosis::calc_Lstar(L, m, p);
}

//' @rdname crossover
// // [[Rcpp::export]]
Rcpp::NumericVector crossover(const double L,
                              const int m,
                              const double p,
                              const bool obligate_chiasma,
                              const double Lstar
                              )
{
  return Meiosis::crossover<Rcpp::NumericVector>(L, m, p, obligate_chiasma, Lstar, engine);
}

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
  const auto& xlocations = Meiosis::crossover<std::vector<double>>(L, m, p, obligate_chiasma, Lstar, engine);
  return Meiosis::meiosis_geno<Rcpp::IntegerVector, std::vector<double>, Rcpp::NumericVector>(
        patalle, matalle, xlocations, pos, engine);
}

//' @rdname meiosis_geno
// // [[Rcpp::export]]
Rcpp::List meiosis_geno(const Rcpp::List& individual,
                        const Rcpp::List& positions,
                        const Rcpp::List& xoparam)
{
  const auto& L = Rcpp::as<Rcpp::NumericVector>(xoparam[0]);
  const auto& m = Rcpp::as<int>(xoparam[1]);
  const auto& p = Rcpp::as<double>(xoparam[2]);
  const auto& obligate_chiasma = Rcpp::as<bool>(xoparam[3]);
  const auto& Lstar = Rcpp::as<Rcpp::NumericVector>(xoparam[4]);
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

Rcpp::List meiosis_xo_(const Rcpp::IntegerVector& patalle,
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
  const auto& xlocations = Meiosis::crossover<Rcpp::NumericVector>(L, m, p, obligate_chiasma, Lstar, engine);
  const auto& ret = Meiosis::meiosis_xodat<Rcpp::IntegerVector,
                                           Rcpp::NumericVector,
                                           Rcpp::NumericVector>(
       patalle, patloc, matalle, matloc, xlocations, engine);
  return Rcpp::List::create(Rcpp::Named("alleles") = ret.first,
                            Rcpp::Named("locations") = ret.second);
}



//' @rdname meiosis_xo
// // [[Rcpp::export]]
Rcpp::List meiosis_xo(const Rcpp::List& individual,
                      const Rcpp::List& xoparam)
{
  const auto& L = Rcpp::as<Rcpp::NumericVector>(xoparam[0]);
  const auto& m = Rcpp::as<int>(xoparam[1]);
  const auto& p = Rcpp::as<double>(xoparam[2]);
  const auto& obligate_chiasma = Rcpp::as<bool>(xoparam[3]);
  const auto& Lstar = Rcpp::as<Rcpp::NumericVector>(xoparam[4]);
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

//' @rdname cross_geno
// // [[Rcpp::export]]
Rcpp::List cross_geno(const Rcpp::List& father,
                      const Rcpp::List& mother,
                      const Rcpp::List& positions,
                      const Rcpp::List& xoparam)
{
  const auto& pat = meiosis_geno(father, positions, xoparam);
  const auto& mat = meiosis_geno(mother, positions, xoparam);
  return Rcpp::List::create(Rcpp::Named("paternal") = pat,
                            Rcpp::Named("maternal") = mat);
}


//' @rdname cross_xo
// // [[Rcpp::export]]
Rcpp::List cross_xo(const Rcpp::List& father,
                       const Rcpp::List& mother,
                       const Rcpp::List& xoparam)
{
  const auto& pat = meiosis_xo(father, xoparam);
  const auto& mat = meiosis_xo(mother, xoparam);
  return Rcpp::List::create(Rcpp::Named("paternal") = pat,
                            Rcpp::Named("maternal") = mat);
}

//' @rdname self_geno
// // [[Rcpp::export]]
Rcpp::List self_geno(const Rcpp::List& individual,
                     const Rcpp::List& positions,
                     const Rcpp::List& xoparam)
{
  return cross_geno(individual, individual, positions, xoparam);
}


//' @rdname self_xo
// // [[Rcpp::export]]
Rcpp::List self_xo(const Rcpp::List& individual,
                   const Rcpp::List& xoparam)
{
  return cross_xo(individual, individual, xoparam);
}


//' @rdname dh_geno
// // [[Rcpp::export]]
Rcpp::List dh_geno(const Rcpp::List& individual,
                   const Rcpp::List& positions,
                   const Rcpp::List& xoparam)
{
  const auto& gam = meiosis_geno(individual, positions, xoparam);
  return Rcpp::List::create(Rcpp::Named("paternal") = gam,
                            Rcpp::Named("maternal") = gam);
}


//' @rdname dh_xo
// // [[Rcpp::export]]
Rcpp::List dh_xo(const Rcpp::List& individual,
                    const Rcpp::List& xoparam)
{
  const auto& gam = meiosis_xo(individual, xoparam);
  return Rcpp::List::create(Rcpp::Named("paternal") = gam,
                            Rcpp::Named("maternal") = gam);
}

double realized_coancestry_self(const Rcpp::List& individual) {

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

//' @rdname realized_heter
// // [[Rcpp::export]]
double realized_heter(const Rcpp::List& individual) {
  return 2 * (1 - realized_coancestry_self(individual));
}

//' @rdname realized_coancestry
// // [[Rcpp::export]]
double realized_coancestry(const Rcpp::List& individual_1,
                           const Rcpp::Nullable<Rcpp::List>& individual_2 = R_NilValue) {

  if (individual_2.isNull())
    return realized_coancestry_self(individual_1);

  const auto& individual_2_ = Rcpp::as<Rcpp::List>(individual_2);

  const Rcpp::List& pat_1 = individual_1[0];
  const Rcpp::List& mat_1 = individual_1[1];
  const Rcpp::List& pat_2 = individual_2_[0];
  const Rcpp::List& mat_2 = individual_2_[1];
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
     tot += f(pat_1_i[0], pat_1_i[1], pat_2_i[0], pat_2_i[1]);
     tot += f(pat_1_i[0], pat_1_i[1], mat_2_i[0], mat_2_i[1]);
     tot += f(mat_1_i[0], mat_1_i[1], pat_2_i[0], pat_2_i[1]);
     tot += f(mat_1_i[0], mat_1_i[1], mat_2_i[0], mat_2_i[1]);

   }
   return tot / (4.0 * len);
}

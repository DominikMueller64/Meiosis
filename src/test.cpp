// [[Rcpp::plugins(C++11)]]
#include <Rcpp.h>
#include <vector>
#include <tuple>
#include <random>
#include <Meiosis.h>
#include <algorithm>

std::random_device rdev55;
std::mt19937 engine55;

// std::pair<Rcpp::IntegerVector,
//           Rcpp::NumericVector> meiosis_xodat_test(const Rcpp::IntegerVector& patalle,
std::pair<std::vector<int>,
          std::vector<double>> meiosis_xodat_test(const std::vector<int>& patalle,
                                                  const std::vector<double>& patloc,
                                                  const std::vector<int>& matalle,
                                                  const std::vector<double>& matloc,
                                                  const std::vector<double>& xlocations,
                                                  std::mt19937& engine)
{
  // Draw the parent where to start from.
  std::uniform_int_distribution<int> dist(0, 1);
  int cur_allele = dist(engine);

  // If no crossover occurs, early exit.
  if (xlocations.size() == 0) {
    if (cur_allele == 0) {
      return std::make_pair(matalle, matloc);
    }
    else {
      return std::make_pair(patalle, patloc);
    }
  }

  const std::size_t matloc_size = matloc.size();
  const std::size_t patloc_size = patloc.size();
  const std::size_t product_size = xlocations.size() + 1;
  std::vector<double> product(product_size);
  product[0] = -1.0;
  std::copy(xlocations.begin(), xlocations.end(), product.begin() + 1);

  const std::size_t biggest_length = product_size + matloc_size + patloc_size;
  std::vector<int> alle(biggest_length);
  std::vector<double> loc(biggest_length);
  // Rcpp::IntegerVector alle(biggest_length);
  // Rcpp::NumericVector loc(biggest_length);

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
  // copy over to short vectors. This should be faster than resizing.
  std::vector<int> alle_result(alle.begin(), alle.begin() + curpos);
  std::vector<double> loc_result(loc.begin(), loc.begin() + curpos);
  return std::make_pair(alle_result, loc_result);
}


std::pair<std::vector<int>,
          std::vector<double>> meiosis_xodat_test2(const std::vector<int>& patalle,
                                                  const std::vector<double>& patloc,
                                                  const std::vector<int>& matalle,
                                                  const std::vector<double>& matloc,
                                                  const std::vector<double>& xlocations,
                                                  std::mt19937& engine)
{
  // Draw the parent where to start from.
  std::uniform_int_distribution<int> dist(0, 1);
  bool cur_allele = (dist(engine) != 0);

  // If no crossover occurs, early exit.
  if (xlocations.size() == 0) {
    if (cur_allele == 0) {
      return std::make_pair(matalle, matloc);
    }
    else {
      return std::make_pair(patalle, patloc);
    }
  }

  const std::size_t matloc_size = matloc.size();
  const std::size_t patloc_size = patloc.size();
  const std::size_t product_size = xlocations.size() + 1;
  std::vector<double> product(product_size);
  product[0] = -1.0;
  std::copy(xlocations.begin(), xlocations.end(), product.begin() + 1);

  const std::size_t biggest_length = product_size + matloc_size + patloc_size;
  std::vector<int> alle(biggest_length);
  std::vector<double> loc(biggest_length);
  auto italle = alle.begin();
  auto itloc = loc.begin();
  // Rcpp::IntegerVector alle(biggest_length);
  // Rcpp::NumericVector loc(biggest_length);

  int curpos = 0;
  std::size_t i, j;
  // std::size_t jmat = 0;
  // std::size_t jpat = 0;
  auto itmatloc = matloc.begin();
  auto itmatalle = matalle.begin();
  auto itpatloc = patloc.begin();
  auto itpatalle = patalle.begin();
  double prod0, prod1;
  for(i = 1; i < product_size; ++i) {
    prod0 = product[i - 1];
    prod1 = product[i];
    if(cur_allele == 0) { // mat chr
      itmatloc = std::lower_bound(itmatloc, matloc.end(), prod0);
      itmatalle += itmatloc - matloc.begin();
      while(*itmatloc < prod1) {
        *itloc++ = *itmatloc++;
        *italle++ = *itmatalle++;
      }
      *itloc++ = prod1;
      *italle++ = *itmatalle;
    }
    else { // pat chr
      itpatloc = std::lower_bound(itpatloc, patloc.end(), prod0);
      itpatalle += itpatloc - patloc.begin();
      while(*itpatloc < prod1) {
        *itloc++ = *itpatloc++;
        *italle++ = *itpatalle++;
      }
      *itloc++ = prod1;
      *italle++ = *itpatalle;
    }
    cur_allele = !cur_allele;
  }

  // const double lastxo = product[product_size - 1];

  if (cur_allele == 0) { // mat chr
    while (itmatloc != matloc.end()) {
      *itloc++ = *itmatloc++;
      *italle++ = *itmatalle++;
    }
    // for(j = 0; j < matloc_size; ++j) {
    //   if(matloc[j] > lastxo) {
    //     loc[curpos] = matloc[j];
    //     alle[curpos] = matalle[j];
    //     ++curpos;
    //   }
    // }
  }
  else { // pat chr
    while (itpatloc != patloc.end()) {
      *itloc++ = *itpatloc++;
      *italle++ = *itpatalle++;
    }
    // for(j = 0; j < patloc_size; ++j) {
    //   if(patloc[j] > lastxo) {
    //     loc[curpos] = patloc[j];
    //     alle[curpos] = patalle[j];
    //     ++curpos;
    //   }
    // }
  }

  curpos = itloc - loc.begin();

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
  // copy over to short vectors. This should be faster than resizing.
  std::vector<int> alle_result(alle.begin(), alle.begin() + curpos);
  std::vector<double> loc_result(loc.begin(), loc.begin() + curpos);
  return std::make_pair(alle_result, loc_result);

  // alle.resize(curpos);
  // loc.resize(curpos);
  // return std::make_pair(alle, loc);
}

// [[Rcpp::export]]
Rcpp::List meiosis_xodat_test_R(const std::vector<int>& patalle,
                                const std::vector<double>& patloc,
                                const std::vector<int>& matalle,
                                const std::vector<double>& matloc,
                                const double L,
                                const int m,
                                const double p,
                                const bool obligate_chiasma,
                                const double Lstar
                                )
{
  engine55.seed(rdev55()); // costs about 3 microseconds
  const auto& xlocations = Meiosis::crossover<std::vector<double>>(L, m, p, obligate_chiasma, Lstar, engine55);

  const auto& ret = meiosis_xodat_test(patalle, patloc, matalle, matloc, xlocations, engine55);
  return Rcpp::List::create(Rcpp::Named("alleles") = ret.first,
                            Rcpp::Named("locations") = ret.second);
}




// [[Rcpp::export]]
Rcpp::List meiosis_xodat_test_R2(const std::vector<int>& patalle,
                                const std::vector<double>& patloc,
                                const std::vector<int>& matalle,
                                const std::vector<double>& matloc,
                                const double L,
                                const int m,
                                const double p,
                                const bool obligate_chiasma,
                                const double Lstar
                                )
{
  engine55.seed(rdev55()); // costs about 3 microseconds
  const auto& xlocations = Meiosis::crossover<std::vector<double>>(L, m, p, obligate_chiasma, Lstar, engine55);

  const auto& ret = meiosis_xodat_test2(patalle, patloc, matalle, matloc, xlocations, engine55);
  return Rcpp::List::create(Rcpp::Named("alleles") = ret.first,
                            Rcpp::Named("locations") = ret.second);
}

// [[Rcpp::export]]
int find_lb(Rcpp::NumericVector x, double val) {
  return std::lower_bound(x.begin(), x.end(), val) - x.begin();
}


// [[Rcpp::export]]
int find_bf(Rcpp::NumericVector x, double val) {
  int i = 0;
  while(x[i] < val) ++i;
  return i;
}


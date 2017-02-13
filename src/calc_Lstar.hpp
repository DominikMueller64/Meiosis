#include <limits>

#ifndef CALC_LSTAR_H
#define CALC_LSTAR_H

double calc_Lstar(double L,
                  int m,
                  double p,
                  double epsilon = std::numeric_limits<double>::epsilon());

#endif

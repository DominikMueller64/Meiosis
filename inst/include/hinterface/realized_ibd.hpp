// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

#ifndef realized_ibd__hpp
#define realized_ibd__hpp

namespace realized_ibd_ns
{
  template<typename t_alleles, typename t_locations>
  inline double realized_ibd(const t_alleles& patalle,
                             const t_locations& patloc,
                             const t_alleles& matalle,
                             const t_locations& matloc)
  {
    double tot = 0.0;
    int ipatmax = patloc.size();
    int imatmax = matloc.size();
    int imat = 0;
    int ipat = 0;
    int cur_alle = (matloc[ipat] <= patloc[imat]);
    double left = 0;
    double right;
    while(true) {
      if (cur_alle) {
        do
        {
          right = matloc[imat];
          if (patalle[ipat] == matalle[imat]) {
            tot += (right - left);
          }
          imat++;
          left = right;
        } while (patloc[ipat] >= matloc[imat] && imat < imatmax);
        cur_alle = 1 - cur_alle;
      }
      else {
        do
        {
          right = patloc[ipat];
          if (patalle[ipat] == matalle[imat]) {
            tot += (right - left);
          }
          ipat++;
          left = right;
        } while (patloc[ipat] < matloc[imat] && ipat < ipatmax);
        cur_alle = 1 - cur_alle;
      }
      if (ipat == ipatmax || imat == imatmax ) {
        break;
      }
    }
    return tot;
  }

} // realized_ibd_ns
#endif

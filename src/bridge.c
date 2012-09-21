/*
Copyright 2010-2011 Gabriele Sales <gabriele.sales@unipd.it>


This file is part of parmigene.

knnmi is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License 
version 3 as published by the Free Software Foundation.

knnmi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with parmigene. If not, see <http://www.gnu.org/licenses/>.
*/

#include "mi.h"
#include "points.h"
#include <math.h>


unsigned int gen_seed(const double* const cs, const int n, const int k) {
  return n * k * ((int)cs[n/2]*100);
}

//compute gcc for all pairs of genes
void c_cor_all(const int* const corIndexp, double* const xs, int* const xsix, const int* const lp, const int* const np, const int* const kp, const double* const noisep, double* res) {
  const int corIndex = *corIndexp;
  const int l = *lp;
  const int n = *np;
  const int k = *kp;
  const double noise = *noisep;
  int xnormed[l];

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    int i, j;
    unsigned int seed = gen_seed(xs, l*n, k);

#ifdef _OPENMP
    #pragma omp for nowait
#endif
    for (i = 0; i < l; i++) {
      double* const p = xs+(i*n);
      //xnormed[i] = normalize(p, n);
      xnormed[i] = 1.0;
      //add_noise(p, n, noise, &seed);
    }

    #pragma omp for
    for (i = 0; i < l; i++)
      res[i*l+i] = 0.0;

    mi_t mi;
    make_mi(&mi, n, k);

#ifdef _OPENMP
    #pragma omp for schedule(dynamic)
#endif
    for (i = 1; i < l; i++){
    //  res[i*l+j] = 1.0;
      for (j = 0; j < i; j++){
       	  if( corIndex == 1 )  res[i*l+j] = res[j*l+i] = c_gcc(&mi, xs+(i*n), xs+(j*n), xsix+(i*n), xsix+(j*n) );  //GCC
          else if(corIndex == 2 ) res[i*l+j] = res[j*l+i] = c_pcc(&mi, xs+(i*n), xs+(j*n));                        //PCC
          else if(corIndex == 3 ) res[i*l+j] = res[j*l+i] = c_scc(&mi, xs+(i*n), xs+(j*n), xsix+(i*n), xsix+(j*n) ); //SCC
          else if( corIndex == 4 ) res[i*l+j] = res[j*l+i] = c_kcc(&mi, xs+(i*n), xs+(j*n));                         //KCC
          else                     res[i*l+j] = res[j*l+i] = NAN; 
      }
    }

    destroy_mi(&mi);
  }
}




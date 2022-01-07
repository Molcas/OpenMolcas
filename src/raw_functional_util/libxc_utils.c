/*
************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2022 Molecular Sciences Software Institute             *
************************************************************************
*/

#include <xc.h>
#include <assert.h>

/* Check whether the functional is a mixed functional */
int libxc_num_aux_funcs(int id) {
  /* Initialize functional */
  xc_func_type func;
  int naux;
  assert(xc_func_init(&func, id, XC_UNPOLARIZED) == 0);
  naux = func.n_func_aux;
  xc_func_end(&func);
  return naux;
}

/* Get the functional IDs and coefficients */
void libxc_get_aux_funcs(int id, int *ids, double *weights) {
  /* Initialize functional */
  xc_func_type func;
  assert(xc_func_init(&func, id, XC_UNPOLARIZED) == 0);
  for (int i = 0; i < func.n_func_aux; i++) {
    ids[i] = func.func_aux[i]->info->number;
    weights[i] = func.mix_coef[i];
  }
  xc_func_end(&func);
}

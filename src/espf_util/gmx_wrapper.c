/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
#ifdef _GROMACS_

#include "gmx_wrapper.h"
#include <stdio.h>

// Define wrappers for all Gromacs functions used in Molcas

FILE *fp ;

struct t_commrec *init_commrec_(void) {
  return init_commrec() ;
}

gmx_mmslave_t mmslave_init_(const struct t_commrec *cr, const char *log) {
  fp = fopen(log,"a") ;
  return mmslave_init(cr) ;
}

void mmslave_done_(gmx_mmslave_t gms) {
  fclose(fp) ;
  mmslave_done(gms) ;
}

int mmslave_read_tpr_(const char *tpr, gmx_mmslave_t gms) {
  return mmslave_read_tpr(tpr, fp, gms) ;
}

int mmslave_ngroups_(gmx_mmslave_t gms) {
  return mmslave_ngroups(gms) ;
}

int mmslave_natoms_(gmx_mmslave_t gms) {
  return mmslave_natoms(gms) ;
}

int mmslave_copyx_(gmx_mmslave_t gms, int natoms, rvec *x) {
  return mmslave_copyX(gms, natoms, x) ;
}

int mmslave_copyv_(gmx_mmslave_t gms, int natoms, rvec *v) {
  return mmslave_copyV(gms, natoms, v) ;
}

int mmslave_copyf_(gmx_mmslave_t gms, int natoms, rvec *f) {
  return mmslave_copyF(gms, natoms, f) ;
}

void mmslave_clean_(gmx_mmslave_t gms) {
  mmslave_clean(gms) ;
}

int mmslave_set_q_(gmx_mmslave_t gms, atom_id id, double q) {
  return mmslave_set_q(gms, id, q) ;
}

int mmslave_get_q_(gmx_mmslave_t gms, atom_id id, double *q) {
  return mmslave_get_q(gms, id, q) ;
}

int mmslave_get_atomnumber_(gmx_mmslave_t gms, atom_id id) {
  return mmslave_get_atomnumber(gms, id) ;
}

int mmslave_get_group_id_(gmx_mmslave_t gms, atom_id id) {
  return mmslave_get_group_id(gms, id) ;
}

int mmslave_calc_energy_(gmx_mmslave_t gms, const rvec *x, rvec *f, rvec *A, real *phi, double *energy) {
  return mmslave_calc_energy(gms, fp, x, f, A, phi, energy) ;
}

#else
typedef int make_iso_compilers_happy;
#endif

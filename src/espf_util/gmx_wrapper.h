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

#include "gromacs/mmslave.h"

struct t_commrec *init_commrec_(void) ;

gmx_mmslave_t mmslave_init_(const struct t_commrec *cr, const char *log) ;

void mmslave_done_(gmx_mmslave_t gms) ;

int mmslave_read_tpr_(const char *tpr, gmx_mmslave_t gms) ;

int mmslave_ngroups_(gmx_mmslave_t gms) ;

int mmslave_natoms_(gmx_mmslave_t gms) ;

int mmslave_copyx_(gmx_mmslave_t gms, int natoms, rvec *x) ;

int mmslave_copyv_(gmx_mmslave_t gms, int natoms, rvec *v) ;

int mmslave_copyf_(gmx_mmslave_t gms, int natoms, rvec *f) ;

void mmslave_clean_(gmx_mmslave_t gms) ;

int mmslave_set_q_(gmx_mmslave_t gms, atom_id id, double q) ;

int mmslave_get_q_(gmx_mmslave_t gms, atom_id id, double *q) ;

int mmslave_get_atomnumber_(gmx_mmslave_t gms, atom_id id) ;

int mmslave_get_group_id_(gmx_mmslave_t gms, atom_id id) ;

int mmslave_calc_energy_(gmx_mmslave_t gms, const rvec *x, rvec *f, rvec *A, real *phi, double *energy) ;

#endif

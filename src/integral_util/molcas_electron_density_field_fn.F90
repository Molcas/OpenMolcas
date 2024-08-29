!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************

#include "compiler_features.h"
#ifdef _EFP_

function Molcas_ELECTRON_DENSITY_FIELD_FN(n_pt,xyz,field,user_field)
!***********************************************************************
!                                                                      *
!     This is a callback routine for the computation of the electric   *
!     field (variable field), at a set of points in space (n_pt of     *
!     them, coordinates stored in the variable xyz).                   *
!                                                                      *
!***********************************************************************

use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_double, c_ptr
use EFP, only: EFP_RESULT_SUCCESS
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(c_int) :: Molcas_ELECTRON_DENSITY_FIELD_FN
integer(c_size_t), intent(in) :: n_pt
real(c_double), intent(in) :: xyz(3,n_pt)
real(c_double), intent(out) :: field(3,n_pt)
type(c_ptr), intent(in) :: user_field
integer(kind=iwp) :: nDens
integer(c_size_t) :: nCmp, nOrdOp
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: D1ao(:)

#include "macros.fh"
unused_var(user_field)

! Pick up the one-particle density matrix.

call Qpg_dArray('D1ao',Found,nDens)
if (Found .and. (nDens /= 0)) then
  call mma_allocate(D1ao,nDens)
  call Get_dArray_chk('D1ao',D1ao,nDens)
else
  write(u6,*) 'D1ao not found!'
  call abend()
end if

! Compute the electric field at the points xyz

nCmp = 3
nOrdOp = 1
call Drv1_Pot(D1ao,xyz,field,n_pt,ncmp,nordop)

call mma_deallocate(D1ao)

Molcas_ELECTRON_DENSITY_FIELD_FN = EFP_RESULT_SUCCESS

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(Molcas_ELECTRON_DENSITY_FIELD_FN)

#endif

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

#ifdef _EFP_
function Molcas_ELECTRON_DENSITY_FIELD_FN(n_pt,xyz,field,&user_field)
!***********************************************************************
!                                                                      *
!     This is a callback routine for the computation of the electric   *
!     field (variable field), at a set of points in space (n_pt of     *
!     them, coordinates stored in the variable xyz).                   *
!                                                                      *
!***********************************************************************

use EFP
use iso_c_binding, only: c_int, c_size_t, c_double, c_ptr
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer(c_int) :: Molcas_ELECTRON_DENSITY_FIELD_FN
integer(c_size_t), intent(in) :: n_pt
real(c_double), intent(in) :: xyz(3,n_pt)
real(c_double), intent(out) :: field(3,n_pt)
type(c_ptr) :: user_field
type(c_ptr) :: Dummy
real*8, allocatable :: D1ao(:)
logical Found
integer(c_size_t) :: nCmp, nOrdOp
integer :: nDens

Dummy = user_field

! Pick up the one-particle density matrix.

call Qpg_dArray('D1ao',Found,nDens)
if (Found .and. (nDens /= 0)) then
  call mma_allocate(D1ao,nDens)
  call Get_dArray_chk('D1ao',D1ao,nDens)
else
  write(6,*) 'D1ao not found!'
  call abend()
end if

! Compute the electric field at the points xyz

nCmp = 3
nOrdOp = 1
call Drv1_Pot(D1ao,xyz,field,n_pt,ncmp,nordop)

call mma_deallocate(D1ao)

Molcas_ELECTRON_DENSITY_FIELD_FN = EFP_RESULT_SUCCESS
#else
function Molcas_ELECTRON_DENSITY_FIELD_FN()

implicit none
integer Molcas_ELECTRON_DENSITY_FIELD_FN

Molcas_ELECTRON_DENSITY_FIELD_FN = 1
#endif

return

end function Molcas_ELECTRON_DENSITY_FIELD_FN

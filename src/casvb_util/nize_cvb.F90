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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine nize_cvb(c,nnrm,s,n,metr,ierr)
! Normalizes NNRM vectors in C.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nnrm, n, metr, ierr
real(kind=wp) :: c(n,nnrm), s(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i1
real(kind=wp) :: cnrm
logical(kind=iwp) :: safe
real(kind=wp), parameter :: thresh = 1.0e-8_wp
integer(kind=iwp), external :: mstackr_cvb
real(kind=wp), external :: ddot_, dnrm2_

if (metr /= 0) i1 = mstackr_cvb(n)
safe = ierr /= 0
do i=1,nnrm
  if (metr == 0) then
    cnrm = dnrm2_(n,c(1,i),1)
  else
    call saoon_cvb(c(1,i),work(i1),1,s,n,metr)
    cnrm = sqrt(ddot_(n,c(1,i),1,work(i1),1))
  end if
  if (safe .and. (cnrm < thresh)) then
    ierr = ierr+1
  else
    call dscal_(n,One/cnrm,c(1,i),1)
  end if
end do
if (metr /= 0) call mfreer_cvb(i1)

return

end subroutine nize_cvb

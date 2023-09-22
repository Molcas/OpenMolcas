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

subroutine schmidt2_cvb(c,sxc,cnrm,nvec,sao,n,metr)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nvec, n, metr
real(kind=wp) :: c(n,nvec), sxc(n,nvec), cnrm(nvec), sao(*)
integer(kind=iwp) :: i, j
real(kind=wp), parameter :: thresh = 1.0e-20_wp
real(kind=wp), external :: ddot_

do i=1,nvec
  do j=1,i-1
    if (cnrm(j) > thresh) call daxpy_(n,-ddot_(n,c(1,i),1,sxc(1,j),1)/cnrm(j),c(1,j),1,c(1,i),1)
  end do
  if (metr /= 0) call saoon_cvb(c(1,i),sxc(1,i),1,sao,n,metr)
  cnrm(i) = ddot_(n,c(1,i),1,sxc(1,i),1)
end do

return

end subroutine schmidt2_cvb

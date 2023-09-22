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

subroutine schmidtt2_cvb(c,sxc,nvec,t,nt,sao,n,metr)

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nvec, nt, n, metr
real(kind=wp) :: c(n,nvec), sxc(n,nvec), t(nt,nvec), sao(*)
real(kind=wp) :: thresh = 1.0e-20_wp
integer(kind=iwp) :: i, j
real(kind=wp) :: cnrm, fac
real(kind=wp), external :: ddot_

do i=1,nvec
  do j=1,i-1
    fac = -ddot_(n,c(1,i),1,sxc(1,j),1)
    call daxpy_(n,fac,c(1,j),1,c(1,i),1)
    call daxpy_(nt,fac,t(1,j),1,t(1,i),1)
  end do
  if (metr /= 0) call saoon_cvb(c(1,i),sxc(1,i),1,sao,n,metr)
  cnrm = ddot_(n,c(1,i),1,sxc(1,i),1)
  if (cnrm < thresh) then
    write(u6,*) ' Warning : near-singularity in orthonormalization.'
    write(u6,*) ' Vector norm :',cnrm
  end if
  fac = One/sqrt(cnrm)
  call dscal_(n,fac,c(1,i),1)
  if (metr /= 0) call dscal_(n,fac,sxc(1,i),1)
  call dscal_(nt,fac,t(1,i),1)
end do

return

end subroutine schmidtt2_cvb

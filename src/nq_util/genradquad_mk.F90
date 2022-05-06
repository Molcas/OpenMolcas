!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine GenRadQuad_MK(R,nR,nR_Eff,m,Alpha,iNQ)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "debug.fh"
real*8 R(2,nR-1), Alpha, m

! Last point at infinity is eliminated

if (Debug) then
  write(6,*) 'Log3 Algorithm (Mura-Knowles)'
  write(6,*) 'Alpha,m=',Alpha,m
  write(6,*) 'nR=',nR
end if
do iR=1,nR-1
  x = dble(iR)/dble(nR)
  R(1,iR) = -Alpha*log(One-x**m)
  R(2,iR) = R(1,iR)**2*Alpha*m*x**(m-One)/(One-x**m)/dble(nR)
end do
nR_Eff = nR-1

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iNQ)

end subroutine GenRadQuad_MK

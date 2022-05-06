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

subroutine GenRadQuad_TA(R,nR,nR_Eff,Alpha)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "debug.fh"
real*8 R(2,nR-1), Alpha

! Last point at infinity is eliminated

if (Debug) then
  write(6,*) 'Treutler-Ahlrichs Algorithm'
  write(6,*) 'Alpha=',Alpha
  write(6,*) 'nR=',nR
end if
Fact = Alpha/log(Two)
F6 = 0.6d00
do iR=1,nR-1
  x = Two*dble(iR)/dble(nR)-One
  F1 = (One+x)
  F2 = log(Two/(One-x))
  R(1,iR) = Fact*F1**F6*F2
  R(2,iR) = R(1,iR)**2*Fact*(F6*F1**(F6-One)*F2+F1**F6*One/(One-x))*Two/dble(nR)
end do
nR_Eff = nR-1

return

end subroutine GenRadQuad_TA

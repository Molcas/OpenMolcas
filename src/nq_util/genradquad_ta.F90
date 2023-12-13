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

use Constants, only: One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nR
real(kind=wp), intent(out) :: R(2,nR-1)
integer(kind=iwp), intent(out) :: nR_Eff
real(kind=wp), intent(in) :: Alpha
integer(kind=iwp) :: iR
real(kind=wp) :: F1, F2, F6, Fact, x

! Last point at infinity is eliminated

#ifdef _DEBUGPRINT_
write(u6,*) 'Treutler-Ahlrichs Algorithm'
write(u6,*) 'Alpha=',Alpha
write(u6,*) 'nR=',nR
#endif
Fact = Alpha/log(Two)
F6 = 0.6_wp
do iR=1,nR-1
  x = Two*real(iR,kind=wp)/real(nR,kind=wp)-One
  F1 = (One+x)
  F2 = log(Two/(One-x))
  R(1,iR) = Fact*F1**F6*F2
  R(2,iR) = R(1,iR)**2*Fact*(F6*F1**(F6-One)*F2+F1**F6*One/(One-x))*Two/real(nR,kind=wp)
end do
nR_Eff = nR-1

return

end subroutine GenRadQuad_TA

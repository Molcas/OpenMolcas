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

subroutine GenRadQuad_MK(R,nR,nR_Eff,m,Alpha)

use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nR
real(kind=wp), intent(out) :: R(2,nR-1)
integer(kind=iwp), intent(out) :: nR_Eff
real(kind=wp), intent(in) :: m, Alpha
integer(kind=iwp) :: iR
real(kind=wp) :: x

! Last point at infinity is eliminated

#ifdef _DEBUGPRINT_
write(u6,*) 'Log3 Algorithm (Mura-Knowles)'
write(u6,*) 'Alpha,m=',Alpha,m
write(u6,*) 'nR=',nR
#endif
do iR=1,nR-1
  x = real(iR,kind=wp)/real(nR,kind=wp)
  R(1,iR) = -Alpha*log(One-x**m)
  R(2,iR) = R(1,iR)**2*Alpha*m*x**(m-One)/(One-x**m)/real(nR,kind=wp)
end do
nR_Eff = nR-1

return

end subroutine GenRadQuad_MK

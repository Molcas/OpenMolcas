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

function SelfEn(nChg)
! Compute the self interaction energy of input point charges

use external_centers, only: XF
use Constants, only: Zero, Three
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: SelfEn
integer(kind=iwp), intent(in) :: nChg
integer(kind=iwp) :: iChg, iPL, jChg
real(kind=wp) :: E, Qi, Qj, R, R2, R3, R5, X, XMui, XMuj, Y, YMui, YMuj, Z, ZMui, ZMuj
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()
E = Zero
do iChg=2,nChg
  Qi = XF(4,iChg)
  XMui = XF(5,iChg)
  YMui = XF(6,iChg)
  ZMui = XF(7,iChg)
  do jChg=1,iChg-1
    X = XF(1,iChg)-XF(1,jChg)
    Y = XF(2,iChg)-XF(2,jChg)
    Z = XF(3,iChg)-XF(3,jChg)
    Qj = XF(4,jChg)
    XMuj = XF(5,jChg)
    YMuj = XF(6,jChg)
    ZMuj = XF(7,jChg)
    R2 = X*X+Y*Y+Z*Z
    R = sqrt(R2)
    R3 = R*R2
    R5 = R2*R3
    ! q_i interacts with q_j and mu_j
    if (Qi /= Zero) E = E+Qi*(Qj/R-XMuj*X/R3-YMuj*Y/R3-ZMuj*Z/R3)
    ! mu_x_i interacts with q_j and mu_j
    if (XMui /= Zero) E = E+XMui*(-Qj*X/R3+XMuj*(Three*X*X-R2)/R5+YMuj*Three*X*Y/R5+ZMuj*Three*X*Z/R5)
    ! mu_y_i interacts with q_j and mu_j
    if (YMui /= Zero) E = E+YMui*(-Qj*Y/R3+YMuj*(Three*Y*Y-R2)/R5+XMuj*Three*X*Y/R5+ZMuj*Three*Y*Z/R5)
    ! mu_z_i interacts with q_j and mu_j
    if (ZMui /= Zero) E = E+ZMui*(-Qj*Z/R3+ZMuj*(Three*Z*Z-R2)/R5+XMuj*Three*X*Z/R5+YMuj*Three*Y*Z/R5)
  end do
end do
SelfEn = E
if ((nChg > 0) .and. (iPL >= 2)) write(u6,'(A,F16.10)') ' (For info only) Self Energy of the charges =',SelfEn

return

end function SelfEn

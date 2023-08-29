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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine FindAM(n,X,XMax,IMax)
!-----------------------------------------------------------------------
! Function : Find Absolute Maximum in X(n)    XMax = X(IMax)
!-----------------------------------------------------------------------

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: X(n)
real(kind=wp), intent(out) :: XMax
integer(kind=iwp), intent(out) :: IMax
integer(kind=iwp) :: I
real(kind=wp) :: XVal

IMax = 1
XMax = Zero

do I=1,n
  XVal = abs(X(I))
  if (XMax < XVal) then
    XMax = XVal
    IMax = I
  end if
end do

return

end subroutine FindAM

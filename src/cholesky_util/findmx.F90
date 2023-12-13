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

function FindMx(n,X)
!-----------------------------------------------------------------------
! Function : Find Maximum norm in X(n)
!-----------------------------------------------------------------------

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: FindMx
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: X(n)
integer(kind=iwp) :: I
real(kind=wp) :: XMax

XMax = Zero
do I=1,n
  XMax = max(XMax,abs(X(I)))
end do
FindMx = XMax

return

end function FindMx

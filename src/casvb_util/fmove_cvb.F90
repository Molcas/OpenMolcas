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

!IFG trivial
subroutine FMOVE_CVB(IA,IB,N)

use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp) :: N
real(kind=wp) :: IA(N*RtoI), IB(N*RtoI)
integer(kind=iwp) :: I

!call DCOPY_(N*RtoI,IA,1,IB,1)
do I=1,N*RtoI
  IB(I) = IA(I)
end do

return

end subroutine FMOVE_CVB

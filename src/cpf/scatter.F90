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

subroutine SCATTER(N,A,IND,B)

#include "intent.fh"

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, IND(N)
real(kind=wp), intent(_OUT_) :: A(*)
real(kind=wp), intent(in) :: B(N)
integer(kind=iwp) :: I

do I=1,N
  A(IND(I)) = B(I)
end do

return

end subroutine SCATTER

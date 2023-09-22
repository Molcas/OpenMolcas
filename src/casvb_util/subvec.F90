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
subroutine SUBVEC(A,B,C,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: N
real(kind=wp) :: A(N), B(N), C(N)
integer(kind=iwp) :: I

do I=1,N
  A(I) = B(I)-C(I)
end do

return

end subroutine SUBVEC

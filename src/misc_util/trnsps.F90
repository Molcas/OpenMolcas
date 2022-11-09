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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************

subroutine TRNSPS(N,M,A,B)
! Return B = Transpose of A.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, M
real(kind=wp), intent(in) :: A(N,M)
real(kind=wp), intent(out) :: B(M,N)
integer(kind=iwp) :: I

do I=1,N
  B(:,I) = A(I,:)
end do

return

end subroutine TRNSPS

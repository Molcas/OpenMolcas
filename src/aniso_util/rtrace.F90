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

subroutine rtrace(N,A,B)
! removes the trace of a Real array

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N ! size of the array
real(kind=wp), intent(in) :: A(N) ! input
real(kind=wp), intent(out) :: B(N) ! output
real(kind=wp) :: AS

! compute the equal-weighted average
AS = sum(A(:))/real(N,kind=wp)
! translate each of the elements
B(:) = A(:)-AS

return

end subroutine rtrace

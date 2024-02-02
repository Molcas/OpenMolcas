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

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: N ! size of the array
real(kind=8), intent(in) :: A(N) ! input
real(kind=8), intent(out) :: B(N) ! output
! local variables
integer :: i
real(kind=8) :: AS

AS = 0.0_wp
call dcopy_(N,[0.0_wp],0,B,1)
! compute the equal-weighted average
do i=1,N
  AS = AS+A(i)/dble(N)
end do
! translate each of the elements
do i=1,N
  B(i) = A(i)-AS
end do

return

end subroutine rtrace

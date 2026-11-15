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

!> Fix the phase of the vector A: change its sign if its projection on a
!> standard vector with random structure is negative. The standard vector
!> has the components x_i/(2**31-1) - 1/2, x_i from the minimal standard
!> generator x_i = 16807 x_(i-1) mod (2**31-1) with the same start value on
!> every call, so the phase depends only on A (not on earlier calls or on
!> MOLCAS_RANDOM).
subroutine VecPhase(A,nA)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp
use, intrinsic :: iso_fortran_env, only: int64

implicit none
integer(kind=iwp), intent(in) :: nA
real(kind=wp), intent(inout) :: A(nA)
integer(kind=iwp) :: i
integer(kind=int64) :: iRan
real(kind=wp) :: Phase
integer(kind=int64), parameter :: iMult = 16807_int64, iMod = 2147483647_int64, iStart = 314159265_int64

Phase = Zero
iRan = iStart
do i=1,nA
  iRan = mod(iMult*iRan,iMod)
  Phase = Phase+A(i)*(real(iRan,kind=wp)/real(iMod,kind=wp)-Half)
end do
if (Phase < Zero) A(:) = -A(:)

end subroutine VecPhase

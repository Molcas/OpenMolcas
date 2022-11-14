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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Chk_Unitary(irc,U,n,Thr)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: check that U is unitary.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: U(n,n), Thr
integer(kind=iwp) :: i, n2
real(kind=wp) :: RMS, x2
real(kind=wp), allocatable :: UTU(:,:)
real(kind=wp), external :: ddot_

if (n < 1) then
  irc = 0
  return
end if

n2 = n**2
call mma_allocate(UTU,n,n,label='UTU')

call DGEMM_('T','N',n,n,n,One,U,n,U,n,Zero,UTU,n)
do i=1,n
  UTU(i,i) = UTU(i,i)-One
end do

x2 = real(n2,kind=wp)
RMS = sqrt(dDot_(n2,UTU,1,UTU,1)/x2)
if (RMS > Thr) then
  irc = 1
else
  irc = 0
end if

call mma_deallocate(UTU)

end subroutine Chk_Unitary

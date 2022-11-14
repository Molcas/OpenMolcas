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
! Copyright (C) 2009, Per Ake Malmqvist                                *
!***********************************************************************

subroutine NewGS(N,S,C,Temp,M)
! PAM 2009, New Gram-Schmidt 090216

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
integer(kind=iwp), intent(out) :: M
real(kind=wp), intent(in) :: S(N,N)
real(kind=wp), intent(inout) :: C(N,N)
real(kind=wp), intent(out) :: Temp(N)
integer(kind=iwp) :: i, isfail, k, Loop
real(kind=wp) :: ovl, X, xn2, Y
real(kind=wp), external :: DDOT_

M = 0
outer: do i=1,N
  X = S(i,i)
  if (X < 1.0e-6_wp) cycle
  Y = One/sqrt(X)
  C(:,M+1) = Zero
  C(i,M+1) = Y
  Temp(:) = Y*S(:,i)

  Loop = 0
  do
    Loop = Loop+1
    do k=1,M
      ovl = DDOT_(N,Temp,1,C(1,k),1)
      C(:,M+1) = C(:,M+1)-ovl*C(:,k)
    end do
    call dGeMV_('N',N,N,One,S,N,C(1,M+1),1,Zero,Temp,1)
    xn2 = DDOT_(N,Temp,1,C(1,M+1),1)

    if (xn2 < 1.0e-6_wp) cycle outer

    Y = One/sqrt(xn2)
    C(:,M+1) = Y*C(:,M+1)
    call dGeMV_('N',N,N,One,S,N,C(1,M+1),1,Zero,Temp,1)
    if ((Loop /= 1) .or. (Y <= 100.0_wp)) exit
  end do
  ! MGD issues with many states : to be very safe, test the result
  isfail = 0
  do k=1,M
    ovl = DDOT_(N,Temp,1,C(1,k),1)
    if (abs(ovl) > 1.0e-4_wp) isfail = 1
  end do
  if (isfail == 0) M = M+1

end do outer

return

end subroutine NewGS

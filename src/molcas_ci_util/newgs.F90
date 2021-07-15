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
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp), intent(in) :: N
integer(kind=iwp), intent(out) :: M
real(kind=wp), intent(in) :: S(N,N)
real(kind=wp), intent(inout) :: C(N,N)
real(kind=wp), intent(out) :: Temp(N)
integer(kind=iwp) :: i, isfail, k, Loop
real(kind=wp) :: ovl, X, xn2, Y
real(kind=r8), external :: DDOT_

M = 0
do i=1,N
  X = S(i,i)
  if (X < 1.0e-6_wp) goto 90
  Y = One/sqrt(X)
  call dcopy_(N,[Zero],0,C(1,M+1),1)
  C(i,M+1) = Y
  call dcopy_(N,S(1,i),1,Temp,1)
  call DSCAL_(N,Y,Temp,1)

  Loop = 0
10 continue
  Loop = Loop+1
  do k=1,M
    ovl = DDOT_(N,Temp,1,C(1,k),1)
    call daxpy_(N,-ovl,C(1,k),1,C(1,M+1),1)
  end do
  call dGeMV_('N',N,N,One,S,N,C(1,M+1),1,Zero,Temp,1)
  xn2 = DDOT_(N,Temp,1,C(1,M+1),1)

  if (xn2 < 1.0e-6_wp) goto 90

  Y = One/sqrt(xn2)
  call DSCAL_(N,Y,C(1,M+1),1)
  call dGeMV_('N',N,N,One,S,N,C(1,M+1),1,Zero,Temp,1)
  if ((Loop == 1) .and. (Y > 100.0_wp)) goto 10
  ! MGD issues with many states : to be very safe, test the result
  isfail = 0
  do k=1,M
    ovl = DDOT_(N,Temp,1,C(1,k),1)
    if (abs(ovl) > 1.0e-4_wp) isfail = 1
  end do
  if (isfail == 0) M = M+1

90 continue
end do

return

end subroutine NewGS

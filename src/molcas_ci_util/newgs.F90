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
! Copyright (C) 2016, Per Ake Malmqvist                                *
!***********************************************************************

subroutine NewGS(N,S,C,Temp,M)

implicit real*8(A-H,O-Z)
dimension S(N,N), C(N,N), Temp(N)
! PAM 2009, New Gram-Schmidt 090216
#include "warnings.h"

M = 0
do i=1,N
  X = S(i,i)
  if (X < 1.0D-6) goto 90
  Y = 1.0d0/sqrt(X)
  call dcopy_(N,[0.0d0],0,C(1,M+1),1)
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
  call dGeMV_('N',N,N,1.0d0,S,N,C(1,M+1),1,0.0d0,Temp,1)
  xn2 = DDOT_(N,Temp,1,C(1,M+1),1)

  if (xn2 < 1.0D-6) goto 90

  Y = 1.0d0/sqrt(xn2)
  call DSCAL_(N,Y,C(1,M+1),1)
  call dGeMV_('N',N,N,1.0d0,S,N,C(1,M+1),1,0.0d0,Temp,1)
  if ((Loop == 1) .and. (Y > 100.0d0)) goto 10
  ! MGD issues with many states : to be very safe, test the result
  isfail = 0
  do k=1,M
    ovl = DDOT_(N,Temp,1,C(1,k),1)
    if (abs(ovl) > 1.0d-4) isfail = 1
  end do
  if (isfail == 0) M = M+1

90 continue
end do

return

end subroutine NewGS

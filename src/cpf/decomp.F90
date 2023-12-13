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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine DECOMP(NN,A)

use cpf_global, only: IPS
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NN
real(kind=wp), intent(inout) :: A(NN,NN)
integer(kind=iwp) :: I, IDXPIV, IP, J, K, KP, KP1, N, NM1
real(kind=wp) :: BIG, EM, PIVOT, ROWNRM, SCALES(200), RSIZE

N = NN
IDXPIV = 0 ! dummy initialize
! INITIALIZE IPS AND SCALES
do I=1,N
  IPS(I) = I
  ROWNRM = Zero
  do J=1,N
    if (ROWNRM < abs(A(I,J))) ROWNRM = abs(A(I,J))
  end do
  if (ROWNRM == Zero) then
    call SING(1)
    SCALES(I) = Zero
  else
    SCALES(I) = One/ROWNRM
  end if
end do
! GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
NM1 = N-1
do K=1,NM1
  BIG = Zero
  do I=K,N
    IP = IPS(I)
    RSIZE = abs(A(IP,K))*SCALES(IP)
    if (RSIZE > BIG) then
      BIG = RSIZE
      IDXPIV = I
    end if
  end do
  if (BIG == 0) then
    call SING(2)
  else
    if (IDXPIV /= K) then
      J = IPS(K)
      IPS(K) = IPS(IDXPIV)
      IPS(IDXPIV) = J
    end if
    KP = IPS(K)
    PIVOT = A(KP,K)
    KP1 = K+1
    do I=KP1,N
      IP = IPS(I)
      EM = -A(IP,K)/PIVOT
      A(IP,K) = -EM
      do J=KP1,N
        A(IP,J) = A(IP,J)+EM*A(KP,J)
      end do
    end do
  end if
end do
KP = IPS(N)
if (A(KP,N) == 0) call SING(2)

return

end subroutine DECOMP

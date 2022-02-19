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

subroutine DECOMP(NN,A,UL)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NN
real(kind=wp) :: A(NN,NN), UL(NN,NN)
#include "ips.fh"
integer(kind=iwp) :: I, IDXPIV, IP, J, K, KP, KP1, N, NM1
real(kind=wp) :: BIG, EM, PIVOT, ROWNRM, SCALES(200), RSIZE

N = NN
IDXPIV = 0 ! dummy initialize
! INITIALIZE IPS, UL AND SCALES
do I=1,N
  IPS(I) = I
  ROWNRM = Zero
  do J=1,N
    UL(I,J) = A(I,J)
    if (ROWNRM-abs(UL(I,J)) >= 0) goto 2
    ROWNRM = abs(UL(I,J))
2   continue
  end do
  if (ROWNRM == 0) goto 4
  SCALES(I) = One/ROWNRM
  GO TO 5
4 call SING(1)
  SCALES(I) = Zero
5 continue
end do
! GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
NM1 = N-1
do K=1,NM1
  BIG = Zero
  do I=K,N
    IP = IPS(I)
    RSIZE = abs(UL(IP,K))*SCALES(IP)
    if (RSIZE-BIG <= 0) goto 11
    BIG = RSIZE
    IDXPIV = I
11  continue
  end do
  if (BIG /= 0) goto 13
  call SING(2)
  goto 17
13 if (IDXPIV == K) goto 15
  J = IPS(K)
  IPS(K) = IPS(IDXPIV)
  IPS(IDXPIV) = J
15 KP = IPS(K)
  PIVOT = UL(KP,K)
  KP1 = K+1
  do I=KP1,N
    IP = IPS(I)
    EM = -UL(IP,K)/PIVOT
    UL(IP,K) = -EM
    do J=KP1,N
      UL(IP,J) = UL(IP,J)+EM*UL(KP,J)
    end do
  end do
17 continue
end do
KP = IPS(N)
if (UL(KP,N) /= 0) goto 19
call SING(2)

19 return

end subroutine DECOMP

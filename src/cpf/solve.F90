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

subroutine SOLVE(NN,UL,B,X)

implicit real*8(A-H,O-Z)
dimension UL(NN,NN), B(*), X(*)

#include "ips.fh"

N = NN
NP1 = N+1

IP = IPS(1)
X(1) = B(IP)
do I=2,N
  IP = IPS(I)
  IM1 = I-1
  SUM = 0.0d00
  do J=1,IM1
    SUM = SUM+UL(IP,J)*X(J)
  end do
  X(I) = B(IP)-SUM
end do
IP = IPS(N)
X(N) = X(N)/UL(IP,N)
do IBACK=2,N
  I = NP1-IBACK
  ! I GOES (N-1),...,1
  IP = IPS(I)
  IP1 = I+1
  SUM = 0.0d00
  do J=IP1,N
    SUM = SUM+UL(IP,J)*X(J)
  end do
  X(I) = (X(I)-SUM)/UL(IP,I)
end do

return

end subroutine SOLVE

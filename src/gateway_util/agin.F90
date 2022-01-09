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

subroutine AGIN

implicit real*8(A-H,O-Z)
!...  auxiliar constant pool:       ready only up to g-valence/g-core
#include "const.fh"

DFAC(1) = 1.d0
DFAC(2) = 1.d0
do I=3,lp12
  DFAC(I) = DFAC(I-2)*dble(I-1)
end do
do J=1,lp13
  do I=1,lp1
    RCA(I,J) = 0.d0
  end do
end do

! RCA(i,j) = c(k) (la,0;lb,0) * sqrt((2*la+1)(2*lb+1))
!   j identifies la,lb:  ss,ps,pp,ds,dp,dd,... 1,2,3,4,5,6,...
!   i identifies k: i=1 for the lowest k whose c(k) is non-zero,
!                   i=2 for the next   k  "     "    "   ", etc.

RCA(1,1) = 1.d0
RCA(1,2) = 1.d0/3.d0
RCA(1,3) = 1.d0/3.d0
RCA(2,3) = 2.d0/15.d0
RCA(1,4) = 1.d0/5.d0
RCA(1,5) = 2.d0/15.d0
RCA(2,5) = 3.d0/35.d0
RCA(1,6) = 1.d0/5.d0
RCA(2,6) = 2.d0/35.d0
RCA(3,6) = 2.d0/35.d0
RCA(1,7) = 1.d0/7.d0
RCA(1,8) = 3.d0/35.d0
RCA(2,8) = 4.d0/63.d0
RCA(1,9) = 3.d0/35.d0
RCA(2,9) = 4.d0/105.d0
RCA(3,9) = 10.d0/231.d0
RCA(1,10) = 1.d0/7.d0
RCA(2,10) = 4.d0/105.d0
RCA(3,10) = 2.d0/77.d0
RCA(4,10) = 100.d0/(13.d0*231.d0)
RCA(1,11) = 1.d0/9.d0
RCA(1,12) = 4.d0/63.d0
RCA(2,12) = 10.d0/198.d0
RCA(1,13) = 2.d0/35.d0
RCA(2,13) = 20.d0/693.d0
RCA(3,13) = 10.d0/286.d0
RCA(1,14) = 4.d0/63.d0
RCA(2,14) = 2.d0/77.d0
RCA(3,14) = 20.d0/1001.d0
RCA(4,14) = 70.d0/2574.d0
RCA(1,15) = 1.d0/9.d0
RCA(2,15) = 20.d0/693.d0
RCA(3,15) = 162.d0/9009.d0
RCA(4,15) = 20.d0/1287.d0
RCA(5,15) = 490.d0/21879.d0

IJ = 0
do I=1,lp1
  do J=1,I
    IJ = IJ+1
    KOSUU(IJ) = J
  end do
end do

ICOL = 0
do L=1,lp1
  do I=1,L
    ICOL = ICOL+1
    IVAL = L-I-2
    do IROW=1,I
      IVAL = IVAL+2
      NYU(IROW,ICOL) = IVAL
    end do
  end do
end do

return

end subroutine AGIN

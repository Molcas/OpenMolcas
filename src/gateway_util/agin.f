************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE AGIN
      IMPLICIT REAL*8 (A-H,O-Z)
C...  auxiliar constant pool:       ready only up to g-valence/g-core
      PARAMETER (lp1=5,lp12=lp1*lp1,lp13=(lp1*lp1+lp1)/2)
      COMMON/CONST/RCA(lp1,lp13),DFAC(lp12),KOSUU(lp13),NYU(lp1,lp13)
C
      DFAC(1)=1.D0
      DFAC(2)=1.D0
      DO 10 I=3,lp12
        DFAC(I)=DFAC(I-2)*DBLE(I-1)
10    CONTINUE
      DO 20 J=1,lp13
      DO 20 I=1,lp1
      RCA(I,J)=0.D0
   20 CONTINUE
C
C...  RCA(i,j) = c(k) (la,0;lb,0) * sqrt((2*la+1)(2*lb+1))
C       j identifies la,lb:  ss,ps,pp,ds,dp,dd,... 1,2,3,4,5,6,...
C       i identifies k: i=1 for the lowest k whose c(k) is non-zero,
C                       i=2 for the next   k  "     "    "   ", etc.
C
      RCA(1,1)=1.D0
      RCA(1,2)=1.D0/3.D0
      RCA(1,3)=1.D0/3.D0
      RCA(2,3)=2.D0/15.D0
      RCA(1,4)=1.D0/5.D0
      RCA(1,5)=2.D0/15.D0
      RCA(2,5)=3.D0/35.D0
      RCA(1,6)=1.D0/5.D0
      RCA(2,6)=2.D0/35.D0
      RCA(3,6)=2.D0/35.D0
      RCA(1,7)=1.D0/7.D0
      RCA(1,8)=3.D0/35.D0
      RCA(2,8)=4.D0/63.D0
      RCA(1,9)=3.D0/35.D0
      RCA(2,9)=4.D0/105.D0
      RCA(3,9)=10.D0/231.D0
      RCA(1,10)=1.D0/7.D0
      RCA(2,10)=4.D0/105.D0
      RCA(3,10)=2.D0/77.D0
      RCA(4,10)=100.D0/(13.D0*231.D0)
      RCA(1,11)=1.D0/9.D0
      RCA(1,12)=4.D0/63.D0
      RCA(2,12)=10.D0/198.D0
      RCA(1,13)=2.D0/35.D0
      RCA(2,13)=20.D0/693.D0
      RCA(3,13)=10.D0/286.D0
      RCA(1,14)=4.D0/63.D0
      RCA(2,14)=2.D0/77.D0
      RCA(3,14)=20.D0/1001.D0
      RCA(4,14)=70.D0/2574.D0
      RCA(1,15)=1.D0/9.D0
      RCA(2,15)=20.D0/693.D0
      RCA(3,15)=162.D0/9009.D0
      RCA(4,15)=20.D0/1287.D0
      RCA(5,15)=490.D0/21879.D0
C
      IJ=0
      DO 30 I=1,lp1
        DO 30 J=1,I
          IJ=IJ+1
          KOSUU(IJ)=J
30    CONTINUE
C
      ICOL=0
      DO 40 L=1,lp1
        DO 40 I=1,L
          ICOL=ICOL+1
          IVAL=L-I-2
          DO 40 IROW=1,I
            IVAL=IVAL+2
            NYU(IROW,ICOL)=IVAL
40    CONTINUE
C
      RETURN
      END

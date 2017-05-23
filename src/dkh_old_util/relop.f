************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Bernd Artur Hess                                 *
************************************************************************
      SUBROUTINE RELOP
      IMPLICIT REAL*8(A-H,O-Z)
C
C     SUBROUTINE RELOP INITIALIZES THE COMMON BLOCK USED BY
C     THE RELOP PACKAGE
C     V 1.0 - 12.3.86 - BERND HESS
C
#include "RelLight.fh"
      COMMON /CRELOP/ PI,ZWP,SQPI,VELIT,PREA,CSQ,ZWPH32,FAK(26),
     *                ZWPH12,BCO(210),GA(20),IMAX
      SAVE /CRELOP/
c      WRITE (6,100)
c100   FORMAT(/,
c     X' ****** RELATIVISTIC OPERATORS V 1.0 - BERND HESS ******'
c     X//)
      PI=4.D0*ATAN(1.D0)
      ZWP=2.D0*PI
      ZWPH32=ZWP**1.5D0
      ZWPH12=sqrt(ZWP)
      SQPI=sqrt(PI)
      VELIT=CLight
      PREA=1.D0/(VELIT*VELIT)
      CSQ=VELIT*VELIT
      FAK(1)=1.D0
C     GHALB(1)=SQPI
      DO 1 I=2,26
C     GHALB(I)=GHALB(I-1)*(DBLE(I)-1.5D0)
1     FAK(I)=FAK(I-1)*DBLE(I-1)
C
C     BINOMIALKOEFFIZIENTEN
C
      IMAX=20
      BCO(1)=1.D0
      IBIAS=1
      JBIAS=1
      K=IMAX-1
      DO 53 I=1,K
      ADD=0.D0
      DO 52 J=1,I
      JBIAS=JBIAS+1
      BCO(JBIAS)=ADD+BCO(IBIAS)
      ADD=BCO(IBIAS)
      IBIAS=IBIAS+1
52    CONTINUE
      JBIAS=JBIAS+1
      BCO(JBIAS)=1.D0
53    CONTINUE
C
      DO 10 N=1,20
10    GA(N)=GAM(N-1)
      RETURN
      END
      SUBROUTINE ADDMA(N,S,OVE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S(N),OVE(N)
      DO 1 I=1,N
      OVE(I)=OVE(I)+S(I)
1     CONTINUE
      RETURN
      END
C     COMPILER (XM=3)
*                                                                      *
************************************************************************
*                                                                      *

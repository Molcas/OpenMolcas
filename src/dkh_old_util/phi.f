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
* Copyright (C) 1984,1986, Bernd Artur Hess                            *
************************************************************************
      REAL*8 FUNCTION PHI(M,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CRELOP/ PI,ZWP,SQPI,VELIT,PREA,CSQ,ZWPH32,FAK(26),
     *ZWPH12,BCO(210),GAM(20),IMAX
      IF (MOD(N,2).EQ.1.OR.MOD(M,2).EQ.1) GOTO 10
      PHI=2.D0*GAM(M+1)*GAM(N+1)/GAM(M+N+2)
      RETURN
10    PHI=0.D0
      RETURN
      END
C *** PVP V 1.0 - 19.1.84 - BERND HESS
C *** PVP V 1.1 -  5.6.86 - BERND HESS
      REAL*8 FUNCTION
     *PVP(CHARGE,AL,BE,L1,M1,N1,L2,M2,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CRELOP/ PI,ZWP,SQPI,VELIT,PREA,CSQ,ZWPH32,FAK(26),
     *ZWPH12,BCO(210),GAM(20),IMAX
      INTEGER IS1(3),IS2(3)
C
C    CALCULATE ANGULAR AND RADIAL PART
C
      II=L1+L2
      JJ=M1+M2
      KK=N1+N2
CBS   LAMBDA set to zero, as it is not initialized.   98/02/17
      LAMBDA=0
      IMAX=II+JJ+KK+LAMBDA+3
      IF (IMAX.LE.20) GOTO 2
C
C    ERROR BRANCH: ANGULAR MOMENTUM  > MAXIMUM GIVEN BY ARRAY GAM
C
      WRITE (6,1002) L1,M1,N1,L2,M2,N2,LAMBDA
1002  FORMAT(' ILLEGAL ANGULAR MOMENTUM (PVP)'/,
     *       ' L1,M1,N1,L2,M2,N2,LAMBDA PRINTED'/,1X,7I5)
      Call Abend
C
C    COMPUTE INTEGRAL OVER DERIVATIVE OF THE FUNCTIONS
C
2     IS1(1)=L1
      IS1(2)=M1
      IS1(3)=N1
      IS2(1)=L2
      IS2(2)=M2
      IS2(3)=N2
      SUM=DER(1,IS1,IS2,AL,BE)+DER(2,IS1,IS2,AL,BE)+
     \    DER(3,IS1,IS2,AL,BE)
      SUM=SUM*CHARGE
      PVP=SUM
      RETURN
      END

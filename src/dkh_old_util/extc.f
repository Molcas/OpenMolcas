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
      real*8 FUNCTION EXTC(LAMBDA,AL,BE,L1,M1,N1,L2,M2,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CRELOP/ PI,ZWP,SQPI,VELIT,PREA,CSQ,ZWPH32,FAK(26),
     *ZWPH12,BCO(210),GA(20),IMAX
      INTEGER IS1(3),IS2(3)
C
C    CALCULATE ANGULAR AND RADIAL PART
C
      II=L1+L2
      JJ=M1+M2
      KK=N1+N2
      IMAX=II+JJ+KK+3
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
C
C     NORMALIZATION
C
      II=L1+L1
      JJ=M1+M1
      KK=N1+N1
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      EX=-0.5D0*DBLE(II+JJ+KK+3)
      OV1=0.5D0*ANG*GA(II+JJ+KK+3)*((AL+AL)**EX)
C@    WRITE (6,*) ' EXTC OV1',L1,M1,N1,AL,ANG,sqrt(1/OV1)
      II=L2+L2
      JJ=M2+M2
      KK=N2+N2
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      EX=-0.5D0*DBLE(II+JJ+KK+3)
      OV2=0.5D0*ANG*GA(II+JJ+KK+3)*((BE+BE)**EX)
C@    WRITE (6,*) ' EXTC OV1',L2,M2,N2,BE,ANG,sqrt(1/OV2)
      EXTC=SUM/sqrt(OV1*OV2)
      RETURN
      END

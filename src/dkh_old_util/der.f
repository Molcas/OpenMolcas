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
      REAL*8 FUNCTION DER(IDER,IS1,IS2,AL,BE)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C    CALCULATE INTEGRAL OVER DERIVATIVE OF THE FUNCTIONS
C
      INTEGER IS1(3),IS2(3),I1(2,3),I2(2,3)
      REAL*8 F(2),G(2)
      COMMON /CRELOP/ PI,ZWP,SQPI,VELIT,PREA,CSQ,ZWPH32,FAK(26),
     *ZWPH12,BCO(210),GAM(20),IMAX
      DO 1 I=1,3
      DO 2 J=1,2
      I1(J,I)=IS1(I)
      I2(J,I)=IS2(I)
2     CONTINUE
1     CONTINUE
      I1(1,IDER)=I1(1,IDER)+1
      I1(2,IDER)=I1(2,IDER)-1
      I2(1,IDER)=I2(1,IDER)+1
      I2(2,IDER)=I2(2,IDER)-1
      L1=IS1(IDER)+1
      GOTO (10,11,12,13,14),L1
101   WRITE (6,100) IDER,IS1,IS2,AL,BE
100   FORMAT(' ILLEGAL ANGULAR MOMENTUM (DER)'/,
     \       ' IDER,IS1,IS2,AL,BE PRINTED'/,1X,7I5,3X,2D20.8)
      Call Abend
10    F(1)=-2.D0*AL
      J1=1
      GOTO 19
C
11    F(2)=1.D0
      F(1)=-2.D0*AL
      J1=2
      GOTO 19
C
12    F(2)=2.D0
      F(1)=-2.D0*AL
      J1=2
      GOTO 19
C
13    F(2)=3.D0
      F(1)=-2.D0*AL
      J1=2
      GOTO 19
C
14    F(2)=4.0D0
      F(1)=-2.0D0*AL
      J1=2
C
19    L2=IS2(IDER)+1
      GOTO (20,21,22,23,24),L2
      GOTO 101
C
20    G(1)=-2.D0*BE
      J2=1
      GOTO 29
C
21    G(2)=1.D0
      G(1)=-2.D0*BE
      J2=2
      GOTO 29
C
22    G(2)=2.D0
      G(1)=-2.D0*BE
      J2=2
      GOTO 29
C
23    G(2)=3.D0
      G(1)=-2.D0*BE
      J2=2
      GOTO 29
C
24    G(2)=4.0D0
      G(1)=-2.0D0*BE
      J2=2
C
29    SUM=0.D0
      DO 30 I=1,J1
      DO 31 J=1,J2
      II=I1(I,1)+I2(J,1)
      JJ=I1(I,2)+I2(J,2)
      KK=I1(I,3)+I2(J,3)
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      IF (ANG.EQ.0.D0) GOTO 31
      EX=-DBLE(II+JJ+KK+2)*0.5D0
      SUM=SUM+F(I)*G(J)*0.5D0*ANG*GAM(II+JJ+KK+2)*
     \        (AL+BE)**EX
31    CONTINUE
30    CONTINUE
      DER=SUM
      RETURN
      END

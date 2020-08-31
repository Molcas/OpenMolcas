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
      REAL*8 FUNCTION GAM(M)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CRELOP/ PI,ZWP,SQPI,VELIT,PREA,CSQ,ZWPH32,FAK(26),
     *ZWPH12,BCO(210),GA(20),IMAX
*rl   COMMON /CRELOP/ PI,ZWP,SQPI
      SAVE /CRELOP/
      IF (MOD(M,2).EQ.0) GOTO 10
      MA=(M+1)/2
      G=1.D0
      IF (MA.EQ.1) GOTO 11
      DO 1 I=2,MA
      G=G*DBLE(I-1)
1     CONTINUE
      GAM=G
      RETURN
10    MA=M
      G=SQPI
      IF (MA.EQ.0) GOTO 11
      DO 2 I=1,MA,2
      G=G*0.5D0*DBLE(I)
2     CONTINUE
11    GAM=G
      RETURN
      END

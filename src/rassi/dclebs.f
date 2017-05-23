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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
      REAL*8 FUNCTION DCLEBS(XJ1,XJ2,XJ3,XM1,XM2,XM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXJ=10, MAXF=3*MAXJ+1)
      DIMENSION DFACT(0:MAXF)
      SAVE ICALL,DFACT
      DATA ICALL / 0 /
C
C DCLEBS: REAL*8 Clebsch-Gordan coefficients. From a
C modification of Racah''s formula. Coded: Malmqvist 1998.
C
C Note carefully: The input values XJ1..XM3 are REAL*8, not
C integers. Half-integer spins are allowed. Half-integers
C are assumed exactly represented.

      IF(ICALL.EQ.0) THEN
        ICALL=ICALL+1
        DF=1.0D0
        DFACT(0)=DF
        DO I=1,MAXF
          DF=DBLE(I)*DF
          DFACT(I)=DF
        END DO
      END IF

      DCLEBS=0.0D0

      XJSUM=XJ1+XJ2+XJ3
      JSUM=NINT(XJSUM)
      IF(XJSUM.NE.DBLE(JSUM)) RETURN
      IF(XM1+XM2.NE.XM3) RETURN

      IA1=NINT(XJ1+XM1)
      IF(IA1.LT.0) RETURN
      IB1=NINT(XJ1-XM1)
      IF(IB1.LT.0) RETURN
      IA2=NINT(XJ2+XM2)
      IF(IA2.LT.0) RETURN
      IB2=NINT(XJ2-XM2)
      IF(IB2.LT.0) RETURN
      IA3=NINT(XJ3-XM3)
      IF(IA3.LT.0) RETURN
      IB3=NINT(XJ3+XM3)
      IF(IB3.LT.0) RETURN
      IF(JSUM-IA1-IB1.LT.0) RETURN
      IF(JSUM-IA2-IB2.LT.0) RETURN
      IF(JSUM-IA3-IB3.LT.0) RETURN

      PRE2=DBLE(1+IA3+IB3)*DFACT(JSUM-IA1-IB1)
     &          *DFACT(JSUM-IA2-IB2)*DFACT(JSUM-IA3-IB3)
     &          *DFACT(IA1)*DFACT(IA2)*DFACT(IA3)
     &          *DFACT(IB1)*DFACT(IB2)*DFACT(IB3)
     &          /DFACT(JSUM+1)
      PRE=SQRT(PRE2)

      IY0=(JSUM-IA3-IB3)
      IX1=(IA2+IB1-JSUM)+IB2
      IX2=(IA2+IB1-JSUM)+IA1
      IX=MAX(0,IX1,IX2)
      IY=MIN(IY0,IB1,IA2)

      SUMMA=0.0D0
      DO I=IX,IY
        DEN=DFACT(I)*DFACT(I-IX1)*DFACT(I-IX2)*DFACT(IY0-I)*
     &      DFACT(IB1-I)*DFACT(IA2-I)
        TERM=1.0D0/DEN
        SUMMA=SUMMA+DBLE((-1)**I)*TERM
      END DO

      DCLEBS=PRE*SUMMA

      RETURN
      END

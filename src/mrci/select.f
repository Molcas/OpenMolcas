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
      SUBROUTINE CI_SELECT_MRCI(NREF,AREF,PLEN,NSEL,CISEL,NRROOT,IROOT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AREF(NREF,NREF),CISEL(NREF,*),IROOT(NRROOT)
      DIMENSION PLEN(NREF)
      IF(NSEL.EQ.0) RETURN
C SELECTION BY PROJECTION ONTO SPACE SPANNED BY CISEL VECTORS. IROOT()
C IS SET TO SELECT THE NRROOT VECTORS WITH MAX PROJECTED LENGTH.
      DO 140 J=1,NREF
        SUM=0.0D00
        DO 130 ISEL=1,NSEL
          SUM1=0.0D00
          DO 120 I=1,NREF
            SUM1=SUM1+AREF(I,J)*CISEL(I,ISEL)
120       CONTINUE
          SUM=SUM+SUM1**2
130     CONTINUE
        PLEN(J)=SUM+J*1.0D-12
140   CONTINUE
C SELECT BY MAGNITUDE OF PLEN:
      DO 160 J=1,NRROOT
        PMAX=PLEN(1)
        JMAX=1
        DO 150 JJ=2,NREF
          IF(PMAX.GE.PLEN(JJ)) GOTO 150
          PMAX=PLEN(JJ)
          JMAX=JJ
150     CONTINUE
        PLEN(JMAX)=-PMAX
160   CONTINUE
      I=0
      DO 170 IR=1,NREF
        PL=PLEN(IR)
        IF(PL.LT.0.0D00) THEN
          I=I+1
          IROOT(I)=IR
          PL=-PL
        END IF
        PLEN(IR)=PL-IR*1.0D-12
170   CONTINUE
      RETURN
      END

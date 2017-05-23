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
      SUBROUTINE WRH(LU,NSYM,NBAS,NORB,CMO,OCC,LOCC,TITLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NBAS(NSYM),NORB(NSYM),CMO(*),OCC(*)
      CHARACTER*(*) TITLE
      CHARACTER FMT*40
      FMT='(4E20.12)'
*     REWIND (LU)
      KCMO  = 0
      NDIV  = 4
      IF(TITLE(1:1).NE.'*') TITLE='*'//TITLE
      If (locc.ne.2) Then
      DO 100 ISYM=1,NSYM
         DO 110 IORB=1,NORB(ISYM)
            WRITE(LU,'(A,I5)') '* Column    ',IORB
            DO 111 IBAS=1,NBAS(ISYM),NDIV
               WRITE(LU,FMT) (CMO(I+KCMO),I=IBAS,MIN(IBAS+3,NBAS(ISYM)))
111         CONTINUE
            KCMO=KCMO+NBAS(ISYM)
110      CONTINUE
100   CONTINUE
      End If
      IF(LOCC.EQ.0) RETURN
      WRITE(LU,'(A)') Title
      KOCC=0
      DO 200 ISYM=1,NSYM
         DO 210 IORB=1,NORB(ISYM),NDIV
            WRITE(LU,FMT) (OCC(I+KOCC),I=IORB,MIN(IORB+3,NORB(ISYM)))
210      CONTINUE
         KOCC=KOCC+NORB(ISYM)
200   CONTINUE
      RETURN
      END

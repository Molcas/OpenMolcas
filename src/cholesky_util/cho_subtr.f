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
      SUBROUTINE CHO_SUBTR(XINT,WRK,LWRK,ISYM)
C
C     Purpose: driver for subtracting contributions from previous vectors
C              from the qualified integrals (in XINT).
C
#include "implicit.fh"
      DIMENSION XINT(*), WRK(LWRK)
#include "cholesky.fh"

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_SUBTR')

      LOGICAL LOCDBG, FXDMEM
      PARAMETER (LOCDBG = .FALSE.)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      SubRoutine Cho_VecBuf_Subtr(xInt,Wrk,lWrk,iSym,DoTime,DoStat)
      Real*8, Target::  xInt(*), Wrk(lWrk)
      Logical DoTime, DoStat
      End SubRoutine Cho_VecBuf_Subtr
      End Interface
*                                                                      *
************************************************************************
*                                                                      *

C     Return if nothing to do.
C     ------------------------

      IF (NUMCHO(ISYM) .LT. 1) THEN ! no prev. vectors.
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': nothing done because NUMCHO = ',
     &                     NUMCHO(ISYM),' (sym. ',ISYM,')'
         END IF
         return
      ELSE IF (NNBSTR(ISYM,2) .LT. 1) THEN ! nothing to do (this sym.)
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': nothing done because NNBSTR = ',
     &                     NNBSTR(ISYM,2),' (sym. ',ISYM,')'
         END IF
         return
      ELSE IF (NQUAL(ISYM) .LT. 1) THEN ! no qualifieds in this sym.
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': nothing done because NQUAL  = ',
     &                     NQUAL(ISYM),' (sym. ',ISYM,')'
         END IF
         return
      END IF

C     Debug: read original diagonal and check that these elements are
C            included in the integrals
C     ---------------------------------------------------------------

      IF (CHO_DIACHK .OR. LOCDBG) THEN
         KDIAG = 1
         KEND  = KDIAG + NNBSTRT(1)
         LWRK  = LWRK  - KEND + 1
         IF (LWRK .LT. 0) THEN
            WRITE(LUPRI,*) SECNAM,': diagonal/integral check skipped ',
     &                     'due to insufficient memory'
         ELSE
            TOL  = TOL_DIACHK
            NERR = 0
            CALL CHO_CHKINTO(XINT,WRK(KDIAG),ISYM,NERR,TOL,.TRUE.)
            IF (NERR .NE. 0) THEN
               WRITE(LUPRI,*) SECNAM,': ',NERR,' diagonal errors found!'
               WRITE(LUPRI,*) '          #tests: ',NQUAL(ISYM)
c              WRITE(LUPRI,*) '          Printing integrals:'
c              CALL CHO_OUTPUT(XINT,1,NNBSTR(ISYM,2),1,NQUAL(ISYM),
c    &                         NNBSTR(ISYM,2),NQUAL(ISYM),1,LUPRI)
               CALL CHO_QUIT('Diagonal errors in '//SECNAM,104)
            ELSE
               WRITE(LUPRI,*) SECNAM,': comparison of qual. integrals ',
     &                     'and original diagonal: no errors !'
            END IF
         END IF
      END IF

C     Subtract contributions for vectors in buffer.
C     (Returns immediately if nothing to do.)
C     ---------------------------------------------

      CALL CHO_VECBUF_SUBTR(XINT,WRK,LWRK,ISYM,.TRUE.,.TRUE.)

C     Subtract contributions for vectors on disk.
C     -------------------------------------------

      IF (CHO_IOVEC.EQ.3 .OR. CHO_IOVEC.EQ.4) THEN
         FXDMEM = CHO_IOVEC .EQ. 4
         CALL CHO_SUBTR1(XINT,WRK,LWRK,ISYM,FXDMEM)
      ELSE
         CALL CHO_SUBTR0(XINT,WRK,LWRK,ISYM)
      END IF
      END

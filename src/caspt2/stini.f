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
      SUBROUTINE STINI
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
      CHARACTER(50)  STLNE2
C     timers
      REAL*8 CPU0,CPU1,CPU,
     &       TIO0,TIO1,TIO
C     indices
      INTEGER I,J,IFTEST
      ! INTEGER IDCI

      CALL QENTER('STINI')
      Write(STLNE2,'(A,I4)')
     &                ' Compute H0 matrices for state ',MSTATE(JSTATE)
      Call StatusLine('CASPT2: ',trim(STLNE2))
      IF(IPRGLB.GE.USUAL) THEN
        WRITE(6,'(20A4)')('****',I=1,20)
        WRITE(6,'(A,I4)')
     &   ' Compute H0 matrices for state ',MSTATE(JSTATE)
        WRITE(6,'(20A4)')('----',I=1,20)
        CALL XFlush(6)
      END IF

* WITH NEW CMOS, TRANSFORM ONE- AND TWO-ELECTRON INTEGRALS.
* LUSOLV, LUSBT and LUDMAT will be reused in TRACTL.

* WITH NEW CI, RECOMPUTE 1- AND 2-DENSITY FOR ROOT STATE JSTATE:
    !   CALL GETMEM('LCI','ALLO','REAL',LCI,NCONF)
    !   IF(.NOT.DoCumulant.AND.ISCF.EQ.0) THEN
    !     IDCI=IDTCEX
    !     DO J=1,JSTATE-1
    !       CALL DDAFILE(LUCIEX,0,WORK(LCI),NCONF,IDCI)
    !     END DO
    !     CALL DDAFILE(LUCIEX,2,WORK(LCI),NCONF,IDCI)
    !     IF(IPRGLB.GE.VERBOSE) THEN
    !       WRITE(6,*)
    !       IF(NSTATE.GT.1) THEN
    !         WRITE(6,'(A,I4)')
    !  &      ' With new orbitals, the CI array of state ',MSTATE(JSTATE)
    !       ELSE
    !         WRITE(6,*)' With new orbitals, the CI array is:'
    !       END IF
    !       CALL PRWF_CP2(LSYM,NCONF,WORK(LCI),CITHR)
    !     ENDIF
    !   ELSE
    !     WORK(LCI)=1.0D0
    !   END IF

      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' STINI calling POLY3...'
      END IF
      CALL TIMING(CPU0,CPU,TIO0,TIO)
      CALL POLY3(1)
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUFG3=CPU1-CPU0
      TIOFG3=TIO1-TIO0
      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' STINI back from POLY3.'
      END IF

      ! IF(IPRGLB.GE.DEBUG) THEN
      !   WRITE(6,*)' STINI calling POLY2...'
      ! END IF
      ! CALL POLY2(WORK(LCI))
* GETDPREF: Restructure GAMMA1 and GAMMA2, as DREF and PREF arrays.
      CALL GETDPREF(WORK(LDREF),WORK(LPREF))
      ! IF(IPRGLB.GE.DEBUG) THEN
      !   WRITE(6,*)' STINI back from POLY2.'
      ! END IF

      IFTEST = 0
      IF ( IFTEST.NE.0 ) THEN
        WRITE(6,*)' DREF for state nr. ',MSTATE(JSTATE)
        DO I=1,NASHT
          WRITE(6,'(1x,14f10.6)')(WORK(LDREF+(I*(I-1))/2+J-1),J=1,I)
        END DO
        WRITE(6,*)
      END IF


      EREF=REFENE(JSTATE)
* With new DREF, recompute EASUM:
      EASUM=0.0D0
      DO I=1,NASHT
        EASUM=EASUM+EPSA(I)*WORK(LDREF-1+(I*(I+1))/2)
      END DO


      IF(IPRGLB.GE.USUAL) THEN
       WRITE(6,'(20A4)')('----',I=1,20)
       WRITE(6,'(A)')' H0 matrices have been computed.'
       WRITE(6,*)
      ENDIF
      CALL QEXIT('STINI')

      RETURN
      END

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
* Copyright (C) 1991,1997, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE INTIM
*
* Interface to external integrals
*
* If NOINT .ne. 0, only pointers are constructed
* Jeppe Olsen, Winter of 1991
*
* Version : Fall 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "crun.fh"
#include "glbbas.fh"
#include "clunit.fh"
#include "lucinp.fh"
#include "csm.fh"
#include "orbinp.fh"
#include "cintfo.fh"
#include "oper.fh"
#include "cecore.fh"
*
*. : Pointers for symmetry blocks of integrals
*
      CALL INTPNT(iWORK(KPINT1),iWORK(KLSM1),
     &            iWORK(KPINT2),iWORK(KLSM2))
*
*. Pointer for orbital indeces for symmetry blocked matrices
      CALL ORBINH1(iWORK(KINH1),iWORK(KINH1_NOCCSYM),NTOOBS,NTOOB,NSMOB)
*
*. Change one-electron integrals to inactive fock matrix
      IF(NOINT.EQ.0) THEN
C?      WRITE(6,*) ' INTIM : IUSE_PH', IUSE_PH
        CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
c        IF(IUSE_PH.EQ.1) THEN
c           IF (ENVIRO(1:6) .EQ. 'RASSCF') THEN
c              ECORE_HEX = 0.0D0
cc           ELSE
cc              CALL FI(WORK(KINT1),ECORE_HEX,1)
c           END IF
c        ELSE
           ECORE_HEX = 0.0D0
c        END IF
      END IF
      ECORE_ORIG = ECORE
      ECORE = ECORE + ECORE_HEX
c      IF (ENVIRO(1:6) .NE. 'RASSCF')
c     &   WRITE(6,*) ' Updated core energy ',ECORE
*
C?    WRITE(6,*) ' IDMPIN ', IDMPIN
c      IF (IDMPIN.EQ.1 ) THEN
c        WRITE(6,*)
c     &   ' Integrals written formatted (E22.15) on unit 90'
c        LU90 = 90
c        REWIND LU90
c*.1 : One-electron integrals
c        WRITE(LU90,'(E22.15)')
c     &   (WORK(KINT1O-1+INT1),INT1=1,NINT1)
c*.2 : Two-electron integrals
c        WRITE(LU90,'(E22.15)')
c     &   (WORK(KINT2-1+INT2),INT2=1,NINT2)
c*.3. Core energy
c        WRITE(LU90,'(E22.15)')ECORE_ORIG
c*.4  Rewind to empty buffer
c        REWIND LU90
c*.   Symmetry info etc two LU91
c        LU91 = 91
c        CALL DUMP_1EL_INFO(LU91)
c      END IF
*
c      IF (ENVIRO(1:6) .NE. 'RASSCF' .OR. IPRNT .GE. 20)
c     &       WRITE(6,*) ' INTIM : First integrals in WORK(KINT1) '
      LLL = MIN(10,NINT1)
      LLL = NINT1
c      IF (ENVIRO(1:6) .NE. 'RASSCF' .OR. IPRNT .GE. 20) THEN
c         WRITE(6,*) ' NINT1 = ',NINT1
c         CALL WRTMAT(WORK(KINT1),1,LLL,1,LLL)
c         WRITE(6,*) ' INTIM : First integrals in WORK(KINT2) '
c      ENDIF
      LLL = MIN(10,NINT2)
      LLL = NINT2
c      IF (ENVIRO(1:6) .NE. 'RASSCF' .OR. IPRNT .GE. 20) THEN
c         WRITE(6,*) ' NINT2 = ',NINT2
c         CALL WRTMAT(WORK(KINT2),1,LLL,1,LLL)
c      ENDIF

C!    stop ' Jeppe forced my to stop in INTIM '
      RETURN
      END

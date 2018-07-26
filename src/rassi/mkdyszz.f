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

*****************************************************************
*  SUBROUTINE MKDYSZZ
*  PURPOSE: CALCULATE DYSON ORBITAL COEFFICIENTS FOR CI EXPANSIONS IN
*  BASIS FUNCTION BASE BASE Z,
*  IN ANALOGUE TO MKTDZZ FOR TRANSITION DENSITY MATRIX.
*****************************************************************

      SUBROUTINE MKDYSZZ(CMOA,DYSAB,DYSZZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DYSAB(*),DYSZZ(*)
      DIMENSION CMOA(NCMO)
      DIMENSION ISTCMO(8)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"

C *** Symmetry is likely not handled correctly, since the effect
C *** of the annihilated electron is not accountd for.

      IST=1
      ISTTA=1
      ISTCA=1
      ISTTZ=1
      DO 20 ISY1=1,NSYM
        NO1=NOSH(ISY1)
        NB1=NBASF(ISY1)
        IF(NB1.EQ.0) GOTO 15
        CALL FZERO(DYSZZ(ISTTZ),NB1)
        CALL DGEMM_('N','T',1,NB1,NO1,1.0D0,
     &              DYSAB(ISTTA),1,CMOA(ISTCA),NB1,
     &       0.0D0,  DYSZZ(ISTTZ),1)
        ISTTA=ISTTA+NO1
15      CONTINUE
        ISTCA=ISTCA+NB1*NO1
        ISTTZ=ISTTZ+NB1
20    CONTINUE

C *** Helpful debugging printouts

!      WRITE(*,*)
!      WRITE(*,*)'--- MKDYSZZ ---'
!      WRITE(*,*)'NCMO:',NCMO
!      WRITE(*,*)'DYSAB:'
!      DO I=1,NO1
!       WRITE(*,'(F6.2)', advance="no"),DYSAB(I)
!      END DO
!      WRITE(*,*)
!      WRITE(*,*)'CMOA:'
!      DO I=1,NO1
!       DO J=1,NB1
!        IST=NB1*(I-1)+J
!        WRITE(*,'(F10.2)', advance="no"),CMOA(IST)
!       END DO
!       WRITE(*,*)
!      END DO
!      WRITE(*,*)'DYSZZ:'
!      DO I=1,NB1
!       WRITE(*,'(F6.2)', advance="no"),DYSZZ(I)
!      END DO

      RETURN
      END



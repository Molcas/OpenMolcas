!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE CIDIA_CI_UTIL(NORB,NCONF,IREFSM,                       &
     &                 CSFDIA,G,TUVX,LUDAVID)
!
!     PURPOSE: - COMPUTE DIAGONAL ELEMENTS OF THE CI-MATRIX
!                THE DETERMINANT BASIS
!              - TRANSLATE FORM DET => CSF BASIS
!
!     CALLING PARAMETERS:
!     NORB  :  NO. OF ACTIVE ORBITALS
!     NCONF :  NO. OF CSF
!     IREFSM:  REFERENCE SYMMETRY
!     CSFDIA:  DIAGONAL OF CI MATRIX IN CSF BASIS
!     G     :  MODIFIED ONE ELECTRON HAMILTONIAN INCLUDING CORE ELECTR.
!     TUVX  :  TWO-ELECTRON INTEGRALS
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
#include "ciinfo.fh"
#include "spinfo.fh"
#include "detbas.fh"
#include "csfbas.fh"
#include "strnum.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "rasscf_lucia.fh"
#include "output_ras.fh"
!
      DIMENSION CSFDIA(*)
      DIMENSION G(*)
      DIMENSION TUVX(*)
      DIMENSION Dummy(1)
!
      Call Timing(Tissot_1,Swatch,Swatch,Swatch)
      IPRLEV=IPRLOC(3)
!
!     ALLOCATE LOCAL MEMORY
!
      CALL GETMEM('XA','ALLO','REAL',KX,NORB)
      CALL GETMEM('SCR','ALLO','REAL',KSCR,2*NORB)
      CALL GETMEM('H1DIA','ALLO','REAL',KH1DIA,NORB)
!
!     SELECT DIAGONAL ONEBODY INTEGRALS
!
      II=0
      DO 100 I = 1,NORB
        II=II+I
        Work(KH1DIA-1+I) = G(II)
100   CONTINUE
!
!     COMPUTE CI DIAGONAL IN DETERMINANT BASIS
!
!
      CALL Lucia_Util('Diag',iDummy,iDummy,Dummy)
!
      CALL GETMEM('DETDIA','ALLO','REAL',KDDIA,NDET)
      CALL get_diag(work(kddia),ndet)
!
!     TRANSFORM CI DIAGONAL FROM DET TO CSF BASIS
!
      IPRINT=0
      IF(IPRLEV.EQ.INSANE) IPRINT=40
      CALL CSDIAG_CI_UTIL(CSFDIA,Work(KDDIA),NCNFTP(1,IREFSM),          &
     &     NTYP,iWork(KICTS(1)),NDTFTP,NCSFTP,IPRINT)
      One       = 1.0d0
      eCore_Hex = Get_eCore()
      CALL DAXPY_(NCONF,eCore_Hex,[One],0,CSFDIA,1)
!
!     DEALLOCATE LOCAL MEMORY
!
      CALL GETMEM('XA','FREE','REAL',KX,NORB)
      CALL GETMEM('SCR','FREE','REAL',KSCR,2*NORB)
      CALL GETMEM('H1DIA','FREE','REAL',KH1DIA,NORB)
      CALL GETMEM('DETDIA','FREE','REAL',KDDIA,NDET)
!
!     PRINT CI-DIAGONAL
!
      IF ( IPRLEV.GE.DEBUG ) THEN
        IPRL=NCONF
        IF( IPRLEV.LT.INSANE ) IPRL=MIN(IPRL,200)
        Call dVcPrt('CI-DIAGONAL (max.200 elemwnts)',' ',CSFDIA,NCONF)
      END IF
!
!     SAVE THE CI_DIAGONAL ON TAPE
!
      Call Save_H_diag(nConf,CSFDIA,LUDAVID)
!
      Call Timing(Tissot_2,Swatch,Swatch,Swatch)
      Tissot_2 = Tissot_2 - Tissot_1
      Tissot_3 = Tissot_3 + Tissot_2
!
      RETURN
! Avoid unused argument warnings
      If (.False.) Call Unused_real_array(TUVX)
      END

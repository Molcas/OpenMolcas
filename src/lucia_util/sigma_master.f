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
      SUBROUTINE sigma_master()
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS
*
* Controls the calculation of the sigma vector, when Lucia is called
* from Molcas Rasscf.
*
      implicit real*8 (a-h,o-z)
#include "mxpdim.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "clunit.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "cecore.fh"
#include "crun.fh"
#include "rasscf_lucia.fh"
      Integer, Allocatable:: lVec(:)
*
* Put CI-vector from RASSCF on luc and get h0 from Molcas enviroment.
*
      IF (INI_H0 .EQ. 0) THEN
         ECORE = ECORE_ORIG
      ENDIF
      INI_H0 = 0
      NDIM = NTOOB**2
      CALL COPVEC(WORK(KINT1O_POINTER),WORK(KINT1_POINTER),NDIM)
      ECORE_ORIG = ECORE
c      IF (IUSE_PH .EQ. 1) THEN
c         CALL FI(WORK(KINT1_POINTER),ECORE_HEX,1)
c         ECORE = ECORE + ECORE_HEX
c      END IF
      call mma_allocate(lVec,MXNTTS,Label='lVec')
      CALL CPCIVC(LUC, MXNTTS, IREFSM, 1, lVec)
      call mma_deallocate(lVec)
*
* Calculate the sigma vector:
*
      Call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
      CALL MV7(WORK(KCI_POINTER), WORK(KSIGMA_POINTER), LUC, LUSC34)
      Call mma_deallocate(VEC3)
*
* Export lusc34 to RASSCF
*
      call mma_allocate(lVec,MXNTTS,Label='lVec')
      CALL CPCIVC(LUSC34, MXNTTS, IREFSM, 2, lVec)
      call mma_deallocate(lVec)
*
      RETURN
      END
******************************
*                            *
*   Now for the CASVB case   *
*                            *
******************************
      SUBROUTINE SIGMA_MASTER_CVB(IREFSM_CASVB)
      use GLBBAS
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mxpdim.fh"
#include "cands.fh"
!#include "cicisp.fh"
      COMMON/CICISP/IDUMMY,NICISP,
     &              NELCI(MXPICI),
     &              XISPSM(MXPCSM,MXPICI),
     &              ISMOST(MXPCSM,MXPCSM),MXSB,MXSOOB,
     &              NBLKIC(MXPCSM,MXPICI),LCOLIC(MXPCSM,MXPICI),
     &              MXNTTS,MXSOOB_AS
#include "cstate.fh"
#include "clunit.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "cecore.fh"
#include "crun.fh"
#include "rasscf_lucia.fh"
      Integer, Allocatable:: lVec(:)
*
* Set ICSM and ISSM (from cands.fh) to the correct symmetry for this call
*
      ICSM  = IREFSM_CASVB
      ISSM  = IREFSM_CASVB
*
* Get h0 from Molcas enviroment.
*
      IF (INI_H0 .EQ. 0) THEN
         ECORE = ECORE_ORIG
      ENDIF
      INI_H0 = 0
      NDIM = NTOOB**2
      CALL COPVEC(WORK(KINT1O_POINTER),WORK(KINT1_POINTER),NDIM)
      ECORE_ORIG = ECORE
c      IF (IUSE_PH .EQ. 1) THEN
c         CALL FI(WORK(KINT1_POINTER),ECORE_HEX,1)
c         ECORE = ECORE + ECORE_HEX
c      END IF
*
* Write CI-vector to disc
*
      call mma_allocate(lVec,MXNTTS,Label='lVec')
      CALL CPCIVC(LUC, MXNTTS, ISSM, 1,lVec)
      call mma_deallocate(lVec)
*
* Calculate the sigma vector:
*
      CALL DIAG_MASTER()
      Call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
      CALL MV7(WORK(KCI_POINTER), WORK(KSIGMA_POINTER), LUC, LUSC34)
      Call mma_deallocate(VEC3)
*
* Export lusc34 to RASSCF
*
      ISIGMA_ON_DISK = 1
*
* Set ICSM and ISSM (from cands.fh) back to IREFSM
*
      ICSM  = IREFSM
      ISSM  = IREFSM
*
      RETURN
      END

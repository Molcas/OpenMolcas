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
      SUBROUTINE sigma_master(CIVEC,nCIVEC)
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS
      use rasscf_lucia, only: INI_H0, KVEC3_LENGTH
*
* Controls the calculation of the sigma vector, when Lucia is called
* from Molcas Rasscf.
*
      implicit real*8 (a-h,o-z)
#include "mxpdim.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "clunit.fh"
#include "orbinp.fh"
#include "cecore.fh"
#include "crun.fh"
#include "spinfo_lucia.fh"
      Integer nCIVEC
      Real*8 CIVEC(nCIVEC)

      Integer, Allocatable:: lVec(:)
      Integer nSD
*
* Put CI-vector from RASSCF on luc and get h0 from Molcas enviroment.
*
      nSD = NSD_PER_SYM(IREFSM)
      IF (INI_H0 .EQ. 0) THEN
         ECORE = ECORE_ORIG
      ENDIF
      INI_H0 = 0
      INT1(:)=INT1O(:)
      ECORE_ORIG = ECORE
c      IF (IUSE_PH .EQ. 1) THEN
c         CALL FI(INT1,ECORE_HEX,1)
c         ECORE = ECORE + ECORE_HEX
c      END IF
      call mma_allocate(lVec,MXNTTS,Label='lVec')
      CALL CPCIVC(CIVEC,nSD,LUC, MXNTTS, IREFSM, 1, lVec)
      call mma_deallocate(lVec)
*
* Calculate the sigma vector:
*
      Call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
      CALL MV7(CI_Vec, SIGMA_Vec, LUC, LUSC34)
      Call mma_deallocate(VEC3)
*
* Export lusc34 to RASSCF
*
      call mma_allocate(lVec,MXNTTS,Label='lVec')
      CALL CPCIVC(CIVEC,nSD,LUSC34, MXNTTS, IREFSM, 2, lVec)
      call mma_deallocate(lVec)
*
      END
******************************
*                            *
*   Now for the CASVB case   *
*                            *
******************************
      SUBROUTINE SIGMA_MASTER_CVB(IREFSM_CASVB)
      use GLBBAS
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_lucia
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
#include "orbinp.fh"
#include "cecore.fh"
#include "crun.fh"
      Integer, Allocatable:: lVec(:)
      Integer nSD

      nSD=Size(C_Pointer)
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
      INT1(:)=INT1O(:)
      ECORE_ORIG = ECORE
c      IF (IUSE_PH .EQ. 1) THEN
c         CALL FI(INT1,ECORE_HEX,1)
c         ECORE = ECORE + ECORE_HEX
c      END IF
*
* Write CI-vector to disc
*
      call mma_allocate(lVec,MXNTTS,Label='lVec')
      CALL CPCIVC(C_Pointer,nSD,LUC, MXNTTS, ISSM, 1,lVec)
      call mma_deallocate(lVec)
*
* Calculate the sigma vector:
*
      CALL DIAG_MASTER()
      Call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
      CALL MV7(CI_VEC, SIGMA_Vec, LUC, LUSC34)
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
      END

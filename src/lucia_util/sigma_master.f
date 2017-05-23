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
      SUBROUTINE sigma_master
*
* Controls the calculation of the sigma vector, when Lucia is called
* from Molcas Rasscf.
*
      implicit real*8 (a-h,o-z)
#include "mxpdim.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "clunit.fh"
#include "glbbas.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "cecore.fh"
#include "crun.fh"
#include "rasscf_lucia.fh"
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
      call GetMem('lvec','Allo','inte',ivlrec,MXNTTS)
      CALL CPCIVC(LUC, MXNTTS, IREFSM, 1, iWork(ivlrec))
      call GetMem('lvec','Free','inte',ivlrec,MXNTTS)
*
* Calculate the sigma vector:
*
      CALL GETMEM('KC2   ','ALLO','REAL',KVEC3,KVEC3_LENGTH)
      CALL MV7(WORK(KCI_POINTER), WORK(KSIGMA_POINTER), LUC, LUSC34)
      CALL GETMEM('KC2   ','FREE','REAL',KVEC3,KVEC3_LENGTH)
*
* Export lusc34 to RASSCF
*
      call GetMem('lvec','Allo','inte',ivlrec,MXNTTS)

      CALL CPCIVC(LUSC34, MXNTTS, IREFSM, 2, iWork(ivlrec))
      call GetMem('lvec','Free','inte',ivlrec,MXNTTS)
*
      RETURN
      END
******************************
*                            *
*   Now for the CASVB case   *
*                            *
******************************
      SUBROUTINE SIGMA_MASTER_CVB(IREFSM_CASVB)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mxpdim.fh"
#include "cands.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "clunit.fh"
#include "glbbas.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "cecore.fh"
#include "crun.fh"
#include "rasscf_lucia.fh"
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
      call GetMem('lvec','Allo','inte',ivlrec,MXNTTS)

      CALL CPCIVC(LUC, MXNTTS, ISSM, 1,iWork(ivlrec))
      call GetMem('lvec','Free','inte',ivlrec,MXNTTS)
*
* Calculate the sigma vector:
*
      CALL DIAG_MASTER
      CALL GETMEM('KC2   ','ALLO','REAL',KVEC3,KVEC3_LENGTH)
      CALL MV7(WORK(KCI_POINTER), WORK(KSIGMA_POINTER), LUC, LUSC34)
      CALL GETMEM('KC2   ','FREE','REAL',KVEC3,KVEC3_LENGTH)
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

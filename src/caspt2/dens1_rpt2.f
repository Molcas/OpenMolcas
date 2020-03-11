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
* Copyright (C) 2006, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2006  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE DENS1_RPT2 (CI,SGM1,G1)
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

#include "para_info.fh"
      LOGICAL RSV_TSK

      REAL*8 CI(MXCI),SGM1(MXCI)
      REAL*8 G1(NLEV,NLEV)
#ifdef _ENABLE_CHEMPS2_DMRG_
      REAL*8 G2(NLEV,NLEV,NLEV,NLEV)
#endif

      REAL*8 GTU

      INTEGER ID
      INTEGER IST,ISU,ISTU
      INTEGER IT,IU,LT,LU,LTU

      INTEGER ITASK,LTASK,LTASK2T,LTASK2U,NTASKS

      INTEGER ISSG,NSGM

      REAL*8, EXTERNAL :: DDOT_,DNRM2_

* Purpose: Compute the 1-electron density matrix array G1.

      CALL QENTER('DENS1_RPT2')

      CALL DCOPY_(NG1,[0.0D0],0,G1,1)

* For the special cases, there is no actual CI-routines involved:
* Special code for hi-spin case:
      IF(ISCF.EQ.2) THEN
        DO IT=1,NASHT
          G1(IT,IT)=1.0D00
        END DO
        GOTO 99
      END IF
* Special code for closed-shell:
      IF(ISCF.EQ.1 .AND. NACTEL.GT.0) THEN
        DO IT=1,NASHT
          G1(IT,IT)=2.0D00
        END DO
        GOTO 99
      END IF

* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level inices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in pt2_guga.fh

* SVC20100311: set up a task table with LT,LU
* SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
      nTasks=(nLev**2+nLev)/2

      CALL GETMEM ('Tasks','ALLO','INTE',lTask,2*nTasks)
      lTask2T=lTask
      lTask2U=lTask+nTasks

      iTask=0
      DO LT=1,nLev
        DO LU=1,LT
          iTask=iTask+1
          iWork(lTask2T+iTask-1)=LT
          iWork(lTask2U+iTask-1)=LU
        ENDDO
      ENDDO
      IF (iTask.NE.nTasks) WRITE(6,*) "ERROR nTasks"

      Call Init_Tsk(ID, nTasks)

* SVC20100311: BEGIN SEPARATE TASK EXECUTION
 500  If (.NOT.Rsv_Tsk (ID,iTask)) GOTO 501

* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
      LT=iWork(lTask2T+iTask-1)
        IST=ISM(LT)
        IT=L2ACT(LT)
        LU=iWork(lTask2U+iTask-1)
          LTU=iTask
          ISU=ISM(LU)
          IU=L2ACT(LU)
          ISTU=MUL(IST,ISU)
          ISSG=MUL(ISTU,LSYM)
          NSGM=NCSF(ISSG)
          IF(NSGM.EQ.0) GOTO 500
* GETSGM2 computes E_UT acting on CI and saves it on SGM1
          CALL GETSGM2(LU,LT,LSYM,CI,SGM1)
          IF(ISTU.EQ.1) THEN
            GTU=DDOT_(NSGM,CI,1,SGM1,1)
            G1(IT,IU)=GTU
            G1(IU,IT)=GTU
          END IF

* SVC: The master node now continues to only handle task scheduling,
*      needed to achieve better load balancing. So it exits from the task
*      list. It has to do it here since each process gets at least one
*      task.
#if defined (_MOLCAS_MPP_) && !defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif

      GOTO 500
 501  CONTINUE
      CALL Free_Tsk(ID)

      CALL GETMEM ('Tasks','FREE','INTE',lTask,2*nTasks)

      CALL GAdSUM (G1,NG1)

  99  CONTINUE

#ifdef _ENABLE_CHEMPS2_DMRG_
      If (DoCumulant) THEN
      If(NACTEL.GT.1) Then
*QP: At this point, only load 2RDM of one state, JSTATE=1
        Call chemps2_load2pdm( nlev, G2, MSTATE(1) )
        Call two2onerdm_bis( nlev, NACTEL, G2, G1 )
      Else
        write(6,*) "FATAL ERROR: DMRG-CASPT2 with
     & CHEMPS2 does not work with NACTEL=1"
      End If
      End If
#endif

      IF(iPrGlb.GE.DEBUG) THEN
        WRITE(6,'("DEBUG> ",A)')
     &   "DENS1_RPT2: norms of the 1-el density matrix:"
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G1:", DNRM2_(NG1,G1,1)
      ENDIF

      CALL QEXIT('DENS1_RPT2')

      RETURN
      END

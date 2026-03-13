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
      SUBROUTINE DENS1_RPT2 (CI,nCI,SGM1,nSGM1,G1,nLev)
      use Symmetry_Info, only: Mul
      use caspt2_global, only:iPrGlb
      use fciqmc_interface, only: load_fciqmc_g1, DoFCIQMC
#ifdef _DMRG_
      use qcmaquis_interface, only:qcmaquis_interface_get_1rdm_full
      use caspt2_module, only: DMRG
#endif
#ifdef _ENABLE_CHEMPS2_DMRG_
      use caspt2_module, only: DoCumulant
#endif
      use PrintLevel, only: DEBUG
      use gugx, only: SGS, L2ACT, CIS
      use stdalloc, only: mma_allocate, mma_deallocate
      use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk
      use caspt2_module, only: iSCF, jState, nActEl, nAshT, STSym,
     &                         mState
      use pt2_guga, only: nG1
      use constants, only: Zero, One, Two
      use definitions, only: iwp, wp, u6
      IMPLICIT NONE


      Integer(kind=iwp), Intent(In):: nCI, nSGM1, nLev
      REAL(kind=wp), Intent(in):: CI(nCI)
      REAL(kind=wp), Intent(out):: SGM1(nSGM1)
      REAL(kind=wp), Intent(out):: G1(NLEV,NLEV)
#ifdef _ENABLE_CHEMPS2_DMRG_
      REAL(kind=wp) G2(NLEV,NLEV,NLEV,NLEV)
#endif

      REAL(kind=wp) GTU

      INTEGER(kind=iwp) ID
      INTEGER(kind=iwp) IST,ISU,ISTU
      INTEGER(kind=iwp) IT,IU,LT,LU

      INTEGER(kind=iwp) ITASK,NTASKS

      INTEGER(kind=iwp) ISSG,NSGM

      REAL(kind=wp), EXTERNAL :: DDOT_,DNRM2_
      INTEGER(kind=iwp), ALLOCATABLE:: TASK(:,:)

* Purpose: Compute the 1-electron density matrix array G1.

      G1(:,:)=Zero

      if (DoFCIQMC) then
        call load_fciqmc_g1(g1, mstate(jstate),nLev)
        Call End_Stuff()
        Return
      end if

#ifdef _DMRG_
      if (DMRG) then
        if (iPrGlb >= DEBUG) then
            write (u6,*) 'DENS1_RPT2> Calculating 1RDM...'
        end if
        call qcmaquis_interface_get_1rdm_full(G1)
        Call End_Stuff()
        Return
      end if
#endif

* For the special cases, there is no actual CI-routines involved:
* Special code for hi-spin case:
      IF(ISCF.EQ.2) THEN
        DO IT=1,NASHT
          G1(IT,IT)=One
        END DO
        Call End_Stuff()
        Return
      END IF
* Special code for closed-shell:
      IF(ISCF.EQ.1 .AND. NACTEL.GT.0) THEN
        DO IT=1,NASHT
          G1(IT,IT)=Two
        END DO
        Call End_Stuff()
        Return
      END IF

! TODO: skip completely this part if using a DMRG reference!

* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level indices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in pt2_guga.F90

* SVC20100311: set up a task table with LT,LU
* SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
      nTasks=(nLev**2+nLev)/2

      CALL mma_allocate (Task,nTasks,2,Label='TASK')

      iTask=0
      DO LT=1,nLev
        DO LU=1,LT
          iTask=iTask+1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        ENDDO
      ENDDO
      IF (iTask.NE.nTasks) WRITE(u6,*) "ERROR nTasks"

      Call Init_Tsk(ID, nTasks)

* SVC20100311: BEGIN SEPARATE TASK EXECUTION
      Do While (Rsv_Tsk (ID,iTask))

* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
      LT=TASK(iTask,1)
        IST=SGS%ISM(LT)
        IT=L2ACT(LT)
        LU=Task(iTask,2)
          ISU=SGS%ISM(LU)
          IU=L2ACT(LU)
          ISTU=Mul(IST,ISU)
          IF (ISTU/=1) CYCLE
          ISSG=Mul(ISTU,STSYM)
          NSGM=CIS%NCSF(ISSG)
          IF(NSGM.EQ.0) Cycle
* GETSGM2 computes E_UT acting on CI and saves it on SGM1
          CALL GETSGM2(LU,LT,STSYM,CI,nCI,SGM1,NSGM)
          GTU=DDOT_(NSGM,CI,1,SGM1,1)
          G1(IT,IU)=GTU
          G1(IU,IT)=GTU

* SVC: The master node now continues to only handle task scheduling,
*      needed to achieve better load balancing. So it exits from the task
*      list. It has to do it here since each process gets at least one
*      task.

      End Do

      CALL Free_Tsk(ID)
      CALL mma_deallocate(Task)

      CALL GAdGOP (G1,NG1,'+')

      Call End_Stuff()

      Contains

      Subroutine End_Stuff()
#ifdef _ENABLE_CHEMPS2_DMRG_
      If (DoCumulant) THEN
      If(NACTEL.GT.1) Then
*QP: At this point, only load 2RDM of one state, JSTATE=1
        Call chemps2_load2pdm( nlev, G2, MSTATE(1) )
        Call two2onerdm( nlev, NACTEL, G2, G1 )
      Else
        write(6,*) "FATAL ERROR: DMRG-CASPT2 with
     & CHEMPS2 does not work with NACTEL=1"
      End If
      End If
#endif


      ! DEBUG print while developing DMRG-CASPT2
      ! do i = 1,nLev
      !   write (6,'(1x,14f10.6)') (G1(i,j),j=1,nLev)
      ! end do
      ! write(6,*)


      IF(iPrGlb.GE.DEBUG) THEN
        WRITE(u6,'("DEBUG> ",A)')
     &   "DENS1_RPT2: norms of the 1-el density matrix:"
        WRITE(u6,'("DEBUG> ",A,1X,ES21.14)') "G1:", DNRM2_(NG1,G1,1)
      ENDIF
      End Subroutine End_Stuff


      END SUBROUTINE DENS1_RPT2

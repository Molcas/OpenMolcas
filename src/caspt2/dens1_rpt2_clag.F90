!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

      SUBROUTINE DENS1_RPT2_CLag (CI,NCI,SGM1,NSGM1,CLag,nConf,RDMEIG,  &
     &                            nLev)
      use Symmetry_Info, only: Mul
      use sguga, only: SGS, L2ACT, CIS
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: STSym
      use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk
#if defined (_MOLCAS_MPP_) && ! defined (_GA_)
      USE Para_Info, ONLY: Is_Real_Par, King, nProcs
#endif

      IMPLICIT NONE

      integer(kind=iwp), intent(in) :: NCI, NSGM1, nLev, nConf
      real(kind=wp), intent(in) :: CI(NCI), RDMEIG(NLEV,NLEV)
      real(kind=wp), intent(inout) :: SGM1(NSGM1), CLag(nConf)

      integer(kind=iwp), allocatable :: TASK(:,:)

      integer(kind=iwp) :: ID, IST, ISU, ISTU, IT, IU, LT, LU, ITASK,   &
     &                     NTASKS, ISSG, NSGM

! Purpose: Compute the 1-electron density matrix array G1.


! For the general cases, we use actual CI routine calls, and
! have to take account of orbital order.
! We will use level inices LT,LU... in these calls, but produce
! the density matrices with usual active orbital indices.
! Translation tables L2ACT and LEVEL, in caspt2_module.F90

! SVC20100311: set up a task table with LT,LU
! SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
      nTasks = nLev**2
      CALL mma_allocate (Task,nTasks,2,Label='TASK')

      iTask=0
      ! First, IL < JL pairs.
      Do LT = 1, nLev-1
        Do LU = LT+1, nLev
          iTask = iTask + 1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        End Do
      End Do
      ! Then, IL = JL pairs.
      Do LT = 1, nLev
        iTask = iTask + 1
        TASK(iTask,1)=LT
        TASK(iTask,2)=LT
      End Do
      ! Last, IL > JL pairs.
      Do LT = 2, nLev
        Do LU = 1, LT-1
          iTask = iTask + 1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        End Do
      End Do
      IF (iTask /= nTasks) WRITE(u6,*) 'ERROR nTasks'

      Call Init_Tsk(ID, nTasks)

! SVC20100311: BEGIN SEPARATE TASK EXECUTION
      do while (Rsv_Tsk(ID,iTask))

! Compute SGM1 = E_UT acting on CI, with T.ge.U,
! i.e., lowering operations. These are allowed in RAS.
        LT=TASK(iTask,1)
        IST=SGS%ISM(LT)
        IT=L2ACT(LT)
        LU=Task(iTask,2)
        ISU=SGS%ISM(LU)
        IU=L2ACT(LU)
        ISTU=Mul(IST,ISU)
        ISSG=Mul(ISTU,STSYM)
        NSGM=CIS%NCSF(ISSG)
        IF(NSGM == 0) cycle
! GETSGM2 computes E_UT acting on CI and saves it on SGM1
        CALL GETSGM2(LU,LT,STSYM,CI,NCI,SGM1,NSGM)
        IF(ISTU == 1) THEN
          ! Symmetry not yet
!          write(u6,*) 'it,iu = ', it,iu
           CLag(1:NSGM) = CLag(1:NSGM) + RDMEIG(IT,IU)*SGM1(1:NSGM)
        END IF

! SVC: The master node now continues to only handle task scheduling,
!     needed to achieve better load balancing. So it exits from the task
!      list. It has to do it here since each process gets at least one
!      task.
#if defined (_MOLCAS_MPP_) && ! defined (_GA_)
        IF (IS_REAL_PAR() .AND. KING() .AND. NPROCS > 1) exit
#endif
      end do

      CALL Free_Tsk(ID)
      CALL mma_deallocate(Task)

      return
      end subroutine DENS1_RPT2_CLag

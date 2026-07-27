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

subroutine DENS1_RPT2_CLag(CI,NCI,SGM1,NSGM1,CLag,nConf,RDMEIG,nLev)
! Purpose: Compute the 1-electron density matrix array G1.

! For the general cases, we use actual CI routine calls, and
! have to take account of orbital order.
! We will use level inices LT,LU... in these calls, but produce
! the density matrices with usual active orbital indices.
! Translation tables L2ACT and LEVEL, in SGS

use Symmetry_Info, only: Mul
use sguga_states, only: CIS, SGS
use general_data, only: STSym
use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
#if defined (_MOLCAS_MPP_) && ! defined (_GA_)
use Para_Info, only: Is_Real_Par, King, nProcs
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NCI, NSGM1, nConf, nLev
real(kind=wp), intent(in) :: CI(NCI), RDMEIG(NLEV,NLEV)
real(kind=wp), intent(inout) :: SGM1(NSGM1), CLag(nConf)
integer(kind=iwp) :: ID, ISSG, IST, ISTU, ISU, IT, ITASK, IU, LT, LU, NSGM, NTASKS
integer(kind=iwp), allocatable :: TASK(:,:)
integer(kind=iwp), parameter :: istate=1

! SVC20100311: set up a task table with LT,LU
! SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
nTasks = nLev**2
call mma_allocate(Task,nTasks,2,Label='TASK')

iTask = 0
! First, IL < JL pairs.
do LT=1,nLev-1
  TASK(iTask+1:iTask+nLev-LT,1) = LT
  TASK(iTask+1:iTask+nLev-LT,2) = [(LU,LU=LT+1,nLev)]
  iTask = iTask+nLev-LT
end do
! Then, IL = JL pairs.
TASK(iTask+1:iTask+nLev,1) = [(LT,LT=1,nLev)]
TASK(iTask+1:iTask+nLev,2) = [(LT,LT=1,nLev)]
iTask = iTask+nLev
! Last, IL > JL pairs.
do LT=2,nLev
  TASK(iTask+1:iTask+LT-1,1) = LT
  TASK(iTask+1:iTask+LT-1,2) = [(LU,LU=1,LT-1)]
  iTask = iTask+LT-1
end do
if (iTask /= nTasks) write(u6,*) 'ERROR nTasks'

call Init_Tsk(ID,nTasks)

! SVC20100311: BEGIN SEPARATE TASK EXECUTION
do while (Rsv_Tsk(ID,iTask))

  ! Compute SGM1 = E_UT acting on CI, with T >= U,
  ! i.e., lowering operations. These are allowed in RAS.
  LT = TASK(iTask,1)
  IST = SGS(istate)%ISM(LT)
  IT = SGS(istate)%L2ACT(LT)
  LU = Task(iTask,2)
  ISU = SGS(istate)%ISM(LU)
  IU = SGS(istate)%L2ACT(LU)
  ISTU = Mul(IST,ISU)
  ISSG = Mul(ISTU,STSYM)
  NSGM = CIS(istate)%NCSF(ISSG)
  if (NSGM == 0) cycle
  ! GETSGM2 computes E_UT acting on CI and saves it on SGM1
  call GETSGM2(LU,LT,STSYM,CI,NCI,SGM1,NSGM)
  ! Symmetry not yet
  if (ISTU == 1) CLag(1:NSGM) = CLag(1:NSGM)+RDMEIG(IT,IU)*SGM1(1:NSGM)

  ! SVC: The master node now continues to only handle task scheduling,
  !     needed to achieve better load balancing. So it exits from the task
  !      list. It has to do it here since each process gets at least one
  !      task.
# if defined (_MOLCAS_MPP_) && ! defined (_GA_)
  if (IS_REAL_PAR() .and. KING() .and. (NPROCS > 1)) exit
# endif
end do

call Free_Tsk(ID)
call mma_deallocate(Task)

return

end subroutine DENS1_RPT2_CLag

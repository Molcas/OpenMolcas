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

subroutine DENS1T_RPT2_CLag(CI1,CI2,SGM1,CLag1,CLag2,RDMEIG,SCAL,nLev)
! Purpose: Compute the 1-electron density matrix array G1.

! For the general cases, we use actual CI routine calls, and
! have to take account of orbital order.
! We will use level indices LT,LU... in these calls, but produce
! the density matrices with usual active orbital indices.
! Translation tables L2ACT and LEVEL, in caspt2_module

use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
use Symmetry_Info, only: Mul
use sguga, only: SGS, L2ACT, CIS
use caspt2_module, only: MxCI, nConf, STSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nLev
real(kind=wp), intent(in) :: CI1(MXCI), CI2(MXCI), RDMEIG(NLEV,NLEV), SCAL
real(kind=wp), intent(inout) :: SGM1(MXCI), CLag1(nConf), CLag2(nConf)
integer(kind=iwp) :: ID, ISSG, IST, ISTU, ISU, IT, ITASK, IU, LT, LU, NSGM, NTASKS
integer(kind=iwp), allocatable :: TASK(:,:)

! SVC20100311: set up a task table with LT,LU
! SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
nTasks = (nLev**2+nLev)/2
nTasks = nLev**2
call mma_allocate(Task,nTasks,2,Label='TASK')

iTask = 0
! First, IL < JL pairs.
do LT=1,nLev-1
  do LU=LT+1,nLev
    iTask = iTask+1
    TASK(iTask,1) = LT
    TASK(iTask,2) = LU
  end do
end do
! Then, IL = JL pairs.
do LT=1,nLev
  iTask = iTask+1
  TASK(iTask,1) = LT
  TASK(iTask,2) = LT
end do
! Last, IL > JL pairs.
do LT=2,nLev
  do LU=1,LT-1
    iTask = iTask+1
    TASK(iTask,1) = LT
    TASK(iTask,2) = LU
  end do
end do
if (iTask /= nTasks) write(u6,*) 'ERROR nTasks'

call Init_Tsk(ID,nTasks)

! SVC20100311: BEGIN SEPARATE TASK EXECUTION
do while (Rsv_Tsk(ID,iTask))
  ! Compute SGM1 = E_UT acting on CI, with T >= U,
  ! i.e., lowering operations. These are allowed in RAS.
  LT = TASK(iTask,1)
  IST = SGS%ISM(LT)
  IT = L2ACT(LT)
  LU = Task(iTask,2)
  !LTU = iTask
  ISU = SGS%ISM(LU)
  IU = L2ACT(LU)
  ISTU = Mul(IST,ISU)
  ISSG = Mul(ISTU,STSYM)
  NSGM = CIS%NCSF(ISSG)
  if (NSGM == 0) cycle
  ! GETSGM2 computes E_UT acting on CI and saves it on SGM1
  if (ISTU == 1) then
    call GETSGM2(LU,LT,STSYM,CI1,MXCI,SGM1,NSGM)
    CLag2(1:NSGM) = CLag2(1:NSGM)+SCAL*RDMEIG(IT,IU)*SGM1(1:NSGM)
    call GETSGM2(LU,LT,STSYM,CI2,MXCI,SGM1,NSGM)
    CLag1(1:NSGM) = CLag1(1:NSGM)+SCAL*RDMEIG(IT,IU)*SGM1(1:NSGM)
  end if
  ! SVC: The master node now continues to only handle task scheduling,
  !      needed to achieve better load balancing. So it exits from the task
  !      list. It has to do it here since each process gets at least one
  !      task.
end do
call Free_Tsk(ID)
call mma_deallocate(Task)

end subroutine DENS1T_RPT2_CLag

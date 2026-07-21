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

subroutine DENS1T_RPT2(CI1,CI2,SGM1,G1,NLEV)

use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
use Symmetry_Info, only: Mul
use PrintLevel, only: DEBUG
use sguga_states, only: SGS, CIS
use caspt2_global, only: iPrGlb
use caspt2_module, only: iSCF, MxCI, nActEl, nAshT, nG1, STSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CI1(MXCI), CI2(MXCI)
integer(kind=iwp), intent(in) :: nLev
real(kind=wp), intent(out) :: SGM1(MXCI), G1(NLEV,NLEV)
integer(kind=iwp) :: ID, ISSG, IST, ISTU, ISU, IT, ITASK, IU, LT, LU, NSGM, NTASKS
real(kind=wp) :: GTU
integer(kind=iwp), allocatable :: TASK(:,:)
real(kind=wp), external :: ddot_, dnrm2_
integer(kind=iwp), parameter :: istate=1

! Purpose: Compute the 1- and 2-electron density matrix
! arrays G1 and G2.
G1(:,:) = Zero

! For the special cases, there is no actual CI-routines involved:
if (ISCF == 2) then
  ! Special code for hi-spin case:
  do IT=1,NASHT
    G1(IT,IT) = One
  end do
else if ((ISCF == 1) .and. (NACTEL > 0)) then
  ! Special code for closed-shell:
  do IT=1,NASHT
    G1(IT,IT) = Two
  end do
else

  ! For the general cases, we use actual CI routine calls, and
  ! have to take account of orbital order.
  ! We will use level inices LT,LU... in these calls, but produce
  ! the density matrices with usual active orbital indices.
  ! Translation tables L2ACT and LEVEL, in SGS

  !-SVC20100311: set up a task table with LT,LU
  nTasks = nLev**2
  call mma_allocate(Task,nTasks,2,Label='TASK')

  iTask = 0
  do LT=1,nLev
    TASK(iTask+1:iTask+nLev,1) = LT
    TASK(iTask+1:iTask+nLev,2) = [(LU,LU=1,nLev)]
    iTask = iTask+nLev
  end do
  if (iTask /= nTasks) write(u6,*) 'ERROR nTasks'

  call Init_Tsk(ID,nTasks)

  !-SVC20100311: BEGIN SEPARATE TASK EXECUTION
  do while (Rsv_Tsk(ID,iTask))
    ! Compute SGM1 = E_UT acting on CI, with T >= U,
    ! i.e., lowering operations. These are allowed in RAS.
    !LTU = 0
    !do LT=1,NLEV
    LT = TASK(iTask,1)
    IST = SGS(istate)%ISM(LT)
    IT = SGS(istate)%L2ACT(LT)
    !do LU=1,LT
    LU = Task(iTask,2)
    !LTU = LTU+1
    !LTU = iTask
    ISU = SGS(istate)%ISM(LU)
    IU = SGS(istate)%L2ACT(LU)
    ISTU = Mul(IST,ISU)
    ISSG = Mul(ISTU,STSYM)
    NSGM = CIS(istate)%NCSF(ISSG)
    if (NSGM == 0) cycle
    call GETSGM2(LU,LT,STSYM,CI1,MXCI,SGM1,NSGM)
    if (ISTU == 1) then
      GTU = DDOT_(NSGM,CI2,1,SGM1,1)
      G1(IT,IU) = G1(IT,IU)+GTU
    end if

    !SVC: The master node now continues to only handle task scheduling,
    !     needed to achieve better load balancing. So it exits from the task
    !     list.  It has to do it here since each process gets at least one
    !     task.
  end do
  call Free_Tsk(ID)
  call mma_deallocate(Task)
  call GAdGOP(G1,NG1,'+')
end if

if (iPrGlb >= DEBUG) then
  write(u6,'("DEBUG> ",A)') 'DENS1_RPT2: norms of the density matrices:'
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G1:',DNRM2_(NG1,G1,1)
end if

end subroutine DENS1T_RPT2

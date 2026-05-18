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
! Copyright (C) 2006, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 2006  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine DENS1_RPT2(CI,nCI,SGM1,nSGM1,G1,nLev)

use Symmetry_Info, only: Mul
use fciqmc_interface, only: DoFCIQMC, load_fciqmc_g1
use PrintLevel, only: DEBUG
use sguga, only: CIS, L2ACT, SGS
use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
use caspt2_global, only: iPrGlb
use caspt2_module, only: iSCF, jState, mState, nActEl, nAshT, nG1, STSym
#ifdef _DMRG_
use qcmaquis_interface, only: qcmaquis_interface_get_1rdm_full
use caspt2_module, only: DMRG
#endif
#ifdef _ENABLE_CHEMPS2_DMRG_
use caspt2_module, only: DoCumulant
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCI, nSGM1, nLev
real(kind=wp), intent(in) :: CI(nCI)
real(kind=wp), intent(out) :: SGM1(nSGM1), G1(NLEV,NLEV)
integer(kind=iwp) :: ID, ISSG, IST, ISTU, ISU, IT, ITASK, IU, LT, LU, NSGM, NTASKS
real(kind=wp) :: GTU
#ifdef _ENABLE_CHEMPS2_DMRG_
real(kind=wp) :: G2(NLEV,NLEV,NLEV,NLEV)
#endif
integer(kind=iwp), allocatable :: TASK(:,:)
real(kind=wp), external :: DDOT_, DNRM2_

! Purpose: Compute the 1-electron density matrix array G1.

G1(:,:) = Zero

if (DoFCIQMC) then
  call load_fciqmc_g1(g1,mstate(jstate),nLev)
  call End_Stuff()
  return
end if

#ifdef _DMRG_
if (DMRG) then
  if (iPrGlb >= DEBUG) write(u6,*) 'DENS1_RPT2> Calculating 1RDM...'
  call qcmaquis_interface_get_1rdm_full(G1)
  call End_Stuff()
  return
end if
#endif

! For the special cases, there is no actual CI-routines involved:
! Special code for hi-spin case:
if (ISCF == 2) then
  do IT=1,NASHT
    G1(IT,IT) = One
  end do
  call End_Stuff()
  return
end if
! Special code for closed-shell:
if ((ISCF == 1) .and. (NACTEL > 0)) then
  do IT=1,NASHT
    G1(IT,IT) = Two
  end do
  call End_Stuff()
  return
end if

! TODO: skip completely this part if using a DMRG reference!

! For the general cases, we use actual CI routine calls, and
! have to take account of orbital order.
! We will use level indices LT,LU... in these calls, but produce
! the density matrices with usual active orbital indices.
! Translation tables L2ACT and LEVEL, in caspt2_module

! SVC20100311: set up a task table with LT,LU
! SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
nTasks = (nLev**2+nLev)/2

call mma_allocate(Task,nTasks,2,Label='TASK')

iTask = 0
do LT=1,nLev
  do LU=1,LT
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
  ISU = SGS%ISM(LU)
  IU = L2ACT(LU)
  ISTU = Mul(IST,ISU)
  if (ISTU /= 1) cycle
  ISSG = Mul(ISTU,STSYM)
  NSGM = CIS%NCSF(ISSG)
  if (NSGM == 0) cycle
  ! GETSGM2 computes E_UT acting on CI and saves it on SGM1
  call GETSGM2(LU,LT,STSYM,CI,nCI,SGM1,NSGM)
  GTU = DDOT_(NSGM,CI,1,SGM1,1)
  G1(IT,IU) = GTU
  G1(IU,IT) = GTU

  ! SVC: The master node now continues to only handle task scheduling,
  !      needed to achieve better load balancing. So it exits from the task
  !      list. It has to do it here since each process gets at least one
  !      task.

end do

call Free_Tsk(ID)
call mma_deallocate(Task)

call GAdGOP(G1,NG1,'+')

call End_Stuff()

contains

subroutine End_Stuff()

# ifdef _ENABLE_CHEMPS2_DMRG_
  if (DoCumulant) then
    if (NACTEL > 1) then
      !QP: At this point, only load 2RDM of one state, JSTATE=1
      call chemps2_load2pdm(nlev,G2,MSTATE(1))
      call two2onerdm(nlev,NACTEL,G2,G1)
    else
      write(u6,*) 'FATAL ERROR: DMRG-CASPT2 with CHEMPS2 does not work with NACTEL=1'
    end if
  end if
# endif

  ! DEBUG print while developing DMRG-CASPT2
  !do i=1,nLev
  !  write(u6,'(1x,14f10.6)') (G1(i,j),j=1,nLev)
  !end do
  !write (u6,*)

  if (iPrGlb >= DEBUG) then
    write(u6,'("DEBUG> ",A)') 'DENS1_RPT2: norms of the 1-el density matrix:'
    write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G1:',DNRM2_(NG1,G1,1)
  end if

end subroutine End_Stuff

end subroutine DENS1_RPT2

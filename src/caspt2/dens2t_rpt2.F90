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

subroutine DENS2T_RPT2(NLEV,NCONF,MXCI,CI1,CI2,SGM1,SGM2,G1,G2)

use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
use Symmetry_Info, only: Mul
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG
use sguga, only: SGS, L2ACT, CIS
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: iSCF, nActEl, nAshT, STSym
use caspt2_module, only: nG1, nG2
use constants, only: Zero, One, Two, Four
use definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NLEV, NCONF, MXCI
real(kind=wp), intent(in) :: CI1(NCONF), CI2(NCONF)
real(kind=wp), intent(out) :: SGM1(MXCI), SGM2(MXCI), G1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV)
real(kind=wp) :: GTU, GTUVX !! ,GTUXV
integer(kind=iwp) :: ID
integer(kind=iwp) :: IST, ISU, ISV, ISX, ISTU, ISVX
integer(kind=iwp) :: IT, IU, IV, IX, LT, LU, LV, LX, LVX
integer(kind=iwp) :: itu, ivx
integer(kind=iwp) :: ITASK, NTASKS
integer(kind=iwp) :: ISSG1, ISSG2, NSGM1, NSGM2
real(kind=wp), external :: DDOT_, DNRM2_
integer(kind=iwp), allocatable :: Task(:,:)

! Purpose: Compute the 1- and 2-electron density matrix
! arrays G1 and G2.

G1(:,:) = Zero
G2(:,:,:,:) = Zero

! For the special cases, there is no actual CI-routines involved:
! Special code for hi-spin case:
if (ISCF == 2) then
  do IT=1,NASHT
    G1(IT,IT) = One
    do IU=1,IT-1
      G2(IT,IT,IU,IU) = One
      G2(IU,IU,IT,IT) = One
      G2(IT,IU,IU,IT) = -One
      G2(IU,IT,IT,IU) = -One
    end do
  end do
  return
end if
! Special code for closed-shell:
if ((ISCF == 1) .and. (NACTEL > 0)) then
  do IT=1,NASHT
    G1(IT,IT) = Two
    G2(IT,IT,IT,IT) = Two
    do IU=1,IT-1
      G2(IT,IT,IU,IU) = Four
      G2(IU,IU,IT,IT) = Four
      G2(IT,IU,IU,IT) = -Two
      G2(IU,IT,IT,IU) = -Two
    end do
  end do
  return
end if

! For the general cases, we use actual CI routine calls, and
! have to take account of orbital order.
! We will use level inices LT,LU... in these calls, but produce
! the density matrices with usual active orbital indices.
! Translation tables L2ACT and LEVEL, in caspt2_module

!-SVC20100311: set up a task table with LT,LU
nTasks = (nLev**2+nLev)/2
nTasks = nLev**2

call mma_allocate(Task,nTasks,2,Label='Task')

iTask = 0
do LT=1,nLev
  do LU=1,nLev!LT
    iTask = iTask+1
    Task(iTask,1) = LT
    Task(iTask,2) = LU
  end do
end do
if (iTask /= nTasks) write(u6,*) 'ERROR nTasks'

call Init_Tsk(ID,nTasks)

!-SVC20100311: BEGIN SEPARATE TASK EXECUTION
do while (Rsv_Tsk(ID,iTask))

  ! Compute SGM1 = E_UT acting on CI, with T >= U,
  ! i.e., lowering operations. These are allowed in RAS.
  LT = Task(iTask,1)
  IST = SGS%ISM(LT)
  IT = L2ACT(LT)
  LU = Task(iTask,2)
  ISU = SGS%ISM(LU)
  IU = L2ACT(LU)
  ISTU = Mul(IST,ISU)
  ISSG1 = Mul(ISTU,STSYM)
  NSGM1 = CIS%NCSF(ISSG1)
  if (NSGM1 == 0) cycle
  call GETSGM2(LU,LT,STSYM,CI1,NCONF,SGM1,NSGM1)
  if (ISTU == STSYM) then
    GTU = DDOT_(NSGM1,CI2,1,SGM1,1)
    G1(IT,IU) = G1(IT,IU)+GTU
    !G1(IU,IT) = GTU
  end if
  LVX = 0
  do LV=1,NLEV!LT
    ISV = SGS%ISM(LV)
    IV = L2ACT(LV)
    do LX=1,NLEV!LV
      LVX = LVX+1
      ISX = SGS%ISM(LX)
      ISVX = Mul(ISV,ISX)
      !if (ISVX /= ISTU) cycle
      IX = L2ACT(LX)
      ISSG2 = Mul(ISVX,ISSG1)
      NSGM2 = CIS%NCSF(ISSG2)
      if (NSGM2 == 0) cycle
      if (ISSG2 == STSYM) then
        !if (LX == LT) then
        !  ! then actually T=U=V=X.
        !  GTUVX = DDOT_(NSGM,SGM1,1,SGM1,1)
        !else
        call GETSGM2(LX,LV,ISSG1,SGM1,NCONF,SGM2,NSGM2)
        GTUVX = DDOT_(NSGM2,CI2,1,SGM2,1)
        !end if

        !if (LV == LX) then
        !  GTUXV = GTUVX
        !else
        !  if (LVX == LTU) then
        !    GTUXV = DDOT_(NSGM,SGM1,1,SGM1,1)
        !  else
        !    call GETSGM2(LX,LV,STSYM,CI,MXCI,SGM2,NSGM)
        !    GTUXV = DDOT_(NSGM,SGM1,1,SGM2,1)
        !  end if
        !end if
        G2(IT,IU,IV,IX) = G2(IT,IU,IV,IX)+GTUVX
        !G2(IT,IU,IX,IV) = GTUXV
      end if
    end do
  end do

  call GETSGM2(LU,LT,STSYM,CI2,NCONF,SGM1,NSGM1)
  !if (ISTU == 1) then
  if (ISTU == STSYM) then
    GTU = DDOT_(NSGM1,CI1,1,SGM1,1)
    G1(IT,IU) = G1(IT,IU)+GTU
  end if
  LVX = 0
  do LV=1,NLEV
    ISV = SGS%ISM(LV)
    IV = L2ACT(LV)
    do LX=1,NLEV
      LVX = LVX+1
      !if (LVX > LTU) cycle
      ISX = SGS%ISM(LX)
      IX = L2ACT(LX)
      ISVX = Mul(ISV,ISX)
      !if (ISVX /= ISTU) cycle
      ISSG2 = Mul(ISVX,ISSG1)
      NSGM2 = CIS%NCSF(ISSG2)
      if (NSGM2 == 0) cycle
      if (ISSG2 == STSYM) then
        call GETSGM2(LX,LV,ISSG1,SGM1,NCONF,SGM2,NSGM2)
        GTUVX = DDOT_(NSGM2,CI1,1,SGM2,1)
        G2(IT,IU,IV,IX) = G2(IT,IU,IV,IX)+GTUVX
      end if
    end do
  end do

  !SVC: The master node now continues to only handle task scheduling,
  !     needed to achieve better load balancing. So it exits from the task
  !     list.  It has to do it here since each process gets at least one
  !     task.

end do

call Free_Tsk(ID)

call mma_deallocate(Task)

call GAdGOP(G1,NG1,'+')
call GAdGOP(G2,NG2,'+')

!write(u6,*) 'before'
!call sqprt(g2,nlev**2)
do LT=1,NLEV
  IT = L2ACT(LT)
  do LU=1,NLEV
    IU = L2ACT(LU)
    do LV=1,NLEV
      IV = L2ACT(LV)
      G2(IT,IV,IV,IU) = G2(IT,IV,IV,IU)-G1(IT,IU)
    end do
  end do
end do
do it=1,nlev
  do iu=1,nlev
    itu = it+nasht*(iu-1)
    do iv=1,nlev
      do ix=1,nlev
        ivx = iv+nasht*(ix-1)
        if (ivx > itu) g2(iv,ix,it,iu) = g2(it,iu,iv,ix)
      end do
    end do
  end do
end do
!-SVC20100311: serial part: add corrections to G2
!do LT=1,NLEV
!  IT = L2ACT(LT)
!  do LX=1,LT
!    IX = L2ACT(LX)
!    do LU=LX,LT
!      IU = L2ACT(LU)
!      G2(IT,IU,IU,IX) = G2(IT,IU,IU,IX)-G1(IT,IX)
!    end do
!  end do
!end do
!do LT=2,NLEV
!  IT = L2ACT(LT)
!  do LX=2,LT
!    IX = L2ACT(LX)
!    do LU=1,LX-1
!      IU = L2ACT(LU)
!      G2(IT,IU,IU,IX) = G2(IT,IU,IU,IX)-G1(IT,IX)
!    end do
!  end do
!end do
!LTU = 0
!do LT=1,NLEV
!  IT = L2ACT(LT)
!  do LU=1,LT
!    LTU = LTU+1
!    IU = L2ACT(LU)
!    LVX = 0
!    outer: do LV=1,LT
!      IV = L2ACT(LV)
!      do LX=1,LV
!        LVX = LVX+1
!        IX = L2ACT(LX)
!        if (LVX > LTU) exit outer
!        GTUVX = G2(IT,IU,IV,IX)
!        G2(IU,IT,IX,IV) = GTUVX
!        G2(IV,IX,IT,IU) = GTUVX
!        G2(IX,IV,IU,IT) = GTUVX
!        GTUXV = G2(IT,IU,IX,IV)
!        G2(IU,IT,IV,IX) = GTUXV
!        G2(IX,IV,IT,IU) = GTUXV
!        G2(IV,IX,IU,IT) = GTUXV
!      end do
!    end do outer
!  end do
!end do

if (iPrGlb >= DEBUG) then
  write(u6,'("DEBUG> ",A)') 'DENS2_RPT2: norms of the density matrices:'
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G1:',DNRM2_(NG1,G1,1)
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G2:',DNRM2_(NG2,G2,1)
end if

return

end subroutine DENS2T_RPT2

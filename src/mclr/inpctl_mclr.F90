!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine InpCtl_MCLR(iPL)
!***********************************************************************
!                                                                      *
!     Read all relevant input data and display them                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Str_Info, only: DTOC
use negpre, only: nGP
use ipPage, only: W
use gugx, only: SGS, CIS, EXS
use MCLR_Data, only: ipCI
use MCLR_Data, only: SA, ISTATE
use MCLR_Data, only: LuPT2
use input_mclr, only: PT2, iMethod, TimeDep, nCSF, nSym, State_Sym, iMCPD, nDisp, iRoot, iSpin, nActEl, nElec3, nHole1, nRS1, &
                      nRS2, nRS3, Page, nRoots, nConf
use dmrginfo, only: DoDMRG, DoMCLR, nDets_RGLR
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: u6

implicit none
integer iPL
character(len=8) Method
real*8, allocatable :: CIVec(:,:), CITmp(:)
integer i, ii, ipCII, iSSM
integer, external :: ipGet
integer, external :: IsFreeUnit
integer, allocatable :: index_SD(:) ! not final version
real*8, allocatable :: vector_cidmrg(:)
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine RdJobIph_td(CIVec)
    real*8, allocatable :: CIVec(:,:)
  end subroutine RdJobIph_td
  subroutine RdJobIph(CIVec)
    real*8, allocatable :: CIVec(:,:)
  end subroutine RdJobIph
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
! Read in interesting info from RUNFILE and ONEINT
call Rd1Int_MCLR()
call RdAB()   ! Read in orbitals, perturbation type, etc.
!                                                                      *
!***********************************************************************
!                                                                      *
call Rd2Int(iPL) ! Read in 2el header
!                                                                      *
!***********************************************************************
!                                                                      *
call RdInp_MCLR()  ! Read in input
!                                                                      *
!***********************************************************************
!                                                                      *
! Default activate ippage utility

call ipopen(.true.)

PT2 = .false.
call Get_cArray('Relax Method',Method,8)
if (Method == 'CASPT2') then
  PT2 = .true.
  ! Read the states requested by CASPT2
  ! This means that the root(s) specified in &MCLR is usually
  ! ignored for CASPT2 gradient/NAC.
  call check_caspt2(1)
end if

!write(u6,*) 'iMethod:',iMethod,2
if (iMethod == 2) then
  if (TimeDep) then
    call RdJobIph_td(CIVec)
  else
    call RdJobIph(CIVec)
  end if

  !write(u6,*) 'Setup of Determinant tables'
  call DetCtl()   ! set up determinant tables
  ! Read in tables from disk
  call InCsfSD(State_sym,State_sym)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !write(u6,*) 'Transformation of CI vector to symmetric group from GUGA pepresentation'

  ! scratch  ! yma testing
  !if (doDMRG .and. doMCLR) then
  !  call xflush(117)
  !  close(117)
  !end if

  do i=1,nroots
    ! yma
    ! No need to copy,since there are no CI-vectors
    if (doDMRG .and. doMCLR) then
      call mma_allocate(CITmp,ndets_RGLR,Label='CITmp')
    else
      call mma_allocate(CITmp,nconf,Label='CITmp')
      call dcopy_(nconf,CIVec(:,i),1,CITmp,1)
    end if

    ! If doDMRG
    if (.not. (doDMRG .and. doMCLR)) then ! yma
      ! transform to sym. group
      call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,CITmp,1,State_Sym,State_Sym)
      NCSF(1:nSym) = CIS%NCSF(1:nSym)
      NCONF = CIS%NCSF(State_Sym)
      call mkGuga_Free(SGS,CIS,EXS)

    end if

    ! Here should be the position for introducing the CI(SR) coefficients
    !iSSM = 1     ! yma
    !write(u6,*) 'Set ISSM == 1 ',ISSM

    if (doDMRG) then !yma
      call mma_allocate(index_SD,ndets_RGLR,label='index_SD')
      call mma_allocate(vector_cidmrg,ndets_RGLR,label='vector_cidmrg')
      call ci_reconstruct(i,ndets_RGLR,vector_cidmrg,index_SD)
      do ii=1,ndets_RGLR
        if (abs(vector_cidmrg(ii)) < Zero) vector_cidmrg(ii) = Zero
      end do
      call CSDTVC_dmrg(CITmp,vector_cidmrg,2,DTOC,index_SD,ISSM,1)
      call mma_deallocate(index_SD)
      call mma_deallocate(vector_cidmrg)
    end if

    call dcopy_(nconf,CITmp,1,CIVec(:,i),1)
    call mma_deallocate(CITmp)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call ipopen(page)

  ! If we are computing Lagrangian multipliers we pick up all CI
  ! vectors. For Hessian calculations we pick up just one vector.

  !write(u6,*) 'iState,SA,nroots=',iState,SA,nroots
  if (SA .or. iMCPD .or. PT2) then
    ipcii = ipget(nconf*nroots)
    call ipin(ipcii)
    call dcopy_(nconf*nroots,CIVec,1,W(ipcii)%Vec,1)
    nDisp = 1
  else
    ipcii = ipget(nconf)
    call ipin(ipcii)
    call dcopy_(nConf,CIVec(:,iState),1,W(ipcii)%Vec,1)
    if (iRoot(iState) /= 1) then
      write(u6,*) 'McKinley does not support computation of harmonic frequencies of excited states'
      call Abend()
    end if
  end if
  !call ipin(ipcii)
  !call RecPrt('CI vector',' ',W(ipcii)%Vec,1,nConf)
  call mma_deallocate(CIVec)

  ! At this point we change to ipci being the index of the CI
  ! vector in the ipage utility.

  ipci = ipcii
  call ipout(ipci)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (ngp) call rdciv()
  if (PT2) then
    LuPT2 = isFreeUnit(LuPT2)
    call Molcas_Open(LuPT2,'PT2_Lag')
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call InpOne()         ! read in oneham
call PrInp_MCLR(iPL)  ! Print all info
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine InpCtl_MCLR

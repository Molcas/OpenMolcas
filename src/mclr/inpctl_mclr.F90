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

use ipPage, only: ipget, ipin, ipout, ipopen, W
use MCLR_Data, only: ipCI, ISTATE, LuPT2, nGP, SA
use MCLR_procedures, only: RdJobIph, RdJobIph_td
use CI_Input, only: iSpin, nActEl, STSym
use input_mclr, only: CILib, iMCPD, iMethod, iRoot, nConf, nDisp, nRoots, Page, PT2, TimeDep
use dmrginfo, only: DoDMRG, DoMCLR, nDets_RGLR
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

use CI_Operations, only: CI_Lib_Dets_To_Basis, CI_Lib_Faroald, CI_Lib_File_To_Basis, CI_Lib_Initialize, CI_Lib_MCLR_Lucia, &
                         CI_Lib_Select, CI_Lib_SGUGA
use MCLR_Lucia_CI, only: MCLR_Lucia_CI_Setup
implicit none
integer(kind=iwp), intent(in) :: iPL
! the CI state of the run (SG_Setup_MCLR); ISTATE of MCLR_Data is the root
integer(kind=iwp), parameter :: iState_CI = 1
integer(kind=iwp) :: i, ii, ipCII
character(len=8) :: Method
integer(kind=iwp), allocatable :: index_SD(:) ! not final version
real(kind=wp), allocatable :: CITmp(:), CIVec(:,:), vector_cidmrg(:)
integer(kind=iwp), external :: IsFreeUnit

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
call RdInp_MCLR(iPL)  ! Read in input
!                                                                      *
!***********************************************************************
!                                                                      *
! Default activate ippage utility

call ipopen(.true.)
! the CI library of MCLR, through the CI interface (keyword CILIbrary):
! its copy of Lucia (default), SGUGA or Faroald
call MCLR_Lucia_CI_Setup()
if (CILib == 'SGUGA') then
  call CI_Lib_Select(CI_Lib_SGUGA)
else if (CILib == 'FAROALD') then
  call CI_Lib_Select(CI_Lib_Faroald)
else
  call CI_Lib_Select(CI_Lib_MCLR_Lucia)
end if

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
  ! set up the CI library (determinant tables) for the CI state
  call CI_Lib_Initialize(iState_CI)
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
      CITmp(:) = CIVec(:,i)
    end if

    ! If doDMRG
    if (.not. (doDMRG .and. doMCLR)) then ! yma
      ! transform to sym. group

      call CI_Lib_File_To_Basis(iState_CI,STSym,size(CITmp),CITmp)

    end if

    ! Here should be the position for introducing the CI(SR) coefficients
    !iSSM = 1     ! yma

    if (doDMRG) then !yma
      call mma_allocate(index_SD,ndets_RGLR,label='index_SD')
      call mma_allocate(vector_cidmrg,ndets_RGLR,label='vector_cidmrg')
      call ci_reconstruct(i,ndets_RGLR,vector_cidmrg,index_SD)
      do ii=1,ndets_RGLR
        if (abs(vector_cidmrg(ii)) < Zero) vector_cidmrg(ii) = Zero
      end do
      ! (the symmetry of the CI state; ISSM was used uninitialized here)
      call CI_Lib_Dets_To_Basis(iState_CI,STSym,ndets_RGLR,vector_cidmrg,index_SD,size(CITmp),CITmp)
      call mma_deallocate(index_SD)
      call mma_deallocate(vector_cidmrg)
    end if

    CIVec(:,i) = CITmp(:)
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
    W(ipcii)%A(1:nconf*nroots) = pack(CIVec(:,:),.true.)
    nDisp = 1
  else
    ipcii = ipget(nconf)
    call ipin(ipcii)
    W(ipcii)%A(1:nConf) = CIVec(:,iState)
    if (iRoot(iState) /= 1) then
      write(u6,*) 'McKinley does not support computation of harmonic frequencies of excited states'
      call Abend()
    end if
  end if
  !call ipin(ipcii)
  !call RecPrt('CI vector',' ',W(ipcii)%A,1,nConf)
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
else
  ! SCF reference: a CI state without active electrons, so that the
  ! state's values (CI_St) are defined for the whole run as well
  nActEl = 0
  iSpin = 1
  call SG_Setup_MCLR('NONE')
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

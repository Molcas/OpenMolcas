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

subroutine Start_MCLR()
!***********************************************************************
!                                                                      *
!     Precompute whatever can be before starting the response section  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use transform_procedures, only: SetUp_CASPT2_Tra
use MCLR_Data, only: CMO, CMO_Inv, FnMotra, FnQDat, FnTri1, INT1, LuHlf2, LuHlf3, LuMotra, LuQDat, LuTri1, LuTri2, nDens, SA
use input_mclr, only: ERASSCF, kPrint, LuAChoVec, LuChoInt, LuIChoVec, nAsh, nBas, nDel, NewCho, nFro, nIsh, nOrb, nSym, PT2, &
                      StepType, TwoStep
use PCM_grad, only: do_RF, PCM_grad_init, PCM_mod_ERASSCF, PrepPCM
use rctfld_module, only: iCharge_Ref, NonEQ_Ref, PCM
!use DWSol, only: DWSCF_init
use ISRotation, only: InvEne, InvSCF, InvSol
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iComp, Indx, iOff, iOff1, iOff2, iOpt, iRC, iSeed, iSym, iSymLbl, iType, j, lSqrDens, lTriDens, nOrbBas
real(kind=wp) :: BufFrac
character(len=8) :: Label
character(len=5) :: Fname
real(kind=wp), allocatable :: Smat(:), STmat(:)
integer(kind=iwp), external :: IsFreeUnit
logical(kind=iwp), external :: RF_On

!----------------------------------------------------------------------*
!     start                                                            *
!----------------------------------------------------------------------*

!                                                                      *
!***********************************************************************
!                                                                      *
call setup_MCLR(1)

if ((StepType /= 'RUN2') .and. btest(kPrint,2)) write(u6,'(6X,A)') 'Transformation of integrals'
! For the mp2-gradient calculations we want the transformation
! routine to produce all integrals of the occupied and virtual
! orbitals so we tell it that the whole space is inactive and
! that the active is empty

! Use the driver from the CASPT2 code (only none-squared).

! re-direct the transformed integrals to the MOTRA file
! which is preserved at the end of the calculation.
! LuTri1 is deleted.
if (TwoStep) then
  LuTri1 = LuMOTRA
  FnTri1 = FnMOTRA
  call DaName_MF_wa(LuQDAT,FnQDAT)
end if

call DaName_MF_wa(LuTri1,FnTri1)

if (newCho) then

  ! Compute inverse CMO

  lSqrDens = sum(nBas(1:nSym)**2)
  nOrbBas = sum(nOrb(1:nSym)*nBas(1:nSym))
  lTriDens = 0
  do iSym=1,nSym
    lTriDens = lTriDens+nTri_Elem(nBas(iSym))
  end do
  call mma_allocate(STmat,lTriDens,Label='STmat')
  call mma_allocate(Smat,lSqrDens,Label='Smat')

  iSymlbl = 1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  Label = 'Mltpl  0'
  iComp = 1
  call RdOne(irc,iOpt,Label,iComp,STmat,iSymlbl)

  Indx = 1
  iOff = 0
  do iSym=1,nSym
    do i=1,nBas(iSym)
      do j=1,i-1
        Smat(j+nBas(iSym)*(i-1)+iOff) = STmat(Indx)
        Smat(i+nBas(iSym)*(j-1)+iOff) = STmat(Indx)
        Indx = Indx+1
      end do
      Smat(i+nBas(iSym)*(i-1)+iOff) = STmat(Indx)
      Indx = Indx+1
    end do
    ioff = ioff+nBas(iSym)**2
  end do
  call mma_deallocate(STmat)

  call mma_allocate(CMO_inv,nOrbBas,Label='CMO_Inv')
  iOff1 = 1
  iOff2 = 1
  do iSym=1,nSym
    call dGemm_('T','N',nOrb(iSym),nBas(iSym),nBas(iSym),One,CMO(iOff2),nBas(iSym),Smat(iOff1),nBas(iSym),Zero,CMO_Inv(iOff2), &
                nOrb(iSym))

    iOff1 = iOff1+nBas(iSym)**2
    iOff2 = iOff2+nOrb(iSym)*nBas(iSym)
  end do

  call mma_deallocate(Smat)
end if

call SetUp_CASPT2_Tra(nSym,nBas,nOrb,nIsh,nAsh,nFro,nDel,CMO,nDens,LuTri1,LuTri2,LuHlf2,LuHlf3)
iType = 3  ! Means that TraCtl is called by MCLR

if ((.not. newCho) .and. (StepType /= 'RUN2')) call TraCtl_Drv(iType,.true.,1)

! fetch some data from existing file LuTri1
! (from a previous MCLR run)
! and make it available to the module intgrl
! (LuTRI1=LuMOTRA)
if (TwoStep .and. (StepType == 'RUN2')) call put_temp_data_on_intgrl(LuMOTRA,nSym,nOrb)

! Init Cholesky informations
if (NewCho) then
  BufFrac = 0.3_wp
  call Cho_X_Init(irc,BufFrac)
  iSeed = 10
  do i=1,nsym
    LuAChoVec(i) = IsFreeUnit(iSeed)
    iseed = LuAChoVec(i)+1
    write(Fname,'(A4,I1)') 'CHTA',i
    call DANAME_MF_WA(LuAChoVec(i),Fname)
  end do
  do i=1,nsym
    LuIChoVec(i) = IsFreeUnit(iSeed)
    iSeed = LuIChoVec(i)+1
    write(Fname,'(A4,I1)') 'CHTI',i
    call DANAME_MF_WA(LuIChoVec(i),Fname)
  end do
  LuChoInt(1) = IsFreeUnit(iSeed)
  write(Fname,'(A4)') 'CHIN'
  call DANAME_MF_WA(LuChoInt(1),Fname)
  LuChoInt(2) = IsFreeUnit(iSeed)
  write(Fname,'(A4)') 'CHTU'
  call DANAME_MF_WA(LuChoInt(2),Fname)
end if

call DaClos(LuTri2)
call DaClos(LuHlf2)
call DaClos(LuHlf3)

! initialize non-invariance flags (just in case)

InvEne = .true.
InvSCF = .true.
InvSol = .true.

! For dynamically weighted MCSCF
! if (SA .or. PT2) call DWSCF_init(2,nRoots)

! MCSCF with PCM
if (RF_On()) then
  if (NonEq_Ref) then
    call WarningMessage(2,'Error in MCLR')
    write(u6,*) 'NonEq=.true., invalid option'
    call Abend()
  end if
  call Init_RctFld(.false.,iCharge_Ref)
  if ((.not. SA) .and. (.not. PT2)) then ! only for SA-CASSCF grad
    call WarningMessage(2,'Error in MCLR')
    write(u6,*) 'PCM can be combined with SA-CASSCF or CASPT2 gradient only'
    call Abend()
  end if
  if (.not. PCM) then
    call WarningMessage(2,'Error in MCLR')
    write(u6,*) 'PCM-model must be used in RF-Input'
    call Abend()
  end if
  call PCM_grad_init()
  !! Construct INT1 with the correct PCM contributions
  call mma_deallocate(Int1)
  call PrepPCM()
  call InpOne()
  !! Modify ERASSCF so that they are made eigenvalues
  if (do_RF) call PCM_mod_ERASSCF(ERASSCF)
end if

call FckMat()
call StPert()

! With Cholesky there is no other choice than computing some
! integrals used for the preconditioner

if (NewCho) call cho_prec_mclr(CMO,nIsh,nASh,LuAChoVec,LuChoInt)

!----------------------------------------------------------------------*
!     exit                                                             *
!----------------------------------------------------------------------*

end subroutine Start_MCLR

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
! Copyright (C) 2025, Yoshio Nishimoto                                 *
!***********************************************************************
Module PCM_alaska

  use Constants, only: Zero, One, Two
  use definitions, only: iwp,wp,u6
  use pso_stuff, only: nDens
  use rctfld_module, only: nTS, PCM
  use stdalloc, only: mma_allocate, mma_deallocate

  Implicit None

  !! PCM_SQ_ind(1,:) :: nuclear-induced charges
  !! PCM_SQ_ind(2,:) :: electron-induced charges
  ! ASCs polarized by the relaxed density matrix (after Z-vector)
  real(kind=wp), Allocatable :: PCM_SQ_ind(:,:)
  ! CMO obtained in RASSCF (PCM + CMO)
  real(kind=wp), Allocatable :: PCMO(:)
  ! density matrix that polarizes the ASCs during SCF (in AO)
  ! def_solv = 1: state-specific
  ! def_solv = 3: state-averaged
  real(kind=wp), Allocatable :: DSA_AO(:)

  logical(kind=iwp) :: lSA = .false.
  integer(kind=iwp) :: def_solv = 0_iwp

  contains
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_alaska_lSA()

  implicit none

  integer(kind=iwp) :: iGo
  Character(Len=8) Method

  Call Get_cArray('Relax Method',Method,8)
  lSA=.false.
  if (Method=='CASSCFSA' .or. Method=='DMRGSCFS' .or. Method=='GASSCFSA' .or. &
      Method=='RASSCFSA' .or. Method=='CASPT2  ') then
    Call Get_iScalar('SA ready',iGo)
    if (iGO==1) lSA=.true.
  end if

  if (.not.PCM) return
  if (.not.lSA) return

  call PCM_alaska_init()

  return

  End Subroutine PCM_alaska_lSA
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_alaska_init()

  use Sizes_of_Seward, only: S

  implicit none

  integer(kind=iwp) :: mCMO, iPCMRoot

  if (def_solv == 0_iwp) then
    call Get_iScalar('RF CASSCF root',iPCMRoot)
    if (iPCMRoot >   0) def_solv = 1_iwp
    if (iPCMRoot ==  0) def_solv = 3_iwp
    if (iPCMRoot == -1) def_solv = 4_iwp
    if (iPCMRoot == -2) def_solv = 5_iwp
    if (iPCMRoot == -3) def_solv = 6_iwp

    if (def_solv /= 1 .and. def_solv /= 3) then
      write (u6,*) ' RFROOT in &RASSCF must be non-negative for SA-MCSCF/PCM gradients'
      call abend()
    end if
  end if

  call mma_allocate(PCM_SQ_ind,2,nTS,Label='PCM_SQ_ind')
  PCM_SQ_ind = Zero

  mCMO = S%n2Tot
  Call mma_allocate(PCMO,mCMO,Label='CMO')
  Call Get_dArray_chk('Last orbitals',PCMO,mCMO)

  return

  End Subroutine PCM_alaska_init
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_alaska_final()

  implicit none

  if (.not.PCM) return
  if (.not.lSA) return

  if (allocated(PCM_SQ_ind)) call mma_deallocate(PCM_SQ_ind)
  if (allocated(PCMO))       call mma_deallocate(PCMO)
  if (allocated(DSA_AO))     call mma_deallocate(DSA_AO)

  return

  End Subroutine PCM_alaska_final
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_alaska_prep()

  use Basis_Info, only: nBas
  use etwas, only: nAsh
  use PCM_arrays, only: PCM_SQ
  use pso_stuff, only: nDens
  use Symmetry_Info, only: nIrrep

  implicit none

  real(kind=wp), allocatable :: Dtmp(:), h1(:), TwoHam(:)

  integer(kind=iwp) :: iCharge, iIrrep, nAct, nG1
  logical(kind=iwp) :: Dff, First, NonEq
  real(kind=wp) :: RepNuc, Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge

  call Set_Basis_Mode('Valence')
  call Setup_iSD()
  call Get_iArray('nAsh',nAsh,nIrrep)

  ! First, construct the density matrix that polarizes the ASCs during SCF (in AO)

  nAct = 0
  nDens = 0
  Do iIrrep = 0, nIrrep-1
    nAct = nAct + nAsh(iIrrep)
    nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
  End Do
  nG1 = nAct*(nAct+1)/2

  call mma_allocate(DSA_AO,nDens,Label='DSA_AO')
  call mma_allocate(Dtmp,nDens,Label='Dtmp') ! ntBtri length

  ! Read the density matrix that polarized the ASCs during SCF
  ! constructed in mclr/inpone.F90
  DSA_AO(:) = Zero
  Call Get_dArray_chk('D1ao_PCM',DSA_AO,nDens)

  ! Second, obtain ASCs induced by the effective density (after Z-vector)
  ! The induced charges are put in PCM_SQ_ind

  call Get_D1ao_Var(Dtmp,nDens) ! Read from D1aoVar
  call mma_allocate(h1,nDens,Label='h1')
  call mma_allocate(TwoHam,nDens,Label='TwoHam')
  RepNuc = Zero
  First = .true.
  Dff = .false.
  NonEq = .false.
  Call DrvPCM(h1,TwoHam,Dtmp,RepNuc,nDens,First,Dff,NonEq)
  Call Get_dArray('PCM Charges',PCM_SQ_ind,2*nTs)

  ! Third, obtain ASCs induced by the ASC-polarizing density
  ! In the conventional implementation, the ASCs are stored in PCM_SQ.
  ! However, the ASCs have been overwritten during MCLR so the correct ASCs are newly computed here.

  Call DrvPCM(h1,TwoHam,DSA_AO,RepNuc,nDens,First,Dff,NonEq)
  Call Get_dArray('PCM Charges',PCM_SQ,2*nTs)

  call Free_iSD()

  ! tentative patch for MECI calculations with MS-type CASPT2
  ! "Reaction Field" has been overwritten in MCLR, but the correct reaction field is needed
  ! for computing H0 in MS-type calculations for other roots, so recompute it here.
  ! Note that ALASKA should not be called for PTED solvation.
  h1(:) = Zero
  TwoHam(:) = Zero

  call Get_dScalar('Total Nuclear Charge',Tot_Nuc_Charge)
  call Get_dScalar('Total Charge',Tot_El_Charge)
  Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
  iCharge=Int(Tot_Charge)
  Call DrvRF(h1,TwoHam,DSA_AO,RepNuc,nDens,First,Dff,NonEq,iCharge)

  if (allocated(Dtmp))   call mma_deallocate(Dtmp)
  if (allocated(h1))     call mma_deallocate(h1)
  if (allocated(TwoHam)) call mma_deallocate(TwoHam)

  End Subroutine PCM_alaska_prep

End Module PCM_alaska

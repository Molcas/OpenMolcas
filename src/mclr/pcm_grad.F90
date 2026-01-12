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
! Some memo
! "RASSCF energy for state X" is D^SS*V(D^SS)
! State averaged energies, which are used for the optimization in RASSCF, are D^SS*V(D^SS)
!
! electron-electron contributions are symmetric: D1*V(D2) = D2*V(D1),
! whereas electron-nuclear contributions are not always: D1*V(N) /= D2*V(N),
! so D*V/2 + D*V/2 = D*V is valid only if the explicit density used for defining the (SS) energy
! is the same to the ASC density
! For D^SS*V(D^SA), the correct nuclear contribution is (D^SS + D^SA)/2 * V^N
!
! The definition of the PCM density during SCF:
! def_solv = 1     : state-specific
! def_solv = 3,4,5 : state-average
! Only def_solv = 1, 3, and 4 are supported, but def_solv = 1 alone
! is physically meaningful while not mathematically.
! def_solv = 3 is the opposite situation.
! State-averaged PCM is invariant wrt rotations between internal states,
! whereas state-specific PCM is not invariant, so the state rotation contributions have to be
! included during Z-vector.
! The state rotation contribution may be computed at the first cycle unless def_solv = 3 (only?).

module PCM_grad

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use DWSol, only: DWSol_init, DWSol_final, DWSol_wgt, DWSolv, W_SOLV
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none

! RF_On() = .true. is equivalent to lRF = .true.
! PCM is for PCM (must be combined with cond: C-PCM),
! so PCM can be false even if lRF = .true. (not considered, though)
! do_RF below is different: do_RF is whether PCM equation is solved in MCLR or not.
! If RFPERT is true, we use the fixed reaction field and do not consider the response of the reaction field

! iStpPCM = 0 : before initization
! iStpPCM = 1 : after initialization, before Z-vector
! iStpPCM = 2 : during Z-vector
! iStpPCM = 3 : after Z-vector
integer(kind=iwp) :: iStpPCM = 0
! flag for PCM
logical(kind=iwp) :: do_RF = .false.
! definition of the solvation energy (WIP)
! def_solv = 1 : \Delta F = D^SS*V(D^SS)/2
! def_solv = 3 : \Delta F = D^SS * (V^N + V(D^SA)) - D^SA*V(D^SA)/2
! def_solv = 4 : \Delta F = D^SS * V(D^SA)/2 = D^SS * (V^N + V(D^SA)) - D^SS*V(D^SA)/2
! def_solv = 5 : \Delta F = D^SA * V(D^SA)/2
integer(kind=iwp) :: def_solv = 0_iwp
! whether RFpert has been used in CASPT2
logical(kind=iwp) :: PT2_solv = .false.
! RFPERT
logical(kind=iwp) :: RFPERT = .false.

!!!  One-electron integral contirbutions (in AO)
! from the SCF density
! PCMSCFAO(:,1) :: nuclear + electron contributions
! PCMSCFAO(:,2) :: nuclear contributions (1st index of DrvXV)
! PCMSCFAO(:,3) :: electron contributions (2nd index of DrvXV)
real(kind=wp), allocatable :: PCMSCFAO(:,:)
! from the state-specific (RLXROOT) density
! The meaning of the second index is the same above
real(kind=wp), allocatable :: PCMSSAO(:,:)
! from the inactive density
real(kind=wp), allocatable :: PCMIAO(:,:)

!!!  One-electron integral contirbutions (in MO)
! from the SCF density
! The meaning of the second index is the same to AO
real(kind=wp), allocatable :: PCMSCFMO(:,:)
! from the state-specific (RLXROOT) density
real(kind=wp), allocatable :: PCMSSMO(:,:)
! from the state-specific (RLXROOT) density
real(kind=wp), allocatable :: PCMSSMOori(:,:)
! from the unrelaxed PT2 density
real(kind=wp), allocatable :: PCMPT2MO(:,:)

!!! Square active density (in MO)
! used for PCM in SCF
real(kind=wp), allocatable, target :: DSCFMO(:,:)
! state-specific (RLXROOT) density
real(kind=wp), allocatable, target :: DSSMO(:,:)

!!! Square active density (in AO)
! used for PCM in SCF
real(kind=wp), allocatable :: DSCFAO(:)
! state-specific (RLXROOT) density
real(kind=wp), allocatable :: DSSAO(:)

!!! Response densities and integrals
real(kind=wp), allocatable :: DZMO(:)      ! square density in MO (all orbitals)
real(kind=wp), allocatable :: DZACTMO(:,:) ! active response density in MO
real(kind=wp), allocatable :: DZAO(:)      ! square density in AO
real(kind=wp), allocatable :: PCMZMO(:,:)  ! square integral in MO
real(kind=wp), allocatable :: PCMZAO(:,:)  ! square integral in AO

!! E_NN^PCM: nuclear-nuclear interaction energies (< 0)
!! To be precise, interaction energies between nuclei (positive) and
!! ASCs induced by nuclei (negative)
real(kind=wp) :: potnuc_pcm = Zero

integer(kind=iwp) :: iCharge_PCM ! from inpone.F90
integer(kind=iwp) :: IPCMROOT

contains

!-----------------------------------------------------------------------

subroutine PCM_grad_init()

  use input_mclr, only: ERASSCF, nRoots, ntAsh, ntBsqr
  use ISRotation, only: InvEne, InvSCF, InvSol

  integer(kind=iwp) :: iR

  call mma_allocate(PCMSCFAO,ntBsqr,3,Label='PCMSCFAO')
  call mma_allocate(PCMSSAO,ntBsqr,3,Label='PCMSSAO')

  call mma_allocate(PCMSCFMO,ntBsqr,3,Label='PCMSCFMO')
  call mma_allocate(PCMSSMO,ntBsqr,3,Label='PCMSSMO')
  call mma_allocate(PCMSSMOori,ntBsqr,3,Label='PCMSSMOori')
  call mma_allocate(PCMPT2MO,ntBsqr,3,Label='PCMPT2MO')

  call mma_allocate(DSCFMO,ntAsh,ntAsh,Label='DSCFMO')
  call mma_allocate(DSSMO,ntAsh,ntAsh,Label='DSSMO')

  call mma_allocate(DSCFAO,ntBsqr,Label='DSCFAO')
  call mma_allocate(DSSAO,ntBsqr,Label='DSSAO')

  call mma_allocate(DZMO,ntBsqr,Label='DZMO')
  call mma_allocate(DZACTMO,ntAsh,ntAsh,Label='DZACTMO')
  call mma_allocate(DZAO,ntBsqr,Label='DZAO')
  call mma_allocate(PCMZMO,ntBsqr,3,Label='PCMZMO')
  call mma_allocate(PCMZAO,ntBsqr,3,Label='PCMZAO')

  iStpPCM = 1 !! after initialization, before Z-vector
  do_RF = .true.
  if (RFPERT) do_RF = .false.

  if (def_solv == 0) then
    call Get_iScalar('RF CASSCF root',iPCMRoot)
    select case (iPCMRoot)
      case (1:)
        def_solv = 1
      case (0)
        def_solv = 3
      !case (-1)
      !  def_solv = 4
      !case (-2)
      !  def_solv = 5
      !case (-3)
      !  def_solv = 6
    end select

    if (def_solv == 0) then
      write(u6,'(" RFROOT = ",i3," is not considered yet")') iPCMRoot
      call abend()
    end if
    if ((def_solv /= 1) .and. (def_solv /= 3)) then
      write(u6,'(" RFROOT = ",i3," cannot be used")') iPCMRoot
      write(u6,'(" Please use RFROOT >= 0 in &RASSCF")')
      call abend()
    end if
  end if

  !! For weighted solvation
  !! assume that noneq is irrelevant in MCLR
  call DWSol_init(IPCMROOT,nRoots,.false.)
  call DWSol_wgt(2,ERASSCF)
  !write(u6,*) 'weight of the solvation density'
  !do i=1,nroots
  !  write (u6,*) 'weight:',i,W_solv(i)
  !end do

  InvSol = .true.

  if (nRoots == 1) return

  if (def_solv == 1) then
    InvEne = .false.
    InvSCF = .false.
    InvSol = .false.
  end if
  if (def_solv == 3) then
    InvEne = .true.
    InvSCF = .true.
    InvSol = .true.
    do iR=2,nRoots
      if (W_solv(1) /= W_solv(iR)) InvSol = .false.
    end do
  end if

  if (DWSolv%DWZeta /= One) InvSCF = .false.

end subroutine PCM_grad_init

!-----------------------------------------------------------------------

subroutine PCM_grad_final()

  !! nullify the pointer

  call mma_deallocate(PCMSCFAO)
  call mma_deallocate(PCMSSAO)

  call mma_deallocate(PCMSCFMO)
  call mma_deallocate(PCMSSMO)
  call mma_deallocate(PCMSSMOori)
  call mma_deallocate(PCMPT2MO)

  call mma_deallocate(DSCFMO)
  call mma_deallocate(DSSMO)

  call mma_deallocate(DSCFAO)
  call mma_deallocate(DSSAO)

  call mma_deallocate(DZMO)
  call mma_deallocate(DZACTMO)
  call mma_deallocate(DZAO)
  call mma_deallocate(PCMZMO)
  call mma_deallocate(PCMZAO)

  !call mma_deallocate(W_SOLV)
  call DWSol_final()

end subroutine PCM_grad_final

!-----------------------------------------------------------------------

subroutine PrepPCM()

  use MCLR_Data, only: ipCM, isNAC
  use input_mclr, only: iSpin, nBas, nSym, ntBsqr, ntBtri

  integer(kind=iwp) :: ip1, iSym, leng
  real(kind=wp) :: ExFac, PotNuc
  logical(kind=iwp) :: Do_DFT, Dff, First, lRF, NonEq
  real(kind=wp), allocatable :: D1ao(:), Gtmp(:), Htmp(:)

  ! Prepare PCM-related integrals etc

  leng = ntBtri
  call mma_allocate(Htmp,leng,Label='Htmp')
  call mma_allocate(Gtmp,leng,Label='Gtmp')
  call mma_allocate(D1ao,ntBsqr,Label='D1ao')

  !! Compute the density used for solvation during SCF
  call PCM_grad_dens(1) ! SCF
  call PCM_grad_dens2(1,DSCFMO,DSCFAO)
  call fold(nSym,nBas,DSCFAO,D1ao)

  !write(6,*) 'DSCFAO'
  !do isym=1,leng
  !  write(6,'(i3,f20.10)') isym,d1ao(isym)
  !end do

  NonEq = .false.
  First = .true.
  Dff = .false.
  Do_DFT = .true.
  lRF = .true.
  ExFac = Zero
  PotNuc = Zero
  call Get_dScalar('PotNuc',PotNuc)

  !! Htmp: nuclear-electron contributions
  !! Gtmp: electron-electron contributions
  Htmp(:) = Zero
  Gtmp(:) = Zero
  call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF,'SCF',ExFac,iCharge_PCM,iSpin,'1234',Do_DFT)

  ip1 = 1
  do iSym=1,nSym
    call square(Htmp(ip1),PCMSCFAO(ipCM(iSym),2),1,nBas(iSym),nBas(iSym))
    call square(Gtmp(ip1),PCMSCFAO(ipCM(iSym),3),1,nBas(iSym),nBas(iSym))
    ip1 = ip1+nTri_Elem(nBas(iSym))
  end do
  PCMSCFAO(:,1) = PCMSCFAO(:,2)
  PCMSCFAO(:,1) = PCMSCFAO(:,1)+PCMSCFAO(:,3)

  PCMSCFMO(:,:) = PCMSCFAO(:,:)
  call tcmo(PCMSCFMO(1,1),1,1)
  call tcmo(PCMSCFMO(1,2),1,1)
  call tcmo(PCMSCFMO(1,3),1,1)

  !! Compute SS density
  call PCM_grad_dens(2) ! SS
  if (isNAC) then
    call PCM_grad_dens2(2,DSSMO,DSSAO)
  else
    call PCM_grad_dens2(1,DSSMO,DSSAO)
  end if
  call fold(nSym,nBas,DSSAO,D1ao)

  Htmp(:) = Zero
  Gtmp(:) = Zero
  call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF,'SCF',ExFac,iCharge_PCM,iSpin,'1234',Do_DFT)

  ip1 = 1
  do iSym=1,nSym
    call square(Htmp(ip1),PCMSSAO(ipCM(iSym),2),1,nBas(iSym),nBas(iSym))
    call square(Gtmp(ip1),PCMSSAO(ipCM(iSym),3),1,nBas(iSym),nBas(iSym))
    ip1 = ip1+nTri_Elem(nBas(iSym))
  end do
  PCMSSAO(:,1) = PCMSSAO(:,2)
  PCMSSAO(:,1) = PCMSSAO(:,1)+PCMSSAO(:,3)

  PCMSSMO(:,:) = PCMSSAO(:,:)
  call tcmo(PCMSSMO(1,1),1,1)
  call tcmo(PCMSSMO(1,2),1,1)
  call tcmo(PCMSSMO(1,3),1,1)

  PCMSSMOori(:,:) = PCMSSMO(:,:)

  call mma_deallocate(Htmp)
  call mma_deallocate(Gtmp)
  call mma_deallocate(D1ao)

end subroutine PrepPCM

!-----------------------------------------------------------------------

subroutine PrepPCM2(mode,DMO,DAO,PCMAO,PCMMO)

  use input_mclr, only: iSpin, nBas, nSym, ntAsh, ntBsqr, ntBtri

  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: DMO(ntAsh**2)
  real(kind=wp), intent(inout) :: DAO(ntBsqr), PCMAO(ntBsqr,3), PCMMO(ntBsqr,3)
  integer(kind=iwp) :: ip1, iSym, leng
  real(kind=wp) :: ExFac, PotNuc
  logical(kind=iwp) :: Dff, Do_DFT, First, lRF, NonEq
  real(kind=wp), allocatable :: D1ao(:), Gtmp(:), Htmp(:)

  ! Prepare PCM-related integrals etc
  ! mode = 1: gradient
  ! mode = 2: non-adiabatic coupling

  leng = ntBtri
  call mma_allocate(Htmp,leng,Label='Htmp')
  call mma_allocate(Gtmp,leng,Label='Gtmp')
  call mma_allocate(D1ao,leng,Label='D1ao')

  !! Compute SCF density
  call PCM_grad_dens2(mode,DMO,DAO)
  call fold(nSym,nBas,DAO,D1ao)

  NonEq = .false.
  First = .true.
  Dff = .false.
  Do_DFT = .true.
  lRF = .true.
  ExFac = Zero
  PotNuc = Zero
  call Get_dScalar('PotNuc',PotNuc)

  !! Htmp: nuclear-electron contributions
  !! Gtmp: electron-electron contributions
  Htmp(:) = Zero
  Gtmp(:) = Zero
  call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF,'SCF',ExFac,iCharge_PCM,iSpin,'1234',Do_DFT)

  ip1 = 1
  do iSym=1,nSym
    call square(Htmp(ip1),PCMAO(1,2),1,nBas(iSym),nBas(iSym))
    call square(Gtmp(ip1),PCMAO(1,3),1,nBas(iSym),nBas(iSym))
    ip1 = ip1+nTri_Elem(nBas(iSym))
  end do
  PCMAO(:,1) = PCMAO(:,2)
  PCMAO(:,1) = PCMAO(:,1)+PCMAO(:,3)

  PCMMO(:,:) = PCMAO(:,:)
  call tcmo(PCMMO(1,1),1,1)
  call tcmo(PCMMO(1,2),1,1)
  call tcmo(PCMMO(1,3),1,1)

  call mma_deallocate(Htmp)
  call mma_deallocate(Gtmp)
  call mma_deallocate(D1ao)

end subroutine PrepPCM2

!-----------------------------------------------------------------------

subroutine PCM_grad_dens(mode)

  use ipPage, only: W
  use MCLR_Data, only: ipCI, iRlxRoot, isNAC, NSSA, xispsm
  use input_mclr, only: nConf, ncsf, nRoots, ntAsh, State_Sym

  integer(kind=iwp), intent(in) :: mode
  integer(kind=iwp) :: i, j, jR, kR, nConf1, nConfL, nConfR, ng1, ng2
  real(kind=wp), pointer :: G1q(:,:)
  real(kind=wp), allocatable :: CIL(:), CIR(:), G1r(:,:), G2r(:)

  ! Compute active density in MO
  ! mode=1: construct density used for PCM in SCF
  ! mode=2: mostly state-specific density

  nconf1 = ncsf(state_sym)
  nConfL = max(nconf1,nint(xispsm(state_sym,1)))
  nConfR = nConfL

  ng1 = ntash*ntash
  ng2 = nTri_Elem(ng1)

  call mma_allocate(CIL,nConfL,Label='CIL')
  call mma_allocate(CIR,nConfR,Label='CIR')
  call mma_allocate(G1r,ntAsh,ntAsh,Label='G1r')
  call mma_allocate(G2r,ng2,Label='G2r') ! not used?

  select case (mode)
    case (1)
      G1q(1:ntAsh,1:ntAsh) => DSCFMO(:,1:ntAsh)
      G1r(:,:) = Zero
      G1q(:,:) = Zero
      do jR=1,nRoots
        if (W_SOLV(jR) <= 1.0e-10_wp) cycle
        call CSF2SD(W(ipCI)%A(1+(jR-1)*nconf1),CIL,state_sym)
        call CSF2SD(W(ipCI)%A(1+(jR-1)*nconf1),CIR,state_sym)
        G1r(:,:) = Zero
        call Densi2_MCLR(1,G1r,G2r,CIL,CIR,0,0,0,ng1,ng2)
        !! For RDM1
        G1q(:,:) = G1q(:,:)+W_SOLV(jR)*G1r(:,:)
      end do
    case (2)
      G1q(1:ntAsh,1:ntAsh) => DSSMO(:,1:ntAsh)
      G1r(:,:) = Zero
      G1q(:,:) = Zero
      if (isNAC) then
        jR = NSSA(2)
        kR = NSSA(1)
      else
        jR = iRlxRoot
        kR = iRlxRoot
      end if
      call CSF2SD(W(ipCI)%A(1+(jR-1)*nconf1),CIL,state_sym)
      call CSF2SD(W(ipCI)%A(1+(kR-1)*nconf1),CIR,state_sym)
      call Densi2_MCLR(1,G1r,G2r,CIL,CIR,0,0,0,ng1,ng2)
      !! For RDM1
      G1q(:,:) = G1q(:,:)+G1r(:,:)
      do i=1,ntAsh
        do j=1,i-1
          G1q(i,j) = (G1q(i,j)+G1q(j,i))*Half
          G1q(j,i) = G1q(i,j)
        end do
      end do
    case default
      nullify(G1q)
      G1r(:,:) = Zero
  end select

  call mma_deallocate(CIL)
  call mma_deallocate(CIR)
  call mma_deallocate(G1r)
  call mma_deallocate(G2r)

  nConf = ncsf(state_sym)
  !! use weights
  !if (mode == 1) G1q(:,:) = G1q(:,:)/real(nRoots,kind=wp)

end subroutine PCM_grad_dens

!-----------------------------------------------------------------------

subroutine PCM_grad_dens2(mode,DMO,DAO)

  use input_mclr, only: nAsh, nFro, nIsh, nOrb, nSym, ntAsh, ntBsqr
  use MCLR_Data, only: ipMat, nA

  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: DMO(ntAsh,ntAsh)
  real(kind=wp), intent(out) :: DAO(:)
  integer(kind=iwp) :: iB, iiB, ijB, iOrb, iS, jB
  real(kind=wp) :: rd
  real(kind=wp), allocatable :: WRK(:)

  ! Compute AO density using MO density (DMO)

  call mma_allocate(WRK,ntBsqr,Label='WRK')

  !! transform MO to AO
  WRK(:) = Zero
  do iS=1,nSym
    iOrb = nOrb(iS)
    if (mode == 1) then
      !! For NAC, do not include the inactive density
      do iB=1,nFro(is)+nIsh(is)
        WRK(ipmat(iS,iS)+iB-1+iOrb*(iB-1)) = Two
      end do
    end if

    do iB=1,nAsh(iS)
      do jB=1,nAsh(iS)
        iiB = nA(iS)+ib
        ijB = nA(iS)+jb
        rd = (DMO(iiB,ijB)+DMO(ijB,iiB))*Half
        iiB = nFro(iS)+nIsh(iS)+ib
        ijB = nFro(iS)+nIsh(iS)+jb
        WRK(ipmat(iS,iS)+iiB-1+iOrb*(ijB-1)) = rd
      end do
    end do
  end do
  call tcmo(WRK,1,-2)
  DAO(1:ntBsqr) = WRK(1:ntBsqr)

  call mma_deallocate(WRK)

end subroutine PCM_grad_dens2

!-----------------------------------------------------------------------

subroutine PCM_grad_D2V(dens,vintMO,vintAO,first_,Dff_,NonEq_,Do_DFT_)

  use input_mclr, only: iSpin, nBas, nSym, ntBsqr, ntBtri
  use MCLR_Data, only: ipCM

  real(kind=wp), intent(in) :: dens(ntBsqr)
  real(kind=wp), intent(out) :: vintMO(ntBsqr,3), vintAO(ntBsqr,3)
  logical(kind=iwp), intent(in) :: first_, Dff_, NonEq_, Do_DFT_
  integer(kind=iwp) :: ip1, iSym, leng
  real(kind=wp) :: ExFac, PotNuc
  logical(kind=iwp) :: Dff, Do_DFT, first, lRF, NonEq
  real(kind=wp), allocatable :: D1ao(:), Gtmp(:), Htmp(:)

  ! density matrix (D; AO) --> electrostatic potential integral (V; AO)

  vintMO(:,:) = Zero
  vintAO(:,:) = Zero

  leng = ntBtri
  First = First_
  Dff = Dff_
  NonEq = NonEq_
  Do_DFT = Do_DFT_
  lRF = .true.
  ExFac = Zero
  PotNuc = Zero ! not used actually

  call mma_allocate(Htmp,leng,Label='Htmp')
  call mma_allocate(Gtmp,leng,Label='Gtmp')
  call mma_allocate(D1ao,leng,Label='D1ao')

  call fold(nSym,nBas,dens,D1ao)
  Htmp(:) = Zero
  Gtmp(:) = Zero

  ! h1, vint, and dens is triangular in the AO basis
  !call DrvRF(Htmp,Gtmp,D1ao,RepNuc,nh1,First,Dff,NonEq,iCharge)
  call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF,'SCF',ExFac,iCharge_PCM,iSpin,'1234',Do_DFT)

  ip1 = 1
  do iSym=1,nSym
    call square(Htmp(ip1),vintAO(ipCM(iSym),2),1,nBas(iSym),nBas(iSym))
    call square(Gtmp(ip1),vintAO(ipCM(iSym),3),1,nBas(iSym),nBas(iSym))
    ip1 = ip1+nTri_Elem(nBas(iSym))
  end do
  vintAO(:,1) = vintAO(:,2)
  vintAO(:,1) = vintAO(:,1)+vintAO(:,3)

  vintMO(:,:) = vintAO(:,:)
  call tcmo(vintMO(1,1),1,1)
  call tcmo(vintMO(1,2),1,1)
  call tcmo(vintMO(1,3),1,1)

  call mma_deallocate(Htmp)
  call mma_deallocate(Gtmp)
  call mma_deallocate(D1ao)

end subroutine PCM_grad_D2V

!-----------------------------------------------------------------------

subroutine PCM_grad_TimesE2(idSym,rKappa,FockOut,ipCIOut)

  use input_mclr, only: nAsh, nBas, nIsh, nOrb, nRoots, nSym, ntBsqr, State_Sym, weight
  use ISRotation, only: ISR
  use MCLR_Data, only: ipCI, ipCM, ipMat, nA
  use MCLR_procedures, only: CISigma_sa

  integer(kind=iwp), intent(in) :: idSym, ipCIOut
  real(kind=wp), intent(in) :: rKappa(*)
  real(kind=wp), intent(inout) :: FockOut(*)
  integer(kind=iwp) :: iA, ip1, ip2, iS, jA, jS
  real(kind=wp) :: rd, rdum(1)
  real(kind=wp), allocatable :: WRK(:)

  ! Compute implicit contributions
  ! Explicit contributions are computed in existing subroutines
  !
  ! The PCM contribution has a quartic dependence on the CI coefficient.
  ! The energy contribution to the S state is like (with state-averaged PCM):
  !   E^S = sum_pq Dpq^S K_{pq,rs} D_rs^SA
  !       = sum_pq sum_IJ <CI^S|Epq^IJ|CJ^S> K_{pq,rs} sum_T sum_rs w_T <CK^T|E_rs|CL^T>
  ! Second derivatives of the same density matrix can be evaluated by one-center-type
  ! integrals, whereas those of the different density matrix is evaluated by using
  ! response density
  ! The former is like: CId(I) = <CI|Epq^IJ|ZJ>*Hpq(D)
  ! The latter is like: CId(I) = <CI|Epq^IJ|CJ>*Hpq(Z)
  ! The former is evaluated through FIMO (or FAMO), but the latter is newly evaluated
  !
  ! Different from the CI Lagrangian, TimesE2 is not dependent on the choice of the solvation energy,
  ! because Z-vector uses the diagonality of the Hamiltonian, but not the definition of the state-specific energy.

  !! Construct the orbital response density using SCF (state-averaged) density
  !! The Lagrangian (for generalized Brillouin theorem) is defined with the SCF density
  call mma_allocate(WRK,ntBsqr,Label='WRK')
  call OITD(rKappa,1,DZMO,WRK,.true.)
  call mma_deallocate(WRK)

  !! Then, CI response density (DZACTMO), which has been computed in CIDens_sa
  !! through ci_kap (iStpPCM=2) or out_pt2 (iStpPCM=3)
  do iS=1,nSym
    if (nAsh(iS) > 0) then
      jS = Mul(is,idSym)
      do jA=1,nAsh(js)
        !ip1 = nOrb(iS)*(nIsh(is)+iA-1)+ipCM(is)
        ip2 = nIsh(iS)+nOrb(iS)*(nIsh(js)+jA-1)+ipmat(is,js)
        DZMO(ip2:ip2+nAsh(iS)-1) = DZMO(ip2:ip2+nAsh(iS)-1)+DZACTMO(1+nA(iS):nAsh(iS)+nA(iS),jA+nA(jS))
      end do
    end if
  end do

  !! Transform to AO
  DZAO(:) = DZMO(:)
  call TCMO(DZAO,1,-2)

  !! Compute V(D^Z) integrals
  !! integrals used for implicit contributions are just V(D^SA), so we already have them
  !! we only need PCMZMO(:,3) actually
  call PCM_grad_D2V(DZAO,PCMZMO,PCMZAO,.false.,.false.,.false.,.true.)

  !! For kap*Z (= RInt_generic + CI_KAP) contributions
  !! Similar to the fockgen. However, only the implicit contribution in the eigenenergy has to be
  !! evaluated here and the evaluation of the erfx term is not needed
  DZMO = Zero
  do iS=1,nSym
    if (nBas(iS) > 0) then
      jS = Mul(is,idSym)
      do iA=1,nAsh(is)
        do jA=1,nAsh(js)
          rd = DSCFMO(nA(is)+iA,nA(js)+jA) ! solvent density during SCF
          ip1 = nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
          ip2 = nBas(iS)*(nIsh(js)+jA-1)+ipMat(is,js)
          DZMO(ip2:ip2+nBas(iS)-1) = DZMO(ip2:ip2+nBas(iS)-1)+Rd*PCMZMO(ip1:ip1+nBas(iS)-1,3)
        end do
      end do
    end if
  end do

  if (idSym == 1) then
  !if (state_sym == 1) then !?
    do iS=1,nSym
      if (nBas(iS)*nIsh(iS) > 0) &
        DZMO(ipMat(iS,iS):ipMat(iS,iS)+nBas(iS)*nIsh(iS)-1) = DZMO(ipMat(iS,iS):ipMat(iS,iS)+nBas(iS)*nIsh(iS)-1)+ &
                                                              Two*PCMZMO(ipMat(iS,iS):ipMat(iS,iS)+nBas(iS)*nIsh(iS)-1,3)
    end do
  end if

  select case (iStpPCM)
    case (2)
      !! during Z-vector, compute kappa (DZAO)
      DZAO(:) = Zero
      do iS=1,nSym
        jS = Mul(iS,idSym)
        if (nBas(is)*nBas(jS) /= 0) &
          call DGeSub(DZMO(ipMat(iS,jS)),nBas(iS),'N',DZMO(ipMat(jS,iS)),nBas(jS),'T',DZAO(ipMat(iS,jS)),nBas(iS),nBas(iS),nBas(jS))
      end do
      FockOut(1:ntBsqr) = FockOut(1:ntBsqr)+Two*DZAO(1:ntBsqr)
    case (3)
      !! only for the connection term
      FockOut(1:ntBsqr) = FockOut(1:ntBsqr)+Two*DZMO(1:ntBsqr)
      !! No CI response
      return
  end select

  ! For ci*Z (= Kap_CI + Ci_Ci) contributions

  ! First, construct the contribution in the MO basis
  DZMO(:) = Zero
  do iS=1,nSym
    if (nBas(iS) > 0) then
      jS = Mul(iS,idSym)
      do jA=1,nAsh(js)
        !ip1 = nIsh(iS)+1+nBas(jS)*(nIsh(js)+jA-1)!+ipCM(is)
        ip1 = nIsh(iS)+1+nBas(iS)*(nIsh(js)+jA-1)+ipMat(iS,iS)-1
        ip2 = nIsh(iS)+1+nBas(iS)*(nIsh(js)+jA-1)+ipMat(iS,jS)-1
        DZMO(ip2:ip2+nAsh(iS)-1) = DZMO(ip2:ip2+nAsh(iS)-1)+PCMZMO(ip1:ip1+nAsh(iS)-1,3)
      end do
    end if
  end do

  ! Evaluate the CI derivative
  ! Note that ipCIOUT will be weighted with (2)W_SOLV
  call dswap_(nRoots,weight,1,W_SOLV,1)
  Weight(:) = Two*Weight(:)
  call CISigma_sa(0,state_sym,state_sym,DZMO,size(DZMO),rdum,1,rdum,1,ipCI,ipCIOUT,.false.)
  Weight(:) = Half*Weight(:)
  call dswap_(nRoots,weight,1,W_SOLV,1)

  ! consider the weight of the derivative
  if (DWSolv%DWZeta /= Zero) call DWder_MCLR(2,idSym,DZMO,size(DZMO),DZMO,size(DZMO),ISR%Ap)

end subroutine PCM_grad_TimesE2

!-----------------------------------------------------------------------

subroutine PCM_grad_CLag(mode,ipCI,ipCID)

  use MCLR_Data, only: ipMat, isNAC, nDens
  use MCLR_procedures, only: CISigma_sa
  use input_mclr, only: nAsh, nBas, ncsf, nIsh, nRoots, nSym, State_Sym, weight
  use ipPage, only: W
  use ISRotation, only: ISR

  integer(kind=iwp), intent(in) :: mode, ipCI, ipCID
  integer(kind=iwp) :: nconf1, idsym, ip1, is, ja, js
  real(kind=wp) :: rtmp(1)
  real(kind=wp), allocatable :: SCFcont(:)

  ! Evaluate derivatives of the SCF energy with respect to CI coefficients
  ! mode = 1: evaluate derivatives of the state-specific energy
  ! mode = 2: evaluate derivatives of the rotated state

  nconf1 = ncsf(State_Sym)
  rtmp(1) = Zero

  call mma_allocate(SCFcont,nDens,Label='SCFcont')

  ! CI derivative of the Lagrangian (dependent on the density that polarizes ASCs)

  if (PT2_solv .and. (mode == 1)) PCMSSMO(:,3) = PCMSSMO(:,3)+PCMPT2MO(:,3)
  if (mode == 2) PCMSSMO(:,3) = PCMSSMO(:,3)-PCMSSMOori(:,3) ! this has been evaluated with mode = 1

  ! CI derivative of the state-specific energy that is dependent on the definition of the solvation

  SCFcont = Zero !! contributions to the density that polarizes ASCs during SCF

  ! For NAC, we just need the implicit part of D^SS*V(e,SCF) in the eigenenergy
  do iS=1,nSym
    if (nAsh(iS) > 0) then
      jS = Mul(iS,state_sym)
      do jA=1,nAsh(js)
        !ip1 = nIsh(iS)+1+nBas(jS)*(nIsh(js)+jA-1)!+ipCM(is)
        ip1 = nIsh(iS)+1+nBas(iS)*(nIsh(js)+jA-1)+ipMat(iS,jS)-1 !+ipCM(is)
        !! def_solv = 1: subtract implicit D^SS*V(e,SCF)/2
        !! def_solv = 3: implicit V(e) in FIMO
        SCFcont(ip1:ip1+nAsh(iS)-1) = SCFcont(ip1:ip1+nAsh(iS)-1)+PCMSSMO(ip1:ip1+nAsh(iS)-1,3)
        !! def_solv = 1: explicit and implicit D^SS*V(e,SCF)/2
        !! def_solv = 3: explicit + implicit erfx
        if ((.not. isNAC) .and. (mode == 1)) SCFcont(ip1:ip1+nAsh(iS)-1) = SCFcont(ip1:ip1+nAsh(iS)-1)- &
                                                                           PCMSCFMO(ip1:ip1+nAsh(iS)-1,3)
        !! sum of these contributions are usually (RlxRoot = PCMRoot) zero,
        !! because D^SS is used for polarizing ASCs
      end do
    end if
  end do

  if (PT2_solv .and. (mode == 1)) PCMSSMO(:,3) = PCMSSMO(:,3)-PCMPT2MO(:,3)
  if (mode == 2) PCMSSMO(:,3) = PCMSSMO(:,3)+PCMSSMOori(:,3)

  !! state averaged quantities (ipCID is overwritten)
  !! it is not actually state-averaged; it is solvent density, so use W_SOLV
  call dswap_(nRoots,weight,1,W_SOLV,1)
  call CISigma_sa(0,state_sym,state_sym,SCFcont,nDens,rtmp,1,rtmp,1,ipCI,ipCID,.false.)
  call dswap_(nRoots,weight,1,W_SOLV,1)

  !! consider the weight of the derivative
  idsym = state_sym !?
  if (DWSolv%DWZeta /= Zero) call DWder_MCLR(2,idsym,SCFcont,nDens,SCFcont,nDens,ISR%RVec)

  !! State rotations are not needed, because CI vectors are obtained by diagonalization (?)
  W(ipCID)%A(1:nConf1*nRoots) = Two*W(ipCID)%A(1:nConf1*nRoots)

  call mma_deallocate(SCFcont)

end subroutine PCM_grad_CLag

!-----------------------------------------------------------------------

subroutine PCM_mod_ERASSCF(ERASSCF_)

  use input_mclr, only: nAsh, nIsh, nOrb, nRoots, nSym
  use MCLR_Data, only: ipMat, nA

  real(kind=wp), intent(inout) :: ERASSCF_(nRoots)
  integer(kind=iwp) :: iA, idSym, ip2, iRoot, iS, jA, jS
  real(kind=wp) :: ecorr

  ! Subtract (add) the surplus solvation energies from ERASSCF.
  ! ERASSCF at present is not eigenvalues of the Hamiltonian, because of some additional
  ! contributions from the solvent. The purpose of this subroutine is
  ! to subtract the surplus energies so that ERASSCF is an eigenvalue.

  idSym = 1

  do iRoot=1,nRoots
    !write(6,*) 'original erasscf = ',erasscf(iroot),erasscf_(iroot)
    ecorr = Zero

    do iS=1,nSym
      !jS = Mul(iS,state_sym)
      !! inactive
      do iA=1,nIsh(iS)
        ip2 = iA-1+nOrb(iS)*(iA-1)+ipMat(iS,iS)
        !! - D^SS * V(e,SS) / 2 or - D^SS * V(e,SA) / 2
        ecorr = ecorr+Two*pcmscfmo(ip2,3)
      end do
    end do
    ecorr = ecorr*Half

    !! active
    do iS=1,nSym
      !jS = Mul(iS,state_sym)
      jS = Mul(iS,idSym)
      do iA=1,nash(iS)
        do jA=1,nash(jS)
          ip2 = nIsh(iS)+iA-1+nOrb(iS)*(nIsh(jS)+jA-1)+ipMat(iS,jS) ! ipCM(iS)
          !! - D^SS * V(e,SA) / 2
          ecorr = ecorr+Half*DSCFMO(iA+nA(iS),jA+nA(jS))*pcmscfmo(ip2,3)
        end do
      end do
    end do

    ERASSCF_(iRoot) = ERASSCF_(iRoot)+ecorr
    !write(6,*) 'modified erasscf = ',erasscf_(iroot)
  end do

end subroutine PCM_mod_ERASSCF

!-----------------------------------------------------------------------

subroutine PCM_grad_PT2()

  use input_mclr, only: nBas, nSym, ntBsqr, ntBtri

  integer(kind=iwp) :: mData
  logical(kind=iwp) :: Found
  real(kind=wp), external :: ddot_

  ! Add implicit contributions due to the unrelaxed CASPT2 density

  !! check D1aoVar
  !! if it is non-zero, RFpert is on in CASPT2,
  !! so the response contributions should be evaluated.
  !! Otherwise, exit this subroutine
  call qpg_dArray('D1aoVar',Found,mData)
  if ((.not. Found) .or. (mData == 0)) then
    call warningMessage(2,'Gradients without RFpert in &CASPT2 is incorrect.')
    return
  end if

  !! PCMPT2MO here is used as a temporary array
  call Get_dArray('D1aoVar',PCMPT2MO,ntBtri)
  if (ddot_(ntBtri,PCMPT2MO,1,PCMPT2MO,1) <= 1.0e-10_wp) return

  PT2_solv = .true.

  call unfold(PCMPT2MO,ntBtri,DZAO,ntBsqr,nSym,nBas)
  call PCM_grad_D2V(DZAO,PCMPT2MO,PCMZAO,.false.,.false.,.false.,.true.)

end subroutine PCM_grad_PT2

end module PCM_grad

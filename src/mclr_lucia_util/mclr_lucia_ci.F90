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

!> The MCLR copy of Lucia as a CI library of the CI interface
!> (CI_Operations, which calls these routines directly): the CI operations
!> on top of the Lucia routines of the copy, which keep their data in
!> MCLR_Lucia_Data and Str_Info. The CI space comes with MCLR_Lucia_Init
!> (copies L_* in MCLR_Lucia_Data), the integrals with each operation; no
!> data of the programs or of the CI interface is used (iState only
!> identifies the CI state towards the CI interface: CI_ReOrd_Orig).
!> The vectors are CSF vectors of their natural length; the expansion in
!> determinants and its work space stay here.
module MCLR_Lucia_CI

use MCLR_Lucia_Data, only: L_nCSF, L_STSym
use Definitions, only: wp, iwp, u6

implicit none
private

logical(kind=iwp) :: Initialized = .false.

public :: MCLR_Lucia_Close, MCLR_Lucia_Conv, MCLR_Lucia_Dets, MCLR_Lucia_Diag, MCLR_Lucia_Init, MCLR_Lucia_PSpace, &
          MCLR_Lucia_Sigma, MCLR_Lucia_TDens

contains

!> Set up the determinant tables for the CI state (State_Proc of
!> CI_Operations) and publish the dimensions of the determinant space
subroutine MCLR_Lucia_Init(nSym,nOrb,nIsh,nAsh,nRs1,nRs2,nRs3,nHole1,nElec3,nActEl,iSpin,STSym,nCSF,TD,nDet)
  use MCLR_Lucia_Data, only: L_iSpin, L_nActEl, L_nAsh, L_nElec3, L_nHole1, L_nIsh, L_nOrb, L_nRs1, L_nRs2, L_nRs3, L_nSym, &
                             TD_Ints, XISPSM
  ! the CI space: symmetries; orbitals, inactive and active orbitals, RAS1,
  ! RAS2 and RAS3 orbitals per symmetry; holes in RAS1, electrons in RAS3;
  ! active electrons, multiplicity, symmetry of the CI state, CSFs per
  ! symmetry
  integer(kind=iwp), intent(in) :: nSym, nOrb(nSym), nIsh(nSym), nAsh(nSym), nRs1(nSym), nRs2(nSym), nRs3(nSym), nHole1, nElec3, &
                                   nActEl, iSpin, STSym, nCSF(nSym)
  ! TD: the two-body integrals of the run are in the square-pair layout
  logical(kind=iwp), intent(in) :: TD
  ! the dimensions of the determinant spaces per symmetry
  integer(kind=iwp), intent(out) :: nDet(nSym)
  L_nSym = nSym
  L_nOrb(:) = 0
  L_nIsh(:) = 0
  L_nAsh(:) = 0
  L_nRs1(:) = 0
  L_nRs2(:) = 0
  L_nRs3(:) = 0
  L_nCSF(:) = 0
  L_nOrb(1:nSym) = nOrb(:)
  L_nIsh(1:nSym) = nIsh(:)
  L_nAsh(1:nSym) = nAsh(:)
  L_nRs1(1:nSym) = nRs1(:)
  L_nRs2(1:nSym) = nRs2(:)
  L_nRs3(1:nSym) = nRs3(:)
  L_nHole1 = nHole1
  L_nElec3 = nElec3
  L_nActEl = nActEl
  L_iSpin = iSpin
  L_STSym = STSym
  L_nCSF(1:nSym) = nCSF(:)
  TD_Ints = TD
  call DetCtl()
  call InCsfSD(STSym,STSym)
  nDet(:) = nint(XISPSM(1:nSym,1))
  Initialized = .true.
end subroutine MCLR_Lucia_Init

!> Free the determinant tables (State_Proc of CI_Operations)
subroutine MCLR_Lucia_Close(iState)
  use Str_Info, only: CFTP, CNSM, DFTP, DTOC
  use MCLR_Lucia_Data, only: LuCSF2SD, pInt1
  use stdalloc, only: mma_deallocate
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: i
#include "macros.fh"
  unused_var(iState)
  if (Initialized) then
    call DaClos(LuCSF2SD)
    call FreeStr()
    call mma_deallocate(DTOC)
    call mma_deallocate(CFTP)
    call mma_deallocate(DFTP)
    call mma_deallocate(pInt1)
  end if
  do i=1,size(CNSM)
    call mma_deallocate(CNSM(i)%ICONF,safe='*')
    call mma_deallocate(CNSM(i)%ICTS,safe='*')
  end do
  Initialized = .false.
end subroutine MCLR_Lucia_Close

!> The CSF <-> SD tables (InCsfSD) for the symmetries of the vectors of an
!> operation: those of the CI state are always there, one other symmetry
!> is cached
subroutine Use_Syms(iSym1,iSym2)
  integer(kind=iwp), intent(in) :: iSym1, iSym2
  integer(kind=iwp) :: iRef
  iRef = L_STSym
  if ((iSym1 /= iRef) .and. (iSym2 /= iRef) .and. (iSym1 /= iSym2)) then
    write(u6,'(A,3I3)') ' MCLR copy of Lucia: two symmetries other than that of the CI state:',iSym1,iSym2,iRef
    call Abend()
  end if
  if (iSym1 /= iRef) then
    call InCsfSD(iSym1,iRef)
  else if (iSym2 /= iRef) then
    call InCsfSD(iSym2,iRef)
  end if
end subroutine Use_Syms

!> sigma = H C (Sigma_Proc of CI_Operations)
subroutine MCLR_Lucia_Sigma(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
  use Symmetry_Info, only: Mul
  use MCLR_Lucia_Data, only: i12, iRefSM, iST, KAIN1, KINT2, pInt1, Square_ => Square, XISPSM
  use CandS, only: ICSM, ISSM
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nInt1, nInt2, nC, nS
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  logical(kind=iwp), intent(in) :: Have_2_el, Square
  real(kind=wp), intent(in) :: C(nC)
  real(kind=wp), intent(out) :: S(nS)
  integer(kind=iwp) :: iOp, kic(2), nDet
  real(kind=wp), allocatable :: CDet(:), SDet(:)
#include "macros.fh"
  unused_var(iState)

  S(:) = Zero
  call Use_Syms(iCSym,iSSym)
  ! the integrals, as in CISigma
  KAIN1(1:nInt1) => Int1(:)
  KINT2(1:nInt2) => Int2(:)
  irefsm = iCSym
  if (Have_2_el) then
    i12 = 2
  else
    i12 = 1
  end if
  ! symmetries of sigma and C (kic = 1: the symmetry of the CI state)
  iSSM = iSSym
  kic(2) = 2
  if (iSSM == L_STSym) kic(2) = 1
  iCSM = iCSym
  kic(1) = 2
  if (iCSM == L_STSym) kic(1) = 1
  nDet = nint(max(XISPSM(iSSym,1),XISPSM(iCSym,1)))
  nDet = max(nDet,L_nCSF(iCSym),L_nCSF(iSSym))
  if (nDet == 0) then
    nullify(KAIN1,KINT2)
    return
  end if
  ! symmetry of the operator: the one-body blocks
  iOp = Mul(iCSM,iSSM)
  call Int1_Offsets(iOp,pInt1)
  ! singlet/triplet operator, symmetric/square integrals
  ist = iOpSpin+1
  Square_ = Square

  ! C and sigma are expanded in determinants in place: work space of the
  ! determinant length
  call mma_allocate(CDet,nDet,Label='CDet')
  call mma_allocate(SDet,nDet,Label='SDet')
  CDet(1:nC) = C(:)
  CDet(nC+1:) = Zero
  call SigmaVec(CDet,SDet,kic)
  S(:) = SDet(1:nS)
  call mma_deallocate(SDet)
  call mma_deallocate(CDet)
  Square_ = .false.
  nullify(KAIN1,KINT2)

end subroutine MCLR_Lucia_Sigma

!> Transition densities between L and R (TDens_Proc of CI_Operations)
subroutine MCLR_Lucia_TDens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  use CandS, only: ICSM, ISSM
  use MCLR_Lucia_Data, only: NOCSF, XISPSM
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  real(kind=wp), allocatable :: LDet(:), RDet(:)

#include "macros.fh"
  unused_var(iState)
  if (iOpSpin /= 0) then
    write(u6,'(A)') ' MCLR_Lucia_TDens: the MCLR copy of Lucia has spin-free densities only'
    call Abend()
  end if
  call Use_Syms(iLSym,iRSym)
  ! the determinant expansions (CSF2SD of the state of the run); without
  ! CSFs (NOCSF /= 0) the vectors are in determinants already
  call mma_allocate(LDet,max(L_nCSF(iLSym),nint(XISPSM(iLSym,1))),Label='LDet')
  call mma_allocate(RDet,max(L_nCSF(iRSym),nint(XISPSM(iRSym,1))),Label='RDet')
  if (NOCSF == 0) then
    call CSF2SD(L,LDet,iLSym)
    call CSF2SD(R,RDet,iRSym)
  else
    LDet(1:nL) = L(:)
    RDet(1:nR) = R(:)
  end if
  ! Densi2_MCLR takes L as the sigma (ISSM) and R as the C vector (ICSM)
  iCSM = iRSym
  iSSM = iLSym
  call Densi2_MCLR(I12,Rho1,Rho2,LDet,RDet,0,0,0,n1,n2)
  call mma_deallocate(RDet)
  call mma_deallocate(LDet)

end subroutine MCLR_Lucia_TDens

!> The CSF diagonal (CI_Diag of CI_Operations), from the integrals of the
!> call (a totally symmetric operator)
subroutine MCLR_Lucia_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
  use CI_Op_Defs, only: CI_Diag_Cnf_Avg
  use Str_Info, only: CNSM
  use MCLR_Lucia_Data, only: ICISTR, KAIN1, KINT2, NCNFTP, NCSFTP, NDPCNT, NOCSF, NTYP, pInt1, XISPSM
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nInt1, nInt2, n
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  real(kind=wp), intent(in) :: ECore
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp) :: i, iAMCmp, iSM(1), iSPC(1), nSD, nSpc
  real(kind=wp), allocatable :: DSD(:)
#include "macros.fh"
  unused_var(iState)

  if (Mode /= CI_Diag_Cnf_Avg) then
    write(u6,'(A)') ' MCLR_Lucia_Diag: the MCLR copy of Lucia has the configuration-averaged diagonal only'
    call Abend()
  end if
  call Use_Syms(iSym,iSym)
  nSpc = 1
  iSPC(1) = 1
  iSM(1) = iSym
  iAMCmp = 0
  ICISTR = 1
  i = 2
  if (iSym == L_STSym) i = 1
  nSD = max(L_nCSF(iSym),nint(XISPSM(iSym,1)))
  call mma_allocate(DSD,nSD,Label='DSD')
  KAIN1(1:nInt1) => Int1(:)
  KINT2(1:nInt2) => Int2(:)
  call Int1_Offsets(1,pInt1)
  call IntDia(DSD,nSpc,iSPC,iSM,iAMCmp,ECore)
  nullify(KAIN1,KINT2)
  if (NOCSF /= 1) then
    call CSDIAG_MCLR(Diag,DSD,NCNFTP(1,iSym),NTYP,CNSM(i)%ICTS,NDPCNT,NCSFTP)
  else
    Diag(:) = DSD(1:n)
  end if
  call mma_deallocate(DSD)

end subroutine MCLR_Lucia_Diag

!> The explicit P space (PSpace_Proc of CI_Operations): H0MAT_MCLR with
!> a P2 space of at most nPMax CSFs (no P1 and Q spaces)
subroutine MCLR_Lucia_PSpace(iState,iSym,nInt1,Int1,nInt2,Int2,n,Diag,nPMax,nP,IdxP,nHP,HP)
  use Index_Functions, only: nTri_Elem
  use MCLR_Lucia_Data, only: IDC, iRefSM, KAIN1, KINT2, NACOB, NAELCI, NBELCI, NCSFTP, NDPCNT, NOCOB, NTYP, pInt1, PSSIGN, XISPSM
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nInt1, nInt2, n, nPMax, nHP
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  real(kind=wp), intent(in) :: Diag(n)
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  integer(kind=iwp) :: ieaw, lH0Scr, lVec2, MxCSFC, MxDTFC, nActEl
  integer(kind=iwp), allocatable :: SBCnf(:)
  real(kind=wp), allocatable :: H0Scr(:), Vec2(:)

#include "macros.fh"
  unused_var(iState)
  call Use_Syms(iSym,iSym)
  MxCSFC = max(0,maxval(NCSFTP(1:NTYP)))
  MxDTFC = max(0,maxval(NDPCNT(1:NTYP)))
  nActEl = NAELCI(1)+NBELCI(1)
  lH0Scr = max(6*nPMax,4*nPMax+4*NOCOB)
  lVec2 = MxCSFC**2+2*MxDTFC+2*MxDTFC**2+MxDTFC*nActEl
  lVec2 = max(lVec2,nint(XISPSM(iSym,1)))
  ! (iRefSM is the symmetry of the last sigma vector, as before)
  if (iSym == iRefSM) then
    ieaw = 1
  else
    ieaw = 2
  end if
  if (nHP < nTri_Elem(nPMax)) call Abend()
  ! the integrals of the (totally symmetric) Hamiltonian, as for sigma
  KAIN1(1:nInt1) => Int1(:)
  KINT2(1:nInt2) => Int2(:)
  call Int1_Offsets(1,pInt1)
  call mma_allocate(SBCnf,nPMax,Label='SBCnf')
  call mma_allocate(H0Scr,lH0Scr,Label='H0Scr')
  call mma_allocate(Vec2,lVec2,Label='Vec2')
  call H0MAT_MCLR(HP,IdxP,SBCnf,0,nPMax,0,NACOB,nP,iSym,IDC,PSSIGN,Diag,Vec2,H0Scr,ieaw)
  call mma_deallocate(Vec2)
  call mma_deallocate(H0Scr)
  call mma_deallocate(SBCnf)
  nullify(KAIN1,KINT2)

end subroutine MCLR_Lucia_PSpace

!> Determinants (of a DMRG reconstruction) into the CSF basis (Dets_Proc
!> of CI_Operations)
subroutine MCLR_Lucia_Dets(iState,iSym,nD,DetVec,IdxD,n,CI)
  use Str_Info, only: DTOC
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nD, IdxD(nD), n
  real(kind=wp), intent(in) :: DetVec(nD)
  real(kind=wp), intent(out) :: CI(n)
  real(kind=wp), allocatable :: DW(:)
#include "macros.fh"
  unused_var(iState)
  call Use_Syms(iSym,iSym)
  call mma_allocate(DW,nD,Label='DW')
  DW(:) = DetVec(:)
  call CSDTVC_dmrg(CI,DW,DTOC,IdxD,iSym,1)
  call mma_deallocate(DW)
end subroutine MCLR_Lucia_Dets

!> Stored basis <-> the CSF basis of the MCLR copy of Lucia (Conv_Proc
!> of CI_Operations; the reordering of CI_ReOrd_Orig, via sg2symg)
subroutine MCLR_Lucia_Conv(iState,iSym,iMode,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, iMode, n
  real(kind=wp), intent(inout) :: CI(n)

#include "macros.fh"
  unused_var(iState)
  call Use_Syms(iSym,iSym)
  call SG2SymG(CI,n,iMode,iSym)

end subroutine MCLR_Lucia_Conv

!> The offsets of the one-electron integral blocks (iS,Mul(iS,iOp)) of an
!> operator of symmetry iOp: square nOrb(iS) x nOrb(Mul(iS,iOp)) blocks
!> one after the other in the order of iS (in MCLR: ipCM, ipMat)
subroutine Int1_Offsets(iOp,pInt1)
  use Symmetry_Info, only: Mul
  use MCLR_Lucia_Data, only: nOrb => L_nOrb, nSym => L_nSym
  integer(kind=iwp), intent(in) :: iOp
  integer(kind=iwp), intent(out) :: pInt1(:)
  integer(kind=iwp) :: iS
  pInt1(:) = 0
  pInt1(1) = 1
  do iS=2,nSym
    pInt1(iS) = pInt1(iS-1)+nOrb(iS-1)*nOrb(Mul(iS-1,iOp))
  end do
end subroutine Int1_Offsets

end module MCLR_Lucia_CI

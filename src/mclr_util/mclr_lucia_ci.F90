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
!> (CI_Operations): the CI operations of MCLR on top of the Lucia
!> routines of mclr_util, which keep their data in MCLR_Data.
!> MCLR_Lucia_CI_Setup selects this library and registers the routines.
!> The vectors are CSF vectors of their natural length; the expansion in
!> determinants and its work space stay here.
module MCLR_Lucia_CI

use Definitions, only: wp, iwp, u6

implicit none
private

logical(kind=iwp) :: Initialized = .false.

public :: MCLR_Lucia_CI_Setup

contains

!> Select the MCLR copy of Lucia and register its routines
subroutine MCLR_Lucia_CI_Setup()
  use CI_Operations, only: CI_Lib_MCLR_Lucia, CI_Lib_Register, CI_Lib_Select
  call CI_Lib_Select(CI_Lib_MCLR_Lucia)
  call CI_Lib_Register(MCLR_Lucia_Init,MCLR_Lucia_Close,MCLR_Lucia_Sigma,MCLR_Lucia_TDens,MCLR_Lucia_Diag,MCLR_Lucia_PSpace, &
                       MCLR_Lucia_Conv,MCLR_Lucia_Dets)
end subroutine MCLR_Lucia_CI_Setup

!> Set up the determinant tables for the CI state (State_Proc of
!> CI_Operations) and publish the dimensions of the determinant space
subroutine MCLR_Lucia_Init(iState)
  use CI_Info, only: CI_St
  use MCLR_Data, only: XISPSM
  use Orbital_Spaces, only: nSym
  integer(kind=iwp), intent(in) :: iState
  call DetCtl()
  call InCsfSD(CI_St(iState)%STSym,CI_St(iState)%STSym)
  CI_St(iState)%nDet(1:nSym) = nint(XISPSM(1:nSym,1))
  Initialized = .true.
end subroutine MCLR_Lucia_Init

!> Free the determinant tables (State_Proc of CI_Operations)
subroutine MCLR_Lucia_Close(iState)
  use Str_Info, only: CFTP, CNSM, DFTP, DTOC
  use MCLR_Data, only: pInt1
  use stdalloc, only: mma_deallocate
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: i
#include "macros.fh"
  unused_var(iState)
  if (Initialized) then
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
subroutine Use_Syms(iState,iSym1,iSym2)
  use CI_Info, only: CI_St
  integer(kind=iwp), intent(in) :: iState, iSym1, iSym2
  integer(kind=iwp) :: iRef
  iRef = CI_St(iState)%STSym
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
  use CI_Info, only: CI_St
  use Symmetry_Info, only: Mul
  use MCLR_Data, only: i12, ipCM, ipMat, iRefSM, iST, KAIN1, KINT2, pInt1, Square_ => Square, XISPSM
  use CandS, only: ICSM, ISSM
  use Orbital_Spaces, only: nSym
  use input_mclr, only: nCSF
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nInt1, nInt2, nC, nS
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  logical(kind=iwp), intent(in) :: Have_2_el, Square
  real(kind=wp), intent(in) :: C(nC)
  real(kind=wp), intent(out) :: S(nS)
  integer(kind=iwp) :: iOp, iS, kic(2), nDet
  real(kind=wp), allocatable :: CDet(:), SDet(:)

  S(:) = Zero
  call Use_Syms(iState,iCSym,iSSym)
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
  if (iSSM == CI_St(iState)%STSym) kic(2) = 1
  iCSM = iCSym
  kic(1) = 2
  if (iCSM == CI_St(iState)%STSym) kic(1) = 1
  nDet = nint(max(XISPSM(iSSym,1),XISPSM(iCSym,1)))
  nDet = max(nDet,nCSF(iCSym),nCSF(iSSym))
  if (nDet == 0) then
    nullify(KAIN1,KINT2)
    return
  end if
  ! symmetry of the operator: the one-body blocks
  iOp = Mul(iCSM,iSSM)
  if (iOp == 1) then
    pInt1(1:nSym) = ipCM(1:nSym)
  else
    do iS=1,nSym
      pInt1(iS) = ipMat(iS,Mul(iS,iOp))
    end do
  end if
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
  use MCLR_Data, only: NOCSF, XISPSM
  use input_mclr, only: nCSF
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  real(kind=wp), allocatable :: LDet(:), RDet(:)

#include "macros.fh"
  if (iOpSpin /= 0) then
    write(u6,'(A)') ' MCLR_Lucia_TDens: the MCLR copy of Lucia has spin-free densities only'
    call Abend()
  end if
  call Use_Syms(iState,iLSym,iRSym)
  ! the determinant expansions (CSF2SD of the state of the run); without
  ! CSFs (NOCSF /= 0) the vectors are in determinants already
  call mma_allocate(LDet,max(nCSF(iLSym),nint(XISPSM(iLSym,1))),Label='LDet')
  call mma_allocate(RDet,max(nCSF(iRSym),nint(XISPSM(iRSym,1))),Label='RDet')
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

!> The CSF diagonal (Diag_Proc of CI_Operations).
!> IntDia reads the active integrals from MCLR (FIMO and the two-body
!> integrals of MCLR_Data), which are those passed by MCLR; Int1 and Int2
!> are not referenced here.
subroutine MCLR_Lucia_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
  use CI_Info, only: CI_St
  use CI_Operations, only: CI_Diag_Cnf_Avg
  use Str_Info, only: CNSM
  use MCLR_Data, only: ICISTR, NCNFTP, NCSFTP, NDPCNT, NOCSF, NTYP, XISPSM
  use input_mclr, only: nCSF
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nInt1, nInt2, n
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  real(kind=wp), intent(in) :: ECore
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp) :: i, iAMCmp, iSM(1), iSPC(1), nSD, nSpc
  real(kind=wp), allocatable :: DSD(:)

#include "macros.fh"
  unused_var(Int1)
  unused_var(Int2)

  if (Mode /= CI_Diag_Cnf_Avg) then
    write(u6,'(A)') ' MCLR_Lucia_Diag: the MCLR copy of Lucia has the configuration-averaged diagonal only'
    call Abend()
  end if
  call Use_Syms(iState,iSym,iSym)
  nSpc = 1
  iSPC(1) = 1
  iSM(1) = iSym
  iAMCmp = 0
  ICISTR = 1
  i = 2
  if (iSym == CI_St(iState)%STSym) i = 1
  nSD = max(nCSF(iSym),nint(XISPSM(iSym,1)))
  call mma_allocate(DSD,nSD,Label='DSD')
  call IntDia(DSD,nSpc,iSPC,iSM,iAMCmp,ECore)
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
  use MCLR_Data, only: IDC, ipCM, iRefSM, KAIN1, KINT2, NACOB, NAELCI, NBELCI, NCSFTP, NDPCNT, NOCOB, NTYP, pInt1, PSSIGN, &
                       XISPSM
  use Orbital_Spaces, only: nSym
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
  call Use_Syms(iState,iSym,iSym)
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
  pInt1(1:nSym) = ipCM(1:nSym)
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
  call Use_Syms(iState,iSym,iSym)
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
  call Use_Syms(iState,iSym,iSym)
  call SG2SymG(CI,n,iMode,iSym)

end subroutine MCLR_Lucia_Conv

end module MCLR_Lucia_CI

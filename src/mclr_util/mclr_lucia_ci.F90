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
  call CI_Lib_Register(MCLR_Lucia_Init,MCLR_Lucia_Close,MCLR_Lucia_Sigma,MCLR_Lucia_TDens,MCLR_Lucia_Diag,MCLR_Lucia_Conv)
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
  integer(kind=iwp), intent(in) :: iState
#include "macros.fh"
  unused_var(iState)
  if (Initialized) call FreeStr()
  Initialized = .false.
end subroutine MCLR_Lucia_Close

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
  unused_var(iState)

  if (iOpSpin /= 0) then
    write(u6,'(A)') ' MCLR_Lucia_TDens: the MCLR copy of Lucia has spin-free densities only'
    call Abend()
  end if
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

!> Stored basis <-> the CSF basis of the MCLR copy of Lucia (Conv_Proc
!> of CI_Operations; the reordering of CI_ReOrd_Orig, via sg2symg)
subroutine MCLR_Lucia_Conv(iState,iSym,iMode,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, iMode, n
  real(kind=wp), intent(inout) :: CI(n)

#include "macros.fh"
  unused_var(iState)

  call SG2SymG(CI,n,iMode,iSym)

end subroutine MCLR_Lucia_Conv

end module MCLR_Lucia_CI

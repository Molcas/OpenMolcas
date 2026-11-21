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

!> Generic CI operations of the CI interface, for any program, and the
!> choice of the CI library that performs them.
!>
!> Operations (all for the CI state iState, CI_St(iState)):
!>  - CI_Lib_Initialize, CI_Lib_Close: set up and free the selected
!>    library for the CI state (CI_Lib_Close is safe to call in any case).
!>    The library publishes its dimensions in CI_St(iState) (nCSF; nDet: the
!>    dimension of its determinant space per symmetry).
!>  - CI_Sigma_Op: sigma = H C for one CI vector. H is a Hermitian
!>    operator of symmetry Mul(iCSym,iSSym) (C of symmetry iCSym, sigma of
!>    symmetry iSSym); iOpSpin (the spin of the operator, not the
!>    multiplicity of the state) = 0: singlet E_pq = a+_pa a_qa + a+_pb a_qb,
!>    iOpSpin = 1: a+_pa a_qa - a+_pb a_qb (= sqrt(2) T^{1,0}_pq, one-body
!>    only). S is always overwritten (no accumulation).
!>    Int1: one-body integrals over all orbitals, in square symmetry blocks
!>    (block (iS,jS), jS = Mul(iS,iOp), nOrb(iS) x nOrb(jS), column-major,
!>    one after the other in the order of iS); only the active part is
!>    used. Int2: two-body integrals over the active orbitals,
!>    symmetry-ordered active indices (all active orbitals of irrep 1, then
!>    irrep 2, ...; RAS1, RAS2, RAS3 inside each irrep), (ij|kl) at
!>    iTri(iTri(i,j),iTri(k,l)). Have_2_el = .false.: one-body operator
!>    only. Square = .true.: the operator is not symmetric (time-dependent
!>    response): Int1 is a general (non-symmetric) matrix and (ij|kl) is at
!>    iTri(ij,kl) with the square pair index ij = (j-1)*nAct+i.
!>  - CI_Trans_Dens: transition densities between the bra L (symmetry
!>    iLSym) and the ket R (symmetry iRSym); iOpSpin = 0: spin-free,
!>    iOpSpin = 1: spin densities (a+_a a_a - a+_b a_b). rho1(i,j) =
!>    <L|E_ij|R>, full square nAct x nAct; for I12 = 2 also rho2(iTri(ij,kl))
!>    = <L|E_ij E_kl - delta_jk E_il|R> with the square pair index
!>    ij = (j-1)*nAct+i (packed over ij >= kl); symmetry-ordered active
!>    indices. Spin-free normalisation: sum_i rho1(i,i) = N_act*<L|R>.
!>    Rho2 is not referenced for I12 = 1.
!>  - CI_Diag: the diagonal of H + ECore (ECore added to every element) in
!>    the CSF basis of symmetry iSym, for the symmetric active Hamiltonian;
!>    Int1 and Int2 as for the sigma with iOp = 1 and Square = .false..
!>    Mode = CI_Diag_Cnf_Avg: averaged over the CSFs of each configuration
!>    (Lucia's CSDIAG); CI_Diag_Exact: the exact CSF diagonal.
!>  - CI_PSpace: the explicit part of the Hamiltonian for the
!>    preconditioner: a P space of at most nPMax CSFs of symmetry iSym,
!>    selected by the library from the diagonal Diag (lowest first, complete
!>    configurations), its nP CSF indices IdxP and the block of H (plus
!>    ECore on the diagonal, as Diag) packed as a lower triangle, HP(iTri(i,j))
!>    = <IdxP(i)|H|IdxP(j)>; Int1 and Int2 as for CI_Diag.
!>  - CI_Lib_Dets_To_Basis: a vector of nD determinants DetVec (IdxD: their
!>    indices in the determinant space of the library) into the CSF basis
!>    of symmetry iSym (used for CI vectors reconstructed from DMRG).
!>  - CI_Lib_File_To_Basis, CI_Lib_Basis_To_File: a CI vector of symmetry
!>    iSym between the stored basis (JobIph) and the CSF basis of the
!>    selected library, in place.
!> All vectors passed to the gateway are in the CSF basis of the selected
!> library, with their natural length (nCSF of their symmetry); each
!> library keeps its own internal basis, ordering and work space.
!>
!> Libraries, all with every operation except the triplet operator
!> (iOpSpin = 1, the MCLR copy of Lucia only):
!> - the MCLR copy of Lucia (mclr_lucia_util, MCLR_Lucia_CI; its own data
!>   in MCLR_Lucia_Data and Str_Info);
!> - SGUGA (CAS/RAS; the split graph has to be set up with the coupling
!>   tables, CI_Lib_Graph_Mode() = 'FULL'). It takes the integrals in
!>   symmetry order and maps them to the levels of its graph itself; its
!>   densities come in level order and are mapped back here
!>   (SG_Level_Map);
!> - Faroald (CAS/RAS; the string-driven RAS CI in the symmetry of the CI
!>   state, the full determinant space for other symmetries and
!>   non-Hermitian operators).
!> The orbital spaces are taken from Orbital_Spaces (nSym, nOrb, nIsh,
!> nAsh) and CI_Input, the active orbitals in symmetry order; the
!> integrals come with each operation.
module CI_Operations

use CI_Info, only: CI_St
use CI_Op_Defs, only: CI_Diag_Cnf_Avg, CI_Diag_Exact, CI_Lib_Faroald, CI_Lib_Lucia, CI_Lib_MCLR_Lucia, CI_Lib_None, CI_Lib_SGUGA, &
                      CI_To_Basis, CI_To_File
use MCLR_Lucia_CI, only: MCLR_Lucia_Close, MCLR_Lucia_Conv, MCLR_Lucia_Dets, MCLR_Lucia_Diag, MCLR_Lucia_Init, MCLR_Lucia_PSpace, &
                         MCLR_Lucia_Sigma, MCLR_Lucia_TDens
use Definitions, only: wp, iwp, u6

implicit none
private

!> The selected library (changed only by CI_Lib_Select)
integer(kind=iwp), protected :: CI_Lib_Op = CI_Lib_None


! Faroald: its orbitals are the levels of the split graph (by RAS
! partition, within each by symmetry); Faro_Lev(iLev) is the
! symmetry-ordered active orbital of level iLev. Faro_Sym is the symmetry
! of the CI state, the space of the string-driven RAS CI (faroald_ras);
! citrans holds the CSF space of one symmetry at a time, Faro_Cur, set up
! again from the stored description (Faro_OrbSym, ...) when another one
! is needed (Faro_Set_Sym).
integer(kind=iwp) :: Faro_Cur = 0, Faro_MaxEl(3) = 0, Faro_MinEl(3) = 0, Faro_nAct = 0, Faro_nPrt = 0, Faro_Sym = 0
integer(kind=iwp), allocatable :: Faro_Lev(:), Faro_OrbRas(:), Faro_OrbSym(:)

public :: CI_Diag, CI_Diag_Cnf_Avg, CI_Diag_Exact, CI_Lib_Basis_To_File, CI_Lib_Close, CI_Lib_Dets_To_Basis, CI_Lib_Faroald, &
          CI_Lib_File_To_Basis, CI_Lib_Graph_Mode, CI_Lib_Initialize, CI_Lib_Lucia, CI_Lib_MCLR_Lucia, CI_Lib_None, CI_Lib_Op, &
          CI_Lib_Select, CI_Lib_SGUGA, CI_PSpace, CI_Sigma_Op, CI_To_Basis, CI_To_File, CI_Trans_Dens

contains

!> Select the CI library for MCLR (CI_Lib_*)
subroutine CI_Lib_Select(iLib)
  integer(kind=iwp), intent(in) :: iLib
  select case (iLib)
    case (CI_Lib_MCLR_Lucia,CI_Lib_SGUGA,CI_Lib_Faroald)
      CI_Lib_Op = iLib
    case (CI_Lib_Lucia)
      write(u6,'(A,I2,A)') ' CI_Lib_Select: CI library',iLib,' is not available for MCLR yet'
      call Abend()
    case default
      write(u6,'(A,I4)') ' CI_Lib_Select: unknown CI library',iLib
      call Abend()
  end select
end subroutine CI_Lib_Select

!> Set up the selected library for the CI state. TD (default .false.): the
!> two-body integrals of all later operations are in the square-pair
!> layout (time-dependent response of MCLR)
subroutine CI_Lib_Initialize(iState,TD)
  integer(kind=iwp), intent(in) :: iState
  logical(kind=iwp), intent(in), optional :: TD
  logical(kind=iwp) :: TD_
  if (CI_Lib_Op == CI_Lib_SGUGA) then
    ! the split graph is the state of the CI interface (set up by the
    ! program with CI_Lib_Graph_Mode()); no determinants
    CI_St(iState)%nDet(:) = CI_St(iState)%nCSF(:)
    return
  end if
  if (CI_Lib_Op == CI_Lib_Faroald) then
    call Faro_Init(iState)
    return
  end if
  call Check_Lib('CI_Lib_Initialize')
  TD_ = .false.
  if (present(TD)) TD_ = TD
  call MCLR_Lucia_Init(iState,TD_)
end subroutine CI_Lib_Initialize

!> Free the library of the CI state and forget the selection (safe to call
!> in any case)
subroutine CI_Lib_Close(iState)
  use SGUGA_Sigma, only: SG_Sigma_Free
  use SGUGA_Densities, only: SG_Densities_Free
  integer(kind=iwp), intent(in) :: iState
  if (CI_Lib_Op == CI_Lib_MCLR_Lucia) call MCLR_Lucia_Close(iState)
  ! the work space of SGUGA (safe to call in any case, as in CI_Close)
  if (CI_Lib_Op == CI_Lib_SGUGA) then
    call SG_Sigma_Free()
    call SG_Densities_Free()
  end if
  if (CI_Lib_Op == CI_Lib_Faroald) call Faro_Close()
  CI_Lib_Op = CI_Lib_None
end subroutine CI_Lib_Close

subroutine CI_Sigma_Op(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nInt1, nInt2, nC, nS
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  logical(kind=iwp), intent(in) :: Have_2_el, Square
  real(kind=wp), intent(in) :: C(nC)
  real(kind=wp), intent(out) :: S(nS)
  if ((iOpSpin /= 0) .and. (iOpSpin /= 1)) then
    write(u6,'(A,I4)') ' CI_Sigma_Op: invalid operator spin iOpSpin =',iOpSpin
    call Abend()
  end if
  if ((iOpSpin == 1) .and. Have_2_el) then
    write(u6,'(A)') ' CI_Sigma_Op: the operator with iOpSpin = 1 is one-body only'
    call Abend()
  end if
  if (CI_Lib_Op == CI_Lib_SGUGA) then
    if (iOpSpin /= 0) call Not_Yet('CI_Sigma_Op','triplet operators')
    if (Square) then
      call SG_Sigma_Sq(iState,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,nC,C,nS,S)
    else
      call SG_Sigma_Op(iState,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,nC,C,nS,S)
    end if
    return
  end if
  if (CI_Lib_Op == CI_Lib_Faroald) then
    if (iOpSpin /= 0) call Not_Yet('CI_Sigma_Op','triplet operators')
    call Faro_Sigma_Op(iState,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,nC,C,nS,S,Square)
    return
  end if
  call Check_Lib('CI_Sigma_Op')
  call MCLR_Lucia_Sigma(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
end subroutine CI_Sigma_Op

subroutine CI_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  if ((iOpSpin /= 0) .and. (iOpSpin /= 1)) then
    write(u6,'(A,I4)') ' CI_Trans_Dens: invalid operator spin iOpSpin =',iOpSpin
    call Abend()
  end if
  if (CI_Lib_Op == CI_Lib_SGUGA) then
    if (iOpSpin /= 0) call Not_Yet('CI_Trans_Dens','spin densities')
    call SG_Trans_Dens(iState,I12,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
    return
  end if
  if (CI_Lib_Op == CI_Lib_Faroald) then
    if (iOpSpin /= 0) call Not_Yet('CI_Trans_Dens','spin densities')
    call Faro_Trans_Dens(iState,I12,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
    return
  end if
  call Check_Lib('CI_Trans_Dens')
  call MCLR_Lucia_TDens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
end subroutine CI_Trans_Dens

subroutine CI_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nInt1, nInt2, n
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  real(kind=wp), intent(in) :: ECore
  real(kind=wp), intent(out) :: Diag(n)
  if ((Mode /= CI_Diag_Cnf_Avg) .and. (Mode /= CI_Diag_Exact)) then
    write(u6,'(A,I4)') ' CI_Diag: invalid Mode =',Mode
    call Abend()
  end if
  if (CI_Lib_Op == CI_Lib_SGUGA) then
    call SG_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
    return
  end if
  if (CI_Lib_Op == CI_Lib_Faroald) then
    call Faro_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
    return
  end if
  call Check_Lib('CI_Diag')
  call MCLR_Lucia_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
end subroutine CI_Diag

subroutine CI_PSpace(iState,iSym,nInt1,Int1,nInt2,Int2,n,Diag,nPMax,nP,IdxP,nHP,HP)
  integer(kind=iwp), intent(in) :: iState, iSym, nInt1, nInt2, n, nPMax, nHP
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  real(kind=wp), intent(in) :: Diag(n)
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  if (CI_Lib_Op == CI_Lib_SGUGA) then
    call SG_PSpace(iState,iSym,nInt1,Int1,nInt2,Int2,nPMax,nP,IdxP,nHP,HP)
    return
  end if
  if (CI_Lib_Op == CI_Lib_Faroald) then
    call Faro_PSpace(iState,iSym,nInt1,Int1,nInt2,Int2,n,nPMax,nP,IdxP,nHP,HP)
    return
  end if
  call Check_Lib('CI_PSpace')
  call MCLR_Lucia_PSpace(iState,iSym,nInt1,Int1,nInt2,Int2,n,Diag,nPMax,nP,IdxP,nHP,HP)
end subroutine CI_PSpace

subroutine CI_Lib_Dets_To_Basis(iState,iSym,nD,DetVec,IdxD,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, nD, IdxD(nD), n
  real(kind=wp), intent(in) :: DetVec(nD)
  real(kind=wp), intent(out) :: CI(n)
  if ((CI_Lib_Op == CI_Lib_SGUGA) .or. (CI_Lib_Op == CI_Lib_Faroald)) &
    call Not_Yet('CI_Lib_Dets_To_Basis','determinant vectors (DMRG)')
  call Check_Lib('CI_Lib_Dets_To_Basis')
  call MCLR_Lucia_Dets(iState,iSym,nD,DetVec,IdxD,n,CI)
end subroutine CI_Lib_Dets_To_Basis

subroutine CI_Lib_File_To_Basis(iState,iSym,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, n
  real(kind=wp), intent(inout) :: CI(n)
  ! SGUGA, Faroald: the split-graph basis is the stored one
  if ((CI_Lib_Op == CI_Lib_SGUGA) .or. (CI_Lib_Op == CI_Lib_Faroald)) return
  call Check_Lib('CI_Lib_File_To_Basis')
  call MCLR_Lucia_Conv(iState,iSym,CI_To_Basis,n,CI)
end subroutine CI_Lib_File_To_Basis

subroutine CI_Lib_Basis_To_File(iState,iSym,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, n
  real(kind=wp), intent(inout) :: CI(n)
  ! SGUGA, Faroald: the split-graph basis is the stored one
  if ((CI_Lib_Op == CI_Lib_SGUGA) .or. (CI_Lib_Op == CI_Lib_Faroald)) return
  call Check_Lib('CI_Lib_Basis_To_File')
  call MCLR_Lucia_Conv(iState,iSym,CI_To_File,n,CI)
end subroutine CI_Lib_Basis_To_File

!> The setup of the split graph that the selected library needs ('FULL':
!> with the coupling tables, for SGUGA; 'GRAPH': graph and walks only)
function CI_Lib_Graph_Mode() result(Mode)
  character(len=5) :: Mode
  if (CI_Lib_Op == CI_Lib_SGUGA) then
    Mode = 'FULL'
  else
    Mode = 'GRAPH'
  end if
end function CI_Lib_Graph_Mode

subroutine Not_Yet(Caller,What)
  character(len=*), intent(in) :: Caller, What
  write(u6,'(1X,A,A,A)') Caller,': not available with this CI library yet: ',What
  call Abend()
end subroutine Not_Yet

!> The active integrals for SGUGA, for an operator of symmetry iOp: the
!> active part of Int1 (square blocks (iS,Mul(iS,iOp)), nOrb(iS) x
!> nOrb(jS), one after the other in the order of iS) packed triangularly
!> and Int2 (zero without two-body part), both over the symmetry-ordered
!> active orbitals. The SGUGA sigma and diagonal routines
!> map them to the levels themselves (do not map them here).
subroutine SG_Ints(iState,iOp,nInt1,Int1,nInt2,Int2,Have_2_el,nTU,TU,nTUVX,TUVX)
  use Symmetry_Info, only: Mul
  use Orbital_Spaces, only: nAsh, nIsh, nOrb, nSym
  use Index_Functions, only: iTri, nTri_Elem
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iOp, nInt1, nInt2, nTU, nTUVX
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2)
  logical(kind=iwp), intent(in) :: Have_2_el
  real(kind=wp), intent(out) :: TU(nTU), TUVX(nTUVX)
  integer(kind=iwp) :: iOff(8), iS, it, iu, jS, nAOff(8), nAct
#include "macros.fh"
  unused_var(iState)
  TU(:) = Zero
  ! offsets of the Int1 blocks (iS,Mul(iS,iOp)) and of the active orbitals
  nAct = 0
  do iS=1,nSym
    jS = Mul(iS,iOp)
    if (iS == 1) then
      iOff(iS) = 0
    else
      iOff(iS) = iOff(iS-1)+nOrb(iS-1)*nOrb(Mul(iS-1,iOp))
    end if
    nAOff(iS) = nAct
    nAct = nAct+nAsh(iS)
  end do
  if (iOff(nSym)+nOrb(nSym)*nOrb(Mul(nSym,iOp)) > nInt1) call Abend()
  do iS=1,nSym
    jS = Mul(iS,iOp)
    ! element (t,u), t of symmetry iS, u of symmetry jS (Hermitian: once)
    do iu=1,nAsh(jS)
      do it=1,nAsh(iS)
        if (nAOff(iS)+it < nAOff(jS)+iu) cycle
        TU(iTri(nAOff(iS)+it,nAOff(jS)+iu)) = Int1(iOff(iS)+(nIsh(jS)+iu-1)*nOrb(iS)+nIsh(iS)+it)
      end do
    end do
  end do
  if (Have_2_el) then
    if (nInt2 < nTUVX) call Abend()
    TUVX(:) = Int2(1:nTUVX)
  else
    TUVX(:) = Zero
  end if
  if (nTU /= nTri_Elem(nAct)) call Abend()
end subroutine SG_Ints

!> SGUGA: sigma for a singlet operator between vectors of symmetry iSym
subroutine SG_Sigma_Op(iState,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,nC,C,nS,S)
  use Symmetry_Info, only: Mul
  use sguga, only: CIS, EXS, SGS
  use SGUGA_Sigma, only: sg_h_psi, sg_h_psi_op, sg_h_psi_ras
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iCSym, iSSym, nInt1, nInt2, nC, nS
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2), C(nC)
  logical(kind=iwp), intent(in) :: Have_2_el
  real(kind=wp), intent(out) :: S(nS)
  integer(kind=iwp) :: iOp, nTU, nTUVX
  real(kind=wp), allocatable :: TU(:), TUVX(:)
  nTU = nTri_Elem(sum(nAsh(1:nSym)))
  nTUVX = nTri_Elem(nTU)
  call mma_allocate(TU,nTU,Label='TU')
  call mma_allocate(TUVX,nTUVX,Label='TUVX')
  iOp = Mul(iCSym,iSSym)
  call SG_Ints(iState,iOp,nInt1,Int1,nInt2,Int2,Have_2_el,nTU,TU,nTUVX,TUVX)
  if (iOp /= 1) then
    ! operator of non-total symmetry (perturbations)
    call sg_h_psi_op(SGS(iState),CIS(iState),EXS(iState),nC,C,iCSym,nS,S,iOp,TUVX,nTUVX,TU,nTU)
  else if (CI_St(iState)%nRsPrt > 1) then
    if (nS /= nC) call Abend()
    call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),C,nC,iCSym,S,TUVX,nTUVX,TU,nTU)
  else
    if (nS /= nC) call Abend()
    call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),C,nC,iCSym,S,TUVX,nTUVX,TU,nTU)
  end if
  call mma_deallocate(TUVX)
  call mma_deallocate(TU)
end subroutine SG_Sigma_Op

!> SGUGA: sigma of a non-Hermitian operator (Square: the time-dependent
!> response): the active part of Int1 (square blocks (iS,Mul(iS,iOp)))
!> as the full nA x nA matrix h(i,j), and Int2 as it is (g(ij,kl) at
!> iTri(ij,kl), square pair index), over the symmetry-ordered active
!> orbitals
subroutine SG_Sigma_Sq(iState,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,nC,C,nS,S)
  use Symmetry_Info, only: Mul
  use sguga, only: CIS, EXS, SGS
  use SGUGA_Sigma, only: sg_h_psi_sq
  use Orbital_Spaces, only: nAsh, nIsh, nOrb, nSym
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iCSym, iSSym, nInt1, nInt2, nC, nS
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2), C(nC)
  logical(kind=iwp), intent(in) :: Have_2_el
  real(kind=wp), intent(out) :: S(nS)
  integer(kind=iwp) :: iOff(8), iOp, iS, it, iu, jS, n2, nA, nAOff(8)
  real(kind=wp), allocatable :: G2(:), H1(:,:)
  iOp = Mul(iCSym,iSSym)
  nA = 0
  do iS=1,nSym
    if (iS == 1) then
      iOff(iS) = 0
    else
      iOff(iS) = iOff(iS-1)+nOrb(iS-1)*nOrb(Mul(iS-1,iOp))
    end if
    nAOff(iS) = nA
    nA = nA+nAsh(iS)
  end do
  if (iOff(nSym)+nOrb(nSym)*nOrb(Mul(nSym,iOp)) > nInt1) call Abend()
  call mma_allocate(H1,nA,nA,Label='H1')
  H1(:,:) = Zero
  do iS=1,nSym
    jS = Mul(iS,iOp)
    do iu=1,nAsh(jS)
      do it=1,nAsh(iS)
        H1(nAOff(iS)+it,nAOff(jS)+iu) = Int1(iOff(iS)+(nIsh(jS)+iu-1)*nOrb(iS)+nIsh(iS)+it)
      end do
    end do
  end do
  n2 = nTri_Elem(nA**2)
  call mma_allocate(G2,n2,Label='G2')
  if (Have_2_el) then
    if (nInt2 < n2) call Abend()
    G2(:) = Int2(1:n2)
  else
    G2(:) = Zero
  end if
  call sg_h_psi_sq(SGS(iState),CIS(iState),EXS(iState),nC,C,iCSym,nS,S,iOp,nA,H1,n2,G2)
  call mma_deallocate(G2)
  call mma_deallocate(H1)
end subroutine SG_Sigma_Sq

!> SGUGA: spin-free transition densities, mapped from the level order to
!> the symmetry-ordered active orbitals
subroutine SG_Trans_Dens(iState,I12,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  use sguga, only: CIS, EXS, SGS, SG_Level_Map
  use SGUGA_Densities, only: sg_trans_two_pdm
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, I12, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  integer(kind=iwp) :: a, b, c, d, i, ij, ijL, j, k, kl, klL, l_, N
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: R1(:,:), R2(:)
  N = sum(nAsh(1:nSym))
  if (SGS(iState)%nLev /= N) call Abend()
  if (n1 < N**2) call Abend()
  call mma_allocate(L2O,N,Label='L2O')
  call mma_allocate(R1,N,N,Label='R1')
  call mma_allocate(R2,nTri_Elem(N**2),Label='R2')
  call SG_Level_Map(SGS(iState),L2O)
  call sg_trans_two_pdm(SGS(iState),CIS(iState),EXS(iState),nL,L,iLSym,nR,R,iRSym,N,R1,size(R2),R2)
  Rho1(:) = Zero
  do b=1,N
    do a=1,N
      Rho1(L2O(a)+(L2O(b)-1)*N) = R1(a,b)
    end do
  end do
  if (I12 == 2) then
    if (n2 < nTri_Elem(N**2)) call Abend()
    Rho2(:) = Zero
    do b=1,N
      do a=1,N
        ijL = a+(b-1)*N
        i = L2O(a)
        j = L2O(b)
        ij = i+(j-1)*N
        do d=1,N
          do c=1,N
            klL = c+(d-1)*N
            if (klL > ijL) cycle
            k = L2O(c)
            l_ = L2O(d)
            kl = k+(l_-1)*N
            Rho2(iTri(ij,kl)) = R2(iTri(ijL,klL))
          end do
        end do
      end do
    end do
  end if
  call mma_deallocate(R2)
  call mma_deallocate(R1)
  call mma_deallocate(L2O)
end subroutine SG_Trans_Dens

!> SGUGA: the explicit P space (as Mk_IHJ for RASSCF): whole configurations
!> by increasing mean of the exact CSF diagonal until nPMax CSFs, and the
!> block of the active Hamiltonian (no core energy, as for the MCLR copy of
!> Lucia; the program shifts it). The diagonal of the program (Diag) is not
!> needed for the selection.
subroutine SG_PSpace(iState,iSym,nInt1,Int1,nInt2,Int2,nPMax,nP,IdxP,nHP,HP)
  use sguga, only: CIS, EXS, SGCnfList, SGS, SG_Explicit_Block, SG_Free_Cnf_List, SG_Mk_Cnf_List, SG_Select_Explicit
  use SGUGA_Sigma, only: sg_hdiag_fast
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iSym, nInt1, nInt2, nPMax, nHP
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2)
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  integer(kind=iwp) :: nCSF, NPCNF, nTU, nTUVX
  integer(kind=iwp), allocatable :: IPCNF(:)
  real(kind=wp), allocatable :: DX(:), TU(:), TUVX(:)
  type(SGCnfList) :: Cnf
  nTU = nTri_Elem(sum(nAsh(1:nSym)))
  nTUVX = nTri_Elem(nTU)
  call mma_allocate(TU,nTU,Label='TU')
  call mma_allocate(TUVX,nTUVX,Label='TUVX')
  call SG_Ints(iState,1,nInt1,Int1,nInt2,Int2,.true.,nTU,TU,nTUVX,TUVX)
  nCSF = CI_St(iState)%nCSF(iSym)
  call mma_allocate(DX,nCSF,Label='DX')
  call sg_hdiag_fast(SGS(iState),CIS(iState),EXS(iState),nCSF,iSym,TUVX,nTUVX,TU,nTU,DX)
  call SG_Mk_Cnf_List(iState,iSym,Cnf)
  call mma_allocate(IPCNF,max(1,Cnf%nCnf),Label='IPCNF')
  call SG_Select_Explicit(Cnf,DX,nPMax,NPCNF,IPCNF,nP,IdxP)
  if (nHP < nTri_Elem(nP)) call Abend()
  call SG_Explicit_Block(iState,Cnf,NPCNF,IPCNF,Zero,TU,nTU,TUVX,nTUVX,HP)
  call SG_Free_Cnf_List(Cnf)
  call mma_deallocate(IPCNF)
  call mma_deallocate(DX)
  call mma_deallocate(TUVX)
  call mma_deallocate(TU)
end subroutine SG_PSpace

!> SGUGA: the CSF diagonal (exact or configuration average) + ECore
subroutine SG_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
  use sguga, only: CIS, EXS, SGCnfList, SGS, SG_Cnf_Avg_Diag, SG_Free_Cnf_List, SG_Mk_Cnf_List
  use SGUGA_Sigma, only: sg_hdiag_fast
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nInt1, nInt2, n
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2), ECore
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp) :: nTU, nTUVX
  real(kind=wp), allocatable :: TU(:), TUVX(:)
  type(SGCnfList) :: Cnf
  nTU = nTri_Elem(sum(nAsh(1:nSym)))
  nTUVX = nTri_Elem(nTU)
  call mma_allocate(TU,nTU,Label='TU')
  call mma_allocate(TUVX,nTUVX,Label='TUVX')
  call SG_Ints(iState,1,nInt1,Int1,nInt2,Int2,.true.,nTU,TU,nTUVX,TUVX)
  if (Mode == CI_Diag_Exact) then
    call sg_hdiag_fast(SGS(iState),CIS(iState),EXS(iState),n,iSym,TUVX,nTUVX,TU,nTU,Diag)
  else
    call SG_Mk_Cnf_List(iState,iSym,Cnf)
    call SG_Cnf_Avg_Diag(iState,Cnf,TU,nTU,TUVX,nTUVX,Diag)
    call SG_Free_Cnf_List(Cnf)
  end if
  Diag(:) = Diag(:)+ECore
  call mma_deallocate(TUVX)
  call mma_deallocate(TU)
end subroutine SG_Diag

!***********************************************************************
! Faroald (string-driven RAS CI, faroald_ras), in its CSF basis (citrans);
! the vectors of the program stay in the split-graph order of the files,
! sorted to and from citrans for each operation. One CSF space, of the
! symmetry of the CI state.
!***********************************************************************

!> Faroald: the CSF space of the CI state (RAS partitions of CI_Input)
subroutine Faro_Init(iState)
  use CI_Input, only: iSpin, nActEl, nElec3, nHole1, nRS1, nRS2, nRS3
  use CI_States, only: RAS_Partitions
  use Orbital_Spaces, only: nAsh, nSym
  use faroald, only: nela, nelb
  use faroald_ras, only: faroald_ras_init
  use citrans, only: citrans_ncsf
  use Molcas, only: MxSym
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: i, iLev, iSym, k, MaxEl(3), MinEl(3), nAbove, nAct, nBelow, nOff(MxSym), nOrbP(3), nPrt, nRas(MxSym,3), &
                       nRasEl(3)
  integer(kind=iwp), allocatable :: OrbRas(:), OrbSym(:)
  Faro_Sym = CI_St(iState)%STSym
  nAct = sum(nAsh(1:nSym))
  call RAS_Partitions(nSym,nActEl,nRs1,nRs2,nRs3,nHole1,nElec3,.false.,nPrt,nRas,nRasEl)
  call mma_allocate(OrbSym,max(1,nAct),Label='OrbSym')
  call mma_allocate(OrbRas,max(1,nAct),Label='OrbRas')
  call mma_deallocate(Faro_Lev,safe='*')
  call mma_allocate(Faro_Lev,max(1,nAct),Label='Faro_Lev')
  ! the levels: partition by partition, within each by symmetry
  iLev = 0
  do k=1,nPrt
    do iSym=1,nSym
      do i=1,nRas(iSym,k)
        iLev = iLev+1
        OrbSym(iLev) = iSym
        OrbRas(iLev) = k
      end do
    end do
  end do
  if (iLev /= nAct) then
    write(u6,*) 'Faro_Init: the RAS partitions do not hold the active orbitals',iLev,nAct
    call Abend()
  end if
  ! cumulative electron bounds of the partitions 1..k
  nBelow = 0
  do k=1,nPrt
    nOrbP(k) = sum(nRas(1:nSym,k))
    nBelow = nBelow+nOrbP(k)
    nAbove = nAct-nBelow
    MinEl(k) = max(nRasEl(k),0,nActEl-2*nAbove)
    MaxEl(k) = min(nActEl,2*nBelow)
  end do
  MinEl(nPrt) = nActEl
  MaxEl(nPrt) = nActEl
  ! levels -> symmetry-ordered active orbitals
  nOff(1) = 0
  do iSym=2,nSym
    nOff(iSym) = nOff(iSym-1)+nAsh(iSym-1)
  end do
  do iLev=1,nAct
    nOff(OrbSym(iLev)) = nOff(OrbSym(iLev))+1
    Faro_Lev(iLev) = nOff(OrbSym(iLev))
  end do
  call faroald_init(nActEl,nAct,iSpin,.false.)
  call citrans_init(nActEl,nAct,iSpin,nSym,OrbSym,Faro_Sym,nPrt,OrbRas,MinEl(1:nPrt),MaxEl(1:nPrt))
  call faroald_ras_init(nela,nelb,nAct,nPrt,nOrbP(1:nPrt),MinEl(1:nPrt),MaxEl(1:nPrt),nSym,OrbSym(1:nAct),Faro_Sym,iSpin == 1)
  ! the description of the CSF spaces, for the other symmetries
  call mma_deallocate(Faro_OrbSym,safe='*')
  call mma_deallocate(Faro_OrbRas,safe='*')
  call mma_allocate(Faro_OrbSym,max(1,nAct),Label='Faro_OrbSym')
  call mma_allocate(Faro_OrbRas,max(1,nAct),Label='Faro_OrbRas')
  Faro_OrbSym(:) = OrbSym(:)
  Faro_OrbRas(:) = OrbRas(:)
  Faro_nAct = nAct
  Faro_nPrt = nPrt
  Faro_MinEl(:) = MinEl(:)
  Faro_MaxEl(:) = MaxEl(:)
  Faro_Cur = Faro_Sym
  call mma_deallocate(OrbSym)
  call mma_deallocate(OrbRas)
  if (citrans_ncsf() /= CI_St(iState)%nCSF(Faro_Sym)) then
    write(u6,*) 'Faro_Init: CSFs of Faroald and of the split graph',citrans_ncsf(),CI_St(iState)%nCSF(Faro_Sym)
    call Abend()
  end if
  ! no determinants
  CI_St(iState)%nDet(:) = CI_St(iState)%nCSF(:)
end subroutine Faro_Init

subroutine Faro_Close()
  use faroald_ras, only: faroald_ras_free
  use stdalloc, only: mma_deallocate
  call faroald_ras_free()
  call citrans_free()
  call faroald_free()
  call mma_deallocate(Faro_Lev,safe='*')
  call mma_deallocate(Faro_OrbSym,safe='*')
  call mma_deallocate(Faro_OrbRas,safe='*')
  Faro_Sym = 0
  Faro_Cur = 0
end subroutine Faro_Close

!> Faroald: citrans for the CSF space of symmetry iSym (the configuration
!> lists; the graph of the file order is common to all symmetries)
subroutine Faro_Set_Sym(iState,iSym)
  use CI_Input, only: nActEl
  use Orbital_Spaces, only: nSym
  use citrans, only: citrans_ncsf, citrans_sym_init
  integer(kind=iwp), intent(in) :: iState, iSym
  if (iSym == Faro_Cur) return
  call citrans_sym_init(nActEl,Faro_nAct,nSym,Faro_OrbSym,iSym,Faro_nPrt,Faro_OrbRas,Faro_MinEl(1:Faro_nPrt), &
                        Faro_MaxEl(1:Faro_nPrt))
  Faro_Cur = iSym
  if (citrans_ncsf() /= CI_St(iState)%nCSF(iSym)) then
    write(u6,*) 'Faro_Set_Sym: CSFs of Faroald and of the split graph',iSym,citrans_ncsf(),CI_St(iState)%nCSF(iSym)
    call Abend()
  end if
end subroutine Faro_Set_Sym

!> Faroald: the tables of the full determinant space (sigma_update and the
!> transition densities over determinant matrices), for the vectors of
!> other symmetries than the CI state; set up on first use
subroutine Faro_Full_Init()
  use CI_Input, only: iSpin, nActEl
  use faroald, only: faroald_full, ndeta, ndetb
  if (faroald_full) return
  call faroald_free()
  call faroald_init(nActEl,Faro_nAct,iSpin,.true.)
  if (real(ndeta,kind=wp)*real(ndetb,kind=wp) > 2.0e8_wp) &
    call Not_Yet('Faroald','vectors of other symmetries than the CI state in this large a determinant space')
end subroutine Faro_Full_Init

!> Faroald: its integrals htu, gtuvx (over the levels) from Int1/Int2 of
!> an Hermitian operator of symmetry iOp
subroutine Faro_Ints(iState,iOp,nInt1,Int1,nInt2,Int2,Have_2_el)
  use faroald, only: gtuvx, htu, my_norb
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iOp, nInt1, nInt2
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2)
  logical(kind=iwp), intent(in) :: Have_2_el
  integer(kind=iwp) :: it, iu, iv, ix, nTU, nTUVX, otu, ovx
  real(kind=wp), allocatable :: TU(:), TUVX(:)
  nTU = nTri_Elem(my_norb)
  nTUVX = nTri_Elem(nTU)
  call mma_allocate(TU,nTU,Label='TU')
  call mma_allocate(TUVX,nTUVX,Label='TUVX')
  call SG_Ints(iState,iOp,nInt1,Int1,nInt2,Int2,Have_2_el,nTU,TU,nTUVX,TUVX)
  do iu=1,my_norb
    do it=1,my_norb
      otu = iTri(Faro_Lev(it),Faro_Lev(iu))
      htu(it,iu) = TU(otu)
      do ix=1,my_norb
        do iv=1,my_norb
          ovx = iTri(Faro_Lev(iv),Faro_Lev(ix))
          gtuvx(it,iu,iv,ix) = TUVX(iTri(otu,ovx))
        end do
      end do
    end do
  end do
  call mma_deallocate(TUVX)
  call mma_deallocate(TU)
end subroutine Faro_Ints

!> Faroald: htu, gtuvx (over the levels) of a non-Hermitian operator of
!> symmetry iOp (Square): h(t,u) from the blocks (iS,Mul(iS,iOp)) of Int1,
!> g(tu,vx) from Int2 at iTri(tu,vx), square pair index tu = (u-1)*N+t
!> over the symmetry-ordered active orbitals
subroutine Faro_Ints_Sq(iOp,nInt1,Int1,nInt2,Int2,Have_2_el)
  use Symmetry_Info, only: Mul
  use faroald, only: gtuvx, htu, my_norb
  use Orbital_Spaces, only: nAsh, nIsh, nOrb, nSym
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, Half
  integer(kind=iwp), intent(in) :: iOp, nInt1, nInt2
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2)
  logical(kind=iwp), intent(in) :: Have_2_el
  integer(kind=iwp) :: iOff(8), iS, it, iu, iv, ix, jS, N, nAOff(8), otu, ovx
  real(kind=wp), allocatable :: G(:,:,:,:), H1(:,:)
  N = my_norb
  call mma_allocate(H1,N,N,Label='H1')
  H1(:,:) = Zero
  nAOff(1) = 0
  iOff(1) = 0
  do iS=2,nSym
    iOff(iS) = iOff(iS-1)+nOrb(iS-1)*nOrb(Mul(iS-1,iOp))
    nAOff(iS) = nAOff(iS-1)+nAsh(iS-1)
  end do
  if (iOff(nSym)+nOrb(nSym)*nOrb(Mul(nSym,iOp)) > nInt1) call Abend()
  do iS=1,nSym
    jS = Mul(iS,iOp)
    do iu=1,nAsh(jS)
      do it=1,nAsh(iS)
        H1(nAOff(iS)+it,nAOff(jS)+iu) = Int1(iOff(iS)+(nIsh(jS)+iu-1)*nOrb(iS)+nIsh(iS)+it)
      end do
    end do
  end do
  if (Have_2_el .and. (nInt2 < nTri_Elem(N**2))) call Abend()
  ! sigma_update applies sum_tu h(t,u) <J|E_tu|I> to C(J), i.e. E_ut with
  ! h(t,u), and likewise within each pair of g (harmless for Hermitian
  ! operators): pass g'(t,u,v,x) = g(u,t,x,v) and h transposed, and correct
  ! h for the effective one-electron term that sigma_update forms itself,
  ! k(t,u) = h'(t,u)-1/2 sum_v g'(v,t,v,u), which must be
  ! h(u,t)-1/2 sum_v g(u,v,v,t) (g over the levels in G)
  call mma_allocate(G,N,N,N,N,Label='G')
  do ix=1,N
    do iv=1,N
      ovx = Faro_Lev(iv)+(Faro_Lev(ix)-1)*N
      do iu=1,N
        do it=1,N
          otu = Faro_Lev(it)+(Faro_Lev(iu)-1)*N
          if (Have_2_el) then
            G(it,iu,iv,ix) = Int2(iTri(otu,ovx))
          else
            G(it,iu,iv,ix) = Zero
          end if
        end do
      end do
    end do
  end do
  do ix=1,N
    do iv=1,N
      do iu=1,N
        do it=1,N
          gtuvx(it,iu,iv,ix) = G(iu,it,ix,iv)
        end do
      end do
    end do
  end do
  do iu=1,N
    do it=1,N
      htu(it,iu) = H1(Faro_Lev(iu),Faro_Lev(it))
      do iv=1,N
        htu(it,iu) = htu(it,iu)+Half*(gtuvx(iv,it,iv,iu)-G(iu,iv,iv,it))
      end do
    end do
  end do
  call mma_deallocate(G)
  call mma_deallocate(H1)
end subroutine Faro_Ints_Sq

!> Faroald: S = H C in the CSF space of the CI state
subroutine Faro_Sigma_Op(iState,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,nC,C,nS,S,Square)
  use Symmetry_Info, only: Mul
  use faroald, only: gtuvx, htu, ndeta, ndetb, sigma_update
  use faroald_ras, only: n_ras, ras_sigma
  use citrans, only: citrans_csf2ras, citrans_csf2sd, citrans_ras2csf, citrans_sd2csf, citrans_sort
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iCSym, iSSym, nInt1, nInt2, nC, nS
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2), C(nC)
  logical(kind=iwp), intent(in) :: Have_2_el
  real(kind=wp), intent(out) :: S(nS)
  logical(kind=iwp), intent(in), optional :: Square
  logical(kind=iwp) :: Sq
  real(kind=wp), allocatable :: PC(:,:), PS(:,:), V(:), VR(:), W(:), WR(:)
  Sq = .false.
  if (present(Square)) Sq = Square
  if ((iCSym /= Faro_Sym) .or. (iSSym /= Faro_Sym) .or. Sq) then
    ! other symmetries or a non-Hermitian operator: in the full determinant
    ! space, P_S H P_C exactly
    call Faro_Full_Init()
    if (Sq) then
      call Faro_Ints_Sq(Mul(iCSym,iSSym),nInt1,Int1,nInt2,Int2,Have_2_el)
    else
      call Faro_Ints(iState,Mul(iCSym,iSSym),nInt1,Int1,nInt2,Int2,Have_2_el)
    end if
    call mma_allocate(PC,ndeta,ndetb,Label='PC')
    call mma_allocate(PS,ndeta,ndetb,Label='PS')
    call mma_allocate(V,nC,Label='V')
    call mma_allocate(W,nS,Label='W')
    call Faro_Set_Sym(iState,iCSym)
    call citrans_sort('C',C,V)
    call citrans_csf2sd(V,PC)
    PS(:,:) = Zero
    call sigma_update(htu,gtuvx,PS,PC)
    call Faro_Set_Sym(iState,iSSym)
    call citrans_sd2csf(PS,W)
    call citrans_sort('O',W,S)
    call mma_deallocate(W)
    call mma_deallocate(V)
    call mma_deallocate(PS)
    call mma_deallocate(PC)
    return
  end if
  if (nS /= nC) call Abend()
  call Faro_Set_Sym(iState,Faro_Sym)
  call Faro_Ints(iState,1,nInt1,Int1,nInt2,Int2,Have_2_el)
  call mma_allocate(V,nC,Label='V')
  call mma_allocate(W,nC,Label='W')
  call mma_allocate(VR,max(1,n_ras),Label='VR')
  call mma_allocate(WR,max(1,n_ras),Label='WR')
  call citrans_sort('C',C,V)
  call citrans_csf2ras(V,VR)
  call ras_sigma(htu,gtuvx,VR,WR)
  call citrans_ras2csf(WR,W)
  call citrans_sort('O',W,S)
  call mma_deallocate(WR)
  call mma_deallocate(VR)
  call mma_deallocate(W)
  call mma_deallocate(V)
end subroutine Faro_Sigma_Op

!> Faroald: the spin-free transition densities (as SG_Trans_Dens), mapped
!> from the levels to the symmetry-ordered active orbitals
subroutine Faro_Trans_Dens(iState,I12,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  use faroald, only: my_norb, ndeta, ndetb, transition_one_pdm, transition_two_pdm
  use faroald_ras, only: n_ras, ras_trans_one_pdm, ras_trans_two_pdm
  use citrans, only: citrans_csf2ras, citrans_csf2sd, citrans_sort
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, I12, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  integer(kind=iwp) :: ij, it, iu, iv, ix, kl, N
  logical(kind=iwp) :: Full
  real(kind=wp), allocatable :: BR(:), D(:,:), KR(:), P2(:,:,:,:), PL(:,:), PR(:,:), SD(:,:), V(:)
  N = my_norb
  if (n1 < N**2) call Abend()
  ! other symmetries: in the full determinant space
  Full = (iLSym /= Faro_Sym) .or. (iRSym /= Faro_Sym)
  call mma_allocate(V,max(nL,nR),Label='V')
  if (Full) then
    call Faro_Full_Init()
    call mma_allocate(PL,ndeta,ndetb,Label='PL')
    call mma_allocate(PR,ndeta,ndetb,Label='PR')
    call Faro_Set_Sym(iState,iLSym)
    call citrans_sort('C',L,V(1:nL))
    call citrans_csf2sd(V(1:nL),PL)
    call Faro_Set_Sym(iState,iRSym)
    call citrans_sort('C',R,V(1:nR))
    call citrans_csf2sd(V(1:nR),PR)
  else
    call Faro_Set_Sym(iState,Faro_Sym)
    call mma_allocate(BR,max(1,n_ras),Label='BR')
    call mma_allocate(KR,max(1,n_ras),Label='KR')
    call citrans_sort('C',L,V(1:nL))
    call citrans_csf2ras(V(1:nL),BR)
    call citrans_sort('C',R,V(1:nR))
    call citrans_csf2ras(V(1:nR),KR)
  end if
  call mma_allocate(D,N,N,Label='D')
  call mma_allocate(SD,N,N,Label='SD')
  if (Full) then
    call transition_one_pdm(PL,PR,D,SD)
  else
    call ras_trans_one_pdm(BR,KR,D,SD)
  end if
  Rho1(:) = Zero
  do iu=1,N
    do it=1,N
      Rho1(Faro_Lev(it)+(Faro_Lev(iu)-1)*N) = D(it,iu)
    end do
  end do
  call mma_deallocate(SD)
  call mma_deallocate(D)
  if (I12 == 2) then
    if (n2 < nTri_Elem(N**2)) call Abend()
    call mma_allocate(P2,N,N,N,N,Label='P2')
    if (Full) then
      call transition_two_pdm(PL,PR,P2)
    else
      call ras_trans_two_pdm(BR,KR,P2)
    end if
    Rho2(:) = Zero
    do iu=1,N
      do it=1,N
        ij = Faro_Lev(it)+(Faro_Lev(iu)-1)*N
        do ix=1,N
          do iv=1,N
            kl = Faro_Lev(iv)+(Faro_Lev(ix)-1)*N
            if (kl > ij) cycle
            Rho2(iTri(ij,kl)) = P2(it,iu,iv,ix)
          end do
        end do
      end do
    end do
    call mma_deallocate(P2)
  end if
  if (Full) then
    call mma_deallocate(PR)
    call mma_deallocate(PL)
  else
    call mma_deallocate(KR)
    call mma_deallocate(BR)
  end if
  call mma_deallocate(V)
end subroutine Faro_Trans_Dens

!> Faroald: the configuration average of the CSF diagonal (as for RASSCF
!> with Faroald: the average of the determinant energies of each
!> configuration, from its occupations) + ECore
subroutine Faro_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
  use faroald, only: gtuvx, htu, my_norb, nela, nelb
  use citrans, only: citrans_conf_occ, citrans_sort
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Two
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nInt1, nInt2, n
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2), ECore
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp) :: i, ic, j, na, nb, ncf, nd2, no
  integer(kind=iwp), allocatable :: dco(:), kf(:), ncs(:), sco(:)
  real(kind=wp) :: E, ni, nj, P
  real(kind=wp), allocatable :: DC(:), Ph(:)
  if (Mode /= CI_Diag_Cnf_Avg) call Not_Yet('CI_Diag','the exact CSF diagonal')
  call Faro_Set_Sym(iState,iSym)
  call Faro_Ints(iState,1,nInt1,Int1,nInt2,Int2,.true.)
  call mma_allocate(DC,n,Label='DC')
  call mma_allocate(Ph,n,Label='Ph')
  call citrans_conf_occ(ncf,kf,ncs,dco,sco)
  do ic=1,ncf
    nd2 = popcnt(dco(ic))
    no = popcnt(sco(ic))
    na = nela-nd2
    nb = nelb-nd2
    P = Zero
    if (no > 1) P = real(na*(na-1)+nb*(nb-1),kind=wp)/real(no*(no-1),kind=wp)
    E = Zero
    do i=1,my_norb
      if (btest(dco(ic),i-1)) then
        ni = Two
        E = E+gtuvx(i,i,i,i)
      else if (btest(sco(ic),i-1)) then
        ni = One
      else
        cycle
      end if
      E = E+ni*htu(i,i)
      do j=i+1,my_norb
        if (btest(dco(ic),j-1)) then
          nj = Two
        else if (btest(sco(ic),j-1)) then
          nj = One
        else
          cycle
        end if
        E = E+ni*nj*gtuvx(i,i,j,j)
        if ((ni == Two) .and. (nj == Two)) then
          E = E-Two*gtuvx(i,j,j,i)
        else if ((ni == Two) .or. (nj == Two)) then
          E = E-gtuvx(i,j,j,i)
        else
          E = E-P*gtuvx(i,j,j,i)
        end if
      end do
    end do
    DC(kf(ic):kf(ic)+ncs(ic)-1) = E
  end do
  call mma_deallocate(kf)
  call mma_deallocate(ncs)
  call mma_deallocate(dco)
  call mma_deallocate(sco)
  ! to the file order (the sort carries the CSF phases: undo them)
  call citrans_sort('O',DC,Diag)
  DC(:) = One
  call citrans_sort('O',DC,Ph)
  Diag(:) = Diag(:)*Ph(:)+ECore
  call mma_deallocate(Ph)
  call mma_deallocate(DC)
end subroutine Faro_Diag

!> Faroald: the explicit P space: whole configurations by increasing
!> configuration average until nPMax CSFs (as SGUGA and the MCLR copy of
!> Lucia), the block of the active Hamiltonian (no core energy) from the
!> sigma vectors of the unit vectors of the P space
subroutine Faro_PSpace(iState,iSym,nInt1,Int1,nInt2,Int2,n,nPMax,nP,IdxP,nHP,HP)
  use citrans, only: citrans_conf_occ, citrans_sort
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  integer(kind=iwp), intent(in) :: iState, iSym, nInt1, nInt2, n, nPMax, nHP
  real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2)
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  integer(kind=iwp) :: i, ic, icMin, j, k, ncf
  integer(kind=iwp), allocatable :: C2F(:), dco(:), kf(:), ncs(:), sco(:)
  logical(kind=iwp), allocatable :: Used(:)
  real(kind=wp) :: EMin
  real(kind=wp), allocatable :: Diag(:), E(:), X(:), Y(:)
  nP = 0
  IdxP(:) = 0
  HP(:) = Zero
  ! the configuration averages (file order) and the CSFs of each
  ! configuration (citrans order -> file order)
  call mma_allocate(Diag,n,Label='Diag')
  call Faro_Diag(iState,CI_Diag_Cnf_Avg,iSym,nInt1,Int1,nInt2,Int2,Zero,n,Diag)
  call mma_allocate(X,n,Label='X')
  call mma_allocate(Y,n,Label='Y')
  do k=1,n
    X(k) = real(k,kind=wp)
  end do
  call citrans_sort('O',X,Y)
  call mma_allocate(C2F,n,Label='C2F')
  do i=1,n
    C2F(nint(abs(Y(i)))) = i
  end do
  call Faro_Set_Sym(iState,iSym)
  call citrans_conf_occ(ncf,kf,ncs,dco,sco)
  call mma_allocate(E,ncf,Label='E')
  call mma_allocate(Used,ncf,Label='Used')
  do ic=1,ncf
    E(ic) = Diag(C2F(kf(ic)))
  end do
  Used(:) = .false.
  ! the configurations by increasing average, while they fit
  do
    icMin = 0
    EMin = huge(EMin)
    do ic=1,ncf
      if (Used(ic)) cycle
      if (E(ic) < EMin) then
        EMin = E(ic)
        icMin = ic
      end if
    end do
    if (icMin == 0) exit
    if (nP+ncs(icMin) > nPMax) exit
    Used(icMin) = .true.
    do k=kf(icMin),kf(icMin)+ncs(icMin)-1
      nP = nP+1
      IdxP(nP) = C2F(k)
    end do
  end do
  if (nHP < nTri_Elem(nP)) call Abend()
  ! the block: H(i,j) = (H e_j)(IdxP(i))
  do j=1,nP
    X(:) = Zero
    X(IdxP(j)) = One
    call Faro_Sigma_Op(iState,iSym,iSym,nInt1,Int1,nInt2,Int2,.true.,n,X,n,Y)
    do i=j,nP
      HP(iTri(i,j)) = Y(IdxP(i))
    end do
  end do
  call mma_deallocate(Used)
  call mma_deallocate(E)
  call mma_deallocate(kf)
  call mma_deallocate(ncs)
  call mma_deallocate(dco)
  call mma_deallocate(sco)
  call mma_deallocate(C2F)
  call mma_deallocate(Y)
  call mma_deallocate(X)
  call mma_deallocate(Diag)
end subroutine Faro_PSpace

!> The remaining case of the operations is the MCLR copy of Lucia
subroutine Check_Lib(Caller)
  character(len=*), intent(in) :: Caller
  if (CI_Lib_Op /= CI_Lib_MCLR_Lucia) then
    write(u6,'(1X,A,A,I2)') Caller,': no CI library selected, or not available:',CI_Lib_Op
    call Abend()
  end if
end subroutine Check_Lib

end module CI_Operations

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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> Generic CI operations of the CI interface, for any program, and the
!> choice of the CI library that performs them.
!>
!> Operations (all for the CI state iState, CI_St(iState)):
!>  - CI_Lib_Initialize, CI_Lib_Close: set up and free the selected
!>    library for the CI state (CI_Lib_Close is safe to call in any case).
!>    The library publishes its dimensions in CI_St(iState) (nCSF; nDet: the
!>    dimension of its determinant space per symmetry).
!>  - CI_Lib_Sigma: sigma = H C for one CI vector, H an active-space
!>    operator of any symmetry, Hermitian or not (layout of the integrals
!>    in the comment of the routine).
!>  - CI_Lib_Dens: the 1- and 2-particle (transition) densities between
!>    two CI vectors, the same one allowed, spin-free or spin (layout in the
!>    comment of the routine).
!>  - CI_Lib_Diag: the diagonal of H + ECore, configuration average or
!>    exact.
!>  - CI_Lib_nCnf, CI_Lib_Cnf_Layout: the configurations of a symmetry in
!>    the library's own order and their CSFs.
!>  - CI_Lib_Block: the explicit block of H + ECore for configurations
!>    chosen by the caller; CI_Lib_PSpace: the P space of the
!>    preconditioner (the configurations chosen from the diagonal by one
!>    rule for every library), then its block.
!>    (Integrals of all of them: those of CI_Lib_Sigma.)
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
use CI_Op_Defs, only: CI_Diag_Cnf_Avg, CI_Diag_Exact, CI_Lib, CI_Lib_Faroald, CI_Lib_Lucia, CI_Lib_MCLR_Lucia, CI_Lib_None, &
                      CI_Lib_Select, CI_Lib_SGUGA, CI_To_Basis, CI_To_File
use MCLR_Lucia_CI, only: MCLR_Lucia_Block, MCLR_Lucia_Close, MCLR_Lucia_Cnf_Layout, MCLR_Lucia_Conv, MCLR_Lucia_Dets, &
                         MCLR_Lucia_Diag_Act, MCLR_Lucia_Init, MCLR_Lucia_nCnf, MCLR_Lucia_Sigma_Act, MCLR_Lucia_TDens
use Definitions, only: wp, iwp, u6

implicit none
private

!> The selected library is CI_Lib of CI_Op_Defs (set by CI_Lib_Select, one
!> selection for every program; re-exported here)


! Faroald: its orbitals are the levels of the split graph (by RAS
! partition, within each by symmetry); Faro_Lev(iLev) is the
! symmetry-ordered active orbital of level iLev. Faro_Sym is the symmetry
! of the CI state, the space of the string-driven RAS CI (faroald_ras);
! citrans holds the CSF space of one symmetry at a time, Faro_Cur, set up
! again from the stored description (Faro_OrbSym, ...) when another one
! is needed (Faro_Set_Sym).
integer(kind=iwp) :: Faro_Cur = 0, Faro_nAct = 0, Faro_nPrt = 0, Faro_Sym = 0
integer(kind=iwp), allocatable :: Faro_MaxEl(:), Faro_MinEl(:), Faro_OrbRas(:), Faro_OrbSym(:)
!> (public, read-only: the level map for the programs' Faroald paths)
integer(kind=iwp), allocatable, protected :: Faro_Lev(:)

! GAS: SGUGA evaluates products of excitation operators (H, 2-RDM, spin and
! transition densities, diagonal, explicit block) in the relaxed space of
! the state (CI_St(iState)%iRel, see SG_Rel_State), which holds all their
! intermediate states; Rel_Map(LSM)%M(I) is the index there of the CSF I of
! symmetry LSM (SG_Rel_Map).
type Rel_Map_Type
  integer(kind=iwp), allocatable :: M(:)
end type Rel_Map_Type
type(Rel_Map_Type) :: Rel_Map(8)

! RASSCF's Lucia: the determinant scratch of its sigma (Lucia1_Sigma)
real(kind=wp), allocatable :: Lucia1_C(:), Lucia1_S(:)
! Faroald: its string-driven sigma checked once against the full space
logical(kind=iwp) :: Faro_RAS_Verified = .false.

public :: CI_Diag_Cnf_Avg, CI_Diag_Exact, CI_Lib, CI_Lib_Basis_To_File, CI_Lib_Block, CI_Lib_Close, CI_Lib_Cnf_Layout, &
          CI_Lib_Dens, CI_Lib_Dets_To_Basis, CI_Lib_Diag, CI_Lib_Faroald, CI_Lib_File_To_Basis, CI_Lib_Graph_Mode, &
          CI_Lib_Initialize, CI_Lib_Lucia, CI_Lib_MCLR_Lucia, CI_Lib_nCnf, CI_Lib_None, CI_Lib_PSpace, CI_Lib_Select, &
          CI_Lib_SGUGA, CI_Lib_Sigma, CI_To_Basis, CI_To_File, Faro_Lev, Rel_Map, Report_Dev, SG_Rel_Embed, SG_Rel_Free, &
          SG_Rel_Map, SG_Rel_Project

contains

!> Set up the selected library for the CI state
subroutine CI_Lib_Initialize(iState)
  use CI_Input, only: iSpin, nActEl, nElec3, nHole1, nRs1, nRs2, nRs3
  use Orbital_Spaces, only: nAsh, nIsh, nOrb, nSym
  use Molcas, only: MxSym
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: nDet(MxSym)
  if (CI_Lib == CI_Lib_SGUGA) then
    ! the split graph is the state of the CI interface (set up by the
    ! program with CI_Lib_Graph_Mode()); no determinants
    CI_St(iState)%nDet(:) = CI_St(iState)%nCSF(:)
    return
  end if
  if (CI_Lib == CI_Lib_Faroald) then
    call Faro_Init(iState)
    return
  end if
  call Check_Lib('CI_Lib_Initialize')
  ! the CI space as arguments: the library keeps its own copies
  call MCLR_Lucia_Init(nSym,nOrb(1:nSym),nIsh(1:nSym),nAsh(1:nSym),nRs1(1:nSym),nRs2(1:nSym),nRs3(1:nSym),nHole1,nElec3,nActEl, &
                       iSpin,CI_St(iState)%STSym,CI_St(iState)%nCSF(1:nSym),nDet(1:nSym))
  CI_St(iState)%nDet(:) = 0
  CI_St(iState)%nDet(1:nSym) = nDet(1:nSym)
end subroutine CI_Lib_Initialize

!> Free the library of the CI state (safe to call in any case); the
!> selection is kept (CI_Op_Defs: other code of the program may still read it)
subroutine CI_Lib_Close(iState)
  use SGUGA_Sigma, only: SG_Sigma_Free
  use SGUGA_Densities, only: SG_Densities_Free
  integer(kind=iwp), intent(in) :: iState
  if (CI_Lib == CI_Lib_MCLR_Lucia) call MCLR_Lucia_Close(iState)
  ! the work space of SGUGA (safe to call in any case, as in CI_Close)
  if (CI_Lib == CI_Lib_SGUGA) then
    call SG_Sigma_Free()
    call SG_Densities_Free()
  end if
  if (CI_Lib == CI_Lib_Faroald) call Faro_Close()
  if (CI_Lib == CI_Lib_Lucia) call Lucia1_Free()
end subroutine CI_Lib_Close

!> S = H C (S is overwritten, never accumulated) for the CI state iState
!> (its CI space, CI_St(iState)), C of symmetry iCSym (nC = nCSF(iCSym)),
!> S of symmetry iSSym (nS = nCSF(iSSym)), both in the CSF basis of the
!> selected library (CI_Lib). That basis is the library's own and differs
!> between the libraries; only CI_Lib_File_To_Basis and
!> CI_Lib_Basis_To_File map it to the order of the files. H is an
!> active-space operator of symmetry iOp = Mul(iCSym,iSSym), with no core
!> energy:
!>   H = sum_tu h_tu E_tu + 1/2 sum_tuvx g_tuvx (E_tu E_vx - delta_uv E_tx)
!> h_tu are the one-body integrals of H as written (for RASSCF the active
!> block of the inactive Fock matrix); the libraries form
!> k_tu = h_tu - 1/2 sum_v g_tvvu themselves.
!> iOpSpin = 0: singlet E_tu; iOpSpin = 1: a+_ta a_ua - a+_tb a_ub
!> (one-body only, Have_2_el = .false., not with Square).
!> Active orbitals: nAct = sum(nAsh), in symmetry order (all active
!> orbitals of irrep 1, then irrep 2, ...; RAS1, RAS2, RAS3 inside each
!> irrep).
!> Square = .false. (H Hermitian):
!>   TU(iTri(t,u)) = h_tu, the full triangle, nTU = nTri(nAct), not
!>   symmetry-compressed; elements outside the operator symmetry
!>   (Mul(Sym(t),Sym(u)) /= iOp) are set to zero by the caller.
!>   TUVX(iTri(iTri(t,u),iTri(v,x))) = g_tuvx = (tu|vx), the full
!>   triangle, nTUVX = nTri(nTri(nAct)), not symmetry-compressed; zero
!>   outside the operator symmetry likewise.
!> Square = .true. (H not Hermitian; time-dependent response):
!>   TU((u-1)*nAct+t) = h_tu, the full square matrix, nTU = nAct**2;
!>   TUVX(iTri(tu,vx)) = g_tuvx with the square pair index
!>   tu = (u-1)*nAct+t, nTUVX = nTri(nAct**2); g(tu,vx) = g(vx,tu) is
!>   assumed (the iTri packing).
!> Have_2_el = .false.: one-body operator only, TUVX not referenced (a
!> dummy of size 1 may be passed).
!> HEXS/DEXS: projected only if the program set up HEXS/DEXS for iState
!> (RASSCF input); otherwise no projection. (Not yet: RASSCF still calls
!> Mk_H_Psi, which projects.)
!> Scratch: the libraries keep their own work space (freed by
!> CI_Lib_Close).
subroutine CI_Lib_Sigma(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nTU, nTUVX, nC, nS
  logical(kind=iwp), intent(in) :: Square, Have_2_el
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), C(nC)
  real(kind=wp), intent(out) :: S(nS)
  call Check_Op('CI_Lib_Sigma',iOpSpin,Have_2_el,Square)
  select case (CI_Lib)
    case (CI_Lib_SGUGA)
      if (iOpSpin /= 0) call Not_Yet('CI_Lib_Sigma','triplet operators')
      call SG_Sigma_Act(iState,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
    case (CI_Lib_Faroald)
      if (iOpSpin /= 0) call Not_Yet('CI_Lib_Sigma','triplet operators')
      call Faro_Sigma_Act(iState,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
    case (CI_Lib_MCLR_Lucia)
      call MCLR_Lucia_Sigma_Act(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
    case (CI_Lib_Lucia)
      if ((iOpSpin /= 0) .or. Square .or. (.not. Have_2_el) .or. (iCSym /= iSSym) .or. (iCSym /= CI_St(iState)%STSym)) &
        call Not_Yet('CI_Lib_Sigma','RASSCF''s Lucia for other operators than H in the symmetry of the state')
      call Lucia1_Sigma(iState,iCSym,nTU,TU,nTUVX,TUVX,nC,C,S)
    case default
      call Check_Lib('CI_Lib_Sigma')
  end select
end subroutine CI_Lib_Sigma

!> RASSCF's Lucia (CI_Lib_Lucia): S = H C in the symmetry of the state, via
!> its determinant basis (csdtvc); the determinant scratch is kept between
!> calls (Lucia1_Free, CI_Lib_Close)
subroutine Lucia1_Sigma(iState,iSym,nTU,TU,nTUVX,TUVX,n,C,S)
  use lucia_data, only: IREFSM, Sigma_on_disk
  use cands, only: ICSM, ISSM
  use Lucia_Interface, only: Lucia_Util
  use stdalloc, only: mma_allocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, n
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), C(n)
  real(kind=wp), intent(out) :: S(n)
  integer(kind=iwp) :: IREFSM_Keep, nD, nW
  nD = CI_St(iState)%nDet(iSym)
  nW = max(n,nD)
  if (allocated(Lucia1_C)) then
    if (size(Lucia1_C) < nW) call Lucia1_Free()
  end if
  if (.not. allocated(Lucia1_C)) then
    call mma_allocate(Lucia1_C,nW,Label='Lucia1_C')
    call mma_allocate(Lucia1_S,nW,Label='Lucia1_S')
  end if
  IREFSM_Keep = IREFSM
  IREFSM = iSym
  ICSM = iSym
  ISSM = iSym
  ! CSFs -> determinants (Lucia1_S is the scratch of the conversion)
  Lucia1_C(1:n) = C(1:n)
  Lucia1_S(:) = Zero
  call csdtvc(Lucia1_C,Lucia1_S,1,iSym,0)
  Lucia1_C(1:nD) = Lucia1_S(1:nD)
  call Lucia_Util('Sigma',CI_Vector=Lucia1_C,Sigma_Vector=Lucia1_S,nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)
  ! the sigma vector exists on disk (densi_master)
  Sigma_on_disk = .true.
  ! determinants -> CSFs
  call csdtvc(S,Lucia1_S,2,iSym,0)
  IREFSM = IREFSM_Keep
  ICSM = IREFSM
  ISSM = IREFSM
end subroutine Lucia1_Sigma

!> RASSCF's Lucia: free the determinant scratch of Lucia1_Sigma
subroutine Lucia1_Free()
  use stdalloc, only: mma_deallocate
  call mma_deallocate(Lucia1_C,safe='*')
  call mma_deallocate(Lucia1_S,safe='*')
end subroutine Lucia1_Free

!> The checks of the operator arguments of CI_Lib_Sigma
subroutine Check_Op(Caller,iOpSpin,Have_2_el,Square)
  character(len=*), intent(in) :: Caller
  integer(kind=iwp), intent(in) :: iOpSpin
  logical(kind=iwp), intent(in) :: Have_2_el, Square
  if ((iOpSpin /= 0) .and. (iOpSpin /= 1)) then
    write(u6,'(1X,A,A,I4)') Caller,': invalid operator spin iOpSpin =',iOpSpin
    call Abend()
  end if
  if ((iOpSpin == 1) .and. (Have_2_el .or. Square)) then
    write(u6,'(1X,A,A)') Caller,': the operator with iOpSpin = 1 is one-body and Hermitian only'
    call Abend()
  end if
end subroutine Check_Op

!> The (transition) densities between the bra L (symmetry iLSym, nL =
!> nCSF(iLSym)) and the ket R (symmetry iRSym, nR = nCSF(iRSym)) of the CI
!> state iState, both in the CSF basis of the selected library (as in
!> CI_Lib_Sigma). L and R may be the same state (the same array may be
!> passed twice): then they are the densities of that state.
!> iOpSpin = 0: spin-free, Rho1(i+(j-1)*nAct) = <L|E_ij|R>, the full
!>   square nAct x nAct (n1 >= nAct**2), sum_i Rho1(i,i) = N_act <L|R>;
!>   I12 = 2 also Rho2(iTri(ij,kl)) = <L|E_ij E_kl - delta_jk E_il|R>
!>   with the square pair index ij = (j-1)*nAct+i, packed over ij >= kl
!>   (n2 >= nTri(nAct**2)); for I12 = 1, Rho2 is not referenced (a dummy of
!>   size 1 may be passed).
!> iOpSpin = 1: spin density, Rho1(i,j) = <L|a+_ia a_ja - a+_ib a_jb|R>,
!>   one-body only (I12 = 1), for M_S = S (the CI vectors are the high-spin
!>   components; the spin-adapted libraries use that convention).
!> Active orbitals in symmetry order, as in CI_Lib_Sigma. The densities
!> are the elements of the operators as written: not symmetrized, not
!> folded (the callers symmetrize where they need it). Programs that need
!> other layouts (RASSCF: D, SD, P, PA, triangular, in its level order)
!> convert after the call.
subroutine CI_Lib_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  if ((iOpSpin /= 0) .and. (iOpSpin /= 1)) then
    write(u6,'(A,I4)') ' CI_Lib_Dens: invalid operator spin iOpSpin =',iOpSpin
    call Abend()
  end if
  if ((iOpSpin == 1) .and. (I12 /= 1)) then
    write(u6,'(A)') ' CI_Lib_Dens: the spin density is one-body only (I12 = 1)'
    call Abend()
  end if
  if (CI_Lib == CI_Lib_SGUGA) then
    call SG_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
    return
  end if
  if (CI_Lib == CI_Lib_Faroald) then
    call Faro_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
    return
  end if
  call Check_Lib('CI_Lib_Dens')
  call MCLR_Lucia_TDens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
end subroutine CI_Lib_Dens

!> The diagonal of H + ECore (ECore added to every element) in the CSF
!> basis of symmetry iSym (n = nCSF(iSym)) of the selected library, for
!> the active Hamiltonian of total symmetry with the integrals of
!> CI_Lib_Sigma (Square = .false., iOp = 1: TU, TUVX full triangles over
!> the symmetry-ordered active orbitals).
!> Mode = CI_Diag_Cnf_Avg: the average over the CSFs of each configuration
!> (Lucia's CSDIAG; every CSF of a configuration gets the same value);
!> Mode = CI_Diag_Exact: the exact CSF diagonal.
subroutine CI_Lib_Diag(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nTU, nTUVX, n
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(out) :: Diag(n)
  call Check_Diag_Mode('CI_Lib_Diag',Mode)
  select case (CI_Lib)
    case (CI_Lib_SGUGA)
      call SG_Diag(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
    case (CI_Lib_Faroald)
      call Faro_Diag(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
    case (CI_Lib_MCLR_Lucia)
      call MCLR_Lucia_Diag_Act(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
    case default
      call Check_Lib('CI_Lib_Diag')
  end select
end subroutine CI_Lib_Diag

!> The check of the Mode of CI_Lib_Diag
subroutine Check_Diag_Mode(Caller,Mode)
  character(len=*), intent(in) :: Caller
  integer(kind=iwp), intent(in) :: Mode
  if ((Mode /= CI_Diag_Cnf_Avg) .and. (Mode /= CI_Diag_Exact)) then
    write(u6,'(1X,A,A,I4)') Caller,': invalid Mode =',Mode
    call Abend()
  end if
end subroutine Check_Diag_Mode

!> The number of configurations of symmetry iSym of the CI state iState in
!> the library's own configuration order (CI_Lib_Cnf_Layout)
function CI_Lib_nCnf(iState,iSym) result(nCnf)
  use sguga, only: SGCnfList, SG_Free_Cnf_List, SG_Mk_Cnf_List
  use citrans, only: citrans_conf_occ
  use stdalloc, only: mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym
  integer(kind=iwp) :: nCnf
  integer(kind=iwp), allocatable :: dco(:), kf(:), ncs(:), sco(:)
  type(SGCnfList) :: Cnf
  nCnf = 0
  select case (CI_Lib)
    case (CI_Lib_SGUGA)
      call SG_Mk_Cnf_List(iState,iSym,Cnf)
      nCnf = Cnf%nCnf
      call SG_Free_Cnf_List(Cnf)
    case (CI_Lib_Faroald)
      call Faro_Set_Sym(iState,iSym)
      call citrans_conf_occ(nCnf,kf,ncs,dco,sco)
      call mma_deallocate(kf)
      call mma_deallocate(ncs)
      call mma_deallocate(dco)
      call mma_deallocate(sco)
    case (CI_Lib_MCLR_Lucia)
      nCnf = MCLR_Lucia_nCnf(iState,iSym)
    case default
      call Check_Lib('CI_Lib_nCnf')
  end select
end function CI_Lib_nCnf

!> The configurations of symmetry iSym of the CI state iState in the
!> library's own configuration order: nCSFCnf(k) the number of CSFs of
!> configuration k (nCnf = CI_Lib_nCnf(iState,iSym)), iCSFCnf(1:n) their
!> indices in the CSF basis of the library (n = nCSF(iSym)), configuration
!> by configuration
subroutine CI_Lib_Cnf_Layout(iState,iSym,nCnf,nCSFCnf,n,iCSFCnf)
  use sguga, only: SGCnfList, SG_Free_Cnf_List, SG_Mk_Cnf_List
  use citrans, only: citrans_conf_occ
  use stdalloc, only: mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nCnf, n
  integer(kind=iwp), intent(out) :: nCSFCnf(nCnf), iCSFCnf(n)
  integer(kind=iwp) :: i, k, m, ncf
  integer(kind=iwp), allocatable :: dco(:), kf(:), ncs(:), sco(:)
  type(SGCnfList) :: Cnf
  select case (CI_Lib)
    case (CI_Lib_SGUGA)
      call SG_Mk_Cnf_List(iState,iSym,Cnf)
      if ((Cnf%nCnf /= nCnf) .or. (Cnf%iOff(nCnf+1) /= n)) call Abend()
      do k=1,nCnf
        nCSFCnf(k) = Cnf%iOff(k+1)-Cnf%iOff(k)
      end do
      iCSFCnf(1:n) = Cnf%iCSF(1:n)
      call SG_Free_Cnf_List(Cnf)
    case (CI_Lib_Faroald)
      ! (the library basis is the configuration order of citrans)
      call Faro_Set_Sym(iState,iSym)
      call citrans_conf_occ(ncf,kf,ncs,dco,sco)
      if (ncf /= nCnf) call Abend()
      m = 0
      do k=1,nCnf
        nCSFCnf(k) = ncs(k)
        do i=1,ncs(k)
          iCSFCnf(m+i) = kf(k)+i-1
        end do
        m = m+ncs(k)
      end do
      if (m /= n) call Abend()
      call mma_deallocate(kf)
      call mma_deallocate(ncs)
      call mma_deallocate(dco)
      call mma_deallocate(sco)
    case (CI_Lib_MCLR_Lucia)
      call MCLR_Lucia_Cnf_Layout(iState,iSym,nCnf,nCSFCnf,n,iCSFCnf)
    case default
      call Check_Lib('CI_Lib_Cnf_Layout')
  end select
end subroutine CI_Lib_Cnf_Layout

!> The explicit block of H + ECore (total symmetry, the integrals of
!> CI_Lib_Diag) for the configurations IPC(1:nPC) of symmetry iSym, given
!> by their indices in the library's configuration order
!> (CI_Lib_Cnf_Layout): IdxP(1:nP) their CSFs (indices in the CSF basis of
!> the library) in the order of IPC, configuration by configuration
!> (nP <= nPMax), and HP(iTri(i,j)) = <IdxP(i)|H|IdxP(j)> + ECore delta_ij,
!> packed as a lower triangle (nHP >= nTri(nP)).
subroutine CI_Lib_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,IdxP,nHP,HP)
  use sguga, only: SGCnfList, SG_Explicit_Block, SG_Free_Cnf_List, SG_Mk_Cnf_List
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nPC, IPC(nPC), nPMax, nHP
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  integer(kind=iwp) :: i, j, k, m, n, nCnf
  integer(kind=iwp), allocatable :: iCSFCnf(:), kOff(:), nCSFCnf(:)
  real(kind=wp), allocatable :: X(:), Y(:)
  type(SGCnfList) :: Cnf
  n = CI_St(iState)%nCSF(iSym)
  nP = 0
  IdxP(:) = 0
  HP(:) = Zero
  if (CI_Lib == CI_Lib_SGUGA) then
    call SG_Mk_Cnf_List(iState,iSym,Cnf)
    do k=1,nPC
      m = Cnf%iOff(IPC(k)+1)-Cnf%iOff(IPC(k))
      if (nP+m > nPMax) call Abend()
      IdxP(nP+1:nP+m) = Cnf%iCSF(Cnf%iOff(IPC(k))+1:Cnf%iOff(IPC(k)+1))
      nP = nP+m
    end do
    if (nHP < nTri_Elem(nP)) call Abend()
    call SG_Explicit_Block(iState,Cnf,nPC,IPC,ECore,TU,nTU,TUVX,nTUVX,HP)
    call SG_Free_Cnf_List(Cnf)
    return
  end if
  if (CI_Lib == CI_Lib_MCLR_Lucia) then
    call MCLR_Lucia_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,IdxP,nHP,HP)
    return
  end if
  if (CI_Lib /= CI_Lib_Faroald) call Check_Lib('CI_Lib_Block')
  ! Faroald: the columns H e_j of the CSFs of the configurations
  nCnf = CI_Lib_nCnf(iState,iSym)
  call mma_allocate(nCSFCnf,nCnf,Label='nCSFCnf')
  call mma_allocate(kOff,nCnf+1,Label='kOff')
  call mma_allocate(iCSFCnf,n,Label='iCSFCnf')
  call CI_Lib_Cnf_Layout(iState,iSym,nCnf,nCSFCnf,n,iCSFCnf)
  kOff(1) = 0
  do k=1,nCnf
    kOff(k+1) = kOff(k)+nCSFCnf(k)
  end do
  do k=1,nPC
    m = nCSFCnf(IPC(k))
    if (nP+m > nPMax) call Abend()
    IdxP(nP+1:nP+m) = iCSFCnf(kOff(IPC(k))+1:kOff(IPC(k)+1))
    nP = nP+m
  end do
  if (nHP < nTri_Elem(nP)) call Abend()
  call mma_allocate(X,n,Label='X')
  call mma_allocate(Y,n,Label='Y')
  do j=1,nP
    X(:) = Zero
    X(IdxP(j)) = One
    call Faro_Sigma_Act(iState,iSym,iSym,.false.,.true.,nTU,TU,nTUVX,TUVX,n,X,n,Y)
    do i=j,nP
      HP(iTri(i,j)) = Y(IdxP(i))
    end do
    HP(iTri(j,j)) = HP(iTri(j,j))+ECore
  end do
  call mma_deallocate(Y)
  call mma_deallocate(X)
  call mma_deallocate(iCSFCnf)
  call mma_deallocate(kOff)
  call mma_deallocate(nCSFCnf)
end subroutine CI_Lib_Block

!> The P space of the preconditioner: at most nPMax CSFs of symmetry iSym,
!> whole configurations, chosen from the diagonal Diag(1:n) (n =
!> nCSF(iSym), the CSF basis of the library, with the ECore of the caller)
!> by one rule for every library: the configurations by increasing key
!> (the mean of Diag over the CSFs of the configuration) while they fit;
!> when the next one does not fit, the configurations already taken with
!> the same key (within 1.0e-10) are dropped again and the selection
!> stops, so that a degenerate group is taken whole or not at all and the
!> chosen set does not depend on the configuration order of the library.
!> Then the block of H + ECore for them (CI_Lib_Block; IdxP, nP, HP as
!> there).
subroutine CI_Lib_PSpace(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag,nPMax,nP,IdxP,nHP,HP)
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, n, nPMax, nHP
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore, Diag(n)
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  real(kind=wp), parameter :: Deg = 1.0e-10_wp
  integer(kind=iwp) :: k, kMin, m, nCnf, nPC, nSel
  integer(kind=iwp), allocatable :: iCSFCnf(:), IPC(:), nCSFCnf(:)
  logical(kind=iwp), allocatable :: Used(:)
  real(kind=wp) :: KMinV
  real(kind=wp), allocatable :: Key(:)
  nCnf = CI_Lib_nCnf(iState,iSym)
  call mma_allocate(nCSFCnf,nCnf,Label='nCSFCnf')
  call mma_allocate(iCSFCnf,n,Label='iCSFCnf')
  call CI_Lib_Cnf_Layout(iState,iSym,nCnf,nCSFCnf,n,iCSFCnf)
  call mma_allocate(Key,nCnf,Label='Key')
  call mma_allocate(Used,nCnf,Label='Used')
  call mma_allocate(IPC,max(1,nCnf),Label='IPC')
  m = 0
  do k=1,nCnf
    Key(k) = sum(Diag(iCSFCnf(m+1:m+nCSFCnf(k))))/real(nCSFCnf(k),kind=wp)
    m = m+nCSFCnf(k)
  end do
  Used(:) = .false.
  nPC = 0
  nSel = 0
  do
    kMin = 0
    KMinV = huge(KMinV)
    do k=1,nCnf
      if (Used(k)) cycle
      if (Key(k) < KMinV) then
        KMinV = Key(k)
        kMin = k
      end if
    end do
    if (kMin == 0) exit
    if (nSel+nCSFCnf(kMin) > nPMax) then
      ! the degenerate group of the configuration that does not fit
      do while (nPC > 0)
        if (abs(Key(IPC(nPC))-KMinV) > Deg) exit
        nSel = nSel-nCSFCnf(IPC(nPC))
        nPC = nPC-1
      end do
      exit
    end if
    Used(kMin) = .true.
    nPC = nPC+1
    IPC(nPC) = kMin
    nSel = nSel+nCSFCnf(kMin)
  end do
  call CI_Lib_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC(1:max(1,nPC)),nPMax,nP,IdxP,nHP,HP)
  call mma_deallocate(IPC)
  call mma_deallocate(Used)
  call mma_deallocate(Key)
  call mma_deallocate(iCSFCnf)
  call mma_deallocate(nCSFCnf)
end subroutine CI_Lib_PSpace

!> Faroald: a CI vector of symmetry iSym between the order of the files
!> and the configuration order of citrans, in place (Dir = 'C': to citrans,
!> 'O': back; citrans_sort)
subroutine Faro_Sort(iState,iSym,Dir,n,CI)
  use citrans, only: citrans_sort
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, n
  character, intent(in) :: Dir
  real(kind=wp), intent(inout) :: CI(n)
  real(kind=wp), allocatable :: V(:)
  call Faro_Set_Sym(iState,iSym)
  call mma_allocate(V,n,Label='V')
  call citrans_sort(Dir,CI,V)
  CI(:) = V(:)
  call mma_deallocate(V)
end subroutine Faro_Sort

subroutine CI_Lib_Dets_To_Basis(iState,iSym,nD,DetVec,IdxD,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, nD, IdxD(nD), n
  real(kind=wp), intent(in) :: DetVec(nD)
  real(kind=wp), intent(out) :: CI(n)
  if ((CI_Lib == CI_Lib_SGUGA) .or. (CI_Lib == CI_Lib_Faroald)) &
    call Not_Yet('CI_Lib_Dets_To_Basis','determinant vectors (DMRG)')
  call Check_Lib('CI_Lib_Dets_To_Basis')
  call MCLR_Lucia_Dets(iState,iSym,nD,DetVec,IdxD,n,CI)
end subroutine CI_Lib_Dets_To_Basis

subroutine CI_Lib_File_To_Basis(iState,iSym,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, n
  real(kind=wp), intent(inout) :: CI(n)
  ! SGUGA: the split-graph basis is the stored one; Faroald: the
  ! configuration order of citrans (a signed permutation)
  if (CI_Lib == CI_Lib_SGUGA) return
  if (CI_Lib == CI_Lib_Faroald) then
    call Faro_Sort(iState,iSym,'C',n,CI)
    return
  end if
  call Check_Lib('CI_Lib_File_To_Basis')
  call MCLR_Lucia_Conv(iState,iSym,CI_To_Basis,n,CI)
end subroutine CI_Lib_File_To_Basis

subroutine CI_Lib_Basis_To_File(iState,iSym,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, n
  real(kind=wp), intent(inout) :: CI(n)
  ! SGUGA: the split-graph basis is the stored one; Faroald: back from
  ! the configuration order of citrans
  if (CI_Lib == CI_Lib_SGUGA) return
  if (CI_Lib == CI_Lib_Faroald) then
    call Faro_Sort(iState,iSym,'O',n,CI)
    return
  end if
  call Check_Lib('CI_Lib_Basis_To_File')
  call MCLR_Lucia_Conv(iState,iSym,CI_To_File,n,CI)
end subroutine CI_Lib_Basis_To_File

!> The setup of the split graph that the selected library needs ('FULL':
!> with the coupling tables, for SGUGA; 'GRAPH': graph and walks only)
function CI_Lib_Graph_Mode() result(Mode)
  character(len=5) :: Mode
  if (CI_Lib == CI_Lib_SGUGA) then
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

!> SGUGA: sigma for a singlet operator (CI_Lib_Sigma; TU, TUVX in its
!> layout)
subroutine SG_Sigma_Act(iState,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
  use Symmetry_Info, only: Mul
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iCSym, iSSym, nTU, nTUVX, nC, nS
  logical(kind=iwp), intent(in) :: Square, Have_2_el
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), C(nC)
  real(kind=wp), intent(out) :: S(nS)
  integer(kind=iwp) :: iOp, nA, n2
  real(kind=wp), allocatable :: G0(:)
  iOp = Mul(iCSym,iSSym)
  nA = sum(nAsh(1:nSym))
  if (Square) then
    if (nTU /= nA**2) call Abend()
    n2 = nTri_Elem(nA**2)
  else
    if (nTU /= nTri_Elem(nA)) call Abend()
    n2 = nTri_Elem(nTU)
  end if
  if (Have_2_el) then
    if (nTUVX < n2) call Abend()
    call SG_Sigma_Run(iState,iCSym,iOp,Square,nA,nTU,TU,n2,TUVX,nC,C,nS,S)
  else
    call mma_allocate(G0,n2,Label='G0')
    G0(:) = Zero
    call SG_Sigma_Run(iState,iCSym,iOp,Square,nA,nTU,TU,n2,G0,nC,C,nS,S)
    call mma_deallocate(G0)
  end if

end subroutine SG_Sigma_Act

!> SGUGA: the sigma routine for the operator of SG_Sigma_Act (G: the
!> two-body integrals, n2 elements)
subroutine SG_Sigma_Run(iState,iCSym,iOp,Square,nA,nTU,TU,n2,G,nC,C,nS,S)
  use sguga, only: CIS, EXS, SGS
  use SGUGA_Sigma, only: sg_h_psi, sg_h_psi_op, sg_h_psi_ras, sg_h_psi_sq
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iCSym, iOp, nA, nTU, n2, nC, nS
  logical(kind=iwp), intent(in) :: Square
  real(kind=wp), intent(in) :: TU(nTU), G(n2), C(nC)
  real(kind=wp), intent(out) :: S(nS)
  integer(kind=iwp) :: iR, nR
  real(kind=wp), allocatable :: PR(:), SR(:)
  iR = CI_St(iState)%iRel
  if ((iR > 0) .and. (Square .or. (iOp /= 1))) call Not_Yet('CI_Lib_Sigma','GAS with other operators than H')
  if (iR > 0) then
    ! GAS: in the relaxed space of the state (it holds all intermediate
    ! states of the products)
    nR = CIS(iR)%NCSF(iCSym)
    call mma_allocate(PR,nR,Label='PR')
    call mma_allocate(SR,nR,Label='SR')
    call SG_Rel_Embed(iState,iCSym,nC,C,PR)
    call sg_h_psi_ras(SGS(iR),CIS(iR),EXS(iR),PR,nR,iCSym,SR,G,n2,TU,nTU)
    call SG_Rel_Project(iState,iCSym,nS,SR,S)
    call mma_deallocate(SR)
    call mma_deallocate(PR)
  else if (Square) then
    call sg_h_psi_sq(SGS(iState),CIS(iState),EXS(iState),nC,C,iCSym,nS,S,iOp,nA,TU,n2,G)
  else if (iOp /= 1) then
    ! operator of non-total symmetry (perturbations)
    call sg_h_psi_op(SGS(iState),CIS(iState),EXS(iState),nC,C,iCSym,nS,S,iOp,G,n2,TU,nTU)
  else if (CI_St(iState)%nRsPrt > 1) then
    if (nS /= nC) call Abend()
    call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),C,nC,iCSym,S,G,n2,TU,nTU)
  else
    if (nS /= nC) call Abend()
    call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),C,nC,iCSym,S,G,n2,TU,nTU)
  end if
end subroutine SG_Sigma_Run

!> SGUGA: spin-free transition densities, mapped from the level order to
!> the symmetry-ordered active orbitals
subroutine SG_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  use sguga, only: CIS, EXS, SGS, SG_Level_Map
  use SGUGA_Densities, only: sg_trans_pdm, sg_trans_two_pdm
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  integer(kind=iwp) :: a, b, c, d, i, ij, ijL, j, k, kl, klL, l_, N
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: R1(:,:), R2(:), RS(:,:)
  N = sum(nAsh(1:nSym))
  if (SGS(iState)%nLev /= N) call Abend()
  if (n1 < N**2) call Abend()
  call mma_allocate(L2O,N,Label='L2O')
  call mma_allocate(R1,N,N,Label='R1')
  call SG_Level_Map(SGS(iState),L2O)
  if (iOpSpin == 1) then
    ! the spin density (M_S = S, Wigner-Eckart form, sg_trans_pdm), the
    ! states of the same symmetry
    if ((iLSym /= iRSym) .or. (nL /= nR)) call Not_Yet('CI_Lib_Dens','spin densities between symmetries')
    call mma_allocate(RS,N,N,Label='RS')
    call sg_trans_pdm(SGS(iState),CIS(iState),EXS(iState),L,R,nL,iLSym,N,RS,R1,.true.)
    call mma_deallocate(RS)
  else
    call mma_allocate(R2,nTri_Elem(N**2),Label='R2')
    call sg_trans_two_pdm(SGS(iState),CIS(iState),EXS(iState),nL,L,iLSym,nR,R,iRSym,N,R1,size(R2),R2)
  end if
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
  call mma_deallocate(R2,safe='*')
  call mma_deallocate(R1)
  call mma_deallocate(L2O)
end subroutine SG_Trans_Dens

!> SGUGA: the CSF diagonal (exact or configuration average) + ECore
subroutine SG_Diag(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
  use sguga, only: CIS, EXS, SGCnfList, SGS, SG_Cnf_Avg_Diag, SG_Free_Cnf_List, SG_Mk_Cnf_List
  use SGUGA_Sigma, only: sg_hdiag_fast
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nTU, nTUVX, n
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(out) :: Diag(n)
  type(SGCnfList) :: Cnf
  if (Mode == CI_Diag_Exact) then
    call sg_hdiag_fast(SGS(iState),CIS(iState),EXS(iState),n,iSym,TUVX,nTUVX,TU,nTU,Diag)
  else
    call SG_Mk_Cnf_List(iState,iSym,Cnf)
    call SG_Cnf_Avg_Diag(iState,Cnf,TU,nTU,TUVX,nTUVX,Diag)
    call SG_Free_Cnf_List(Cnf)
  end if
  Diag(:) = Diag(:)+ECore
end subroutine SG_Diag

!***********************************************************************
! Faroald (string-driven RAS CI, faroald_ras), in its CSF basis (citrans);
! the vectors of the program stay in the split-graph order of the files,
! sorted to and from citrans for each operation. One CSF space, of the
! symmetry of the CI state.
!***********************************************************************

!> Faroald: the setup for the CI state iState, for every program. The
!> partitions (CAS, RAS, GAS) and their cumulative electron bounds come
!> from the CI space Space_CI when the program has it (RASSCF; MCLR from the
!> CI-space record of the JobIph), otherwise from the RAS input (CI_Input).
!> The orbitals of Faroald are the levels: partition by partition, within
!> each by symmetry; Faro_Lev(iLev) is the symmetry-ordered active orbital
!> of level iLev. The full determinant space is set up here only as the
!> reference of MOLCAS_CI_VERIFY=SIGMA, otherwise later on demand
!> (Faro_Full_Init: other symmetries, non-Hermitian operators).
subroutine Faro_Init(iState)
  use CI_Input, only: iSpin, nActEl, nElec3, nHole1, nRS1, nRS2, nRS3
  use CI_Info, only: Space_CI
  use CI_States, only: RAS_Partitions
  use CI_Verify_Options, only: CI_Verify
  use Orbital_Spaces, only: nAsh, nSym
  use faroald, only: nela, nelb
  use faroald_ras, only: faroald_ras_init, n_ras
  use citrans, only: citrans_ncsf
  use second_quantization, only: binom_coef
  use Molcas, only: MxSym
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: i, iLev, iSym, k, nAbove, nAct, nBelow, nOff(MxSym), nPrt, nRas(MxSym,3), nRasEl(3)
  integer(kind=iwp), allocatable :: MaxEl(:), MinEl(:), nOrbP(:), OrbRas(:), OrbSym(:)
  logical(kind=iwp) :: Full
  Faro_Sym = CI_St(iState)%STSym
  nAct = sum(nAsh(1:nSym))
  call mma_allocate(OrbSym,max(1,nAct),Label='OrbSym')
  call mma_allocate(OrbRas,max(1,nAct),Label='OrbRas')
  if ((Space_CI%nLev == nAct) .and. (Space_CI%nPart >= 1)) then
    ! the CI space of the program (also GAS)
    nPrt = Space_CI%nPart
    OrbSym(1:nAct) = Space_CI%LevSym(1:nAct)
    OrbRas(1:nAct) = Space_CI%LevPart(1:nAct)
    call mma_allocate(MinEl,nPrt,Label='MinEl')
    call mma_allocate(MaxEl,nPrt,Label='MaxEl')
    MinEl(:) = Space_CI%MinEl(1:nPrt)
    MaxEl(:) = Space_CI%MaxEl(1:nPrt)
  else
    ! the RAS input
    call RAS_Partitions(nSym,nActEl,nRs1,nRs2,nRs3,nHole1,nElec3,.false.,nPrt,nRas,nRasEl)
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
    call mma_allocate(MinEl,nPrt,Label='MinEl')
    call mma_allocate(MaxEl,nPrt,Label='MaxEl')
    nBelow = 0
    do k=1,nPrt
      nBelow = nBelow+sum(nRas(1:nSym,k))
      nAbove = nAct-nBelow
      MinEl(k) = max(nRasEl(k),0,nActEl-2*nAbove)
      MaxEl(k) = min(nActEl,2*nBelow)
    end do
    MinEl(nPrt) = nActEl
    MaxEl(nPrt) = nActEl
  end if
  ! levels -> symmetry-ordered active orbitals
  call mma_deallocate(Faro_Lev,safe='*')
  call mma_allocate(Faro_Lev,max(1,nAct),Label='Faro_Lev')
  nOff(1) = 0
  do iSym=2,nSym
    nOff(iSym) = nOff(iSym-1)+nAsh(iSym-1)
  end do
  do iLev=1,nAct
    nOff(OrbSym(iLev)) = nOff(OrbSym(iLev))+1
    Faro_Lev(iLev) = nOff(OrbSym(iLev))
  end do
  ! the tables of the full (C1) determinant space up front only as the
  ! reference of MOLCAS_CI_VERIFY=SIGMA (for CAS always, for RAS/GAS for
  ! small spaces)
  Full = .false.
  if (CI_Verify('SIGMA')) then
    Full = (nPrt == 1)
    if (.not. Full) &
      Full = real(binom_coef((nActEl+iSpin-1)/2,nAct),kind=wp)*real(binom_coef((nActEl-iSpin+1)/2,nAct),kind=wp) <= 1.0e6_wp
  end if
  call faroald_init(nActEl,nAct,iSpin,Full)
  call citrans_init(nActEl,nAct,iSpin,nSym,OrbSym,Faro_Sym,nPrt,OrbRas,MinEl,MaxEl)
  ! the string-driven CI (faroald_ras): the determinant blocks of the CI
  ! symmetry, for singlets with the spin-flip symmetry
  if (nAct > 0) then
    call mma_allocate(nOrbP,nPrt,Label='nOrbP')
    do k=1,nPrt
      nOrbP(k) = count(OrbRas(1:nAct) == k)
    end do
    call faroald_ras_init(nela,nelb,nAct,nPrt,nOrbP,MinEl,MaxEl,nSym,OrbSym(1:nAct),Faro_Sym,iSpin == 1)
    call mma_deallocate(nOrbP)
    write(u6,'(1X,A,I0,A,I0,A)') 'Faroald: string-driven CI, ',nPrt,' partition(s), ',n_ras,' determinants of the CI symmetry'
    if (iSpin == 1) write(u6,'(1X,A)') 'Faroald: singlet, spin-flip symmetry C(Ib,Ia) = C(Ia,Ib) used'
  end if
  ! the description of the CSF spaces, for the other symmetries
  call mma_deallocate(Faro_OrbSym,safe='*')
  call mma_deallocate(Faro_OrbRas,safe='*')
  call mma_deallocate(Faro_MinEl,safe='*')
  call mma_deallocate(Faro_MaxEl,safe='*')
  call mma_allocate(Faro_OrbSym,max(1,nAct),Label='Faro_OrbSym')
  call mma_allocate(Faro_OrbRas,max(1,nAct),Label='Faro_OrbRas')
  call mma_allocate(Faro_MinEl,nPrt,Label='Faro_MinEl')
  call mma_allocate(Faro_MaxEl,nPrt,Label='Faro_MaxEl')
  Faro_OrbSym(:) = OrbSym(:)
  Faro_OrbRas(:) = OrbRas(:)
  Faro_MinEl(:) = MinEl(:)
  Faro_MaxEl(:) = MaxEl(:)
  Faro_nAct = nAct
  Faro_nPrt = nPrt
  Faro_Cur = Faro_Sym
  call mma_deallocate(MaxEl)
  call mma_deallocate(MinEl)
  call mma_deallocate(OrbSym)
  call mma_deallocate(OrbRas)
  if ((citrans_ncsf() /= CI_St(iState)%nCSF(Faro_Sym)) .and. (nActEl > 0)) then
    write(u6,*) 'Faro_Init: CSFs of Faroald and of the split graph',citrans_ncsf(),CI_St(iState)%nCSF(Faro_Sym)
    call Abend()
  end if
  ! no determinant vectors outside the library (a larger dimension of
  ! another library initialized as a reference is kept)
  CI_St(iState)%nDet(:) = max(CI_St(iState)%nDet(:),CI_St(iState)%nCSF(:))
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
  call mma_deallocate(Faro_MinEl,safe='*')
  call mma_deallocate(Faro_MaxEl,safe='*')
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
  call citrans_sym_init(nActEl,Faro_nAct,nSym,Faro_OrbSym,iSym,Faro_nPrt,Faro_OrbRas,Faro_MinEl, &
                        Faro_MaxEl)
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

!> Faroald: htu, gtuvx over the levels from the integrals of CI_Lib_Sigma
!> (Hermitian, full triangles over the symmetry-ordered active orbitals);
!> only the elements allowed by the operator symmetry iOp are taken
subroutine Faro_Ints_Act(iOp,nTU,TU,nTUVX,TUVX,Have_2_el)
  use faroald, only: gtuvx, htu, my_norb
  use Index_Functions, only: iTri, nTri_Elem
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iOp, nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  logical(kind=iwp), intent(in) :: Have_2_el
  integer(kind=iwp) :: it, iu, iv, ix, otu, ovx
  if (nTU /= nTri_Elem(my_norb)) call Abend()
  if (Have_2_el .and. (nTUVX < nTri_Elem(nTU))) call Abend()
  do iu=1,my_norb
    do it=1,my_norb
      otu = iTri(Faro_Lev(it),Faro_Lev(iu))
      htu(it,iu) = TU(otu)
      do ix=1,my_norb
        do iv=1,my_norb
          if (Have_2_el .and. (Faro_Sym4(it,iu,iv,ix) == iOp)) then
            ovx = iTri(Faro_Lev(iv),Faro_Lev(ix))
            gtuvx(it,iu,iv,ix) = TUVX(iTri(otu,ovx))
          else
            gtuvx(it,iu,iv,ix) = Zero
          end if
        end do
      end do
    end do
  end do
end subroutine Faro_Ints_Act

!> Faroald: htu, gtuvx (over the levels) of a non-Hermitian operator of
!> symmetry iOp (Square): H1 = h(t,u), the full nA x nA matrix, G2 =
!> g(tu,vx) at iTri(tu,vx), square pair index tu = (u-1)*N+t, both over
!> the symmetry-ordered active orbitals (the layout of CI_Lib_Sigma)
subroutine Faro_Ints_Sq_Act(iOp,N,H1,nG2,G2,Have_2_el)
  use faroald, only: gtuvx, htu, my_norb
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, Half
  integer(kind=iwp), intent(in) :: iOp, N, nG2
  real(kind=wp), intent(in) :: H1(N,N), G2(nG2)
  logical(kind=iwp), intent(in) :: Have_2_el
  integer(kind=iwp) :: it, iu, iv, ix, otu, ovx
  real(kind=wp), allocatable :: G(:,:,:,:)
  if (N /= my_norb) call Abend()
  if (Have_2_el .and. (nG2 < nTri_Elem(N**2))) call Abend()
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
          if (Have_2_el .and. (Faro_Sym4(it,iu,iv,ix) == iOp)) then
            G(it,iu,iv,ix) = G2(iTri(otu,ovx))
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
end subroutine Faro_Ints_Sq_Act

!> Faroald: the symmetry of the product of the orbitals of levels t, u,
!> v, x
function Faro_Sym4(it,iu,iv,ix) result(iS)
  use Symmetry_Info, only: Mul
  integer(kind=iwp), intent(in) :: it, iu, iv, ix
  integer(kind=iwp) :: iS
  iS = Mul(Mul(Faro_OrbSym(it),Faro_OrbSym(iu)),Mul(Faro_OrbSym(iv),Faro_OrbSym(ix)))
end function Faro_Sym4

!> Faroald: S = H C in the CSF space of the CI state (CI_Lib_Sigma; TU,
!> TUVX in its layout)
subroutine Faro_Sigma_Act(iState,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
  use Symmetry_Info, only: Mul
  use faroald, only: faroald_full, gtuvx, htu, my_norb, ndeta, ndetb, sigma_update
  use CI_Verify_Options, only: CI_Verify
  use faroald_ras, only: n_ras, ras_sigma
  use citrans, only: citrans_csf2ras, citrans_csf2sd, citrans_ras2csf, citrans_sd2csf
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  integer(kind=iwp), intent(in) :: iState, iCSym, iSSym, nTU, nTUVX, nC, nS
  logical(kind=iwp), intent(in) :: Square, Have_2_el
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), C(nC)
  real(kind=wp), intent(out) :: S(nS)
  real(kind=wp), allocatable :: PC(:,:), PS(:,:), VR(:), WR(:)
  if ((iCSym /= Faro_Sym) .or. (iSSym /= Faro_Sym) .or. Square) then
    ! other symmetries or a non-Hermitian operator: in the full determinant
    ! space, P_S H P_C exactly
    call Faro_Full_Init()
    if (Square) then
      call Faro_Ints_Sq_Act(Mul(iCSym,iSSym),my_norb,TU,nTUVX,TUVX,Have_2_el)
    else
      call Faro_Ints_Act(Mul(iCSym,iSSym),nTU,TU,nTUVX,TUVX,Have_2_el)
    end if
    call mma_allocate(PC,ndeta,ndetb,Label='PC')
    call mma_allocate(PS,ndeta,ndetb,Label='PS')
    call Faro_Set_Sym(iState,iCSym)
    call citrans_csf2sd(C,PC)
    PS(:,:) = Zero
    call sigma_update(htu,gtuvx,PS,PC)
    call Faro_Set_Sym(iState,iSSym)
    call citrans_sd2csf(PS,S)
    call mma_deallocate(PS)
    call mma_deallocate(PC)
    return
  end if
  if (nS /= nC) call Abend()
  call Faro_Set_Sym(iState,Faro_Sym)
  call Faro_Ints_Act(1,nTU,TU,nTUVX,TUVX,Have_2_el)
  call mma_allocate(VR,max(1,n_ras),Label='VR')
  call mma_allocate(WR,max(1,n_ras),Label='WR')
  call citrans_csf2ras(C,VR)
  call ras_sigma(htu,gtuvx,VR,WR)
  call citrans_ras2csf(WR,S)
  call mma_deallocate(WR)
  call mma_deallocate(VR)
  ! CI development (MOLCAS_CI_VERIFY=SIGMA): once against the sigma of the
  ! full determinant space (P H P), if its tables are set up
  if ((.not. Faro_RAS_Verified) .and. CI_Verify('SIGMA') .and. faroald_full) then
    call mma_allocate(PC,ndeta,ndetb,Label='PC')
    call mma_allocate(PS,ndeta,ndetb,Label='PS')
    call mma_allocate(VR,nC,Label='SF')
    call citrans_csf2sd(C,PC)
    PS(:,:) = Zero
    call sigma_update(htu,gtuvx,PS,PC)
    call citrans_sd2csf(PS,VR)
    call Report_Dev('Faroald: string-driven sigma vs full determinant space',maxval(abs(VR(:)-S(:))), &
                    1.0e-10_wp*max(One,maxval(abs(VR(:)))),Faro_RAS_Verified)
    call mma_deallocate(VR)
    call mma_deallocate(PS)
    call mma_deallocate(PC)
  end if
end subroutine Faro_Sigma_Act

!> Faroald: the spin-free transition densities (as SG_Trans_Dens), mapped
!> from the levels to the symmetry-ordered active orbitals
subroutine Faro_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  use faroald, only: my_norb, ndeta, ndetb, transition_one_pdm, transition_two_pdm
  use faroald_ras, only: n_ras, ras_trans_one_pdm, ras_trans_two_pdm
  use citrans, only: citrans_csf2ras, citrans_csf2sd
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  integer(kind=iwp) :: ij, it, iu, iv, ix, kl, N
  logical(kind=iwp) :: Full
  real(kind=wp), allocatable :: BR(:), D(:,:), KR(:), P2(:,:,:,:), PL(:,:), PR(:,:), SD(:,:)
  N = my_norb
  if (n1 < N**2) call Abend()
  ! other symmetries: in the full determinant space
  Full = (iLSym /= Faro_Sym) .or. (iRSym /= Faro_Sym)
  if (Full) then
    call Faro_Full_Init()
    call mma_allocate(PL,ndeta,ndetb,Label='PL')
    call mma_allocate(PR,ndeta,ndetb,Label='PR')
    call Faro_Set_Sym(iState,iLSym)
    call citrans_csf2sd(L,PL)
    call Faro_Set_Sym(iState,iRSym)
    call citrans_csf2sd(R,PR)
  else
    call Faro_Set_Sym(iState,Faro_Sym)
    call mma_allocate(BR,max(1,n_ras),Label='BR')
    call mma_allocate(KR,max(1,n_ras),Label='KR')
    call citrans_csf2ras(L,BR)
    call citrans_csf2ras(R,KR)
  end if
  call mma_allocate(D,N,N,Label='D')
  call mma_allocate(SD,N,N,Label='SD')
  if (Full) then
    call transition_one_pdm(PL,PR,D,SD)
  else
    call ras_trans_one_pdm(BR,KR,D,SD)
  end if
  ! spin-free D = da+db or the spin density SD = da-db (M_S = S)
  Rho1(:) = Zero
  do iu=1,N
    do it=1,N
      if (iOpSpin == 1) then
        Rho1(Faro_Lev(it)+(Faro_Lev(iu)-1)*N) = SD(it,iu)
      else
        Rho1(Faro_Lev(it)+(Faro_Lev(iu)-1)*N) = D(it,iu)
      end if
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
end subroutine Faro_Trans_Dens

!> Faroald: the configuration average of the CSF diagonal (as for RASSCF
!> with Faroald: the average of the determinant energies of each
!> configuration, from its occupations) + ECore
subroutine Faro_Diag(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
  use faroald, only: gtuvx, htu, my_norb, nela, nelb
  use citrans, only: citrans_conf_occ
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Two
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nTU, nTUVX, n
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp) :: i, ic, j, na, nb, ncf, nd2, no
  integer(kind=iwp), allocatable :: dco(:), kf(:), ncs(:), sco(:)
  real(kind=wp) :: E, ni, nj, P
  real(kind=wp), allocatable :: DC(:)
  if (Mode /= CI_Diag_Cnf_Avg) call Not_Yet('CI_Lib_Diag','the exact CSF diagonal')
  call Faro_Set_Sym(iState,iSym)
  call Faro_Ints_Act(1,nTU,TU,nTUVX,TUVX,.true.)
  call mma_allocate(DC,n,Label='DC')
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
  Diag(:) = DC(:)+ECore
  call mma_deallocate(DC)
end subroutine Faro_Diag

!> GAS: the CSF map of symmetry LSM of state iState into its relaxed space.
subroutine SG_Rel_Map(iState,LSM)
use sguga, only: CIS, SG_Map_CSFs
use stdalloc, only: mma_allocate
integer(kind=iwp), intent(in) :: iState, LSM
integer(kind=iwp) :: n
if (allocated(Rel_Map(LSM)%M)) return
n = CIS(iState)%NCSF(LSM)
call mma_allocate(Rel_Map(LSM)%M,max(1,n),Label='Rel_Map')
Rel_Map(LSM)%M(:) = 0
call SG_Map_CSFs(iState,CI_St(iState)%iRel,LSM,n,Rel_Map(LSM)%M)
if (any(Rel_Map(LSM)%M(1:n) == 0)) then
  write(u6,*) 'SG_Rel_Map: a CSF of the GAS space is missing in the relaxed space'
  call Abend()
end if
end subroutine SG_Rel_Map

subroutine SG_Rel_Free()
use stdalloc, only: mma_deallocate
integer(kind=iwp) :: k
do k=1,size(Rel_Map)
  call mma_deallocate(Rel_Map(k)%M,safe='*')
end do
end subroutine SG_Rel_Free

!> GAS: X (n CSFs of symmetry LSM of iState) into the relaxed space (XR,
!> zero elsewhere), and back (the components of the GAS CSFs).
subroutine SG_Rel_Embed(iState,LSM,n,X,XR)
use Constants, only: Zero
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: X(n)
real(kind=wp), intent(out) :: XR(:)
call SG_Rel_Map(iState,LSM)
XR(:) = Zero
XR(Rel_Map(LSM)%M(1:n)) = X(1:n)
end subroutine SG_Rel_Embed

subroutine SG_Rel_Project(iState,LSM,n,XR,X)
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: XR(:)
real(kind=wp), intent(out) :: X(n)
call SG_Rel_Map(iState,LSM)
X(1:n) = XR(Rel_Map(LSM)%M(1:n))
end subroutine SG_Rel_Project

!> Report (once, Done) and check a maximal deviation of a verification.
subroutine Report_Dev(Text,Dev,Thr,Done)
character(len=*), intent(in) :: Text
real(kind=wp), intent(in) :: Dev, Thr
logical(kind=iwp), intent(inout) :: Done
if (.not. Done) write(u6,'(1X,2A,ES10.3)') Text,', max. abs. difference',Dev
Done = .true.
if (Dev > Thr) then
  write(u6,'(1X,3A,ES10.3)') 'Verification failed: ',Text,', max. abs. difference',Dev
  call Abend()
end if
end subroutine Report_Dev


!> The remaining case of the operations is the MCLR copy of Lucia (also the
!> stop for RASSCF's Lucia, CI_Lib_Lucia, not available for MCLR)
subroutine Check_Lib(Caller)
  character(len=*), intent(in) :: Caller
  if (CI_Lib /= CI_Lib_MCLR_Lucia) then
    write(u6,'(1X,A,A,I2)') Caller,': no CI library selected, or not available:',CI_Lib
    call Abend()
  end if
end subroutine Check_Lib

end module CI_Operations

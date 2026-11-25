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
use Constants, only: Zero
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
! Lucia (RASSCF's) has been initialized (not for SGCI and FARO, except as the
! reference of MOLCAS_CI_VERIFY=SIGMA)
logical(kind=iwp) :: Lucia_Active = .true.
! the CI state set up by CI_Lib_Initialize from the CI input (0: none); only
! that state is torn down by CI_Lib_Close (a program that sets up its state
! itself, as MCLR, frees it itself)
integer(kind=iwp) :: iState_Input = 0
! SGEX given where it does not apply: warned once
logical(kind=iwp) :: SGEx_Warned = .false.
! HEXS/DEXS (RASSCF keywords, CI_Info Space_Inp): the configurations with
! the maximum occupation (HEXS) or at least the maximum - 1 (DEXS) in the
! listed RAS partitions are removed as in Lucia: their CSFs stay in the
! vectors, the CI works with P H P (zero sigma, zero diagonal there). The
! maximum occupation of partition g is maxa(g) + maxb(g) (Space_Max_Occ), the
! largest numbers of alpha and beta electrons in g in the determinants of the
! space (Lucia's rule for the determinant blocks, I_AM_SO_EXCITED).
! HEXS_Keep: the kept CSFs in the CI basis (SGCI, FARO; Lucia does it
! itself); HEXS_Keep_SG: in the split-graph basis (the SGUGA references).
logical(kind=iwp), allocatable :: HEXS_Keep(:), HEXS_Keep_SG(:)

type Rel_Map_Type
  integer(kind=iwp), allocatable :: M(:)
end type Rel_Map_Type
type(Rel_Map_Type) :: Rel_Map(8)

! Faroald: its string-driven sigma checked once against the full space
logical(kind=iwp) :: Faro_RAS_Verified = .false.

public :: CI_Lib_Input_State, CI_Lucia_Active, HEXS_Keep, HEXS_Keep_SG, HEXS_Proj
public :: CI_Diag_Cnf_Avg, CI_Diag_Exact, CI_Lib, CI_Lib_Basis_To_File, CI_Lib_Block, CI_Lib_Close, CI_Lib_Cnf_Layout, &
          CI_Lib_Dens, CI_Lib_Dets_To_Basis, CI_Lib_Diag, CI_Lib_Faroald, CI_Lib_File_To_Basis, CI_Lib_Graph_Mode, &
          CI_Lib_Initialize, CI_Lib_Level_Map, CI_Lib_Lucia, CI_Lib_MCLR_Lucia, CI_Lib_nCnf, CI_Lib_None, CI_Lib_PSel, &
          CI_Lib_PSpace, CI_Lib_Select, CI_Lib_SGEX, CI_Lib_SGUGA, CI_Lib_Sigma, CI_Lib_Sigma_Det, CI_Lib_Transform, CI_To_Basis, &
          CI_To_File, Cnf_Exc, Faro_Ints_Act, Faro_Lev, Faroald_HDet, Lucia1_Block, Lucia1_Block_SG, Rel_Map, Report_Dev, &
          SG_Rel_Embed, SG_Rel_Free, SG_Rel_Map, SG_Rel_Project

contains

!> Set up the selected library for the CI state. If the program has set up
!> the state itself (CI_Setup_SG/CI_Setup_RAS before the call, CI_St%Active:
!> MCLR), only the library; otherwise (RASSCF) the whole state from the CI
!> input (CI_Lib_Setup_Input: space, graph, Lucia, self-tests, library,
!> HEXS/DEXS masks).
subroutine CI_Lib_Initialize(iState)
  use CI_Input, only: iSpin, nActEl, nElec3, nHole1, nRs1, nRs2, nRs3
  use Orbital_Spaces, only: nAsh, nIsh, nOrb, nSym
  use Molcas, only: MxSym
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: nDet(MxSym)
  if (.not. CI_St(iState)%Active) then
    ! no state set up by the program: the state from the CI input
    call CI_Lib_Setup_Input(iState)
    return
  end if
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
  use sguga, only: SG_Free
  use Lucia_Interface, only: Lucia_Util
  use CI_Space_Def, only: Space_Free
  use CI_Info, only: Space_CI
  use stdalloc, only: mma_deallocate
  integer(kind=iwp), intent(in) :: iState
  if (CI_Lib == CI_Lib_MCLR_Lucia) call MCLR_Lucia_Close(iState)
  if (CI_Lib == CI_Lib_Faroald) call Faro_Close()
  ! the work space of SGUGA (safe to call in any case; also used as the
  ! reference of MOLCAS_CI_VERIFY=SIGMA with the other libraries)
  if ((CI_Lib == CI_Lib_SGUGA) .or. (iState == iState_Input)) then
    call SG_Sigma_Free()
    call SG_Densities_Free()
  end if
  if (iState /= iState_Input) return
  ! the state from the CI input (CI_Lib_Setup_Input): Lucia, the graphs, the
  ! CI space and the HEXS/DEXS masks
  if (Lucia_Active) call Lucia_Util('Close')
  if (CI_St(iState)%iRel > 0) then
    call SG_Free(CI_St(iState)%iRel)
    CI_St(iState)%iRel = 0
  end if
  call SG_Rel_Free()
  call SG_Free(iState)
  CI_St(iState)%Active = .false.
  call mma_deallocate(HEXS_Keep,safe='*')
  call mma_deallocate(HEXS_Keep_SG,safe='*')
  call Space_Free(Space_CI)
  iState_Input = 0
end subroutine CI_Lib_Close

!> The CI state iState was set up by CI_Lib_Initialize from the CI input
!> (not by the program itself): the checks of MOLCAS_CI_VERIFY apply to it.
function CI_Lib_Input_State(iState) result(Is)
  integer(kind=iwp), intent(in) :: iState
  logical(kind=iwp) :: Is
  Is = (iState == iState_Input)
end function CI_Lib_Input_State

!> Lucia (RASSCF's) is initialized (CI_Lib_Initialize): its tables
!> (determinant counts, configuration types) may be used.
function CI_Lucia_Active() result(Active)
  logical(kind=iwp) :: Active
  Active = Lucia_Active
end function CI_Lucia_Active

!> The CI state iState from the CI input (CI_Input, Orbital_Spaces), when the
!> program has not set it up itself: the CI space (Space_CI, with the HEXS/DEXS
!> removals), the split graph (for GAS also the relaxed space), Lucia (as the
!> library or as the reference of MOLCAS_CI_VERIFY=SIGMA), the self-tests of
!> SGUGA on request, the selected library and the HEXS/DEXS masks.
subroutine CI_Lib_Setup_Input(iState)
  use CI_Input, only: iDoGAS, STSym
  use Lucia_Interface, only: Lucia_Util
  use sguga, only: CIS, EXS, SG_Test_EpqErs_Ref, SGS
  use SGUGA_Sigma, only: SG_Check_HPsi_Op, sg_h_psi, sg_h_psi_ras
  use SGUGA_Densities, only: SG_Check_Trans_Two_PDM
  use citrans, only: citrans_ncsf
  use CI_Verify_Options, only: CI_Verify
  use Orbital_Spaces, only: nSym
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: nCSF_Lucia, nFail_TT

  iState_Input = iState
  nCSF_Lucia = 0
  ! SGCI and FARO do not need Lucia; it is set up only as the reference of
  ! MOLCAS_CI_VERIFY=SIGMA.
  Lucia_Active = .not. (((CI_Lib == CI_Lib_SGUGA) .or. (CI_Lib == CI_Lib_Faroald)) .and. (.not. CI_Verify('SIGMA')))

  ! the CI space in the common description of the libraries
  call CI_Space_Setup()

  if (Lucia_Active) then
    ! the determinant tables (CI_St%nCSF, %nDet set by lucia2molcas; Lucia's
    ! CSF count of the CI symmetry kept for the check against the split
    ! graph, which redefines CI_St%nCSF below)
    call Lucia_Util('Ini')
    nCSF_Lucia = CI_St(iState)%nCSF(STSym)
  end if

  if (CI_Verify('SIGMA')) then
    write(u6,'(1X,A)') '**EXPERIMENTAL**'
    write(u6,'(1X,A)') 'CI back end cross-checked against SGUGA (MOLCAS_CI_VERIFY=SIGMA).'
    write(u6,'(1X,A)') '**EXPERIMENTAL**'
  end if
  call CI_Lib_Graph_Setup(iState)
  ! GAS: the split graph against Lucia's CSF count of the CI symmetry (with
  ! GAS Lucia counts the CI symmetry only)
  if (iDoGAS .and. Lucia_Active) then
    if (CIS(iState)%nCSF(STSym) /= nCSF_Lucia) then
      write(u6,'(1X,A,2I10)') 'CI_Lib_Initialize: GAS split graph and Lucia, CSFs of the CI symmetry',CIS(iState)%nCSF(STSym), &
                              nCSF_Lucia
      call Abend()
    end if
    if (CI_Verify('ANY')) write(u6,'(1X,A,I0,A)') 'CI_Lib_Initialize: GAS split graph verified against Lucia, ', &
                                                  CIS(iState)%nCSF(STSym),' CSFs of the CI symmetry'
  end if
  if (.not. Lucia_Active) then
    ! the CSF counts from SGUGA; no determinants (the Davidson scratch for
    ! the determinant expansion of the other libraries is not needed)
    CI_St(iState)%nCSF(:) = 0
    CI_St(iState)%nCSF(1:nSym) = CIS(iState)%nCSF(1:nSym)
    CI_St(iState)%nDet(:) = CI_St(iState)%nCSF(:)
    write(u6,'(1X,A)') 'SGCI/FARO: Lucia is not initialized.'
  end if
  ! Exhaustive verification of the SGUGA coupling-coefficient codes, only on
  ! request (used by the test/csf_sd groups), see CI_Verify. With GAS only
  ! the single operators (E_pq, T_pq), which are exact in any restricted
  ! graph: the products of the self-test sum over the intermediate CSFs of
  ! the graph itself, which is exact for CAS/RAS only (the CI interface
  ! evaluates them for GAS in the relaxed space, checked under SIGMA against
  ! Lucia and Faroald).
  if (CI_Verify('YES') .and. (.not. iDoGAS)) then
    call SG_Test_EpqErs_Ref(iState,STSym,sg_h_psi,sg_h_psi_ras)
    ! the 2-body transition densities (bra /= ket, any symmetries)
    call SG_Check_Trans_Two_PDM(SGS(iState),CIS(iState),EXS(iState),nFail_TT)
    if (nFail_TT > 0) then
      write(u6,'(1X,A,I0)') 'CI_Lib_Initialize: sg_trans_two_pdm check failed: ',nFail_TT
      call Abend()
    end if
    ! sigma for operators of any symmetry
    call SG_Check_HPsi_Op(SGS(iState),CIS(iState),EXS(iState),nFail_TT)
    if (nFail_TT > 0) then
      write(u6,'(1X,A,I0)') 'CI_Lib_Initialize: sg_h_psi_op check failed: ',nFail_TT
      call Abend()
    end if
  else if (CI_Verify('YES') .or. CI_Verify('EPQ')) then
    call SG_Test_EpqErs_Ref(iState,STSym,sg_h_psi,sg_h_psi_ras,Epq_Only=.true.)
  end if

  if (CI_Lib == CI_Lib_Faroald) then
    write(u6,'(1X,A)') '**EXPERIMENTAL**'
    write(u6,'(1X,A)') 'CI backend is FAROALD instead of LUCIA.'
    write(u6,'(1X,A)') '**EXPERIMENTAL**'
    ! the one Faroald setup, from Space_CI
    call Faro_Init(iState)
    if ((citrans_ncsf() /= CI_St(iState)%nCSF(STSym)) .and. (CI_St(iState)%nActEl > 0)) then
      write(u6,*) 'CI_Lib_Initialize: number of CSFs, Faroald and the split graph:',citrans_ncsf(),CI_St(iState)%nCSF(STSym)
      call Abend()
    end if
  end if

  ! HEXS/DEXS: the kept CSFs
  call HEXS_Setup(iState)
end subroutine CI_Lib_Setup_Input

!> The CI space (Space_CI) from the GAS description of the input (NGAS,
!> NGSSH, IGSOCCX; CAS and RAS are converted to it by the input readers) and
!> the HEXS/DEXS lists. With MOLCAS_CI_VERIFY (any level) printed and checked
!> against the RAS description (nRsPrt, nRas, nRasEl, NSM) and the HEXS/DEXS
!> thresholds used by the CI interface.
subroutine CI_Space_Setup()
  use Orbital_Spaces, only: nAsh, nSym
  use CI_Input, only: iDoGAS, IGSOCCX, ISPIN, NACTEL, NGAS, NGSSH, nElec3, nHole1, nRs1, nRs2, nRs3
  use CI_Info, only: Space_CI, Space_Inp
  use CI_Space_Def, only: Space_From_GAS, Space_Max_Occ, Space_Print, Space_Remove, Space_Removals
  use CI_States, only: RAS_Partitions
  use CI_Verify_Options, only: CI_Verify
  use stdalloc, only: mma_allocate, mma_deallocate
  use Molcas, only: MxGAS, MxSym
  integer(kind=iwp) :: g, k, i, iLev, iSym, iPart(MxGAS), nBelow, nAbove, MinR, MaxR, MinEO(MxGAS), nRsPrt, nRas(MxSym,3), &
                       nRasEl(3)
  integer(kind=iwp), allocatable :: NMx(:)
  logical(kind=iwp) :: OK

  call Space_From_GAS(Space_CI,nSym,NACTEL,ISPIN,NGAS,NGSSH(1:NGAS,1:nSym),IGSOCCX(1:NGAS,1:2))
  ! partition of each input space (0: no orbitals, dropped)
  k = 0
  do g=1,NGAS
    iPart(g) = 0
    if (sum(NGSSH(g,1:nSym)) == 0) cycle
    k = k+1
    iPart(g) = k
  end do
  ! HEXS (maximum occupation) and DEXS (maximum - 1) removals
  if (Space_Removals(Space_Inp)) then
    call mma_allocate(NMx,max(1,Space_CI%nPart),Label='NMx')
    call Space_Max_Occ(Space_CI,NMx)
    do i=1,Space_Inp%nHEXS
      g = Space_Inp%iHEXS(i)
      if (iPart(g) > 0) call Space_Remove(Space_CI,iPart(g),NMx(iPart(g)))
    end do
    do i=1,Space_Inp%nDEXS
      g = Space_Inp%iDEXS(i)
      if (iPart(g) > 0) call Space_Remove(Space_CI,iPart(g),NMx(iPart(g))-1)
    end do
    call mma_deallocate(NMx)
  end if
  if (.not. CI_Verify('ANY')) return

  call Space_Print(Space_CI)
  ! GAS: no RAS description to compare with
  if (iDoGAS) return
  ! the RAS description (as CI_Setup_RAS builds it), with its empty
  ! partitions dropped (their lower bound merged into the previous one): the
  ! same partitions, the same lower bounds (after tightening to the physical
  ! ones), no upper bounds, the same levels
  call RAS_Partitions(nSym,NACTEL,nRs1(1:nSym),nRs2(1:nSym),nRs3(1:nSym),nHole1,nElec3,.false.,nRsPrt,nRas,nRasEl)
  OK = (Space_CI%nLev == sum(nAsh(1:nSym)))
  k = 0
  nBelow = 0
  do i=1,nRsPrt
    if (sum(nRas(1:nSym,i)) == 0) then
      if (k > 0) MinEO(k) = max(MinEO(k),nRasEl(i))
      cycle
    end if
    k = k+1
    if (k > Space_CI%nPart) then
      OK = .false.
      exit
    end if
    OK = OK .and. all(Space_CI%nOrb(1:nSym,k) == nRas(1:nSym,i))
    MinEO(k) = nRasEl(i)
  end do
  OK = OK .and. (k == Space_CI%nPart)
  if (OK) then
    do k=1,Space_CI%nPart
      nBelow = nBelow+sum(Space_CI%nOrb(1:nSym,k))
      nAbove = sum(Space_CI%nOrb(1:nSym,k+1:Space_CI%nPart))
      MinR = max(MinEO(k),0,NACTEL-2*nAbove)
      MaxR = min(NACTEL,2*nBelow)
      if (k == Space_CI%nPart) then
        MinR = NACTEL
        MaxR = NACTEL
      end if
      OK = OK .and. (Space_CI%MinEl(k) == MinR) .and. (Space_CI%MaxEl(k) == MaxR)
    end do
    ! the levels: partition by partition, within each by symmetry
    iLev = 0
    do i=1,nRsPrt
      do iSym=1,nSym
        OK = OK .and. all(Space_CI%LevSym(iLev+1:iLev+nRas(iSym,i)) == iSym)
        iLev = iLev+nRas(iSym,i)
      end do
    end do
  end if
  if (.not. OK) then
    write(u6,'(1X,A)') 'CI_Space_Setup: the CI space differs from the RAS description (nRsPrt, nRas, nRasEl)'
    call Abend()
  end if
  write(u6,'(1X,A)') 'CI_Space_Setup: CI space verified against the RAS description'
end subroutine CI_Space_Setup

!> The split graph of the CI state iState from Space_CI (levels partition by
!> partition and within each by symmetry, the identity as orbital -> level
!> map): CAS/RAS as RAS graph (CI_Setup_RAS); GAS with the upper bounds, plus
!> the relaxed space (SG_Rel_State: the bounds widened by one electron at
!> every boundary), in which the CI interface evaluates the products.
subroutine CI_Lib_Graph_Setup(iState)
  use Orbital_Spaces, only: nSym
  use CI_Input, only: iDoGAS, iSpin, nActEl, nElec3, nHole1, nRs1, nRs2, nRs3, STSym
  use CI_States, only: CI_Setup_RAS, CI_Setup_SG
  use CI_Info, only: SG_Rel_State, Space_CI
  use stdalloc, only: mma_allocate, mma_deallocate
  use Molcas, only: MxLev, MxSym
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: k, nP, nLev, iq, Level(MxLev), NSM(MxLev)
  integer(kind=iwp), allocatable :: nRasG(:,:), MinR(:), MaxR(:)

  nLev = Space_CI%nLev
  do iq=1,MxLev
    Level(iq) = iq
  end do
  NSM(:) = 0
  NSM(1:nLev) = Space_CI%LevSym(1:nLev)
  if (iDoGAS) then
    call mma_allocate(nRasG,MxSym,Space_CI%nPart,Label='nRasG')
    nRasG(:,:) = 0
    nRasG(1:nSym,:) = Space_CI%nOrb(1:nSym,1:Space_CI%nPart)
    call CI_Setup_SG(iState,nSym,nActEl,iSpin,STSym,nRasG,Space_CI%MinEl,Space_CI%nPart,Level,Level,nLev,NSM,'FULL', &
                     nRasMx=Space_CI%MaxEl)
    nP = Space_CI%nPart
    call mma_allocate(MinR,nP,Label='MinR')
    call mma_allocate(MaxR,nP,Label='MaxR')
    do k=1,nP
      if (k < nP) then
        MinR(k) = max(0,Space_CI%MinEl(k)-1)
        MaxR(k) = min(nActEl,Space_CI%MaxEl(k)+1)
      else
        MinR(k) = Space_CI%MinEl(k)
        MaxR(k) = Space_CI%MaxEl(k)
      end if
    end do
    call CI_Setup_SG(SG_Rel_State,nSym,nActEl,iSpin,STSym,nRasG,MinR,nP,Level,Level,nLev,NSM,'FULL',nRasMx=MaxR)
    CI_St(iState)%iRel = SG_Rel_State
    call mma_deallocate(MaxR)
    call mma_deallocate(MinR)
    call mma_deallocate(nRasG)
  else
    call CI_Setup_RAS(iState,nSym,nActEl,iSpin,STSym,nRs1(1:nSym),nRs2(1:nSym),nRs3(1:nSym),nHole1,nElec3,.false.,nLev, &
                      Level(1:nLev),Level(1:nLev),NSM(1:nLev),'FULL')
  end if
end subroutine CI_Lib_Graph_Setup

!> HEXS/DEXS: is a configuration with nOcc(k) electrons in partition k of
!> Space_CI removed (the removal rules of the CI space)?
function HEXS_Removed(nOcc) result(Rem)
  use CI_Info, only: Space_CI
  use CI_Space_Def, only: Space_Allowed
  integer(kind=iwp), intent(in) :: nOcc(:)
  logical(kind=iwp) :: Rem
  Rem = .not. Space_Allowed(Space_CI,nOcc(1:Space_CI%nPart))
end function HEXS_Removed

!> HEXS/DEXS: the kept CSFs of symmetry STSym in the split-graph basis, from
!> the SGUGA configurations.
subroutine HEXS_Mask_SG(iState,Keep)
  use sguga, only: SGCnfList, SG_Cnf_Occ, SG_Free_Cnf_List, SG_Mk_Cnf_List
  use CI_Info, only: Space_CI
  use stdalloc, only: mma_allocate, mma_deallocate
  use Molcas, only: MxGAS
  integer(kind=iwp), intent(in) :: iState
  logical(kind=iwp), intent(out) :: Keep(:)
  type(SGCnfList) :: Cnf
  integer(kind=iwp) :: k, iLev, nOcc(MxGAS)
  integer(kind=iwp), allocatable :: Occ(:)
  call SG_Mk_Cnf_List(iState,CI_St(iState)%STSym,Cnf)
  call mma_allocate(Occ,max(1,Cnf%nLev),Label='HEXS_Occ')
  Keep(:) = .true.
  do k=1,Cnf%nCnf
    call SG_Cnf_Occ(Cnf,k,Occ)
    nOcc(:) = 0
    do iLev=1,Cnf%nLev
      nOcc(Space_CI%LevPart(iLev)) = nOcc(Space_CI%LevPart(iLev))+Occ(iLev)
    end do
    if (HEXS_Removed(nOcc)) Keep(Cnf%iCSF(Cnf%iOff(k)+1:Cnf%iOff(k+1))) = .false.
  end do
  call mma_deallocate(Occ)
  call SG_Free_Cnf_List(Cnf)
end subroutine HEXS_Mask_SG

!> HEXS/DEXS: the kept CSFs in the basis of Faroald (citrans), from its
!> configurations (the orbitals of Faroald are the levels).
subroutine HEXS_Mask_Faro(nLev,Keep)
  use citrans, only: citrans_conf_dets
  use CI_Info, only: Space_CI
  use stdalloc, only: mma_deallocate
  use Molcas, only: MxGAS
  integer(kind=iwp), intent(in) :: nLev
  logical(kind=iwp), intent(out) :: Keep(:)
  integer(kind=iwp) :: k, k1, nc, nd, iLev, nOcc(MxGAS)
  integer(kind=iwp), allocatable :: da(:), db(:)
  real(kind=wp), allocatable :: cf(:,:)
  Keep(:) = .true.
  k = 1
  do while (k <= size(Keep))
    call citrans_conf_dets(k,k1,nc,nd,da,db,cf)
    nOcc(:) = 0
    do iLev=1,nLev
      if (btest(da(1),iLev-1)) nOcc(Space_CI%LevPart(iLev)) = nOcc(Space_CI%LevPart(iLev))+1
      if (btest(db(1),iLev-1)) nOcc(Space_CI%LevPart(iLev)) = nOcc(Space_CI%LevPart(iLev))+1
    end do
    if (HEXS_Removed(nOcc)) Keep(k1:k1+nc-1) = .false.
    k = k1+nc
    call mma_deallocate(da)
    call mma_deallocate(db)
    call mma_deallocate(cf)
  end do
end subroutine HEXS_Mask_Faro

!> HEXS/DEXS: the masks of the kept CSFs (end of CI_Lib_Setup_Input). Without
!> Lucia the number of kept CSFs (CI_St%nCSF_Kept) comes from here (the
!> program caps its explicit subspace with it); with Lucia they are compared.
subroutine HEXS_Setup(iState)
  use CI_Input, only: iDoGAS
  use CI_Info, only: Space_CI
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: nKeep, nCSF
  call mma_deallocate(HEXS_Keep,safe='*')
  call mma_deallocate(HEXS_Keep_SG,safe='*')
  ! the removals are part of the CI space (CI_Space_Setup)
  if (Space_CI%nRem <= 0) return
  nCSF = CI_St(iState)%nCSF(CI_St(iState)%STSym)
  ! (SGUGA covers CAS and RAS directly, GAS through the relaxed space)
  if ((.not. iDoGAS) .or. (CI_St(iState)%iRel > 0)) then
    call mma_allocate(HEXS_Keep_SG,nCSF,Label='HEXS_Keep_SG')
    call HEXS_Mask_SG(iState,HEXS_Keep_SG)
  end if
  if (CI_Lib == CI_Lib_SGUGA) then
    call mma_allocate(HEXS_Keep,nCSF,Label='HEXS_Keep')
    HEXS_Keep(:) = HEXS_Keep_SG(:)
  else if (CI_Lib == CI_Lib_Faroald) then
    call mma_allocate(HEXS_Keep,nCSF,Label='HEXS_Keep')
    call HEXS_Mask_Faro(Space_CI%nLev,HEXS_Keep)
  end if
  if (allocated(HEXS_Keep)) then
    nKeep = count(HEXS_Keep)
  else if (allocated(HEXS_Keep_SG)) then
    nKeep = count(HEXS_Keep_SG)
  else
    return
  end if
  write(u6,'(1X,A,I0,A,I0,A)') 'HEXS/DEXS: ',nKeep,' of ',nCSF,' CSFs kept'
  if (.not. Lucia_Active) then
    CI_St(iState)%nCSF_Kept = nKeep
  else if (nKeep /= CI_St(iState)%nCSF_Kept) then
    write(u6,'(1X,A,I0)') 'HEXS/DEXS: Lucia keeps (nCSF_HEXS) ',CI_St(iState)%nCSF_Kept
  end if
end subroutine HEXS_Setup

!> HEXS/DEXS: v := P v (zero on the removed CSFs), if the mask is set and
!> matches v.
subroutine HEXS_Proj(Keep,v)
  logical(kind=iwp), allocatable, intent(in) :: Keep(:)
  real(kind=wp), intent(inout) :: v(:)
  if (.not. allocated(Keep)) return
  if (size(Keep) /= size(v)) return
  where (.not. Keep) v = Zero
end subroutine HEXS_Proj

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
!> HEXS/DEXS: projected (P H P) for the state from the CI input with
!> HEXS/DEXS, in its symmetry; otherwise no projection.
!> Scratch: the libraries keep their own work space (freed by
!> CI_Lib_Close).
subroutine CI_Lib_Sigma(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nTU, nTUVX, nC, nS
  logical(kind=iwp), intent(in) :: Square, Have_2_el
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), C(nC)
  real(kind=wp), intent(out) :: S(nS)
  real(kind=wp), allocatable :: CM(:)
  call Check_Op('CI_Lib_Sigma',iOpSpin,Have_2_el,Square)
  ! HEXS/DEXS (the state from the CI input, its symmetry): P H P, the vector
  ! and the result projected on the kept CSFs (HEXS_Keep; Lucia does it
  ! itself, its mask is not allocated)
  if ((iState == iState_Input) .and. allocated(HEXS_Keep) .and. (iCSym == CI_St(iState)%STSym) .and. &
      (iSSym == CI_St(iState)%STSym)) then
    call mma_allocate(CM,nC,Label='HEXS_CM')
    CM(:) = C(:)
    call HEXS_Proj(HEXS_Keep,CM)
    call Lib_Sigma_Run(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,CM,nS,S)
    call HEXS_Proj(HEXS_Keep,S)
    call mma_deallocate(CM)
  else
    call Lib_Sigma_Run(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
  end if
end subroutine CI_Lib_Sigma

!> The sigma vector S = H C in the determinant basis of RASSCF's Lucia
!> (CASVB), symmetry iSym, nDet determinants; TU/TUVX as for CI_Lib_Sigma.
!> Lucia's determinant sigma is not spin-projected: a C that is not spin
!> adapted gives an S that is not either. Only with RASSCF's Lucia (CASVB
!> selects it). Lucia's sigma is left on disk (Sigma_on_disk) for the
!> densities.
subroutine CI_Lib_Sigma_Det(iState,iSym,nTU,TU,nTUVX,TUVX,nDet,C,S)
  use Lucia_Interface, only: Lucia_Util
  use lucia_data, only: IREFSM, Sigma_on_disk
  use cands, only: ICSM, ISSM
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nDet
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(inout) :: C(nDet)
  real(kind=wp), intent(out) :: S(nDet)
  integer(kind=iwp) :: IREFSM_Save
#include "macros.fh"
  unused_var(iState)
  if (CI_Lib /= CI_Lib_Lucia) call Not_Yet('CI_Lib_Sigma_Det','the sigma vector in the determinant basis (CASVB)')
  IREFSM_Save = IREFSM
  IREFSM = iSym
  ICSM = iSym
  ISSM = iSym
  call Lucia_Util('Sigma',CI_Vector=C(:),Sigma_Vector=S(:),nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)
  ! the sigma vector exists on disk (densi_master)
  Sigma_on_disk = .true.
  IREFSM = IREFSM_Save
  ICSM = IREFSM
  ISSM = IREFSM
end subroutine CI_Lib_Sigma_Det

!> CI_Lib_Sigma for the selected library
subroutine Lib_Sigma_Run(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nTU, nTUVX, nC, nS
  logical(kind=iwp), intent(in) :: Square, Have_2_el
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), C(nC)
  real(kind=wp), intent(out) :: S(nS)
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
end subroutine Lib_Sigma_Run

!> RASSCF's Lucia (CI_Lib_Lucia): S = H C in the symmetry of the state, via
!> its determinant basis (csdtvc), with determinant scratch of this call
subroutine Lucia1_Sigma(iState,iSym,nTU,TU,nTUVX,TUVX,n,C,S)
  use lucia_data, only: IREFSM, Sigma_on_disk
  use cands, only: ICSM, ISSM
  use Lucia_Interface, only: Lucia_Util
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, n
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), C(n)
  real(kind=wp), intent(out) :: S(n)
  integer(kind=iwp) :: IREFSM_Keep, nD, nW
  real(kind=wp), allocatable :: Lucia1_C(:), Lucia1_S(:)
  ! the determinant scratch, only for this call (peak memory as before)
  nD = CI_St(iState)%nDet(iSym)
  nW = max(n,nD)
  call mma_allocate(Lucia1_C,nW,Label='Lucia1_C')
  call mma_allocate(Lucia1_S,nW,Label='Lucia1_S')
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
  ! determinants -> CSFs (csdtvc uses its CSF vector as work of nD
  ! elements: through the scratch)
  call csdtvc(Lucia1_C,Lucia1_S,2,iSym,0)
  S(1:n) = Lucia1_C(1:n)
  call mma_deallocate(Lucia1_C)
  call mma_deallocate(Lucia1_S)
  IREFSM = IREFSM_Keep
  ICSM = IREFSM
  ISSM = IREFSM
end subroutine Lucia1_Sigma

!> RASSCF's Lucia (CI_Lib_Lucia): the densities of one state (L = R, the
!> symmetry of the state) from Lucia's DENSI (packed, in level order),
!> brought to the layout of CI_Lib_Dens. The spin density is Lucia's as it
!> is (for singlets Lucia's spin-combination value).
subroutine Lucia1_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,Pk)
  use lucia_data, only: DStmp, Dtmp, PAtmp, Ptmp
  use Lucia_Interface, only: Lucia_Util
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  logical(kind=iwp), intent(in) :: Pk
  integer(kind=iwp) :: N, nP
  logical(kind=iwp) :: Own
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: CL(:)
  if ((iLSym /= iRSym) .or. (iLSym /= CI_St(iState)%STSym) .or. (nL /= nR)) &
    call Not_Yet('CI_Lib_Dens','RASSCF''s Lucia for densities between symmetries')
  if (.not. all(L(:) == R(:))) then
    ! transition densities (one-body): Lucia's DENSI with two vectors gives
    ! DTMP(p+(q-1)*N) = <L|E_qp|R> over the levels
    if ((I12 /= 1) .or. Pk) call Not_Yet('CI_Lib_Dens','RASSCF''s Lucia for two-body or packed transition densities')
    call Lucia1_TDens(iState,iOpSpin,nL,L,nR,R,n1,Rho1)
    return
  end if
  N = sum(nAsh(1:nSym))
  nP = nTri_Elem(nTri_Elem(N))
  call Lucia1_Work(N,Own)
  call mma_allocate(CL,nL,Label='CL')
  CL(:) = L(:)
  call Lucia_Util('Densi',CI_Vector=CL)
  call mma_deallocate(CL)
  if (Pk) then
    ! Lucia's own layout (only the triangles of DTMP/DSTMP are defined)
    if (iOpSpin == 1) then
      Rho1(1:nTri_Elem(N)) = DSTMP(1:nTri_Elem(N))
    else
      Rho1(1:nTri_Elem(N)) = DTMP(1:nTri_Elem(N))
      if (I12 == 2) then
        Rho2(1:nP) = PTMP(1:nP)
        Rho2(nP+1:2*nP) = PATMP(1:nP)
      end if
    end if
  else
    if (n1 < N**2) call Abend()
    if ((I12 == 2) .and. (n2 < nTri_Elem(N**2))) call Abend()
    call mma_allocate(L2O,N,Label='L2O')
    call CI_Lib_Level_Map(iState,N,L2O)
    if (iOpSpin == 1) then
      call Packed_To_Rho(N,L2O,DSTMP,1,PTMP,PATMP,Rho1,Rho2)
    else
      call Packed_To_Rho(N,L2O,DTMP,I12,PTMP,PATMP,Rho1,Rho2)
    end if
    call mma_deallocate(L2O)
  end if
  if (Own) call Lucia1_Work_Free()
end subroutine Lucia1_Dens

!> RASSCF's Lucia: the one-body transition density (or its spin part)
!> <L|E_ij|R> in the layout of CI_Lib_Dens, from Lucia's DENSI with two
!> vectors (DTMP(p+(q-1)*N) = <L|E_qp|R> over the levels; a permutation)
subroutine Lucia1_TDens(iState,iOpSpin,nL,L,nR,R,n1,Rho1)
  use lucia_data, only: DStmp, Dtmp
  use Lucia_Interface, only: Lucia_Util
  use Orbital_Spaces, only: nAsh, nSym
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iOpSpin, nL, nR, n1
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1)
  integer(kind=iwp) :: N, p, q
  logical(kind=iwp) :: Own
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: CL(:), CR(:)
  N = sum(nAsh(1:nSym))
  if (n1 < N**2) call Abend()
  call Lucia1_Work(N,Own)
  call mma_allocate(CL,nL,Label='CL')
  call mma_allocate(CR,nR,Label='CR')
  CL(:) = L(:)
  CR(:) = R(:)
  call Lucia_Util('Densi',CI_Vector=CR,RVec=CL)
  call mma_deallocate(CR)
  call mma_deallocate(CL)
  call mma_allocate(L2O,N,Label='L2O')
  call CI_Lib_Level_Map(iState,N,L2O)
  do q=1,N
    do p=1,N
      if (iOpSpin == 1) then
        Rho1(L2O(q)+(L2O(p)-1)*N) = DSTMP(p+(q-1)*N)
      else
        Rho1(L2O(q)+(L2O(p)-1)*N) = DTMP(p+(q-1)*N)
      end if
    end do
  end do
  call mma_deallocate(L2O)
  if (Own) call Lucia1_Work_Free()
end subroutine Lucia1_TDens

!> RASSCF's Lucia: the result arrays of DENSI (lucia_data), allocated here
!> unless the program has them already (Own: allocated here)
subroutine Lucia1_Work(N,Own)
  use lucia_data, only: DStmp, Dtmp, PAtmp, Ptmp
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate
  integer(kind=iwp), intent(in) :: N
  logical(kind=iwp), intent(out) :: Own
  integer(kind=iwp) :: nP
  Own = .not. allocated(DTMP)
  if (.not. Own) then
    if (size(DTMP) < N**2) call Abend()
    return
  end if
  nP = nTri_Elem(nTri_Elem(N))
  call mma_allocate(DTMP,N**2,Label='DTMP')
  call mma_allocate(DSTMP,N**2,Label='DSTMP')
  call mma_allocate(PTMP,max(1,nP),Label='PTMP')
  call mma_allocate(PATMP,max(1,nP),Label='PATMP')
end subroutine Lucia1_Work

!> RASSCF's Lucia: free the arrays of Lucia1_Work
subroutine Lucia1_Work_Free()
  use lucia_data, only: DStmp, Dtmp, PAtmp, Ptmp
  use stdalloc, only: mma_deallocate
  call mma_deallocate(DTMP,safe='*')
  call mma_deallocate(DSTMP,safe='*')
  call mma_deallocate(PTMP,safe='*')
  call mma_deallocate(PATMP,safe='*')
end subroutine Lucia1_Work_Free

!> RASSCF's Lucia: the configuration-average diagonal (Lucia's CSDIAG of
!> its determinant diagonal) in the symmetry of the state, plus ECore
subroutine Lucia1_Diag(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
  use Lucia_Interface, only: Lucia_Util
  use lucia_data, only: NCNFTP, NCSFTP, NDTFTP, NTYP, SDREO
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, n
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp) :: nD
  real(kind=wp), allocatable :: DD(:)
  if (iSym /= CI_St(iState)%STSym) call Not_Yet('CI_Lib_Diag','RASSCF''s Lucia in other symmetries')
  nD = CI_St(iState)%nDet(iSym)
  call Lucia_Util('Diag',nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)
  call mma_allocate(DD,max(1,nD),Label='DETDIA')
  call get_diag(DD,nD)
  call CSDIAG(n,nD,Diag,DD,NCNFTP(1,iSym),NTYP,SDREO,NDTFTP,NCSFTP,0)
  call mma_deallocate(DD)
  Diag(:) = Diag(:)+ECore
end subroutine Lucia1_Diag

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
!> folded (the callers symmetrize where they need it).
!> Packed = .true. (optional, default .false.): the densities of one state
!> (L = R, iLSym = iRSym) in the packed layout instead, in the level order
!> of the CI space (CI_Lib_Level_Map): Rho1(iTri(a,b)) = D_ba (or the spin
!> density) for a >= b (n1 >= nTri(nAct)); for I12 = 2 Rho2(1:nP) = P and
!> Rho2(nP+1:2*nP) = PA, nP = nTri(nTri(nAct)), P(tuvx) =
!> 1/2 (p(tuvx)+p(tuxv)) (1/2 p(tuvv) for v = x), PA(tuvx) =
!> 1/2 (p(tuvx)-p(tuxv)) (zero for t = u or v = x), t >= u, v >= x,
!> tu >= vx (n2 >= 2*nP). Libraries that produce it directly (RASSCF's
!> Lucia, SGUGA) fill it without a conversion; the others fold once.
subroutine CI_Lib_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,Packed)
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  logical(kind=iwp), intent(in), optional :: Packed
  integer(kind=iwp) :: N, nP
  logical(kind=iwp) :: Pk
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: R1(:), R2(:)
  Pk = .false.
  if (present(Packed)) Pk = Packed
  if ((iOpSpin /= 0) .and. (iOpSpin /= 1)) then
    write(u6,'(A,I4)') ' CI_Lib_Dens: invalid operator spin iOpSpin =',iOpSpin
    call Abend()
  end if
  if ((iOpSpin == 1) .and. (I12 /= 1)) then
    write(u6,'(A)') ' CI_Lib_Dens: the spin density is one-body only (I12 = 1)'
    call Abend()
  end if
  if (Pk) then
    if ((iLSym /= iRSym) .or. (nL /= nR)) then
      write(u6,'(A)') ' CI_Lib_Dens: the packed layout is for the densities of one state'
      call Abend()
    end if
    N = sum(nAsh(1:nSym))
    nP = nTri_Elem(nTri_Elem(N))
    if ((n1 < nTri_Elem(N)) .or. ((I12 == 2) .and. (n2 < 2*nP))) call Abend()
    ! the libraries with the state densities in this layout
    if (CI_Lib == CI_Lib_Lucia) then
      call Lucia1_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,.true.)
      return
    end if
    if (CI_Lib == CI_Lib_SGUGA) then
      if (all(L(:) == R(:))) then
        call SG_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,.true.)
        return
      end if
    end if
    ! the others: the layout of CI_Lib_Dens, folded once
    call mma_allocate(R1,N**2,Label='R1')
    call mma_allocate(R2,max(1,nTri_Elem(N**2)),Label='R2')
    call Lib_Dens_Rho(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,size(R1),R1,size(R2),R2)
    call mma_allocate(L2O,max(1,N),Label='L2O')
    call CI_Lib_Level_Map(iState,N,L2O)
    if (I12 == 2) then
      call Rho_To_Packed(N,L2O,R1,I12,R2,Rho1,Rho2(1:nP),Rho2(nP+1:2*nP))
    else
      call Rho_To_Packed(N,L2O,R1,I12,R2,Rho1,R2,R2)
    end if
    call mma_deallocate(L2O)
    call mma_deallocate(R2)
    call mma_deallocate(R1)
    return
  end if
  call Lib_Dens_Rho(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
end subroutine CI_Lib_Dens

!> CI_Lib_Dens in its own layout, for the selected library
subroutine Lib_Dens_Rho(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  if (CI_Lib == CI_Lib_SGUGA) then
    call SG_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,.false.)
    return
  end if
  if (CI_Lib == CI_Lib_Faroald) then
    call Faro_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
    return
  end if
  if (CI_Lib == CI_Lib_Lucia) then
    call Lucia1_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,.false.)
    return
  end if
  call Check_Lib('CI_Lib_Dens')
  call MCLR_Lucia_TDens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
end subroutine Lib_Dens_Rho

!> The level order of the CI state iState (the order of the split graph
!> and of the RASSCF densities: partition by partition, within each by
!> symmetry, Space_CI): L2O(iLev) is the symmetry-ordered active orbital of
!> level iLev (n = nAct; the identity without a CI space)
subroutine CI_Lib_Level_Map(iState,n,L2O)
  use CI_Info, only: Space_CI
  use Orbital_Spaces, only: nAsh, nSym
  use Molcas, only: MxSym
  integer(kind=iwp), intent(in) :: iState, n
  integer(kind=iwp), intent(out) :: L2O(n)
  integer(kind=iwp) :: iLev, iSym, nOff(MxSym)
#include "macros.fh"
  unused_var(iState)
  if (Space_CI%nLev /= n) then
    do iLev=1,n
      L2O(iLev) = iLev
    end do
    return
  end if
  nOff(1) = 0
  do iSym=2,nSym
    nOff(iSym) = nOff(iSym-1)+nAsh(iSym-1)
  end do
  do iLev=1,n
    nOff(Space_CI%LevSym(iLev)) = nOff(Space_CI%LevSym(iLev))+1
    L2O(iLev) = nOff(Space_CI%LevSym(iLev))
  end do
end subroutine CI_Lib_Level_Map

!> The densities of a real state from the packed layout of its library
!> (levels, D(iTri(a,b)), P and PA as in fold_two_pdm: P(tuvx) =
!> 1/2 (p(tuvx)+p(tuxv)), 1/2 p(tuvv) for v = x, PA(tuvx) = 1/2
!> (p(tuvx)-p(tuxv)), t >= u, v >= x, tu >= vx) to the layout of
!> CI_Lib_Dens (symmetry-ordered orbitals, L2O(level)). For a real state
!> p(tuvx) = p(vxtu) = p(utxv), which gives every element.
subroutine Packed_To_Rho(n,L2O,D,I12,P,PA,Rho1,Rho2)
  use Index_Functions, only: iTri
  integer(kind=iwp), intent(in) :: n, L2O(n), I12
  real(kind=wp), intent(in) :: D(*), P(*), PA(*)
  real(kind=wp), intent(out) :: Rho1(n,n), Rho2(*)
  integer(kind=iwp) :: ij, kl, t, u, v, x
  do u=1,n
    do t=1,n
      Rho1(L2O(t),L2O(u)) = D(iTri(t,u))
    end do
  end do
  if (I12 /= 2) return
  do x=1,n
    do v=1,n
      kl = L2O(v)+(L2O(x)-1)*n
      do u=1,n
        do t=1,n
          ij = L2O(t)+(L2O(u)-1)*n
          if (kl > ij) cycle
          Rho2(iTri(ij,kl)) = P2_Val(P,PA,t,u,v,x)
        end do
      end do
    end do
  end do
end subroutine Packed_To_Rho

!> The densities of CI_Lib_Dens (symmetry-ordered orbitals: Rho1 square,
!> Rho2 at iTri(ij,kl), square pair index) in the layout of RASSCF, in
!> level order (L2O, CI_Lib_Level_Map): D(iTri(a,b)) = Rho1(b,a) for
!> a >= b; for I12 = 2 also P(tuvx) = 1/2 (p(tuvx)+p(tuxv)) (1/2 p(tuvv)
!> for v = x) and PA(tuvx) = 1/2 (p(tuvx)-p(tuxv)) (zero for t = u or
!> v = x), t >= u, v >= x, tu >= vx (the packing of fold_two_pdm).
subroutine Rho_To_Packed(N,L2O,Rho1,I12,Rho2,D,P,PA)
use Index_Functions, only: iTri
use Constants, only: Zero, Half
integer(kind=iwp), intent(in) :: N, L2O(N), I12
real(kind=wp), intent(in) :: Rho1(N,N), Rho2(*)
real(kind=wp), intent(out) :: D(*)
real(kind=wp), intent(inout) :: P(*), PA(*)
integer(kind=iwp) :: t, tu, tuvx, u, v, vx, x, x_max
real(kind=wp) :: A, B
do t=1,N
  do u=1,t
    D(iTri(t,u)) = Rho1(L2O(u),L2O(t))
  end do
end do
if (I12 /= 2) return
do t=1,N
  do u=1,t
    tu = iTri(t,u)
    do v=1,t
      x_max = v
      if (v == t) x_max = u
      do x=1,x_max
        vx = iTri(v,x)
        tuvx = iTri(tu,vx)
        A = Rho2(iTri(L2O(t)+(L2O(u)-1)*N,L2O(v)+(L2O(x)-1)*N))
        if (v == x) then
          P(tuvx) = Half*A
          PA(tuvx) = Zero
        else
          B = Rho2(iTri(L2O(t)+(L2O(u)-1)*N,L2O(x)+(L2O(v)-1)*N))
          P(tuvx) = Half*(A+B)
          PA(tuvx) = Zero
          if (t /= u) PA(tuvx) = Half*(A-B)
        end if
      end do
    end do
  end do
end do

end subroutine Rho_To_Packed

!> Packed_To_Rho: p(abcd) of a real state from the stored element of
!> (t,u,v,x) = (a,b,max(c,d),min(c,d)) with t >= u and tu >= vx, using
!> p(abcd) = p(badc) = p(cdab)
function P2_Val(P,PA,t0,u0,v0,x0) result(Val)
  use Index_Functions, only: iTri
  use Constants, only: Two
  real(kind=wp), intent(in) :: P(*), PA(*)
  integer(kind=iwp), intent(in) :: t0, u0, v0, x0
  real(kind=wp) :: Val
  integer(kind=iwp) :: a, b, c, d, tuvx
  a = t0
  b = u0
  c = v0
  d = x0
  call First_Pair(a,b,c,d)
  if (iTri(a,b) < iTri(max(c,d),min(c,d))) then
    call Swap_Int(a,c)
    call Swap_Int(b,d)
    call First_Pair(a,b,c,d)
  end if
  tuvx = iTri(iTri(a,b),iTri(max(c,d),min(c,d)))
  if (c == d) then
    Val = Two*P(tuvx)
  else if (c > d) then
    Val = P(tuvx)+PA(tuvx)
  else
    Val = P(tuvx)-PA(tuvx)
  end if
end function P2_Val

!> p(abcd) = p(badc): the first pair with a >= b
subroutine First_Pair(a,b,c,d)
  integer(kind=iwp), intent(inout) :: a, b, c, d
  if (a < b) then
    call Swap_Int(a,b)
    call Swap_Int(c,d)
  end if
end subroutine First_Pair

subroutine Swap_Int(i,j)
  integer(kind=iwp), intent(inout) :: i, j
  integer(kind=iwp) :: k
  k = i
  i = j
  j = k
end subroutine Swap_Int

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
    case (CI_Lib_Lucia)
      if (Mode /= CI_Diag_Cnf_Avg) call Not_Yet('CI_Lib_Diag','the exact diagonal with RASSCF''s Lucia')
      call Lucia1_Diag(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
    case (CI_Lib_MCLR_Lucia)
      call MCLR_Lucia_Diag_Act(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
    case default
      call Check_Lib('CI_Lib_Diag')
  end select
  ! HEXS/DEXS (the state from the CI input): zero for the removed CSFs, as
  ! Lucia's (which does it itself)
  if ((iState == iState_Input) .and. (iSym == CI_St(iState)%STSym) .and. &
      ((CI_Lib == CI_Lib_SGUGA) .or. (CI_Lib == CI_Lib_Faroald))) call HEXS_Proj(HEXS_Keep,Diag)
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
  use lucia_data, only: NCNFTP, NTYP
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
    case (CI_Lib_Lucia)
      ! RASSCF's Lucia: its configurations of all types (NCNFTP)
      nCnf = sum(NCNFTP(1:NTYP,iSym))
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
  use lucia_data, only: NCNFTP, NCSFTP, NTYP
  use stdalloc, only: mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nCnf, n
  integer(kind=iwp), intent(out) :: nCSFCnf(nCnf), iCSFCnf(n)
  integer(kind=iwp) :: i, iTyp, k, m, ncf
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
    case (CI_Lib_Lucia)
      ! RASSCF's Lucia: type by type, the CSFs of a configuration contiguous
      m = 0
      k = 0
      do iTyp=1,NTYP
        do i=1,NCNFTP(iTyp,iSym)
          k = k+1
          nCSFCnf(k) = NCSFTP(iTyp)
        end do
      end do
      if (k /= nCnf) call Abend()
      do i=1,n
        iCSFCnf(i) = i
      end do
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
  use CI_Input, only: DoSGEx
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nPC, IPC(nPC), nPMax, nHP
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  integer(kind=iwp) :: k, m, n, nCnf
  integer(kind=iwp), allocatable :: iCSFCnf(:), kOff(:), nCSFCnf(:), nSelCnf(:)
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
    if (CI_St(iState)%iRel > 0) then
      ! GAS: the matrix elements in the relaxed space of the state
      call SG_Rel_Map(iState,iSym)
      call SG_Explicit_Block(iState,Cnf,nPC,IPC,ECore,TU,nTU,TUVX,nTUVX,HP,iRel=CI_St(iState)%iRel,RelMap=Rel_Map(iSym)%M)
    else
      call SG_Explicit_Block(iState,Cnf,nPC,IPC,ECore,TU,nTU,TUVX,nTUVX,HP)
    end if
    call SG_Free_Cnf_List(Cnf)
    return
  end if
  if (CI_Lib == CI_Lib_MCLR_Lucia) then
    call MCLR_Lucia_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,IdxP,nHP,HP)
    return
  end if
  ! Faroald and RASSCF's Lucia: the CSFs of the configurations from the
  ! layout
  nCnf = CI_Lib_nCnf(iState,iSym)
  call mma_allocate(nCSFCnf,max(1,nCnf),Label='nCSFCnf')
  call mma_allocate(kOff,nCnf+1,Label='kOff')
  call mma_allocate(iCSFCnf,max(1,n),Label='iCSFCnf')
  call CI_Lib_Cnf_Layout(iState,iSym,nCnf,nCSFCnf,n,iCSFCnf)
  kOff(1) = 0
  do k=1,nCnf
    kOff(k+1) = kOff(k)+nCSFCnf(k)
  end do
  call mma_allocate(nSelCnf,max(1,nPC),Label='nSelCnf')
  do k=1,nPC
    m = nCSFCnf(IPC(k))
    if (nP+m > nPMax) call Abend()
    IdxP(nP+1:nP+m) = iCSFCnf(kOff(IPC(k))+1:kOff(IPC(k)+1))
    nSelCnf(k) = m
    nP = nP+m
  end do
  if (nHP < nTri_Elem(nP)) call Abend()
  ! SGEX (RASSCF keyword, the state from the CI input): the block of SGUGA
  ! for RASSCF's Lucia; not with GAS without the relaxed space, nor Faroald
  if (DoSGEx .and. (iState == iState_Input) .and. (.not. SGEx_Warned)) then
    if (.not. SG_Block_Avail(iState)) then
      write(u6,'(1X,A)') 'CI_Lib_Block: SGEX ignored, SGUGA is not available with GAS; using the CI back end'
      SGEx_Warned = .true.
    else if (CI_Lib == CI_Lib_Faroald) then
      write(u6,'(1X,A)') 'CI_Lib_Block: SGEX ignored with FARO; the explicit block is computed by Faroald'
      SGEx_Warned = .true.
    end if
  end if
  select case (CI_Lib)
    case (CI_Lib_Faroald)
      call Faro_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,nSelCnf,nP,IdxP,HP)
    case (CI_Lib_Lucia)
      if (CI_Lib_SGEX(iState)) then
        call Lucia1_Block_SG(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,nSelCnf,nP,IdxP,HP)
      else
        call Lucia1_Block(iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nP,HP)
      end if
    case default
      call Check_Lib('CI_Lib_Block')
  end select
  call mma_deallocate(nSelCnf)
  call mma_deallocate(iCSFCnf)
  call mma_deallocate(kOff)
  call mma_deallocate(nCSFCnf)
end subroutine CI_Lib_Block

!> Faroald: the explicit block for the configurations of the CSFs
!> IdxP(1:nP) (citrans order, nSelCnf(k) CSFs of configuration k): for each
!> pair of configurations L, R the determinant block <D_L|H|D_R>
!> (Faroald_HDet, Slater-Condon) transformed to the CSFs, C_L^T H C_R;
!> pairs more than a double excitation apart are skipped (Cnf_Exc)
subroutine Faro_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,nSelCnf,nP,IdxP,H)
  use citrans, only: citrans_conf_dets
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nPC, nSelCnf(nPC), nP, IdxP(nP)
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(inout) :: H(*)
  integer(kind=iwp) :: I, IIL, IILB, IIR, IIRB, IIRMAX, ICNL, ICNR, J, k1, kCnf, nc, nCTot, nd, nDTot, ndL, ndR, nL, nR
  integer(kind=iwp), allocatable :: DA(:), da1(:), DB(:), db1(:), DOcc(:), iCOff(:), iDOff(:), nDCnf(:), SOcc(:)
  real(kind=wp), allocatable :: Blk(:,:), C(:), cf(:,:), HDet(:,:), T(:,:)
  call Faro_Set_Sym(iState,iSym)
  call Faro_Ints_Act(1,nTU,TU,nTUVX,TUVX,.true.)
  ! determinants and CSF coefficients of the configurations
  call mma_allocate(nDCnf,nPC,Label='nDCnf')
  call mma_allocate(iDOff,nPC,Label='iDOff')
  call mma_allocate(iCOff,nPC,Label='iCOff')
  call mma_allocate(DOcc,nPC,Label='DOcc')
  call mma_allocate(SOcc,nPC,Label='SOcc')
  nDTot = 0
  nCTot = 0
  IILB = 0
  do kCnf=1,nPC
    call citrans_conf_dets(IdxP(IILB+1),k1,nc,nd,da1,db1,cf)
    if (nc /= nSelCnf(kCnf)) call Abend()
    nDCnf(kCnf) = nd
    ! doubly and singly occupied orbitals of the configuration
    DOcc(kCnf) = iand(da1(1),db1(1))
    SOcc(kCnf) = ieor(da1(1),db1(1))
    iDOff(kCnf) = nDTot
    iCOff(kCnf) = nCTot
    nDTot = nDTot+nd
    nCTot = nCTot+nd*nc
    IILB = IILB+nc
    call mma_deallocate(da1)
    call mma_deallocate(db1)
    call mma_deallocate(cf)
  end do
  call mma_allocate(DA,nDTot,Label='DA')
  call mma_allocate(DB,nDTot,Label='DB')
  call mma_allocate(C,nCTot,Label='C')
  IILB = 0
  do kCnf=1,nPC
    call citrans_conf_dets(IdxP(IILB+1),k1,nc,nd,da1,db1,cf)
    DA(iDOff(kCnf)+1:iDOff(kCnf)+nd) = da1(:)
    DB(iDOff(kCnf)+1:iDOff(kCnf)+nd) = db1(:)
    do I=1,nc
      J = IdxP(IILB+I)-k1+1
      if ((J < 1) .or. (J > nc)) call Abend()
      C(iCOff(kCnf)+(I-1)*nd+1:iCOff(kCnf)+I*nd) = cf(:,J)
    end do
    IILB = IILB+nc
    call mma_deallocate(da1)
    call mma_deallocate(db1)
    call mma_deallocate(cf)
  end do
  ! blocks C_L^T <D_L|H|D_R> C_R
  nd = maxval(nDCnf)
  nc = maxval(nSelCnf)
  call mma_allocate(HDet,nd,nd,Label='HDet')
  call mma_allocate(T,nd,nc,Label='T')
  call mma_allocate(Blk,nc,nc,Label='Blk')
  IILB = 1
  do ICNL=1,nPC
    nL = nSelCnf(ICNL)
    ndL = nDCnf(ICNL)
    IIRB = 1
    do ICNR=1,ICNL
      nR = nSelCnf(ICNR)
      ndR = nDCnf(ICNR)
      if (Cnf_Exc(DOcc(ICNL),SOcc(ICNL),DOcc(ICNR),SOcc(ICNR)) > 2) then
        Blk(1:nL,1:nR) = Zero
      else
        do J=1,ndR
          do I=1,ndL
            HDet(I,J) = Faroald_HDet(DA(iDOff(ICNL)+I),DB(iDOff(ICNL)+I),DA(iDOff(ICNR)+J),DB(iDOff(ICNR)+J))
          end do
        end do
        call DGEMM_('N','N',ndL,nR,ndR,One,HDet,nd,C(iCOff(ICNR)+1:iCOff(ICNR)+ndR*nR),ndR,Zero,T,nd)
        call DGEMM_('T','N',nL,nR,ndL,One,C(iCOff(ICNL)+1:iCOff(ICNL)+ndL*nL),ndL,T,nd,Zero,Blk,nc)
      end if
      do IIL=1,nL
        IIRMAX = merge(IIL,nR,IILB == IIRB)
        do IIR=1,IIRMAX
          H(nTri_Elem(IILB-1+IIL-1)+IIRB-1+IIR) = Blk(IIL,IIR)
        end do
      end do
      IIRB = IIRB+nR
    end do
    IILB = IILB+nL
  end do
  do I=1,nP
    H(nTri_Elem(I)) = H(nTri_Elem(I))+ECore
  end do
  call mma_deallocate(Blk)
  call mma_deallocate(T)
  call mma_deallocate(HDet)
  call mma_deallocate(C)
  call mma_deallocate(DB)
  call mma_deallocate(DA)
  call mma_deallocate(SOcc)
  call mma_deallocate(DOcc)
  call mma_deallocate(iCOff)
  call mma_deallocate(iDOff)
  call mma_deallocate(nDCnf)
end subroutine Faro_Block

!> RASSCF's Lucia: the explicit block for its configurations IPC(1:nPC)
!> (Lucia's configuration order of symmetry iSym), CNHCN for each pair
!> (pairs more than a double excitation apart skipped)
subroutine Lucia1_Block(iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nP,H)
  use lucia_data, only: Conf_Occ, DFTP, DTOC, IREOTS, MINOP, NCSFTP, NDTFTP, NTYP
  use csfbas, only: NAEL, NBEL
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  integer(kind=iwp), intent(in) :: iSym, nTU, nTUVX, nPC, IPC(nPC), nP
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(inout) :: H(*)
  integer(kind=iwp) :: DL, DR, ICNL, ICNR, IIL, IILB, IIR, IIRB, IIRMAX, ILTYP, IRTYP, MXCS, MXDT, NAC, NCSFL, NCSFR, NEL, &
                       nScr, NTEST, SL, SR
  integer(kind=iwp), allocatable :: ICONFL(:), ICONFR(:)
  real(kind=wp), allocatable :: CNHCNM(:), HONE(:,:), Scr(:)
#include "macros.fh"
  unused_var(nP)
  NAC = sum(nAsh(1:nSym))
  NEL = NAEL+NBEL
  NTEST = 0
  MXDT = maxval(NDTFTP(1:NTYP))
  MXCS = maxval(NCSFTP(1:NTYP))
  ! one-body matrix in square storage, as CNHCN wants it
  call mma_allocate(HONE,NAC,NAC,Label='HONE')
  call SQUARE(TU,HONE,NAC,1,NAC)
  ! scratch of CNHCN
  nScr = 2*MXDT*(NEL+1)+NEL*(MXDT+1)+MXDT*MXDT+4*NAC+MXDT+MXCS*MXDT+64
  call mma_allocate(Scr,nScr,Label='IHJ_Scr')
  call mma_allocate(CNHCNM,MXCS*MXCS,Label='CNHCNM')
  call mma_allocate(ICONFL,NEL,Label='ICONFL')
  call mma_allocate(ICONFR,NEL,Label='ICONFR')
  IILB = 1
  do ICNL=1,nPC
    call GETCNF(ICONFL,ILTYP,IPC(ICNL),Conf_Occ(iSym)%A,iSym,NEL)
    NCSFL = NCSFTP(ILTYP)
    call Occ_Bits(ICONFL,ILTYP-1+MINOP,NEL,DL,SL)
    IIRB = 1
    do ICNR=1,ICNL
      call GETCNF(ICONFR,IRTYP,IPC(ICNR),Conf_Occ(iSym)%A,iSym,NEL)
      NCSFR = NCSFTP(IRTYP)
      call Occ_Bits(ICONFR,IRTYP-1+MINOP,NEL,DR,SR)
      if (Cnf_Exc(DL,SL,DR,SR) > 2) then
        CNHCNM(1:NCSFL*NCSFR) = Zero
      else
        call CNHCN(ICONFL,ILTYP,ICONFR,IRTYP,CNHCNM,Scr,NAEL,NBEL,ECore,HONE,DFTP,DTOC,NAC,TUVX,NTEST,One,IREOTS)
      end if
      do IIL=1,NCSFL
        IIRMAX = merge(IIL,NCSFR,IILB == IIRB)
        do IIR=1,IIRMAX
          H(nTri_Elem(IILB-1+IIL-1)+IIRB-1+IIR) = CNHCNM((IIR-1)*NCSFL+IIL)
        end do
      end do
      IIRB = IIRB+NCSFR
    end do
    IILB = IILB+NCSFL
  end do
  call mma_deallocate(ICONFR)
  call mma_deallocate(ICONFL)
  call mma_deallocate(CNHCNM)
  call mma_deallocate(Scr)
  call mma_deallocate(HONE)
end subroutine Lucia1_Block

!> SGEX (RASSCF keyword): the explicit block of RASSCF's Lucia for the
!> state from the CI input comes from SGUGA (Lucia1_Block_SG)
function CI_Lib_SGEX(iState) result(Use_SG)
  use CI_Input, only: DoSGEx
  integer(kind=iwp), intent(in) :: iState
  logical(kind=iwp) :: Use_SG
  Use_SG = DoSGEx .and. (CI_Lib == CI_Lib_Lucia) .and. (iState == iState_Input) .and. SG_Block_Avail(iState)
end function CI_Lib_SGEX

!> SGUGA can evaluate matrix elements of the CI space of iState: CAS and RAS
!> directly, GAS through the relaxed space
function SG_Block_Avail(iState) result(Avail)
  use CI_Input, only: iDoGAS
  integer(kind=iwp), intent(in) :: iState
  logical(kind=iwp) :: Avail
  Avail = (.not. iDoGAS) .or. (CI_St(iState)%iRel > 0)
end function SG_Block_Avail

!> RASSCF's Lucia with SGEX: <I|H|J> = P_I P_J <I_G|H|J_G> + delta_IJ ECore
!> over the CSFs IdxP(1:nP) (Lucia order, nSelCnf(k) of configuration k),
!> with I_G the GUGA CSF of the Lucia CSF I and P_I the phase of SG_ReOrd
!> (mode 0: c_G(I_G) = P_I c(I)); the block <I_G|H|J_G> from SG_Block_H
!> (pairs of configurations more than a double excitation apart skipped;
!> GAS: in the relaxed space, the CSFs by their walks)
subroutine Lucia1_Block_SG(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,nSelCnf,nP,IdxP,H)
  use SG_Lucia_Bridge, only: SG_Lucia_Walk
  use sguga, only: EXS, SG_Block_H, SG_MatEl_Init, SG_NUM, SGS
  use Molcas, only: MxAct
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nPC, nSelCnf(nPC), nP, IdxP(nP)
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(inout) :: H(*)
  integer(kind=iwp) :: I, iS, J
  integer(kind=iwp) :: IWALK(MxAct)
  integer(kind=iwp), allocatable :: iG(:), iPG(:)
  iS = iState
  if (CI_St(iState)%iRel > 0) iS = CI_St(iState)%iRel
  call SG_MatEl_Init(iS,iSym)
  call mma_allocate(iG,nP,Label='iG')
  call mma_allocate(iPG,nP,Label='iPG')
  do J=1,nP
    call SG_Lucia_Walk(iState,iSym,IdxP(J),IWALK,iPG(J))
    iG(J) = SG_NUM(SGS(iS),EXS(iS),IWALK)
  end do
  call SG_Block_H(iSym,nPC,nSelCnf,nP,iG,TU,nTU,TUVX,nTUVX,H)
  do I=1,nP
    do J=1,I
      H(nTri_Elem(I-1)+J) = iPG(I)*iPG(J)*H(nTri_Elem(I-1)+J)
    end do
    H(nTri_Elem(I)) = H(nTri_Elem(I))+ECore
  end do
  call mma_deallocate(iPG)
  call mma_deallocate(iG)
end subroutine Lucia1_Block_SG

!> Doubly (D) and singly (S) occupied orbitals, as bit strings (orbital t
!> in bit t-1), of the configuration ICONF with nOp open shells (GETCNF
!> form: the closed orbitals first, then the open ones), nEl electrons
subroutine Occ_Bits(ICONF,nOp,nEl,D,S)
  integer(kind=iwp), intent(in) :: ICONF(*), nOp, nEl
  integer(kind=iwp), intent(out) :: D, S
  integer(kind=iwp) :: I, ICL
  ICL = (nEl-nOp)/2
  D = 0
  S = 0
  do I=1,ICL
    D = ibset(D,ICONF(I)-1)
  end do
  do I=ICL+1,ICL+nOp
    S = ibset(S,ICONF(I)-1)
  end do
end subroutine Occ_Bits

!> The configurations of the P space of the preconditioner: at most nPMax
!> CSFs of symmetry iSym, whole configurations, chosen from the diagonal
!> Diag(1:n) (n = nCSF(iSym), the CSF basis of the library) by one rule for
!> every library and program: the configurations by increasing key while
!> they fit, the key being the diagonal value of the configuration if its
!> CSFs share it (configuration average) and otherwise the mean over its
!> CSFs (ties within 1.0e-13 in the configuration order); when the next one
!> does not fit, the configurations already taken with the same key
!> (within 1.0e-10) are dropped again and the selection stops, so that a
!> degenerate group is taken whole or not at all; only if that would empty
!> the P space (the first group alone exceeds nPMax) the group is kept as
!> far as it fits, i.e. split, with a warning that names the size needed
!> to take it whole. nSel <= nPMax holds in every case (checked).
!> IPC(1:nPC) (intent(out), dimension nCnf = CI_Lib_nCnf) are the chosen
!> configurations in the library's configuration order
!> (CI_Lib_Cnf_Layout), nSel their number of CSFs.
subroutine CI_Lib_PSel(iState,iSym,n,Diag,nPMax,nCnf,nPC,IPC,nSel)
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  integer(kind=iwp), intent(in) :: iState, iSym, n, nPMax, nCnf
  real(kind=wp), intent(in) :: Diag(n)
  integer(kind=iwp), intent(out) :: nPC, IPC(nCnf), nSel
  real(kind=wp), parameter :: Acc = 1.0e-13_wp, Deg = 1.0e-10_wp
  integer(kind=iwp) :: k, kMin, m, nKeep, nSelKeep
  integer(kind=iwp), allocatable :: iCSFCnf(:), nCSFCnf(:)
  logical(kind=iwp), allocatable :: Used(:)
  real(kind=wp) :: KMinV, XMax, XMn, XMx
  real(kind=wp), allocatable :: Key(:)
  call mma_allocate(nCSFCnf,max(1,nCnf),Label='nCSFCnf')
  call mma_allocate(iCSFCnf,max(1,n),Label='iCSFCnf')
  call CI_Lib_Cnf_Layout(iState,iSym,nCnf,nCSFCnf,n,iCSFCnf)
  call mma_allocate(Key,max(1,nCnf),Label='Key')
  call mma_allocate(Used,max(1,nCnf),Label='Used')
  m = 0
  do k=1,nCnf
    ! (all equal: that value exactly, so that ties are not decided by rounding)
    XMn = minval(Diag(iCSFCnf(m+1:m+nCSFCnf(k))))
    XMx = maxval(Diag(iCSFCnf(m+1:m+nCSFCnf(k))))
    if (XMn == XMx) then
      Key(k) = XMn
    else
      Key(k) = sum(Diag(iCSFCnf(m+1:m+nCSFCnf(k))))/real(nCSFCnf(k),kind=wp)
    end if
    m = m+nCSFCnf(k)
  end do
  XMax = Zero
  if (nCnf > 0) XMax = maxval(Key(1:nCnf))
  Used(:) = .false.
  nPC = 0
  nSel = 0
  do
    kMin = 0
    KMinV = XMax+One
    do k=1,nCnf
      if (Used(k)) cycle
      if (Key(k)+Acc < KMinV) then
        KMinV = Key(k)
        kMin = k
      end if
    end do
    if (kMin == 0) exit
    if (nSel+nCSFCnf(kMin) > nPMax) then
      ! the degenerate group of the configuration that does not fit
      nKeep = nPC
      nSelKeep = nSel
      do while (nPC > 0)
        if (abs(Key(IPC(nPC))-KMinV) > Deg) exit
        nSel = nSel-nCSFCnf(IPC(nPC))
        nPC = nPC-1
      end do
      if (nPC == 0) then
        nPC = nKeep
        nSel = nSelKeep
        ! the first degenerate group alone exceeds nPMax: it is split
        m = nSel
        do k=1,nCnf
          if ((.not. Used(k)) .and. (abs(Key(k)-KMinV) <= Deg)) m = m+nCSFCnf(k)
        end do
        call WarningMessage(1,'The P space of the CI preconditioner splits a degenerate group of configurations.')
        write(u6,'(1X,A,I0,A,I0,A)') 'The explicit subspace holds at most ',nPMax,' CSFs, the lowest degenerate group has ', &
                                     m,' CSFs;'
        write(u6,'(1X,A,I0,A)') 'increase the size of the explicit subspace (RASSCF: keyword SDAV) to at least ',m, &
                                ' to take it whole.'
      end if
      exit
    end if
    Used(kMin) = .true.
    nPC = nPC+1
    IPC(nPC) = kMin
    nSel = nSel+nCSFCnf(kMin)
  end do
  if (nSel > nPMax) then
    write(u6,'(A,2I10)') ' CI_Lib_PSel: P space larger than allowed',nSel,nPMax
    call Abend()
  end if
  call mma_deallocate(Used)
  call mma_deallocate(Key)
  call mma_deallocate(iCSFCnf)
  call mma_deallocate(nCSFCnf)
end subroutine CI_Lib_PSel

!> The P space of the preconditioner: the configurations of CI_Lib_PSel and
!> the block of H + ECore for them (CI_Lib_Block; IdxP, nP, HP as there).
subroutine CI_Lib_PSpace(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag,nPMax,nP,IdxP,nHP,HP)
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, n, nPMax, nHP
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore, Diag(n)
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  integer(kind=iwp) :: nCnf, nPC, nSel
  integer(kind=iwp), allocatable :: IPC(:)
  nCnf = CI_Lib_nCnf(iState,iSym)
  call mma_allocate(IPC,max(1,nCnf),Label='IPC')
  call CI_Lib_PSel(iState,iSym,n,Diag,nPMax,nCnf,nPC,IPC,nSel)
  call CI_Lib_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC(1:max(1,nPC)),nPMax,nP,IdxP,nHP,HP)
  call mma_deallocate(IPC)
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

!> RASSCF's Lucia: a CI vector between the split-graph order of the files
!> and Lucia's CSF order (SG_ReOrd: iWay = 1 file -> Lucia, 0 back)
subroutine Lucia1_Sort(iState,iSym,iWay,n,CI)
  use SG_Lucia_Bridge, only: SG_ReOrd
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, iWay, n
  real(kind=wp), intent(inout) :: CI(n)
  real(kind=wp), allocatable :: V(:)
  call mma_allocate(V,n,Label='V')
  call SG_ReOrd(iState,iSym,iWay,n,CI,V)
  CI(:) = V(:)
  call mma_deallocate(V)
end subroutine Lucia1_Sort

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
  if (CI_Lib == CI_Lib_Lucia) then
    call Lucia1_Sort(iState,iSym,1,n,CI)
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
  if (CI_Lib == CI_Lib_Lucia) then
    call Lucia1_Sort(iState,iSym,0,n,CI)
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
subroutine SG_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,Pk)
  use sguga, only: CIS
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  logical(kind=iwp), intent(in) :: Pk
  integer(kind=iwp) :: iR, nLR, nRR
  logical(kind=iwp) :: Same
  real(kind=wp), allocatable :: LR(:), RR(:)
  ! the densities of one state (L = R) in the packed layout: the
  ! state-density routines
  Same = Pk .and. (iLSym == iRSym) .and. (nL == nR)
  if (Same) Same = all(L(:) == R(:))
  iR = CI_St(iState)%iRel
  if (iR > 0) then
    ! GAS: in the relaxed space of the state (it holds all intermediate
    ! states of the products)
    nLR = CIS(iR)%NCSF(iLSym)
    nRR = CIS(iR)%NCSF(iRSym)
    call mma_allocate(LR,nLR,Label='LR')
    call mma_allocate(RR,nRR,Label='RR')
    call SG_Rel_Embed(iState,iLSym,nL,L,LR)
    call SG_Rel_Embed(iState,iRSym,nR,R,RR)
    call SG_Dens_Core(iR,Same,Pk,I12,iOpSpin,iLSym,iRSym,nLR,LR,nRR,RR,n1,Rho1,n2,Rho2)
    call mma_deallocate(RR)
    call mma_deallocate(LR)
  else
    call SG_Dens_Core(iState,Same,Pk,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  end if
end subroutine SG_Trans_Dens

!> SGUGA: the densities of SG_Trans_Dens in the graph of the state jS (the
!> state or its relaxed space); Same: L and R are the same state
subroutine SG_Dens_Core(iState,Same,Pk,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  use sguga, only: CIS, EXS, SGS, SG_Level_Map
  use SGUGA_Densities, only: sg_one_pdm, sg_spin_pdm, sg_trans_pdm, sg_trans_two_pdm, sg_two_pdm
  use Orbital_Spaces, only: nAsh, nSym
  use Index_Functions, only: iTri, nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  logical(kind=iwp), intent(in) :: Same, Pk
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  integer(kind=iwp) :: a, b, c, d, i, ij, ijL, j, k, kl, klL, l_, N
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: DT(:), PA(:), PT(:), R1(:,:), R2(:), RS(:,:), SDT(:)
  N = sum(nAsh(1:nSym))
  if (SGS(iState)%nLev /= N) call Abend()
  if ((.not. Pk) .and. (n1 < N**2)) call Abend()
  call mma_allocate(L2O,N,Label='L2O')
  call SG_Level_Map(SGS(iState),L2O)
  if (Same) then
    ! one state: the state densities (packed, levels) and their layout
    call mma_allocate(DT,nTri_Elem(N),Label='DT')
    call sg_one_pdm(SGS(iState),CIS(iState),EXS(iState),L,nL,iLSym,DT,size(DT))
    if (iOpSpin == 1) then
      call mma_allocate(SDT,nTri_Elem(N),Label='SDT')
      call sg_spin_pdm(SGS(iState),CIS(iState),EXS(iState),L,nL,iLSym,DT,SDT,size(DT))
      if (Pk) then
        Rho1(1:size(SDT)) = SDT(:)
      else
        call Packed_To_Rho(N,L2O,SDT,1,SDT,SDT,Rho1,Rho2)
      end if
      call mma_deallocate(SDT)
    else if (I12 == 2) then
      call mma_allocate(PT,nTri_Elem(nTri_Elem(N)),Label='PT')
      call mma_allocate(PA,nTri_Elem(nTri_Elem(N)),Label='PA')
      call sg_two_pdm(SGS(iState),CIS(iState),EXS(iState),L,nL,iLSym,PT,PA,size(PT))
      if (Pk) then
        Rho1(1:size(DT)) = DT(:)
        Rho2(1:size(PT)) = PT(:)
        Rho2(size(PT)+1:2*size(PT)) = PA(:)
      else
        if (n2 < nTri_Elem(N**2)) call Abend()
        call Packed_To_Rho(N,L2O,DT,2,PT,PA,Rho1,Rho2)
      end if
      call mma_deallocate(PA)
      call mma_deallocate(PT)
    else if (Pk) then
      Rho1(1:size(DT)) = DT(:)
    else
      call Packed_To_Rho(N,L2O,DT,1,DT,DT,Rho1,Rho2)
    end if
    call mma_deallocate(DT)
    call mma_deallocate(L2O)
    return
  end if
  call mma_allocate(R1,N,N,Label='R1')
  if (iOpSpin == 1) then
    ! the spin density (M_S = S, Wigner-Eckart form, sg_trans_pdm), the
    ! states of the same symmetry
    if ((iLSym /= iRSym) .or. (nL /= nR)) call Not_Yet('CI_Lib_Dens','spin densities between symmetries')
    call mma_allocate(RS,N,N,Label='RS')
    call sg_trans_pdm(SGS(iState),CIS(iState),EXS(iState),L,R,nL,iLSym,N,RS,R1,.true.)
    call mma_deallocate(RS)
  else if ((I12 == 1) .and. (iLSym == iRSym) .and. (nL == nR)) then
    ! one-body, one symmetry: sg_trans_pdm (no two-body part)
    call mma_allocate(RS,N,N,Label='RS')
    call sg_trans_pdm(SGS(iState),CIS(iState),EXS(iState),L,R,nL,iLSym,N,R1,RS,.false.)
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
end subroutine SG_Dens_Core

!> SGUGA: the CSF diagonal (exact or configuration average) + ECore
subroutine SG_Diag(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
  use sguga, only: CIS, EXS, SGCnfList, SGS, SG_Cnf_Avg_Diag, SG_Free_Cnf_List, SG_Mk_Cnf_List
  use SGUGA_Sigma, only: sg_hdiag_fast
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nTU, nTUVX, n
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp) :: iR, nR
  real(kind=wp), allocatable :: DR(:)
  type(SGCnfList) :: Cnf
  if ((Mode == CI_Diag_Exact) .and. (CI_St(iState)%iRel > 0)) then
    ! GAS: the exact diagonal in the relaxed space of the state
    iR = CI_St(iState)%iRel
    nR = CIS(iR)%NCSF(iSym)
    call mma_allocate(DR,nR,Label='DR')
    call sg_hdiag_fast(SGS(iR),CIS(iR),EXS(iR),nR,iSym,TUVX,nTUVX,TU,nTU,DR)
    call SG_Rel_Project(iState,iSym,n,DR,Diag)
    call mma_deallocate(DR)
  else if (Mode == CI_Diag_Exact) then
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
    call faroald_ras_init(nela,nelb,nAct,nPrt,nOrbP,MinEl,MaxEl,nSym,OrbSym(1:nAct),Faro_Sym,iSpin == 1,.not. Faro_Reproducible())
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

!> Faroald: inside the verification (MOLCAS_TEST set, as for the seed of
!> GetSeed) its alpha-beta sigma path comes from the cost estimate, so that
!> the results do not depend on timings; otherwise the fastest is timed
function Faro_Reproducible() result(Repro)
  logical(kind=iwp) :: Repro
  character(len=16) :: Env
  Env = ' '
  call get_environment_variable('MOLCAS_TEST',Env)
  Repro = Env /= ' '
end function Faro_Reproducible

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


!> Excitation level between two configurations given by their doubly (D)
!> and singly (S) occupied orbitals as bit strings: the number of electrons
!> to move, Sum_t max(0,n_t(L)-n_t(R)). H couples only configurations with
!> an excitation level of at most 2.
function Cnf_Exc(DL,SL,DR,SR) result(nExc)
integer(kind=iwp) :: nExc
integer(kind=iwp), intent(in) :: DL, SL, DR, SR

! n_t(L)-n_t(R) = 2 (doubly -> empty), 1 (doubly -> singly, singly -> empty)
nExc = 2*popcnt(iand(DL,not(ior(DR,SR))))+popcnt(iand(DL,SR))+popcnt(iand(SL,not(ior(DR,SR))))

end function Cnf_Exc

!> <aL bL|H|aR bR> between Faroald determinants (alpha and beta strings as
!> bit strings, orbital t in bit t-1, alpha string before the beta string)
!> by the Slater-Condon rules, with the integrals htu and gtuvx (chemists'
!> notation) of Faroald, H = sum_tu h_tu E_tu + 1/2 sum_tuvx g_tuvx
!> (E_tu E_vx - delta_uv E_tx), as in sigma_update.
function Faroald_HDet(aL,bL,aR,bR) result(E)
use faroald, only: gtuvx, htu, my_norb
use Constants, only: Zero, Half

real(kind=wp) :: E
integer(kind=iwp), intent(in) :: aL, bL, aR, bR
integer(kind=iwp) :: dA, dB, nA, nB, i, j, p, q, r, s, aI

dA = ieor(aL,aR)
dB = ieor(bL,bR)
nA = popcnt(dA)/2
nB = popcnt(dB)/2
E = Zero
if (nA+nB > 2) return

if (nA+nB == 0) then
  do i=1,my_norb
    if (btest(aL,i-1)) then
      E = E+htu(i,i)
      do j=1,my_norb
        if (btest(aL,j-1)) E = E+Half*(gtuvx(i,i,j,j)-gtuvx(i,j,j,i))
        if (btest(bL,j-1)) E = E+gtuvx(i,i,j,j)
      end do
    end if
    if (btest(bL,i-1)) then
      E = E+htu(i,i)
      do j=1,my_norb
        if (btest(bL,j-1)) E = E+Half*(gtuvx(i,i,j,j)-gtuvx(i,j,j,i))
      end do
    end if
  end do
else if (nA == 1 .and. nB == 0) then
  call Single(aL,aR,bR)
else if (nA == 0 .and. nB == 1) then
  call Single(bL,bR,aR)
else if (nA == 2) then
  call Double(aL,aR)
else if (nB == 2) then
  call Double(bL,bR)
else
  ! one alpha (p <- q) and one beta (r <- s) excitation
  p = Orb(iand(dA,aL),1)
  q = Orb(iand(dA,aR),1)
  r = Orb(iand(dB,bL),1)
  s = Orb(iand(dB,bR),1)
  E = Par(aR,p,q)*Par(bR,r,s)*gtuvx(p,q,r,s)
end if

contains

  ! single excitation p <- q in the string xL/xR of one spin, yR the string
  ! of the other spin
  subroutine Single(xL,xR,yR)
    integer(kind=iwp), intent(in) :: xL, xR, yR
    integer(kind=iwp) :: dX, iC, j

    dX = ieor(xL,xR)
    p = Orb(iand(dX,xL),1)
    q = Orb(iand(dX,xR),1)
    iC = iand(xL,xR)
    E = htu(p,q)
    do j=1,my_norb
      if (btest(iC,j-1)) E = E+gtuvx(p,q,j,j)-gtuvx(p,j,j,q)
      if (btest(yR,j-1)) E = E+gtuvx(p,q,j,j)
    end do
    E = Par(xR,p,q)*E

  end subroutine Single

  ! double excitation p,r <- q,s in the string xL/xR of one spin
  subroutine Double(xL,xR)
    integer(kind=iwp), intent(in) :: xL, xR
    integer(kind=iwp) :: dX

    dX = ieor(xL,xR)
    p = Orb(iand(dX,xL),1)
    r = Orb(iand(dX,xL),2)
    q = Orb(iand(dX,xR),1)
    s = Orb(iand(dX,xR),2)
    aI = ibset(ibclr(xR,q-1),p-1)
    E = Par(xR,p,q)*Par(aI,r,s)*(gtuvx(p,q,r,s)-gtuvx(p,s,r,q))

  end subroutine Double

  ! the n-th orbital (bit) set in x
  function Orb(x,n) result(k)
    integer(kind=iwp) :: k
    integer(kind=iwp), intent(in) :: x, n
    integer(kind=iwp) :: m

    m = 0
    do k=1,my_norb
      if (btest(x,k-1)) then
        m = m+1
        if (m == n) return
      end if
    end do
    k = 0

  end function Orb

  ! sign of a+_t a_u acting on the string x: the parity of the number of
  ! orbitals occupied in x strictly between t and u
  function Par(x,t,u) result(iSgn)
    integer(kind=iwp) :: iSgn
    integer(kind=iwp), intent(in) :: x, t, u
    integer(kind=iwp) :: lo, hi

    lo = min(t,u)
    hi = max(t,u)
    iSgn = 1-2*poppar(ibits(x,lo,hi-lo-1))

  end function Par

end function Faroald_HDet

!> The CI vector CI (n CSFs of symmetry iSym, in the library basis) for
!> orbitals transformed by T: the active orbital transformation of each
!> symmetry as a square block (nAsh(iSym) x nAsh(iSym), column-major), the
!> blocks in symmetry order (sum nAsh**2 <= nTU), as in the CI
!> transformation of Malmqvist (sequence of single-orbital steps).
subroutine CI_Lib_Transform(iState,iSym,nTU,TU,n,CI)
  use SGUGA, only: SG_TRACI
  use Lucia_Interface, only: Lucia_Util
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, n
  real(kind=wp), intent(in) :: TU(nTU)
  real(kind=wp), intent(inout) :: CI(n)
  select case (CI_Lib)
    case (CI_Lib_SGUGA)
      call SG_TRACI(iState,iSym,TU,nTU,CI,n)
    case (CI_Lib_Faroald)
      call Faro_Transform(iState,iSym,nTU,TU,n,CI)
    case (CI_Lib_Lucia)
      ! (the vector goes through Lucia's determinant basis inside)
      call Lucia_Util('TRACI',CI_Vector=CI,TU=TU,nTU=nTU)
    case (CI_Lib_MCLR_Lucia)
      call Not_Yet('CI_Lib_Transform','the CI transformation with MCLR''s Lucia')
    case default
      call Check_Lib('CI_Lib_Transform')
  end select
end subroutine CI_Lib_Transform

!> Faroald: the CI transformation (CI_Lib_Transform) through the
!> determinants, with the transformation matrix over its orbitals (levels)
subroutine Faro_Transform(iState,iSym,nTU,TU,n,CI)
  use faroald, only: ndeta, ndetb, traci_faroald
  use faroald_ras, only: n_ras, ras_active, ras_traci
  use citrans, only: citrans_csf2ras, citrans_csf2sd, citrans_ras2csf, citrans_sd2csf
  use Molcas, only: MxSym
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, n
  real(kind=wp), intent(in) :: TU(nTU)
  real(kind=wp), intent(inout) :: CI(n)
  integer(kind=iwp) :: nA, iLev, jLev, iO, jO, iS, nAshS(MxSym), iOffO(MxSym), iOffT(MxSym)
  real(kind=wp), allocatable :: TF(:), PsiR(:), PsiD(:,:)
  nA = Faro_nAct
  ! the symmetry and the offsets of the orbitals and of the blocks of T
  nAshS(:) = 0
  do iLev=1,nA
    nAshS(Faro_OrbSym(iLev)) = nAshS(Faro_OrbSym(iLev))+1
  end do
  iOffO(1) = 0
  iOffT(1) = 0
  do iS=2,MxSym
    iOffO(iS) = iOffO(iS-1)+nAshS(iS-1)
    iOffT(iS) = iOffT(iS-1)+nAshS(iS-1)**2
  end do
  if (iOffT(MxSym)+nAshS(MxSym)**2 > nTU) then
    write(u6,*) 'Faro_Transform: T too short',nTU
    call Abend()
  end if
  call mma_allocate(TF,max(1,nA*nA),Label='TF')
  TF(:) = Zero
  do jLev=1,nA
    jO = Faro_Lev(jLev)
    iS = Faro_OrbSym(jLev)
    do iLev=1,nA
      if (Faro_OrbSym(iLev) /= iS) cycle
      iO = Faro_Lev(iLev)
      TF(iLev+(jLev-1)*nA) = TU(iOffT(iS)+iO-iOffO(iS)+(jO-iOffO(iS)-1)*nAshS(iS))
    end do
  end do
  call Faro_Set_Sym(iState,iSym)
  if (ras_active .and. (iSym == Faro_Sym)) then
    ! in the string classes of the CI symmetry (faroald_ras)
    call mma_allocate(PsiR,max(1,n_ras),Label='PsiR')
    call citrans_csf2ras(CI,PsiR)
    call ras_traci(TF,nA*nA,PsiR)
    call citrans_ras2csf(PsiR,CI)
    call mma_deallocate(PsiR)
  else
    call Faro_Full_Init()
    call mma_allocate(PsiD,ndeta,ndetb,Label='PsiD')
    PsiD(:,:) = Zero
    call citrans_csf2sd(CI,PsiD)
    call traci_faroald(TF,nA*nA,PsiD)
    call citrans_sd2csf(PsiD,CI)
    call mma_deallocate(PsiD)
  end if
  call mma_deallocate(TF)
end subroutine Faro_Transform

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

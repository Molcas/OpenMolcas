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
!>  - CI_Lib_File_To_Basis, CI_Lib_Basis_To_File: a CI vector of symmetry
!>    iSym between the stored basis (JobIph) and the CSF basis of the
!>    selected library, in place.
!> All vectors passed to the gateway are in the CSF basis of the selected
!> library, with their natural length (nCSF of their symmetry); each
!> library keeps its own internal basis, ordering and work space.
!>
!> The MCLR copy of Lucia still keeps its data in MCLR (MCLR_Data), which
!> the gateway must not use; MCLR registers its routines with
!> CI_Lib_Register (procedure pointers) until that library is
!> self-contained. The other libraries are connected behind the same
!> interfaces.
module CI_Operations

use Definitions, only: wp, iwp, u6

implicit none
private

!> The CI libraries (Lucia, SGUGA and Faroald with the values of RASSCF's
!> CI_Basis(); always compare against these names)
integer(kind=iwp), parameter :: CI_Lib_None = 0, CI_Lib_Lucia = 1, CI_Lib_SGUGA = 2, CI_Lib_Faroald = 3, CI_Lib_MCLR_Lucia = 4

!> Kinds of diagonal (CI_Diag)
integer(kind=iwp), parameter :: CI_Diag_Cnf_Avg = 1, CI_Diag_Exact = 2

!> Directions of the basis conversion (Conv_Proc)
integer(kind=iwp), parameter :: CI_To_File = 0, CI_To_Basis = 1

!> The selected library (changed only by CI_Lib_Select)
integer(kind=iwp), protected :: CI_Lib_Op = CI_Lib_None

abstract interface
  subroutine Sigma_Proc(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
    import :: wp, iwp
    integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nInt1, nInt2, nC, nS
    real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
    logical(kind=iwp), intent(in) :: Have_2_el, Square
    real(kind=wp), intent(in) :: C(nC)
    real(kind=wp), intent(out) :: S(nS)
  end subroutine Sigma_Proc
  subroutine TDens_Proc(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
    import :: wp, iwp
    integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
    real(kind=wp), intent(in) :: L(nL), R(nR)
    real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  end subroutine TDens_Proc
  subroutine Diag_Proc(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
    import :: wp, iwp
    integer(kind=iwp), intent(in) :: iState, Mode, iSym, nInt1, nInt2, n
    real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
    real(kind=wp), intent(in) :: ECore
    real(kind=wp), intent(out) :: Diag(n)
  end subroutine Diag_Proc
  subroutine State_Proc(iState)
    import :: iwp
    integer(kind=iwp), intent(in) :: iState
  end subroutine State_Proc
  !> iMode = CI_To_Basis: stored basis -> library basis; CI_To_File: the reverse
  subroutine Conv_Proc(iState,iSym,iMode,n,CI)
    import :: wp, iwp
    integer(kind=iwp), intent(in) :: iState, iSym, iMode, n
    real(kind=wp), intent(inout) :: CI(n)
  end subroutine Conv_Proc
end interface

procedure(Sigma_Proc), pointer :: MCLR_Lucia_Sigma => null()
procedure(TDens_Proc), pointer :: MCLR_Lucia_TDens => null()
procedure(Diag_Proc), pointer :: MCLR_Lucia_Diag => null()
procedure(Conv_Proc), pointer :: MCLR_Lucia_Conv => null()
procedure(State_Proc), pointer :: MCLR_Lucia_Init => null(), MCLR_Lucia_Close => null()

public :: CI_Diag, CI_Diag_Cnf_Avg, CI_Diag_Exact, CI_Lib_Basis_To_File, CI_Lib_Close, CI_Lib_Faroald, CI_Lib_File_To_Basis, &
          CI_Lib_Initialize, CI_Lib_Lucia, CI_Lib_MCLR_Lucia, CI_Lib_None, CI_Lib_Op, CI_Lib_Register, CI_Lib_Select, &
          CI_Lib_SGUGA, CI_Sigma_Op, CI_To_Basis, CI_To_File, CI_Trans_Dens, Conv_Proc, Diag_Proc, Sigma_Proc, State_Proc, &
          TDens_Proc

contains

!> Select the CI library for MCLR (CI_Lib_*)
subroutine CI_Lib_Select(iLib)
  integer(kind=iwp), intent(in) :: iLib
  select case (iLib)
    case (CI_Lib_MCLR_Lucia)
      CI_Lib_Op = iLib
    case (CI_Lib_Lucia,CI_Lib_SGUGA,CI_Lib_Faroald)
      write(u6,'(A,I2,A)') ' CI_Lib_Select: CI library',iLib,' is not available for MCLR yet'
      call Abend()
    case default
      write(u6,'(A,I4)') ' CI_Lib_Select: unknown CI library',iLib
      call Abend()
  end select
end subroutine CI_Lib_Select

!> MCLR registers the routines of its copy of Lucia
subroutine CI_Lib_Register(Init,Close,Sigma,TDens,Diag,Conv)
  procedure(State_Proc) :: Init, Close
  procedure(Sigma_Proc) :: Sigma
  procedure(TDens_Proc) :: TDens
  procedure(Diag_Proc) :: Diag
  procedure(Conv_Proc) :: Conv
  MCLR_Lucia_Sigma => Sigma
  MCLR_Lucia_TDens => TDens
  MCLR_Lucia_Diag => Diag
  MCLR_Lucia_Conv => Conv
  MCLR_Lucia_Init => Init
  MCLR_Lucia_Close => Close
end subroutine CI_Lib_Register

!> Set up the selected library for the CI state
subroutine CI_Lib_Initialize(iState)
  integer(kind=iwp), intent(in) :: iState
  call Check_Lib('CI_Lib_Initialize',associated(MCLR_Lucia_Init))
  call MCLR_Lucia_Init(iState)
end subroutine CI_Lib_Initialize

!> Free the library of the CI state, forget the registered routines and
!> the selection (safe to call in any case)
subroutine CI_Lib_Close(iState)
  integer(kind=iwp), intent(in) :: iState
  if ((CI_Lib_Op == CI_Lib_MCLR_Lucia) .and. associated(MCLR_Lucia_Close)) call MCLR_Lucia_Close(iState)
  MCLR_Lucia_Sigma => null()
  MCLR_Lucia_TDens => null()
  MCLR_Lucia_Diag => null()
  MCLR_Lucia_Conv => null()
  MCLR_Lucia_Init => null()
  MCLR_Lucia_Close => null()
  CI_Lib_Op = CI_Lib_None
end subroutine CI_Lib_Close

subroutine CI_Sigma_Op(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nInt1, nInt2, nC, nS
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  logical(kind=iwp), intent(in) :: Have_2_el, Square
  real(kind=wp), intent(in) :: C(nC)
  real(kind=wp), intent(out) :: S(nS)
  call Check_Lib('CI_Sigma_Op',associated(MCLR_Lucia_Sigma))
  if ((iOpSpin /= 0) .and. (iOpSpin /= 1)) then
    write(u6,'(A,I4)') ' CI_Sigma_Op: invalid operator spin iOpSpin =',iOpSpin
    call Abend()
  end if
  if ((iOpSpin == 1) .and. Have_2_el) then
    write(u6,'(A)') ' CI_Sigma_Op: the operator with iOpSpin = 1 is one-body only'
    call Abend()
  end if
  call MCLR_Lucia_Sigma(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
end subroutine CI_Sigma_Op

subroutine CI_Trans_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  call Check_Lib('CI_Trans_Dens',associated(MCLR_Lucia_TDens))
  if ((iOpSpin /= 0) .and. (iOpSpin /= 1)) then
    write(u6,'(A,I4)') ' CI_Trans_Dens: invalid operator spin iOpSpin =',iOpSpin
    call Abend()
  end if
  call MCLR_Lucia_TDens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
end subroutine CI_Trans_Dens

subroutine CI_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nInt1, nInt2, n
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  real(kind=wp), intent(in) :: ECore
  real(kind=wp), intent(out) :: Diag(n)
  call Check_Lib('CI_Diag',associated(MCLR_Lucia_Diag))
  if ((Mode /= CI_Diag_Cnf_Avg) .and. (Mode /= CI_Diag_Exact)) then
    write(u6,'(A,I4)') ' CI_Diag: invalid Mode =',Mode
    call Abend()
  end if
  call MCLR_Lucia_Diag(iState,Mode,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
end subroutine CI_Diag

subroutine CI_Lib_File_To_Basis(iState,iSym,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, n
  real(kind=wp), intent(inout) :: CI(n)
  call Check_Lib('CI_Lib_File_To_Basis',associated(MCLR_Lucia_Conv))
  call MCLR_Lucia_Conv(iState,iSym,CI_To_Basis,n,CI)
end subroutine CI_Lib_File_To_Basis

subroutine CI_Lib_Basis_To_File(iState,iSym,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, n
  real(kind=wp), intent(inout) :: CI(n)
  call Check_Lib('CI_Lib_Basis_To_File',associated(MCLR_Lucia_Conv))
  call MCLR_Lucia_Conv(iState,iSym,CI_To_File,n,CI)
end subroutine CI_Lib_Basis_To_File

!> Only the MCLR copy of Lucia is connected so far
subroutine Check_Lib(Caller,Registered)
  character(len=*), intent(in) :: Caller
  logical(kind=iwp), intent(in) :: Registered
  if (CI_Lib_Op /= CI_Lib_MCLR_Lucia) then
    write(u6,'(1X,A,A,I2)') Caller,': no CI library selected for MCLR, or not connected yet:',CI_Lib_Op
    call Abend()
  end if
  if (.not. Registered) then
    write(u6,'(1X,A,A)') Caller,': the routines of the MCLR copy of Lucia are not registered'
    call Abend()
  end if
end subroutine Check_Lib

end module CI_Operations

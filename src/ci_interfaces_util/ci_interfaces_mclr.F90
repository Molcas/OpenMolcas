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

!> The CI gateway of MCLR: the CI operations of the linear response
!> equations and the choice of the CI library that performs them.
!>
!> Operations (all for the CI state iState, CI_St(iState)):
!>  - CI_MCLR_Sigma: sigma = H C for one CI vector. H is a Hermitian
!>    operator of symmetry Mul(iCSym,iSSym) (C of symmetry iCSym, sigma of
!>    symmetry iSSym); iOpSpin (the spin of the operator, not the
!>    multiplicity of the state) = 0: singlet E_pq = a+_pa a_qa + a+_pb a_qb,
!>    iOpSpin = 1: a+_pa a_qa - a+_pb a_qb (= sqrt(2) T^{1,0}_pq, one-body
!>    only). S is always overwritten (no accumulation).
!>    Int1: one-body integrals over all orbitals, in MCLR's square
!>    symmetry blocks (block (iS,jS), jS = Mul(iS,iOp), nOrb(iS) x nOrb(jS),
!>    column-major, one after the other in the order of iS); only the
!>    active part is used. Int2: two-body integrals over the active
!>    orbitals, symmetry-ordered active indices (all active orbitals of
!>    irrep 1, then irrep 2, ...; RAS1, RAS2, RAS3 inside each irrep),
!>    (ij|kl) at iTri(iTri(i,j),iTri(k,l)). Have_2_el = .false.: one-body
!>    operator only. Square = .true.: the operator is not symmetric (time-
!>    dependent response): Int1 is a general (non-symmetric) matrix and
!>    (ij|kl) is at iTri(ij,kl) with the square pair index ij = (j-1)*nAct+i.
!>  - CI_MCLR_TDens: transition densities between the bra L (symmetry
!>    iLSym) and the ket R (symmetry iRSym). rho1(i,j) = <L|E_ij|R>, full
!>    square nAct x nAct; for I12 = 2 also rho2(iTri(ij,kl)) =
!>    <L|E_ij E_kl - delta_jk E_il|R> with the square pair index
!>    ij = (j-1)*nAct+i (packed over ij >= kl); symmetry-ordered active
!>    indices. Spin-free normalisation: sum_i rho1(i,i) = N_act*<L|R>.
!>    Rho2 is not referenced for I12 = 1.
!>  - CI_MCLR_Diag: the diagonal <I|H|I> + ECore (ECore added to every
!>    element) in the CSF basis of symmetry iSym, for the symmetric
!>    active Hamiltonian; Int1 and Int2 as for the sigma with iOp = 1 and
!>    Square = .false.. (The MCLR copy of Lucia reads them from MCLR, see
!>    MCLR_Lucia_Diag.)
!>  - CI_MCLR_File_To_Basis, CI_MCLR_Basis_To_File: a CI vector of
!>    symmetry iSym between the stored basis (JobIph) and the CSF basis of
!>    the selected library, in place.
!> All vectors passed to the gateway are in the CSF basis of the selected
!> library, with their natural length (nCSF of their symmetry); each
!> library keeps its own internal basis, ordering and work space.
!>
!> The MCLR copy of Lucia keeps its data in MCLR (MCLR_Data), which the
!> gateway must not use. MCLR therefore registers its routines with
!> CI_MCLR_Register when it sets up the CI state; the gateway calls them
!> through procedure pointers. The other libraries are connected behind
!> the same interfaces.
module CI_Interfaces_MCLR

use Definitions, only: wp, iwp, u6

implicit none
private

!> The CI libraries (Lucia, SGUGA and Faroald with the values of RASSCF's
!> CI_Basis(); always compare against these names)
integer(kind=iwp), parameter :: CI_Lib_None = 0, CI_Lib_Lucia = 1, CI_Lib_SGUGA = 2, CI_Lib_Faroald = 3, CI_Lib_MCLR_Lucia = 4

!> Directions of the basis conversion (Conv_Proc)
integer(kind=iwp), parameter :: CI_To_File = 0, CI_To_Basis = 1

!> The selected library (changed only by CI_MCLR_Select)
integer(kind=iwp), protected :: CI_Lib_MCLR = CI_Lib_None

abstract interface
  subroutine Sigma_Proc(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
    import :: wp, iwp
    integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nInt1, nInt2, nC, nS
    real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
    logical(kind=iwp), intent(in) :: Have_2_el, Square
    real(kind=wp), intent(in) :: C(nC)
    real(kind=wp), intent(out) :: S(nS)
  end subroutine Sigma_Proc
  subroutine TDens_Proc(iState,I12,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
    import :: wp, iwp
    integer(kind=iwp), intent(in) :: iState, I12, iLSym, iRSym, nL, nR, n1, n2
    real(kind=wp), intent(in) :: L(nL), R(nR)
    real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  end subroutine TDens_Proc
  subroutine Diag_Proc(iState,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
    import :: wp, iwp
    integer(kind=iwp), intent(in) :: iState, iSym, nInt1, nInt2, n
    real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
    real(kind=wp), intent(in) :: ECore
    real(kind=wp), intent(out) :: Diag(n)
  end subroutine Diag_Proc
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

public :: CI_Lib_Faroald, CI_Lib_Lucia, CI_Lib_MCLR, CI_Lib_MCLR_Lucia, CI_Lib_None, CI_Lib_SGUGA, CI_MCLR_Basis_To_File, &
          CI_MCLR_Diag, CI_MCLR_File_To_Basis, CI_MCLR_Free, CI_MCLR_Register, CI_MCLR_Select, CI_MCLR_Sigma, CI_MCLR_TDens, &
          CI_To_Basis, CI_To_File, Conv_Proc, Diag_Proc, Sigma_Proc, TDens_Proc

contains

!> Select the CI library for MCLR (CI_Lib_*)
subroutine CI_MCLR_Select(iLib)
  integer(kind=iwp), intent(in) :: iLib
  select case (iLib)
    case (CI_Lib_MCLR_Lucia)
      CI_Lib_MCLR = iLib
    case (CI_Lib_Lucia,CI_Lib_SGUGA,CI_Lib_Faroald)
      write(u6,'(A,I2,A)') ' CI_MCLR_Select: CI library',iLib,' is not available for MCLR yet'
      call Abend()
    case default
      write(u6,'(A,I4)') ' CI_MCLR_Select: unknown CI library',iLib
      call Abend()
  end select
end subroutine CI_MCLR_Select

!> MCLR registers the routines of its copy of Lucia
subroutine CI_MCLR_Register(Sigma,TDens,Diag,Conv)
  procedure(Sigma_Proc) :: Sigma
  procedure(TDens_Proc) :: TDens
  procedure(Diag_Proc) :: Diag
  procedure(Conv_Proc) :: Conv
  MCLR_Lucia_Sigma => Sigma
  MCLR_Lucia_TDens => TDens
  MCLR_Lucia_Diag => Diag
  MCLR_Lucia_Conv => Conv
end subroutine CI_MCLR_Register

!> Forget the registered routines and the selection (end of MCLR)
subroutine CI_MCLR_Free()
  MCLR_Lucia_Sigma => null()
  MCLR_Lucia_TDens => null()
  MCLR_Lucia_Diag => null()
  MCLR_Lucia_Conv => null()
  CI_Lib_MCLR = CI_Lib_None
end subroutine CI_MCLR_Free

subroutine CI_MCLR_Sigma(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nInt1, nInt2, nC, nS
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  logical(kind=iwp), intent(in) :: Have_2_el, Square
  real(kind=wp), intent(in) :: C(nC)
  real(kind=wp), intent(out) :: S(nS)
  call Check_Lib('CI_MCLR_Sigma',associated(MCLR_Lucia_Sigma))
  if ((iOpSpin /= 0) .and. (iOpSpin /= 1)) then
    write(u6,'(A,I4)') ' CI_MCLR_Sigma: invalid operator spin iOpSpin =',iOpSpin
    call Abend()
  end if
  if ((iOpSpin == 1) .and. Have_2_el) then
    write(u6,'(A)') ' CI_MCLR_Sigma: the operator with iOpSpin = 1 is one-body only'
    call Abend()
  end if
  call MCLR_Lucia_Sigma(iState,iOpSpin,iCSym,iSSym,nInt1,Int1,nInt2,Int2,Have_2_el,Square,nC,C,nS,S)
end subroutine CI_MCLR_Sigma

subroutine CI_MCLR_TDens(iState,I12,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
  integer(kind=iwp), intent(in) :: iState, I12, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  call Check_Lib('CI_MCLR_TDens',associated(MCLR_Lucia_TDens))
  call MCLR_Lucia_TDens(iState,I12,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2)
end subroutine CI_MCLR_TDens

subroutine CI_MCLR_Diag(iState,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
  integer(kind=iwp), intent(in) :: iState, iSym, nInt1, nInt2, n
  real(kind=wp), target, intent(in) :: Int1(nInt1), Int2(nInt2)
  real(kind=wp), intent(in) :: ECore
  real(kind=wp), intent(out) :: Diag(n)
  call Check_Lib('CI_MCLR_Diag',associated(MCLR_Lucia_Diag))
  call MCLR_Lucia_Diag(iState,iSym,nInt1,Int1,nInt2,Int2,ECore,n,Diag)
end subroutine CI_MCLR_Diag

subroutine CI_MCLR_File_To_Basis(iState,iSym,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, n
  real(kind=wp), intent(inout) :: CI(n)
  call Check_Lib('CI_MCLR_File_To_Basis',associated(MCLR_Lucia_Conv))
  call MCLR_Lucia_Conv(iState,iSym,CI_To_Basis,n,CI)
end subroutine CI_MCLR_File_To_Basis

subroutine CI_MCLR_Basis_To_File(iState,iSym,n,CI)
  integer(kind=iwp), intent(in) :: iState, iSym, n
  real(kind=wp), intent(inout) :: CI(n)
  call Check_Lib('CI_MCLR_Basis_To_File',associated(MCLR_Lucia_Conv))
  call MCLR_Lucia_Conv(iState,iSym,CI_To_File,n,CI)
end subroutine CI_MCLR_Basis_To_File

!> Only the MCLR copy of Lucia is connected so far
subroutine Check_Lib(Caller,Registered)
  character(len=*), intent(in) :: Caller
  logical(kind=iwp), intent(in) :: Registered
  if (CI_Lib_MCLR /= CI_Lib_MCLR_Lucia) then
    write(u6,'(1X,A,A,I2)') Caller,': no CI library selected for MCLR, or not connected yet:',CI_Lib_MCLR
    call Abend()
  end if
  if (.not. Registered) then
    write(u6,'(1X,A,A)') Caller,': the routines of the MCLR copy of Lucia are not registered'
    call Abend()
  end if
end subroutine Check_Lib

end module CI_Interfaces_MCLR

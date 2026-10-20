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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************

!***********************************************************************
! History:                                                             *
! Jie J. Bao on May 09, 2022, created this file                        *
! Matthew R. Hennefarth Apr 05, 2023, upgraded to module               *
!***********************************************************************

module mspdft_util

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
private

public :: print_effective_ham, print_final_energies, print_mspdft_vectors, replace_diag

contains

! Prints the Final MS-PDFT Energies
!
! Args:
!   e_mspdft: ndarray of len nroots
!     Final MS-PDFT energies
!
!   nroots: integer
!     Number of roots in the calculation
subroutine print_final_energies(e_mspdft,nroots)

  integer(kind=iwp), intent(in) :: nroots
  real(kind=wp), intent(in) :: e_mspdft(nroots)
  integer(kind=iwp) :: state

  do state=1,nroots
    call printresult(u6,'(6X,A,I3,A,F16.8)','MSPDFT root number',state,' Total energy:',e_mspdft(state),1)
  end do
  write(u6,*)

end subroutine print_final_energies

! Prints the final mspdft eigenvectors in the intermediate state
! basis and reference state basis
!
! Args:
!   si_pdft: ndarray of length nroots*nroots
!     Contains the orthonormal eigenvectors in the intermediate
!     state basis.
!
!   nroots: integer
!     number of roots (lroots) or dimension of u and eigenvectors
subroutine print_mspdft_vectors(si_pdft,nroots)

  integer(kind=iwp), intent(in) :: nroots
  real(kind=wp), intent(in) :: si_pdft(nroots,nroots)
  integer(kind=iwp) :: root
  real(kind=wp) :: eig_vecs_in_ref(nroots,nroots), reference_vectors(nroots,nroots)
  logical(kind=iwp) :: refbas
  character(len=30) :: mspdftfmt
  character(Len=18) :: MatInfo
  character(len=9) :: StatVec, VecStat(nroots)

  refbas = .false.

  do root=1,nroots
    write(statvec,'(A5,I4)') 'Root ',root
    vecstat(root) = statvec
  end do

  write(u6,*)

  write(u6,'(6X,A)') 'Intermediate-state Basis'
  write(mspdftfmt,'(A4,I5,A9)') '(6X,',nroots,'(A10,5X))'
  write(u6,mspdftfmt) ((VecStat(root)),root=1,nroots)
  call RecPrt(' ','(7X,10(F9.6,6X))',si_pdft,size(si_pdft,1),size(si_pdft,2))
  write(u6,*)

  call f_inquire('ROT_VEC',refbas)
  if (refbas) then
    ! Generate reference state basis
    eig_vecs_in_ref = Zero
    call readmat2('ROT_VEC',MatInfo,reference_vectors,nroots,nroots,7,18,'T')
    call dgemm_('n','n',nroots,nroots,nroots,One,reference_vectors,nroots,si_pdft,nroots,Zero,eig_vecs_in_ref,nroots)
    write(u6,'(6X,A)') 'Reference-state Basis'
    write(u6,mspdftfmt) ((VecStat(root)),root=1,nroots)
    call RecPrt(' ','(7X,10(F9.6,6X))',eig_vecs_in_ref,nroots,nroots)
    call PrintMat2('FIN_VEC',MatInfo,eig_vecs_in_ref,nroots,nroots,7,18,'T')
  end if

end subroutine print_mspdft_vectors

! Prints the effective Hamiltonian
!
! Args:
!   heff: ndarray of len nroots*nroots
!     Effective MS-PDFT Hamiltonian
!
!   nroots: integer
!     Dimension of heff, or number of roots
!
!   digit: integer
!     Threshold value to shift diagonal elements by when printed
subroutine print_effective_ham(heff,nroots,digit)

  use PrintLevel, only: SILENT
  use mcpdft_output, only: iPrLoc

  integer(kind=iwp), intent(in) :: nroots, digit
  real(kind=wp), intent(in) :: heff(nroots,nroots)
  integer(kind=iwp) :: root
  real(kind=wp) :: shift, shifted_heff(nroots,nroots)

  if (iPrLoc(1) == SILENT) return

  shifted_heff = heff

  call should_shift_diag(heff,nroots,digit,shift)

  if (shift /= Zero) then
    write(u6,'(6X,A,F9.2,A)') '(diagonal values increased by',-shift,' hartree)'
    do root=1,nroots
      shifted_heff(root,root) = shifted_heff(root,root)-shift
    end do
  end if

  call recprt(' ','(7X,10(F9.6,1X))',shifted_heff,nroots,nroots)
  write(u6,*)

end subroutine print_effective_ham

! Determines if diagonal elements should be shifted, and if so the number
! of decimal places to shift.
!
! Args:
!   heff: ndarray of len nroots*nroots
!     Effective MS-PDFT Hamiltonian
!
!   nroots: integer
!     Dimension of heff, or number of roots
!
!   digit: integer
!     Threshold value to shift diagonal elements by when printed
!
! Returns:
!   shift: real
!     Amount to shift diagonal elements by
subroutine should_shift_diag(heff,nroots,digit,shift)

  integer(kind=iwp), intent(in) :: nroots, digit
  real(kind=wp), intent(in) :: heff(nroots,nroots)
  real(kind=wp), intent(out) :: shift
  integer(kind=iwp) :: i, ishift
  real(kind=wp) :: maxelem, rdiag(nroots)

  shift = Zero

  do i=1,nroots
    rdiag(i) = heff(i,i)
  end do

  maxelem = maxval(rdiag)

  if (abs(maxelem) < real(digit,8)) return

  ishift = int(maxelem,8)/digit*digit
  shift = real(ishift,8)

end subroutine should_shift_diag

subroutine replace_diag(mat,diag,ndim)

  integer(kind=iwp), intent(in) :: ndim
  real(kind=wp), intent(inout) :: mat(ndim,ndim)
  real(kind=wp), intent(in) :: diag(ndim)
  integer(kind=iwp) :: i

  do i=1,ndim
    mat(i,i) = diag(i)
  end do

end subroutine replace_diag

end module mspdft_util

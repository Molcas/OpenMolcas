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
  use constants, only: zero
  use definitions, only: wp
  use mcpdft_output, only: lf
  implicit none
  private

  public :: print_final_energies, print_mspdft_vectors, print_effective_ham
  public :: replace_diag

  contains
  subroutine print_final_energies(e_mspdft, nroots, method)
    ! Prints the Final MS-PDFT Energies
    !
    ! Args:
    !   e_mspdft: ndarray of len nroots
    !     Final MS-PDFT energies
    !
    !   nroots: integer
    !     Number of roots in the calculation
    !
    !   method: character(len=8)
    !     MS-PDFT method string

    use hybridpdft, only: do_hybrid

    integer, intent(in) :: nroots
    real(kind=wp), dimension(nroots), intent(in) :: e_mspdft
    character(len=8), intent(in) :: method

    integer :: root

    if(.not.do_hybrid) then
      write(lf,'(6X,2A)') method, ' Energies:'
      do root=1, nroots
        write(lf, '(6X,3A,1X,I4,3X,A13,F18.8)') '::    ', method, &
                ' Root', root, 'Total energy:', e_mspdft(root)
      end do
    else
      write(lf ,'(6X,3A)') 'Hybrid ', method, ' Energies:'
      do root=1, nroots
          write(lf,'(6X,4A,1X,I4,3X,A13,F18.8)') '::    ', 'Hybrid ', &
                  method, ' Root', root, 'Total energy:', e_mspdft(root)
      end do
    end if
  end subroutine print_final_energies

  subroutine print_mspdft_vectors(si_pdft, nroots)
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
    !

    integer, intent(in) :: nroots
    real(kind=wp), dimension(nroots**2), intent(in) :: si_pdft

    logical :: refbas = .false.
    character(len=9), dimension(nroots) :: VecStat
    character(len=9) :: StatVec
    character(len=30)::mspdftfmt
    character(Len=18) :: MatInfo

    integer :: root
    real(kind=wp), dimension(nroots**2) :: reference_vectors, eig_vecs_in_ref

    do root=1, nroots
      write(statvec, '(A5,I4)') 'Root ', root
      vecstat(root) = statvec
    end do

    write(lf, *)

    write(lf,'(7X,A)')'Intermediate-state Basis'
    write(mspdftfmt, '(A4,I5,A9)') '(6X,',nroots,'(A10,5X))'
    write(lf, mspdftfmt) ((VecStat(root)),root=1,nroots)
    Call RecPrt(' ','(7X,10(F9.6,6X))', si_pdft, nroots, nroots)

    call f_inquire('ROT_VEC', refbas)
    if (refbas) then
      ! Generate reference state basis
      call fzero(eig_vecs_in_ref, nroots**2)
      call readmat2('ROT_VEC', MatInfo, reference_vectors, nroots, nroots, 7, 18, 'T')
      call dgemm_('n', 'n', nroots, nroots, nroots, 1.0d0, reference_vectors,&
              nroots, si_pdft, nroots, 0.0d0, eig_vecs_in_ref, nroots)
      write(lf,'(7X,A)')'Reference-state Basis'
      write(lf,mspdftfmt)((VecStat(root)),root=1,nroots)
      call RecPrt(' ','(7X,10(F9.6,6X))', eig_vecs_in_ref, nroots, nroots)
      call PrintMat2('FIN_VEC',MatInfo,eig_vecs_in_ref, nroots, nroots,7,18,'T')
    end if

  end subroutine print_mspdft_vectors

  subroutine print_effective_ham(heff, nroots, digit)
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

    integer, intent(in) :: nroots, digit
    real(kind=wp), dimension(nroots**2), intent(in) :: heff

    real(kind=wp), dimension(nroots**2) :: shifted_heff
    real(kind=wp) :: shift
    integer :: root

    shifted_heff = heff

    call should_shift_diag(heff, nroots, digit, shift)

    if(shift /= zero) then
      write(lf,'(6X,A,F9.2,A)') '(diagonal values increased by', -shift, ' hartree)'
      do root=1, nroots
        shifted_heff((root-1)*nroots + root) = shifted_heff((root-1)*nroots + root) - shift
      end do
    end if

    call recprt(' ', '(7X,10(F9.6,1X))', shifted_heff, nroots, nroots)
    write(lf, *)
  end subroutine print_effective_ham

  subroutine should_shift_diag(heff, nroots, digit, shift)
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

    integer, intent(in) :: nroots, digit
    real(kind=wp), dimension(nroots**2), intent(in) :: heff

    real(kind=wp), intent(out) :: shift

    real(kind=wp), dimension(nroots) :: rdiag
    real(kind=wp) :: maxelem
    integer :: i, ishift

    shift = zero

    do i=1, nroots
      rdiag(i) = heff((i-1)*nroots + i)
    end do

    maxelem = maxval(rdiag)

    if (abs(maxelem) < real(digit, 8)) then
      return
    end if

    ishift = int(maxelem, 8)/digit * digit
    shift = real(ishift, 8)

  end subroutine should_shift_diag

  subroutine replace_diag(mat, diag, ndim)
    integer, intent(in) :: ndim
    real(kind=wp), dimension(ndim), intent(in) :: diag
    real(kind=wp), dimension(ndim**2), intent(inout) :: mat

    integer :: i

    do i=1, ndim
      mat(ndim*(i-1)+i) = diag(i)
    end do

  end subroutine replace_diag

end module mspdft_util

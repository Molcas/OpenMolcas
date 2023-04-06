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

  logical :: lshiftdiag = .false.
  real(kind=wp) :: mspdftshift = zero

  public :: print_final_energies, print_mspdft_vectors, print_effective_ham

  contains
  subroutine print_final_energies(energies, ndim)
    ! Prints the Final MS-PDFT Energies
    !
    ! Args:
    !   energies: ndarray of len ndim
    !     Final MS-PDFT energies
    !
    !   ndim: integer
    !     Number of roots in the calculation

    use hybridpdft, only: do_hybrid
    use mspdft, only: mspdftmethod

    integer, intent(in) :: ndim
    real(kind=wp), dimension(ndim), intent(in) :: energies

    integer :: root

    if(.not.do_hybrid) then
      write(lf,'(6X,2A)') MSPDFTMethod, ' Energies:'
      do root=1, ndim
        write(lf, '(6X,3A,1X,I4,3X,A13,F18.8)') '::    ', MSPDFTMethod, &
                ' Root', root, 'Total energy:', energies(root)
      end do
    else
      write(lf ,'(6X,3A)') 'Hybrid ', MSPDFTMethod, ' Energies:'
      do root=1, ndim
          write(lf,'(6X,4A,1X,I4,3X,A13,F18.8)') '::    ', 'Hybrid ', &
                  MSPDFTMethod, ' Root', root, 'Total energy:', energies(root)
      end do
    end if
  end subroutine print_final_energies

  subroutine print_mspdft_vectors(u, ndim)
    ! Prints the final mspdft eigenvectors in the intermediate state
    ! basis and reference state basis
    !
    ! Args:
    !   u: ndarray of length ndim*ndim
    !     Contains the orthonormal eigenvectors in the intermediate
    !     state basis.
    !
    !   ndim: integer
    !     number of roots (lroots) or dimension of u and eigenvectors
    !
    !   matinfo: character(len=18)
    !     info regarding itermediate state basis, or what type of
    !     MS-PDFT method we are doing.

    use hybridpdft, only: do_hybrid
    use mspdft, only: mspdftmethod

    integer, intent(in) :: ndim
    real(kind=wp), dimension(ndim**2), intent(in) :: u

    logical :: refbas = .false.
    character(len=9), dimension(ndim) :: VecStat
    character(len=9) :: StatVec
    character(len=30)::mspdftfmt
    character(Len=18) :: MatInfo

    integer :: root
    real(kind=wp), dimension(ndim**2) :: reference_vectors, eig_vecs_in_ref

    do root=1, ndim
      write(statvec, '(A5,I4)') 'Root ', root
      vecstat(root) = statvec
    end do

    write(lf, *)
    if(do_hybrid) then
      write(lf, '(6X,3A)') 'Hybrid ',MSPDFTMethod,' Eigenvectors:'
    else
        write(lf,'(6X,2A)') MSPDFTMethod,' Eigenvectors:'
    end if

    write(lf,'(7X,A)')'Intermediate-state Basis'
    write(mspdftfmt, '(A4,I5,A9)') '(6X,',ndim,'(A10,5X))'
    write(lf, mspdftfmt) ((VecStat(root)),root=1,ndim)
    Call RecPrt(' ','(7X,10(F9.6,6X))', u, ndim, ndim)

    call f_inquire('ROT_VEC', refbas)
    if (refbas) then
      ! Generate reference state basis
      call fzero(eig_vecs_in_ref, ndim**2)
      call readmat2('ROT_VEC', MatInfo, reference_vectors, ndim, ndim, 7, 18, 'T')
      call dgemm_('n', 'n', ndim, ndim, ndim, 1.0d0, reference_vectors,&
              ndim, u, ndim, 0.0d0, eig_vecs_in_ref, ndim)
      write(lf,'(7X,A)')'Reference-state Basis'
      write(lf,mspdftfmt)((VecStat(root)),root=1,ndim)
      call RecPrt(' ','(7X,10(F9.6,6X))', eig_vecs_in_ref, ndim, ndim)
      call PrintMat2('FIN_VEC',MatInfo,eig_vecs_in_ref, ndim, ndim,7,18,'T')
    end if

  end subroutine print_mspdft_vectors

  subroutine print_effective_ham(mat, ndim, digit)
    ! Prints the effective Hamiltonian
    !
    ! Args:
    !   mat: ndarray of len ndim*ndim
    !     Effective MS-PDFT Hamiltonian
    !
    !   ndim: integer
    !     Dimension of mat, or number of roots
    !
    !   digit: integer
    !     Threshold value to shift diagonal elements by when printed

    use hybridpdft, only: do_hybrid
    use mspdft, only: mspdftmethod

    integer, intent(in) :: ndim, digit
    real(kind=wp), dimension(ndim**2), intent(in) :: mat

    real(kind=wp), dimension(ndim**2) :: shifted_mat
    integer :: root

    shifted_mat = mat

    call should_shift_diag(mat, ndim, digit)

    if(.not. do_hybrid) then
      write(lf, '(6X,2A)') mspdftmethod, ' Effective Hamiltonian'
    else
      write(lf,'(6X,3A)') 'Hybrid ', MSPDFTMethod, ' Effective Hamiltonian'
    end if

    if(lshiftdiag) then
      write(lf,'(6X,A,F9.2,A)') '(diagonal values increased by', -MSPDFTShift, ' hartree)'
      do root=1, ndim
        shifted_mat((root-1)*ndim + root) = shifted_mat((root-1)*ndim + root) - mspdftshift
      end do
    end if

    call recprt(' ', '(7X,10(F9.6,1X))', shifted_mat, ndim, ndim)
    write(lf, *)
  end subroutine print_effective_ham

  subroutine should_shift_diag(mat, ndim, digit)
    ! Determines if diagonal elements should be shifted, and if so the number
    ! of decimal places to shift.
    !
    ! Args:
    !   mat: ndarray of len ndim*ndim
    !     Effective MS-PDFT Hamiltonian
    !
    !   ndim: integer
    !     Dimension of mat, or number of roots
    !
    !   digit: integer
    !     Threshold value to shift diagonal elements by when printed

    integer, intent(in) :: ndim, digit
    real(kind=wp), dimension(ndim**2), intent(in) :: mat

    real(kind=wp), dimension(ndim) :: rdiag
    real(kind=wp) :: maxelem
    integer :: i, ishift

    do i=1, ndim
      rdiag(i) = mat((i-1)*ndim + i)
    end do

    maxelem = maxval(rdiag)

    lshiftdiag = (abs(maxelem) > real(digit, 8))

    if (.not. lshiftdiag) then
      return
    end if

    ishift = int(maxelem, 8)/digit * digit
    mspdftshift = real(ishift, 8)

  end subroutine should_shift_diag

end module mspdft_util

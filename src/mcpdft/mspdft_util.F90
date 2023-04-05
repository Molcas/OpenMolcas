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
  implicit none
  private

  logical :: lshiftdiag = .false.
  real(kind=wp) :: mspdftshift = zero

  public :: print_final_energies, print_effective_ham

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
    use mcpdft_output, only: lf

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
    use mcpdft_output, only: lf

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

    use definitions, only: wp

    integer, intent(in) :: ndim, digit
    real(kind=wp), dimension(ndim**2), intent(in) :: mat

    real(kind=wp), dimension(ndim) :: rdiag
    real(kind=wp) :: maxelem
    integer :: i, ishift

    do i=1, ndim
      rdiag(i) = mat((i-1)*ndim + i)
    end do

    maxelem = maxval(rdiag)

    lshiftdiag = (abs(maxelem) < real(digit, 8))

    if (.not. lshiftdiag) then
      return
    end if

    ishift = int(maxelem, 8)/digit * digit
    mspdftshift = real(ishift, 8)

  end subroutine should_shift_diag

end module mspdft_util

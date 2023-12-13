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
subroutine print_drt(nvert, drt, down)
  ! Prints the DRT Table
  use mcpdft_output, only: lf

  implicit none

  integer, intent(in) :: nvert
  integer, dimension(nvert, 5), intent(in) :: drt
  integer, dimension(nvert, 0:3), intent(in) :: down

  integer :: i, v, s ! Dummy variables

  write(lf, *)
  write(lf, *) ' VERT      L  N    A  B  C      CHAINING INDICES.'

  do v=1, nvert
    write(lf, '(1X,I4,5X,2I3,2X,3I3,5X,4I4)') V,(DRT(V,I),I=1,5),(DOWN(V,S),S=0,3)
  end do
  write(lf, *)
end subroutine print_drt
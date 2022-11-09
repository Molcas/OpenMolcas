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

subroutine VEIG(N,A,V)
! THE DIAGONAL ELEMENTS OF THE LOWER TRIANGULAR MATRIX, A, STORED
! IN UPPER PACKED STORAGE MODE ARE COPIED TO THE FIELD V.

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(nTri_Elem(N))
real(kind=wp), intent(out) :: V(N)
integer(kind=iwp) :: I

do I=1,N
  V(I) = A(nTri_Elem(I))
end do

return

end subroutine VEIG

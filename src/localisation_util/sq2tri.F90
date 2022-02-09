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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Sq2Tri
!
!> @brief
!>   Convert from square to lower triangular storage
!> @author Thomas Bondo Pedersen
!>
!> @param[in]  Sq  Square storage array
!> @param[out] Tri Lower triangular storage array
!> @param[in]  n   Dimension
!>
!> @details
!> Perform the extraction
!>
!> \code
!> Tri(i*(i-1)/2+j) = Sq(i,j)
!> \endcode
!>
!> where \c i &ge; \c j.
!***********************************************************************

subroutine Sq2Tri(Sq,Tri,n)

use Index_Functions, only: iTri, nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: Sq(n,n)
real(kind=wp), intent(out) :: Tri(nTri_Elem(n))
integer(kind=iwp) :: i, j

do j=1,n
  do i=j,n
    Tri(iTri(i,j)) = Sq(i,j)
  end do
end do

end subroutine Sq2Tri

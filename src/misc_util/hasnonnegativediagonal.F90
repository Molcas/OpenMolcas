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
!  hasNonnegativeDiagonal
!
!> @brief
!>   Return ``.True.`` if all diagonal elements of \p M are non-negative
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Test if all diagonal elements of matrix \p M are non-negative.
!> If \f$ \forall i \; M_{ii} \ge 0 \f$ return ``.True.``.
!>
!> @param[in] M \p n &times; \p n square matrix
!> @param[in] n Dimension of \p M
!>
!> @return ``.True.`` if all diagonal elements of \p M are non-negative.
!***********************************************************************

function hasNonnegativeDiagonal(M,n)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: hasNonnegativeDiagonal
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: M(n,n)
integer(kind=iwp) :: i

hasNonnegativeDiagonal = .true.
do i=1,n
  if (M(i,i) < Zero) then
    hasNonnegativeDiagonal = .false.
    exit
  end if
end do

return

end function hasNonnegativeDiagonal

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
!  isSymmetric
!
!> @brief
!>   Return ``.True.`` if \p M is symmetric within tolerance \p Tol
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Test if matrix \p M is symmetric. If \f$ \forall i,j \; |M_{ij}-M_{ji}| \le \text{Tol} \f$
!> return ``.True.``
!>
!> @param[in] M   \p n &times; \p n square matrix to test
!> @param[in] n   Dimension of \p M
!> @param[in] Tol Tolerance
!>
!> @return ``.True.`` if \p M is symmetric within tolerance \p Tol
!***********************************************************************

function isSymmetric(M,n,Tol)

use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: isSymmetric
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: M(n,n), Tol
integer(kind=iwp) :: i, j

isSymmetric = .true.
outer: do j=1,n
  do i=j+1,n
    if (abs(M(i,j)-M(j,i)) > Tol) then
      isSymmetric = .false.
      exit outer
    end if
  end do
end do outer

return

end function isSymmetric

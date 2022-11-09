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
!  obeysCauchySchwarz
!
!> @brief
!>   Return ``.True.`` if \p M obeys the Cauchy--Schwarz inequality within tolerance \p Tol
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Test if matrix \p M obeys the Cauchy--Schwarz inequality.
!> If \f$ \exists i \gt j \; |M_{ij}^2 - M_{ii}M_{jj}| \gt \text{Tol} \f$
!> return ``.False.``.
!>
!> @param[in] M   \p n &times; \p n square matrix to test
!> @param[in] n   Dimension of \p M
!> @param[in] Tol Tolerance
!>
!> @return ``.True.`` if \p M obeys the Cauchy--Schwarz inequality within tolerance \p Tol
!***********************************************************************

function obeysCauchySchwarz(M,n,Tol)

use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: obeysCauchySchwarz
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: M(n,n), Tol
integer(kind=iwp) :: i, j

obeysCauchySchwarz = .true.
outer: do j=1,n
  do i=j+1,n
    if ((M(i,j)**2 > M(i,i)*M(j,j)) .and. (abs(M(i,j)**2-M(i,i)*M(j,j)) > Tol)) then
      obeysCauchySchwarz = .false.
      exit
    end if
  end do
end do outer

return

end function obeysCauchySchwarz

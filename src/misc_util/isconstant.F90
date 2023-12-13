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
!  isConstant
!
!> @brief
!>   Return ``.True.`` if all elements of \p M are identical to \p Const within tolerance \p Tol
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Test if array \p M is constant. If \f$ \forall i \; |M_i - \text{Const}| \le \text{Tol} \f$
!> return ``.True.``.
!>
!> @param[in] M     Array to test
!> @param[in] n     Dimension of \p M
!> @param[in] Const Constant
!> @param[in] Tol   Tolerance
!>
!> @return ``.True.`` if all elements of \p M are identical to \p Const within tolerance \p Tol
!***********************************************************************

function isConstant(M,n,Const,Tol)

use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: isConstant
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: M(n), Const, Tol
integer(kind=iwp) :: i

isConstant = .true.
i = 0
do while ((i < n) .and. isConstant)
  i = i+1
  isConstant = abs(M(i)-Const) <= Tol
end do

end function isConstant

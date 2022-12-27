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
! Copyright (C) 2014, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Sp_MV
!
!> @ingroup Sparse
!> @brief
!>   Compute a matrix-vector product \f$ y \leftarrow \alpha A x + \beta y \f$, with a sparse matrix
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Equivalent to ::DGeMV or ::DSyMV, with a sparse matrix \p A.
!>   \f[ y \leftarrow \alpha A x + \beta y \f]
!>
!> @param[in]     n     Size of the system
!> @param[in]     alpha Factor for the multiplication
!> @param[in]     A     Input matrix in sparse format
!> @param[in]     ija   Index vector of matrix \p A
!> @param[in]     x     Vector to multiply
!> @param[in]     beta  Factor for the initial vector
!> @param[in,out] y     Result vector
!***********************************************************************

subroutine Sp_MV(n,alpha,A,ija,x,beta,y)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, ija(*)
real(kind=wp), intent(in) :: alpha, A(*), x(n), beta
real(kind=wp), intent(inout) :: y(n)
integer(kind=iwp) :: i, j, k

! Very simple routine, but split in different cases
! to gain efficiency
if (A(n+1) > Zero) then
  if ((beta == Zero) .and. (alpha == One)) then
    do i=1,n
      y(i) = A(i)*x(i)
      do k=ija(i),ija(i+1)-1
        j = ija(k)
        y(i) = y(i)+A(k)*x(j)
        y(j) = y(j)+A(k)*x(i)
      end do
    end do
  else
    do i=1,n
      y(i) = beta*y(i)+alpha*A(i)*x(i)
      do k=ija(i),ija(i+1)-1
        j = ija(k)
        y(i) = y(i)+alpha*A(k)*x(j)
        y(j) = y(j)+alpha*A(k)*x(i)
      end do
    end do
  end if
else
  if ((beta == Zero) .and. (alpha == One)) then
    do i=1,n
      y(i) = A(i)*x(i)
      do k=ija(i),ija(i+1)-1
        y(i) = y(i)+A(k)*x(ija(k))
      end do
    end do
  else
    do i=1,n
      y(i) = beta*y(i)+alpha*A(i)*x(i)
      do k=ija(i),ija(i+1)-1
        y(i) = y(i)+alpha*A(k)*x(ija(k))
      end do
    end do
  end if
end if

end subroutine Sp_MV

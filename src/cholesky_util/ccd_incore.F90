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
!  CCD_InCore
!
!> @brief
!>   Complete Cholesky decomposition of the symmetric positive definite matrix \p X
!> @author Thomas Bondo Pedersen
!>
!> @details
!> The \p n &times; \p n matrix \p X is Cholesky decomposed and the resulting
!> Cholesky vectors are returned in the lower triangle of \p X.
!> A non-zero return code signals
!> that an error has occured (\p X might e.g. be non-positive
!> definite) and the output is ill-defined.
!>
!> @param[in,out] X   Matrix to be Cholesky decomposed;
!>                    on exit, lower triangle contains the vectors
!> @param[in]     n   Linear dimension of \p X
!> @param[out]    irc Return code
!***********************************************************************

subroutine CCD_InCore(X,n,irc)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: X(n,n)
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: i, j, k
real(kind=wp) :: Fac

irc = 0
if (n < 1) return ! return (nothing to do)

do j=1,n
  ! Check for negative diagonal
  if (X(j,j) > Zero) then
    Fac = One/sqrt(X(j,j))
  else
    irc = 1
    return
  end if
  ! Compute vector j
  do i=j,n
    X(i,j) = Fac*X(i,j)
  end do
  ! Subtract from remaining columns
  do k=j+1,n
    Fac = X(k,j)
    do i=k,n
      X(i,k) = X(i,k)-X(i,j)*Fac
    end do
  end do
end do

end subroutine CCD_InCore

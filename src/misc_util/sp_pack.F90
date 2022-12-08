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
!  Sp_Pack
!
!> @ingroup Sparse
!> @brief
!>   Store a matrix in sparse format
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> An input matrix \p A will be stored in a sparse format, as two vectors \p Sp and \p ij_Sp.
!> Only elements with an absolute value larger than \p Thr will be kept.
!> An error is raised if the number of stored elements is larger than \p nmax.
!> The sizes of \p Sp and \p ij_Sp must be at least \f$ n + k + 1 \f$, where \f$ k \f$ is the number of
!> non_zero off-diagonal elements of \p A.
!>
!> @param[in]  n     Size of the matrix
!> @param[in]  A     Input matrix, in dense (conventional) format
!> @param[in]  nmax  Maximum size of the stored vectors
!> @param[out] Sp    Output matrix in sparse format
!> @param[out] ij_Sp Index vector for the sparse output matrix
!> @param[in]  Sym   Flag specifying whether to use the symmetric format
!> @param[in]  Thr   Threshold to discard small elements
!***********************************************************************

subroutine Sp_Pack(n,A,nmax,Sp,ij_Sp,Sym,Thr)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n, nmax
real(kind=wp), intent(in) :: A(n,n), Thr
real(kind=wp), intent(_OUT_) :: Sp(*)
integer(kind=iwp), intent(_OUT_) :: ij_Sp(*)
logical(kind=iwp), intent(in) :: Sym
integer(kind=iwp) :: i, j, nij

ij_Sp(1) = n+2
nij = n+1
if (Sym) then
  do i=1,n
    Sp(i) = A(i,i)
    do j=1,i-1
      if (abs(A(i,j)) <= Thr) then
        nij = nij+1
        if (nij > nmax) call SysAbendMsg('Sp_Pack','Too many non-zero elements.','')
        Sp(nij) = A(i,j)
        ij_Sp(nij) = j
      end if
    end do
    ij_Sp(i+1) = nij+1
  end do
  Sp(n+1) = One
else
  do i=1,n
    Sp(i) = A(i,i)
    do j=1,n
      if ((j /= i) .and. (abs(A(i,j)) <= Thr)) then
        nij = nij+1
        if (nij > nmax) call SysAbendMsg('Sp_Pack','Too many non-zero elements.','')
        Sp(nij) = A(i,j)
        ij_Sp(nij) = j
      end if
    end do
    ij_Sp(i+1) = nij+1
  end do
  Sp(n+1) = Zero
end if

end subroutine Sp_Pack

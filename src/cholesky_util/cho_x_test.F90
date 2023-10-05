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
!  Cho_X_Test
!
!> @brief
!>   Check Cholesky decomposition of matrix \p X
!> @author Thomas Bondo Pedersen
!>
!> @details
!> The difference between the matrix \p X and its Cholesky
!> representation is calculated and returned in array \p Y. If the
!> RMS error is less than \p Thr, \p irc = ``0`` is returned, else a positive
!> value is returned (a negative value signals an input error).
!> The input matrix \p X may be stored as a lower triangle
!> (\p Square = ``.false.``) or as a full square matrix (\p Square = ``.true.``).
!> The dimension of \p Y will be the same as that of \p X
!> (i.e. ``n(n+1)/2`` for lower triangular storage or ``n*n`` for full
!> square storage).
!>
!> @param[in]  X      The \p n &times; \p n matrix to be tested
!> @param[in]  n      Dimension of matrix \p X
!> @param[in]  Square Flag for packing of \p X
!> @param[in]  Vec    Cholesky vectors representing \p X
!> @param[in]  nVec   Number of Cholesky vectors
!> @param[in]  xf     Factor for the scaling of the vectors
!> @param[out] Y      ``Y = X - xf*Vec*VecT``
!> @param[in]  lY     Dimension of array \p Y
!> @param[in]  Thr    Threshold allowed for RMS error
!> @param[out] irc    Return code
!***********************************************************************

subroutine Cho_X_Test(X,n,Square,Vec,nVec,xf,Y,lY,Thr,irc)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp, wp

implicit none
integer(kind=iwp), intent(in) :: n, nVec, lY
real(kind=wp), intent(in) :: X(*), Vec(n,nVec), xf, Thr
real(kind=wp), intent(out) :: Y(lY)
logical(kind=iwp), intent(in) :: Square
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: lX
real(kind=wp) :: RMS
real(kind=wp), external :: ddot_

! Check input.
! ------------

if (n < 1) then ! nothing to do
  irc = 0
  return
else if ((nVec < 0) .or. (Thr < Zero)) then ! input error
  irc = -1
  return
end if
if (Square) then
  lX = n**2
else
  lX = nTri_Elem(n)
end if
if (lY < lX) then ! insufficient memory
  irc = -2
  return
end if

! Compute Y(ij) = X(ij) - xf * sum_J L(iJ) * L(jJ).
! -------------------------------------------------

Y(1:lX) = X(1:lX)
if (Square) then
  call DGEMM_('N','T',n,n,nVec,-xf,Vec,n,Vec,n,One,Y,n)
else
  call dGeMM_Tri('N','T',n,n,nVec,-xf,Vec,n,Vec,n,One,Y,n)
end if

! Check RMS error.
! ----------------

RMS = sqrt(dDot_(lX,Y,1,Y,1)/real(lX,kind=wp))
if (RMS > Thr) then
  irc = 1
else
  irc = 0
end if

end subroutine Cho_X_Test

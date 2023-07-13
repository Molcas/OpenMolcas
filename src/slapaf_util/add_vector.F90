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
!  Add_Vector
!
!> @brief
!>   Add (or not) a vector to an orthonormal set
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Adds a given vector to an existing set of orthonormal vectors.
!> The vector is orthogonalized against the existing set and, if the remainder
!> is large enough, it will be normalized and added to the set.
!> The input matrix with the initial set must have size at least \p n (\p m+1).
!> On output, \p m will be increased by 1 if the vector was added, and unchanged
!> otherwise.
!>
!> @param[in]     n   Size of the vectors
!> @param[in,out] m   Number of vectors in the subspace
!> @param[in,out] Sub Subspace of vectors
!> @param[in,out] Vec Vector to be added
!> @param[in]     Thr Threshold for linear dependence check
!***********************************************************************

subroutine Add_Vector(n,m,Sub,Vec,Thr)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(inout) :: m
real(kind=wp), intent(inout) :: Sub(n,m+1), Vec(n)
real(kind=wp), intent(in) :: Thr
integer(kind=iwp) :: i
real(kind=wp) :: Aux
real(kind=wp), external :: DDot_

do i=1,m
  Vec(:) = Vec(:)-DDot_(n,Sub(:,i),1,Vec,1)*Sub(:,i)
end do
Aux = DDot_(n,Vec,1,Vec,1)
if (abs(Aux) > Thr) then
  ! Safety net: orthonormalize again before adding it
  Vec(:) = Vec(:)/sqrt(Aux)
  do i=1,m
    Vec(:) = Vec(:)-DDot_(n,Sub(:,i),1,Vec,1)*Sub(:,i)
  end do
  Aux = DDot_(n,Vec,1,Vec,1)
  m = m+1
  Sub(:,m) = Vec(:)/sqrt(Aux)
end if

end subroutine Add_Vector

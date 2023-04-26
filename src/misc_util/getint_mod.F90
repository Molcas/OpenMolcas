!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module GetInt_mod

use Definitions, only: wp, iwp

implicit none
private

! Variables for computing integrals from Cholesky vectors.
integer(kind=iwp) :: I, ID_IP, LuCVec(2), mNeed, nBas(8), nPQ, nRS, NumCho(8), NumV, nVec, pq1
real(kind=wp), allocatable :: Vec2(:,:)
integer(kind=iwp), allocatable :: Basis_IDs(:,:), hash_table(:), lists(:,:)

public :: Basis_IDs, hash_table, I, ID_IP, lists, LuCVec, mNeed, nBas, nPQ, nRS, NumCho, NumV, nVec, pq1, Vec2

end module GetInt_mod

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

use Definitions, only: iwp

implicit none
private

! Variables for computing integrals from Cholesky vectors.
integer(kind=iwp) :: LuCVec(2), nBas(8), NumCho(8), pq1

public :: LuCVec, nBas, NumCho, pq1

end module GetInt_mod

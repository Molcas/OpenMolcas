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

! auxiliary constant pool: ready only up to g-valence/g-core

module AMatrix

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: lp1 = 5, lp12 = lp1**2, lp13 = lp1*(lp1+1)/2
integer(kind=iwp) :: KOSUU(lp13), NYU(lp1,lp13)
real(kind=wp) :: DFAC(lp12), RCA(lp1,lp13)

public :: DFAC, KOSUU, lp1, lp12, NYU, RCA

end module AMatrix

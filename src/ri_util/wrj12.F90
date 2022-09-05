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

module Wrj12

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: iOffA(4,0:7), Lu_A(0:7), Lu_Q(0:7), nChV(0:7)
integer(kind=iwp), allocatable :: SO2Ind(:)

public :: iOffA, Lu_A, Lu_Q, nChV, SO2Ind

end module Wrj12

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

module TempTime

use Definitions, only: wp

implicit none
private

real(kind=wp) :: ChoGet_CPU, ChoGet_Wall, Drvg1_CPU, Drvg1_Wall, Pget2_CPU, Pget2_Wall, Pget3_CPU, Pget3_Wall, Prepp_CPU, &
                 Prepp_Wall, rMult_CPU, rMult_Wall, Twoel2_CPU, Twoel2_Wall, Twoel3_CPU, Twoel3_Wall

public :: ChoGet_CPU, ChoGet_Wall, Drvg1_CPU, Drvg1_Wall, Pget2_CPU, Pget2_Wall, Pget3_CPU, Pget3_Wall, Prepp_CPU, Prepp_Wall, &
          rMult_CPU, rMult_Wall, Twoel2_CPU, Twoel2_Wall, Twoel3_CPU, Twoel3_Wall

end module TempTime

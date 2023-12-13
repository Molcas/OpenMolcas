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

function IColAt(NumAt)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IColAt
integer(kind=iwp), intent(in) :: NumAt
integer(kind=iwp), parameter :: ICol(0:108) = [ &
  0, &
  1,2, &
  1,2,3,4,5,6,7,8, &
  1,2,3,4,5,6,7,8, &
  1,2, &
    9,9,9,9,9,9,9,9,9,9, &
      3,4,5,6,7,8, &
  1,2, &
    9,9,9,9,9,9,9,9,9,9, &
      3,4,5,6,7,8, &
  1,2, &
    9,9,9,9,9,9,9,9,9,9,9,9,9,9, &
    9,9,9,9,9,9,9,9,9,9, &
      9,9,9,9,9,9, &
  9,9, &
    9,9,9,9,9,9,9,9,9,9,9,9,9,9, &
    9,9,9,9,9,9 &
]

IcolAt = ICol(NumAt)

return

end function IColAt

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

function Occ2Mrs(NO,IARRAY)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Occ2Mrs
integer(kind=iwp) :: NO, IARRAY(NO)
integer(kind=iwp) :: I, POW2

OCC2MRS = 0
POW2 = 1
do I=1,NO
  if (IARRAY(I) /= 0) OCC2MRS = OCC2MRS+POW2
  POW2 = 2*POW2
end do

return

end function Occ2Mrs

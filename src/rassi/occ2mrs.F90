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

integer function Occ2Mrs(NO,IARRAY)

implicit none
integer NO
integer IARRAY(NO)
integer I, POW2

OCC2MRS = 0
POW2 = 1
do I=1,NO
  if (IARRAY(I) /= 0) OCC2MRS = OCC2MRS+POW2
  POW2 = 2*POW2
end do

return

end function Occ2Mrs

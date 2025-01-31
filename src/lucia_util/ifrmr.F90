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

integer function IFRMR(iArray,IROFF,IELMNT)
! An integer array is stored in real array iArray,
! starting from iArray(IROFF). Obtain element
! IELMNT of this array

use lucia_data, only: irat

implicit none
integer iArray(*), IROFF, IELMNT
integer IIOFF

! offset when iArray is integer array
IIOFF = 1+IRAT*(IROFF-1)
IFRMR = iArray(IIOFF-1+IELMNT)

end function IFRMR

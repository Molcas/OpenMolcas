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

function IPHASE(IDRT,IUP,IWALK)
! PURPOSE: THE SYMMETRIC GROUP APPROACH AND THE UNITARY GROUP
!          APPROACH DIFFER IN THE PHASE CONVENTION. FIND THE
!          PHASE FACTOR RELATING THE CSFS IN EITHER BASIS.

implicit real*8(A-H,O-Z)
#include "gugx.fh"
dimension IDRT(NVERT,5), IUP(NVERT,0:3), IWALK(NLEV)

! FIND THE MIDVERTEX AND THE COMBINED WALK SYMMETRY

IPHASE = 1
IVERT = NVERT
do LEV=1,NLEV
  ICASE = IWALK(LEV)
  IVERT = IUP(IVERT,ICASE)
  if ((ICASE == 2) .or. (ICASE == 3)) then
    ISIGN = (-1)**IDRT(IVERT,4)
    IPHASE = IPHASE*ISIGN
  end if
end do

! EXIT

return

end function IPHASE

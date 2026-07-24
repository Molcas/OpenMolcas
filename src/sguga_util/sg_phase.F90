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

function SG_PHASE(SGS,IWALK)
! PURPOSE: THE SYMMETRIC GROUP APPROACH AND THE UNITARY GROUP
!          APPROACH DIFFER IN THE PHASE CONVENTION. FIND THE
!          PHASE FACTOR RELATING THE CSFS IN EITHER BASIS.

use sguga, only: SGStruct
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: SG_PHASE
type (SGStruct), intent(in):: SGS
integer(kind=iwp), intent(in) :: IWALK(SGS%NLEV)
integer(kind=iwp) :: ICASE, ISGN, IVERT, LEV

! FIND THE MIDVERTEX AND THE COMBINED WALK SYMMETRY

SG_PHASE = 1
IVERT = SGS%NVERT
do LEV=1,SGS%NLEV
  ICASE = IWALK(LEV)
  IVERT = SGS%UP(IVERT,ICASE)
  Select Case(iCase)
    Case(2,3)
      ISGN = (-1)**SGS%DRT(IVERT,4)
    Case Default
      ISGN = 1
  End Select
  SG_PHASE = SG_PHASE*ISGN
end do

end function SG_PHASE

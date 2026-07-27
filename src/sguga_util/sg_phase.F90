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

<<<<<<< HEAD:src/molcas_ci_util/iphase.F90
function IPHASE(NLEV,NVERT,IDRT,IUP,IWALK)
=======
function SG_PHASE(SGS,IWALK)
>>>>>>> upstream-openmolcas/master:src/sguga_util/sg_phase.F90
! PURPOSE: THE SYMMETRIC GROUP APPROACH AND THE UNITARY GROUP
!          APPROACH DIFFER IN THE PHASE CONVENTION. FIND THE
!          PHASE FACTOR RELATING THE CSFS IN EITHER BASIS.

use sguga, only: SGStruct
use Definitions, only: iwp

implicit none
<<<<<<< HEAD:src/molcas_ci_util/iphase.F90
integer(kind=iwp) :: IPHASE
integer(kind=iwp), intent(in) :: NLEV, NVERT, IDRT(NVERT,5), IUP(NVERT,0:3), IWALK(NLEV)
=======
integer(kind=iwp) :: SG_PHASE
type (SGStruct), intent(in):: SGS
integer(kind=iwp), intent(in) :: IWALK(SGS%NLEV)
>>>>>>> upstream-openmolcas/master:src/sguga_util/sg_phase.F90
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

<<<<<<< HEAD:src/molcas_ci_util/iphase.F90
! EXIT

end function IPHASE
=======
end function SG_PHASE
>>>>>>> upstream-openmolcas/master:src/sguga_util/sg_phase.F90

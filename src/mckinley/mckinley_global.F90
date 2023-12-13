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

module McKinley_global

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: SCF = 1, RASSCF = 2, &
                                nOneel = 1, nTwoel = 2, nTwoDens = 3, nIntegrals = 4, nScreen = 5, nTrans = 6, nFckAcc = 7, &
                                nMOTrans = 8, nTotal = 9
integer(kind=iwp) :: nFck(0:7), nMethod
real(kind=wp) :: CPUStat(nTotal)
logical(kind=iwp) :: lGrd, lHss, Nona, PreScr, sIrrep
integer(kind=iwp), allocatable :: ipDisp(:), ipDisp2(:), ipDisp3(:), ipMO(:)

public :: CPUStat, ipDisp, ipDisp2, ipDisp3, ipMO, lGrd, lHss, nFck, nFckAcc, nIntegrals, nMethod, nMOTrans, Nona, nOneel, &
          nScreen, nTotal, nTrans, nTwoDens, nTwoel, PreScr, RASSCF, SCF, sIrrep

end module McKinley_global

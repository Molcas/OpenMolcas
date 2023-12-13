!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Roland Lindh                                     *
!               2000, Valera Veryazov                                  *
!               2014, Thomas Dresselhaus                               *
!***********************************************************************

subroutine MOEvalDel(MOValueD,nMOs,nCoor,CCoor,CMOs,nCMO,DoIt)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: nDrv = 1, mAO = 4
integer(kind=iwp), intent(in) :: nMOs, nCoor, nCMO, DoIt(nMOs)
real(kind=wp), intent(out) :: MOValueD(mAO,nCoor,nMOs)
real(kind=wp), intent(in) :: CCoor(3,nCoor), CMOs(nCMO)

call MOEval(MOValueD,nMOs,nCoor,CCoor,CMOs,nCMO,DoIt,nDrv,mAO)

return

end subroutine MOEvalDel

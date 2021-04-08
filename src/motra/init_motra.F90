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

subroutine init_motra()

use motra_global, only: Debug, FnHalf, FnInpOrb, FnJobIph, FnOneMO, FnTwoAO, FnTwoMO, iAutoCut, iOneOnly, iPrint, iRFpert, &
                        iVecTyp, LuHalf, LuInpOrb, LuJobIph, LuOneMO, LuTwoAO, LuTwoMO

implicit none

!----------------------------------------------------------------------*
! Define  file names and unit numbers.                                 *
!----------------------------------------------------------------------*
FnInpOrb = 'INPORB'
FnJobIph = 'JOBIPH'
FnTwoAO = 'ORDINT'
FnOneMO = 'TRAONE'
FnTwoMO = 'TRAINT'
FnHalf = 'TEMP1'
Debug = 0
iPrint = 0
iOneOnly = 0
iVecTyp = 2
iAutoCut = 0
iRFpert = 0
LuInpOrb = 10
LuJobIph = 15
LuTwoAO = 40
LuOneMO = 30
LuTwoMO = 50
LuHalf = 60

return

end subroutine init_motra

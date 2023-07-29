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

function Cho_P_GetmPass(iLoc)

use Cholesky, only: Cho_Real_Par, nnBstRT, nnBstRT_G
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Cho_P_GetmPass
integer(kind=iwp) :: iLoc

if (Cho_Real_Par) then
  Cho_P_GetmPass = nnBstRT_G(iLoc)
else
  Cho_P_GetmPass = nnBstRT(iLoc)
end if

end function Cho_P_GetmPass

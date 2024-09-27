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
! Copyright (C) 1990, IBM                                              *
!***********************************************************************

#define _CHECK_
function iPrmt(jOper,iChct)
!***********************************************************************
!     Returns the phase factor of a basis function under a symmetry    *
!     operation, jOper. iChct contains the information about the       *
!     character of the basis function.                                 *
!***********************************************************************

use Symmetry_Info, only: iOper
use Definitions, only: iwp
#ifdef _CHECK_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: iPrmt
integer(kind=iwp), intent(in) :: jOper, iChct
integer(kind=iwp) :: i, iCom

#ifdef _CHECK_
if (size(iOper) < 1) then
  write(u6,*) 'iPrmt; iOper not defined.'
  call Abend()
end if
#endif
iPrmt = 1
iCom = iand(iOper(jOper),iChct)
do i=1,3
  if (btest(iCom,i-1)) iPrmt = -iPrmt
end do

return

end function iPrmt

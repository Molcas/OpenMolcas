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
!               1991, Roland Lindh                                     *
!***********************************************************************

function iPrmt_(iCom)
!***********************************************************************
!     Returns the phase factor of a basis function under a symmetry    *
!     operation, jOper. iChct contains the information about the       *
!     character of the basis function.                                 *
!***********************************************************************

implicit none
integer iPrmt_
integer iCom
integer i

iPrmt_ = 1
do i=1,3
  if (iand(iCom,2**(i-1)) /= 0) iPrmt_ = iPrmt_*(-1)
end do

return

end function iPrmt_

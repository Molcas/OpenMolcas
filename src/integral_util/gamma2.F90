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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

function Gamma2(m,T)
!***********************************************************************
!                                                                      *
! Object: to compute the auxiliary function in the high argument       *
!         approximation.                                               *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One, Two

implicit none
real*8 Gamma2
integer m
real*8 T
integer i

Gamma2 = sqrt(Two*acos(Zero)/T)/Two
do i=1,m
  Gamma2 = ((Two*dble(i)-One)/(Two*T))*Gamma2
end do

return

end function Gamma2

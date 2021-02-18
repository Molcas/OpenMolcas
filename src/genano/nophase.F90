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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine NOphase(a,n)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: a(n,n)
integer(kind=iwp) :: i
real(kind=wp) :: sgn

sgn = One
do i=1,n
  if (sgn*a(1,i) < Zero) a(:,i) = -a(:,i)
  sgn = -sgn
end do

return

end subroutine NOphase

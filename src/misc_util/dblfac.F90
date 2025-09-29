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
! Copyright (C) 1990,2020, Roland Lindh                                *
!***********************************************************************

function DblFac(n)
!***********************************************************************
!                                                                      *
! Object: to compute the double factorial of n.                        *
!                                                                      *
! Called from: NrmSph                                                  *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: DblFac
integer(kind=iwp), intent(in) :: n
integer(kind=iwp) :: i
integer(kind=iwp), parameter :: cache_len = 10
real(kind=wp), parameter :: cache(cache_len) = [1.0_wp,2.0_wp,3.0_wp,8.0_wp,15.0_wp,48.0_wp,105.0_wp,384.0_wp,945.0_wp,3840.0_wp]

! ideally, if n < 0 this is undefined, but in molcas it is expected to be 1
if (n < 1) then
  DblFac = One
else if (n <= cache_len) then
  ! use cached values whenever possible
  DblFac = cache(n)
else
  DblFac = One
  do i=n,2,-2
    DblFac = DblFac*real(i,kind=wp)
  end do
end if

return

end function DblFac

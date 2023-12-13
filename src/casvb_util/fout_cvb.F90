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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine fout_cvb(f,a1,a2)

use Definitions, only: wp, u6

implicit none
real(kind=wp), intent(in) :: f
character(len=*), intent(in) :: a1, a2
character(len=46) :: b2
character(len=15) :: b1
character(len=12) :: b3
real(kind=wp), parameter :: hge = 1.0e20_wp

b1 = a1
b2 = a2
if (abs(f) /= hge) then
  write(b3,'(es12.4)') f
else
  b3 = '    Disabled'
end if
write(u6,'(1x,3a)') b1,b2,b3

return

end subroutine fout_cvb

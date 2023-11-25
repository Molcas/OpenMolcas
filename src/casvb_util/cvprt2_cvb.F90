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

subroutine cvprt2_cvb(a1,f1,f2,ic)

use casvb_global, only: formcvp
use Definitions, only: wp, iwp, u6

implicit none
character(len=16), intent(in) :: a1
real(kind=wp), intent(in) :: f1, f2
integer(kind=iwp), intent(in) :: ic
real(kind=wp), parameter :: hge = 1.0e20_wp

if (abs(f2) /= hge) then
  if ((ic == 1) .and. (f1 < f2)) then
    write(u6,formcvp) a1,f1,'     smaller than',f2
  else if (ic == 1) then
    write(u6,formcvp) a1,f1,' not smaller than',f2
  else if ((ic == 2) .and. (f1 > f2)) then
    write(u6,formcvp) a1,f1,'     greater than',f2
  else if (ic == 2) then
    write(u6,formcvp) a1,f1,' not greater than',f2
  end if
end if

return

end subroutine cvprt2_cvb

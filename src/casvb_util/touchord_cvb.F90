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

subroutine touchord_cvb(itouch,iorder,n)

implicit real*8(a-h,o-z)
dimension iorder(n)
logical debug
data debug/.false./

if (debug) then
  write(6,*) ' itouch :',itouch
  write(6,*) ' iorder on entry :',iorder
end if
itchord = iorder(itouch)
do i=1,n
  if ((iorder(i) < itchord) .and. (iorder(i) /= 0)) iorder(i) = iorder(i)+1
end do
iorder(itouch) = 1
if (debug) write(6,*) ' iorder on exit  :',iorder

return

end subroutine touchord_cvb

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

subroutine untouch_cvb(chr)

use casvb_global, only: charobj, iprint, mustdeclare, nobj, up2date
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: chr
integer(kind=iwp) :: i, iobj

do
  iobj = 0
  do i=1,nobj
    if (charobj(i) == chr) iobj = i
  end do
  if (iobj /= 0) exit
  if (mustdeclare) then
    write(u6,*) ' Make object not found :',chr
    call abend_cvb()
  end if
  call decl_cvb(chr)
end do
if (up2date(iobj)) return
if (iprint >= 1) write(u6,'(/,a,i3,2a)') ' Untouch object no.',iobj,', name : ',charobj(iobj)
up2date(iobj) = .true.

return

end subroutine untouch_cvb

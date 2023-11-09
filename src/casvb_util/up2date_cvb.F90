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

function up2date_cvb(chr)

use casvb_global, only: charobj, nobj, up2date
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp) :: up2date_cvb
character(len=*), intent(in) :: chr
integer(kind=iwp) :: i, iobj

iobj = 0
do i=1,nobj
  if (charobj(i) == chr) iobj = i
end do
if (iobj == 0) then
  write(u6,*) ' Make object not found :',chr
  call abend_cvb()
end if
up2date_cvb = up2date(iobj)

return

end function up2date_cvb

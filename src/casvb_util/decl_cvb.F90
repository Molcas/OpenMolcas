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

subroutine decl_cvb(chr)

use casvb_global, only: charobj, ioffs, iprint, joffs, mxobj, nobj, up2date
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: chr
integer(kind=iwp) :: i, ii, iobj

iobj = 0
do i=1,nobj
  if (charobj(i) == chr) iobj = i
end do
if (iobj > 0) then
  if (iprint > 1) write(u6,*) ' Make object exists already :',chr
  return
end if
nobj = nobj+1
if (nobj > mxobj) then
  write(u6,*) ' Too many make objects, max :',mxobj
  call abend_cvb()
end if
charobj(nobj) = chr
up2date(nobj) = .false.
ioffs(nobj+1) = ioffs(nobj)
joffs(nobj+1) = joffs(nobj)
if (iprint >= 10) then
  write(u6,*) ' IOFFS :',(ioffs(ii),ii=1,nobj+1)
  write(u6,*) ' JOFFS :',(joffs(ii),ii=1,nobj+1)
end if

return

end subroutine decl_cvb

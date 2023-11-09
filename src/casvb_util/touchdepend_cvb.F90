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

subroutine touchdepend_cvb(chr1,chr2)

use casvb_global, only: charobj, j_dep_on_i, joffs, mxdep, ndep_ji, nobj, up2date
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: chr1, chr2
integer(kind=iwp) :: i, iobj, jobj

call undepend2_cvb(chr1,chr2,2)

iobj = 0
jobj = 0
do i=1,nobj
  if (charobj(i) == chr1) iobj = i
  if (charobj(i) == chr2) jobj = i
end do
if (iobj == 0) then
  write(u6,*) ' Make object not found :',chr1
  call abend_cvb()
end if
if (jobj == 0) then
  write(u6,*) ' Make object not found :',chr2
  call abend_cvb()
end if
ndep_ji = ndep_ji+1
if (ndep_ji > mxdep) then
  write(u6,*) ' Too many make dependencies, max :',mxdep
  call abend_cvb()
end if
j_dep_on_i(joffs(jobj+1)+2:joffs(nobj+1)+1) = j_dep_on_i(joffs(jobj+1)+1:joffs(nobj+1))
j_dep_on_i(joffs(jobj+1)+1) = iobj
joffs(jobj+1:nobj+1) = joffs(jobj+1:nobj+1)+1

if (.not. up2date(jobj)) up2date(iobj) = .false.

return

end subroutine touchdepend_cvb

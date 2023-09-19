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

subroutine touch_cvb(chr)

implicit real*8(a-h,o-z)
character*(*) chr
#include "make_cvb.fh"

do
  iobj = 0
  do i=1,nobj
    if (charobj(i) == chr) iobj = i
  end do
  if (iobj /= 0) exit
  if (mustdeclare) then
    write(6,*) ' Make object not found :',chr
    call abend_cvb()
  end if
  call decl_cvb(chr)
end do
up2date(iobj) = .false.
if (iprint >= 1) write(6,'(/,a,i3,2a)') ' Touch (1) of object no.',iobj,', name : ',charobj(iobj)

! Mark all "child" objects as out-of-date:
do
  n_touched = 0
  do iobj=1,nobj
    if (.not. up2date(iobj)) then
      do i=joffs(iobj)+1,joffs(iobj+1)
        call touchrules_cvb(charobj(j_dep_on_i(i)))
        if (up2date(j_dep_on_i(i))) then
          up2date(j_dep_on_i(i)) = .false.
          if (iprint >= 1) write(6,'(/,a,i3,2a)') ' Touch (2) of object no.',j_dep_on_i(i),', name : ',charobj(j_dep_on_i(i))
          n_touched = n_touched+1
        end if
      end do
    end if
  end do
  if (n_touched == 0) exit
end do

return

end subroutine touch_cvb

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

subroutine make_cvb(chr)

use casvb_global, only: charobj, i_dep_on_j, ioffs, iprint, mustdeclare, nobj, up2date
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: chr
integer(kind=iwp) :: i, iobj, mkobj
logical(kind=iwp) :: done

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

! Make sure all "parent" objects are up-to-date:
do
  mkobj = iobj
  do
    done = .false.
    do i=ioffs(mkobj)+1,ioffs(mkobj+1)
      if (.not. up2date(i_dep_on_j(i))) then
        mkobj = i_dep_on_j(i)
        done = .true.
        exit
      end if
    end do
    if (.not. done) exit
  end do
  if (.not. up2date(mkobj)) then
    if (iprint >= 1) write(u6,'(/,a,i3,2a)') ' Making object no.',mkobj,', name : ',charobj(mkobj)
    call rules_cvb(charobj(mkobj))
    up2date(mkobj) = .true.
  end if
  if (mkobj == iobj) exit
end do

return

end subroutine make_cvb
